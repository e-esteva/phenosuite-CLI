"""
Gemma inference backend for the PhenoSuite mIF Cell Phenotyper.

Generation-agnostic: the architecture class is read from the checkpoint's own
config.architectures and looked up on the transformers module, so Gemma 3,
Gemma 4 (incl. the "unified" variants) and later releases load without a code
change, provided the installed transformers knows the architecture. Whether a
checkpoint is multimodal is decided structurally (does it carry a vision/audio
tower?) rather than by matching a model_type string.

Only text is ever sent — multimodal checkpoints are loaded with their full
architecture (loading one through a text-only class mis-resolves it and leaves
the vision tower unwired) but are prompted text-only.

Two model formats are supported:
  1. Merged checkpoint (primary): --model_path /path/to/merged/
  2. LoRA adapter: --adapter_path /path/to/adapter/ [--base_model google/gemma-3-1b-it]

Called from R via system2(); progress reported to stderr; results written as
newline-delimited JSON to --output_jsonl so the R side can stream_in().
"""
import json
import re
import sys
import argparse

import torch
import transformers
from transformers import (
    AutoConfig,
    AutoTokenizer,
    AutoProcessor,
    AutoModelForCausalLM,
    BitsAndBytesConfig,
)


# Architecture resolution is deliberately generation-agnostic.
#
# Every Gemma generation so far has changed BOTH the config model_type and the
# architecture class name:
#   gemma-3-1b        text-only      AutoModelForCausalLM
#   gemma-3-12b       gemma3         Gemma3ForConditionalGeneration
#   gemma-4-E2B/31B   gemma4         Gemma4ForCausalLM / Gemma4ForConditionalGeneration
#   gemma-4-12B       gemma4_unified Gemma4UnifiedForConditionalGeneration
#
# Hardcoding any of those means the next release needs a code change, so instead
# the class named in the checkpoint's own config.architectures is looked up on
# the transformers module. That works for any generation the installed
# transformers already knows about, including ones released after this file.
def _resolve_config(path_or_id):
    return AutoConfig.from_pretrained(path_or_id)


def _is_multimodal(cfg):
    # Structural test rather than a model_type string match: a checkpoint that
    # carries a vision/audio tower needs the processor + conditional-generation
    # class, whatever the generation happens to call itself.
    return any(hasattr(cfg, attr) for attr in
               ("vision_config", "audio_config", "text_config"))


def _resolve_model_class(cfg, multimodal):
    for name in (getattr(cfg, "architectures", None) or []):
        cls = getattr(transformers, name, None)
        if cls is not None:
            print(f"[infer.py] Using architecture {name} from config.",
                  file=sys.stderr, flush=True)
            return cls
        print(f"[infer.py] config names architecture '{name}', which this "
              f"transformers ({transformers.__version__}) does not provide — "
              f"falling back to an Auto class. Upgrading transformers is the fix "
              f"if loading fails.", file=sys.stderr, flush=True)

    if multimodal:
        # AutoModelForMultimodalLM is newer; older transformers only has the
        # vision-text auto-class. Try in order of specificity.
        for auto_name in ("AutoModelForMultimodalLM",
                          "AutoModelForImageTextToText",
                          "AutoModelForVision2Seq"):
            cls = getattr(transformers, auto_name, None)
            if cls is not None:
                print(f"[infer.py] Falling back to {auto_name}.",
                      file=sys.stderr, flush=True)
                return cls
    return AutoModelForCausalLM


def _build_kwargs(load_in_8bit):
    kwargs = {"device_map": "auto"}
    if load_in_8bit:
        kwargs["quantization_config"] = BitsAndBytesConfig(load_in_8bit=True)
    else:
        kwargs["torch_dtype"] = torch.bfloat16
    return kwargs


def load_merged(model_path, load_in_8bit, multimodal, model_cls):
    kwargs = _build_kwargs(load_in_8bit)
    processor_cls = AutoProcessor if multimodal else AutoTokenizer
    processor = processor_cls.from_pretrained(model_path)
    model = model_cls.from_pretrained(model_path, **kwargs)
    model.eval()
    return processor, model


def load_adapter(adapter_path, base_model, load_in_8bit, multimodal, model_cls):
    from peft import PeftModel
    kwargs = _build_kwargs(load_in_8bit)
    processor_cls = AutoProcessor if multimodal else AutoTokenizer
    try:
        processor = processor_cls.from_pretrained(adapter_path)
    except OSError:
        # Adapter dir didn't ship its own tokenizer/processor files (e.g. no
        # new tokens were added during fine-tuning) — fall back to the base's.
        processor = processor_cls.from_pretrained(base_model)
    base = model_cls.from_pretrained(base_model, **kwargs)
    model = PeftModel.from_pretrained(base, adapter_path)
    model.eval()
    return processor, model


def _as_multimodal_messages(messages):
    # Gemma3Processor's chat template expects each content field to be a list
    # of typed parts; the R side only ever emits plain-string content.
    return [
        {**m, "content": m["content"] if isinstance(m["content"], list)
                          else [{"type": "text", "text": m["content"]}]}
        for m in messages
    ]


def _apply_template(processor, messages, multimodal):
    # enable_thinking is a Gemma 4+ chat-template variable. Reasoning traces
    # would be emitted before the JSON and break parsing, so it is pinned off
    # rather than left to the template's default. Older templates (Gemma 3)
    # don't accept the kwarg, hence the fallback.
    common = dict(add_generation_prompt=True, return_tensors="pt")
    if multimodal:
        msgs = _as_multimodal_messages(messages)
        common.update(tokenize=True, return_dict=True)
    else:
        msgs = messages
    try:
        return processor.apply_chat_template(msgs, enable_thinking=False, **common)
    except TypeError:
        return processor.apply_chat_template(msgs, **common)


def _generate(messages, processor, model, multimodal, max_new_tokens=128,
              temperature=0.0):
    if multimodal:
        inputs = _apply_template(processor, messages, multimodal).to(model.device)
        pad_token_id = processor.tokenizer.eos_token_id
    else:
        inputs = {
            "input_ids": _apply_template(processor, messages, multimodal).to(model.device)
        }
        pad_token_id = processor.eos_token_id

    input_len = inputs["input_ids"].shape[-1]

    gen_kw = dict(max_new_tokens=max_new_tokens, pad_token_id=pad_token_id)
    if temperature and temperature > 0:
        gen_kw.update(do_sample=True, temperature=float(temperature))
    else:
        gen_kw.update(do_sample=False)   # greedy: reproducible

    with torch.no_grad():
        out = model.generate(**inputs, **gen_kw)

    return processor.decode(out[0][input_len:], skip_special_tokens=True).strip()


# Reasoning-trace wrappers to discard before parsing. Gemma 4 emits
# "<|channel>thought ... <channel|>" when thinking is enabled; <think> tags are
# the common convention elsewhere. Belt-and-braces: thinking is already pinned
# off in _apply_template, but a template default change shouldn't silently turn
# every prediction into a parse failure.
_THINK_RE = [
    re.compile(r"<\|?channel\|?>\s*thought\b.*?<\|?/?channel\|?>", re.DOTALL | re.I),
    re.compile(r"<think>.*?</think>", re.DOTALL | re.I),
]


def _extract_json(text):
    """Best-effort JSON object out of a model reply. Returns dict or None."""
    s = text.strip()
    for pat in _THINK_RE:
        s = pat.sub("", s)
    s = re.sub(r"^```[A-Za-z]*\n?|\n?```$", "", s.strip())

    try:
        obj = json.loads(s)
        return obj if isinstance(obj, dict) else None
    except json.JSONDecodeError:
        pass

    # Fall back to the first balanced {...} span, so leading prose ("Sure, here
    # is the JSON:") or a trailing reasoning tail doesn't cost us the answer.
    start = s.find("{")
    while start != -1:
        depth, in_str, esc = 0, False, False
        for i in range(start, len(s)):
            c = s[i]
            if in_str:
                if esc:      esc = False
                elif c == "\\": esc = True
                elif c == '"':  in_str = False
            elif c == '"':   in_str = True
            elif c == "{":   depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    try:
                        obj = json.loads(s[start:i + 1])
                        if isinstance(obj, dict):
                            return obj
                    except json.JSONDecodeError:
                        pass
                    break
        start = s.find("{", start + 1)
    return None


# Kept as a thin alias so the name still means something if called elsewhere.
def _strip_fences(text):
    return re.sub(r"^```[A-Za-z]*\n?|\n?```$", "", text.strip())


def _infer_one(ex, processor, model, multimodal, max_new_tokens=128,
               temperature=0.0):
    text = _generate(ex["messages"], processor, model, multimodal,
                     max_new_tokens=max_new_tokens, temperature=temperature)

    pred = _extract_json(text)
    if pred is None or "cell_type" not in pred:
        pred = {
            "cell_type":  "unknown",
            "confidence": 0.0,
            "error":      "parse_failure",
            "raw":        text,
        }

    pred["cluster_id"] = ex.get("cluster_id")
    pred["cell_index"] = ex.get("cell_index")
    return pred


# Mirrors harmonize_ollama_labels() in gemma-utils.R, but runs in this process
# so the (large) model is loaded exactly once for both passes — a separate
# invocation would pay the full load cost again just to map a handful of labels.
def harmonize_labels(labels, processor, model, multimodal, max_new_tokens=None):
    unique = sorted({l for l in labels if l})
    if len(unique) <= 1:
        return None

    prompt = (
        "Below is a list of cell type annotation labels produced by automated "
        "immunophenotyping. Many entries are near-duplicates that differ only in "
        "pluralisation, spacing, capitalisation, or minor wording (e.g. 'CD4 T cell' "
        "and 'CD4 T cells', or 'NK cell' and 'Natural Killer cell'). Return a JSON "
        "object mapping EVERY label in the input list to its canonical standardised "
        "form. Use Title Case. Do not merge biologically distinct populations. "
        "Return ONLY the JSON object, no markdown fences, no explanation.\n\n"
        "Labels: " + ", ".join(unique)
    )
    messages = [
        {"role": "system",
         "content": "You are an expert immunologist. Return ONLY valid JSON, "
                    "no markdown, no explanation."},
        {"role": "user", "content": prompt},
    ]

    # Scale the budget with the label count — the reply restates every label.
    budget = max_new_tokens or max(256, 48 * len(unique))
    text = _generate(messages, processor, model, multimodal, max_new_tokens=budget)

    mapping = _extract_json(text)
    if not isinstance(mapping, dict) or not mapping:
        print("[infer.py] harmonisation returned no usable mapping; keeping raw labels",
              file=sys.stderr, flush=True)
        return None
    return {k: v for k, v in mapping.items() if isinstance(v, str)}


def run_inference(args):
    print("[infer.py] Loading model…", file=sys.stderr, flush=True)

    source     = args.model_path or args.base_model
    cfg        = _resolve_config(source)
    multimodal = _is_multimodal(cfg)
    model_cls  = _resolve_model_class(cfg, multimodal)
    print(f"[infer.py] model_type={getattr(cfg, 'model_type', '?')} "
          f"multimodal={multimodal}", file=sys.stderr, flush=True)

    if args.model_path:
        processor, model = load_merged(
            args.model_path, args.load_in_8bit, multimodal, model_cls
        )
    else:
        processor, model = load_adapter(
            args.adapter_path, args.base_model, args.load_in_8bit,
            multimodal, model_cls
        )

    print("[infer.py] Model ready. Running inference…", file=sys.stderr, flush=True)

    examples, preds = [], []
    with open(args.input_jsonl) as fin:
        for line in fin:
            line = line.strip()
            if line:
                examples.append(json.loads(line))

    for idx, ex in enumerate(examples):
        pred = _infer_one(ex, processor, model, multimodal,
                          max_new_tokens=args.max_new_tokens)
        preds.append(pred)
        print(f"[infer.py] {idx + 1} done", file=sys.stderr, flush=True)

    print("[infer.py] Inference complete.", file=sys.stderr, flush=True)

    # Retry the ones that produced no usable label. Done here, in-process, while
    # the model is still resident — a separate invocation would pay the full
    # (multi-minute, for a 12B) load again just to redo a handful of prompts.
    for attempt in range(1, args.retry + 1):
        failed = [i for i, p in enumerate(preds) if p.get("error")]
        if not failed:
            break
        print(f"[infer.py] Retry {attempt}/{args.retry}: {len(failed)} failed "
              f"prediction(s)…", file=sys.stderr, flush=True)
        for i in failed:
            preds[i] = _infer_one(examples[i], processor, model, multimodal,
                                  max_new_tokens=args.max_new_tokens)
        recovered = len(failed) - sum(1 for i in failed if preds[i].get("error"))
        print(f"[infer.py]   recovered {recovered}/{len(failed)}",
              file=sys.stderr, flush=True)

    remaining = sum(1 for p in preds if p.get("error"))
    if remaining:
        print(f"[infer.py] {remaining} prediction(s) still failed after retries.",
              file=sys.stderr, flush=True)

    with open(args.output_jsonl, "w") as fout:
        for pred in preds:
            fout.write(json.dumps(pred) + "\n")

    labels = [p.get("cell_type") for p in preds]

    # N-sample vote on clusters that came back "Unclassified". Runs here so the
    # checkpoint is loaded once, exactly like the retry and harmonisation passes.
    #
    # Mirrors vote_unclassified_ollama() in gemma-utils.R, including the rule
    # that a cluster whose prompt says "No markers are elevated" is re-sampled
    # for diagnosis but never overridden: a majority label there would be the
    # model inventing a lineage from nothing, which is what the min_z floor
    # exists to prevent.
    if args.vote_samples and args.vote_samples > 1:
        uncl = [i for i, p in enumerate(preds)
                if "unclassif" in str(p.get("cell_type", "")).lower()]
        if uncl:
            print(f"[infer.py] Voting on {len(uncl)} Unclassified cluster(s), "
                  f"{args.vote_samples} samples at temperature {args.vote_temperature}…",
                  file=sys.stderr, flush=True)
        resolved = 0
        for i in uncl:
            ex = examples[i]
            floor_fired = "No markers are elevated" in ex["messages"][1]["content"]
            votes = []
            for _ in range(args.vote_samples):
                v = _infer_one(ex, processor, model, multimodal,
                               max_new_tokens=args.max_new_tokens,
                               temperature=args.vote_temperature)
                if not v.get("error"):
                    votes.append(str(v.get("cell_type", "")))
            if not votes:
                continue
            norm = lambda x: re.sub(r"\s+", " ",
                                    re.sub(r"[^a-z0-9 ]", " ", x.lower())).strip()
            tally = {}
            for v in votes:
                tally[norm(v)] = tally.get(norm(v), 0) + 1
            top, cnt = max(tally.items(), key=lambda kv: kv[1])
            agree = cnt / len(votes)
            disp = next(v for v in votes if norm(v) == top)

            preds[i]["vote_top"]       = disp
            preds[i]["vote_agreement"] = round(agree, 3)
            preds[i]["vote_n"]         = len(votes)
            preds[i]["vote_basis"]     = ("no_elevated_markers" if floor_fired
                                          else "model_declined")
            preds[i]["vote_detail"]    = "; ".join(f"{k}:{v}" for k, v in
                                                   sorted(tally.items(), key=lambda kv: -kv[1]))
            if ((not floor_fired or args.vote_override_floor)
                    and "unclassif" not in top
                    and agree >= args.vote_min_agreement):
                preds[i]["cell_type"]  = disp
                preds[i]["confidence"] = round(agree, 3)   # empirical, not self-reported
                resolved += 1
        if uncl:
            print(f"[infer.py]   resolved {resolved}/{len(uncl)}",
                  file=sys.stderr, flush=True)

        with open(args.output_jsonl, "w") as fout:
            for pred in preds:
                fout.write(json.dumps(pred) + "\n")
        labels = [p.get("cell_type") for p in preds]

    if args.harmonize_out:
        print("[infer.py] Harmonising labels…", file=sys.stderr, flush=True)
        mapping = harmonize_labels(labels, processor, model, multimodal)
        # Always write the file, even when harmonisation yielded nothing usable:
        # the R side treats an empty object as "no remapping" and carries on with
        # the raw labels, rather than failing on a missing file.
        with open(args.harmonize_out, "w") as fh:
            json.dump(mapping or {}, fh)
        if mapping:
            print(f"[infer.py] Harmonised {len(mapping)} label(s) -> "
                  f"{len(set(mapping.values()))} canonical form(s).",
                  file=sys.stderr, flush=True)


def main():
    p = argparse.ArgumentParser(
        description="Gemma mIF phenotyper — inference backend. Generation-agnostic: "
                     "resolves the architecture from the checkpoint's own config, so "
                     "Gemma 3 and Gemma 4 (incl. the unified 12B) both work given a "
                     "new enough transformers (>=4.50 / >=5.10 respectively)."
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--model_path",   help="Merged checkpoint directory")
    src.add_argument("--adapter_path", help="LoRA adapter directory")

    p.add_argument("--base_model", default="google/gemma-3-1b-it",
                   help="Base model ID or path (used with --adapter_path)")
    p.add_argument("--input_jsonl",  required=True, help="Input JSONL from R")
    p.add_argument("--output_jsonl", required=True, help="Output predictions JSONL")
    p.add_argument("--load_in_8bit", action="store_true",
                   help="Load in 8-bit (requires bitsandbytes)")
    p.add_argument("--retry", type=int, default=1,
                   help="Extra passes over predictions that produced no usable "
                        "label (request or parse failure). Runs in-process so "
                        "the model is not reloaded. 0 disables. [default: 1]")
    p.add_argument("--max_new_tokens", type=int, default=128,
                   help="Generation budget per prediction. 128 is ample for the "
                        "JSON answer; raise it if a model emits reasoning traces "
                        "before the answer (Gemma 4 thinking mode is pinned off, "
                        "so this should not normally be needed).")
    p.add_argument("--vote_samples", type=int, default=0,
                   help="Re-sample each Unclassified cluster this many times and "
                        "use the agreement rate as an empirical confidence. 0/1 "
                        "disables. Runs in-process, so the model is not reloaded.")
    p.add_argument("--vote_temperature", type=float, default=0.7,
                   help="Sampling temperature for the vote. Must be > 0 or every "
                        "sample is identical (the main pass is greedy). [0.7]")
    p.add_argument("--vote_min_agreement", type=float, default=0.6,
                   help="Fraction of votes one label needs before it replaces "
                        "Unclassified. [0.6]")
    p.add_argument("--vote_override_floor", action="store_true",
                   help="Also allow voting to override clusters that had NO "
                        "elevated markers. Off by default: a majority there is "
                        "the model guessing from nothing.")
    p.add_argument("--harmonize_out",
                   help="If set, run a label-harmonisation pass after inference "
                        "(reusing the already-loaded model) and write the raw -> "
                        "canonical mapping as JSON to this path.")

    run_inference(p.parse_args())


if __name__ == "__main__":
    main()
