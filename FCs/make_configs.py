#!/usr/bin/env python3
"""Generate per-chunk MASH configs from a template JSON.

Replaces the old per-run sed/.py scripts: given a template (e.g. feExample2.json
for fixed-effects, reExample2.json for random-effects/site), it writes one config
per chunk of the phenotype columns o{start}..o{end} into <outdir>/<Method>.<Kind>.<i>.json.

The `out` path of each generated config is derived from the template's `out`
directory: <template out dir>/<prefix>.<Method>.<Kind>.<i>. MASH appends .csv/.log.

Usage:
  # write config for a single chunk
  python make_configs.py --template reExample2.json --method AdjHE --kind RE \
      --chunk 208 --total 61776 --prefix pconns --index 5

  # write configs for all chunks (prints the number of chunks)
  python make_configs.py --template reExample2.json --method AdjHE --kind RE \
      --chunk 208 --total 61776 --prefix pconns --all
"""

import argparse
import json
import os


def chunks(total, size):
    return [(s, min(s + size, total)) for s in range(0, total, size)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, help="Template JSON (feExample2.json or reExample2.json)")
    parser.add_argument("--method", required=True, help="Estimation method (AdjHE, GCTA, HEreg, ...)")
    parser.add_argument("--kind", required=True, choices=["FE", "RE"], help="Fixed-effects or random-effects/site run")
    parser.add_argument("--chunk", type=int, required=True, help="Number of phenotypes per config")
    parser.add_argument("--total", type=int, required=True, help="Total number of phenotype columns (o0..o(total-1))")
    parser.add_argument("--prefix", required=True, help="Output basename prefix (e.g. pconns)")
    parser.add_argument("--outdir", default="temp", help="Directory for generated configs")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--index", type=int, help="Write config for this chunk only")
    group.add_argument("--all", action="store_true", help="Write configs for all chunks")
    args = parser.parse_args()

    with open(args.template) as f:
        cfg = json.load(f)

    base = os.path.dirname(cfg["out"])
    os.makedirs(args.outdir, exist_ok=True)
    cs = chunks(args.total, args.chunk)

    def write_one(i, start, end):
        d = dict(cfg)
        d["mpheno"] = [f"o{j}" for j in range(start, end)]
        d["Method"] = args.method
        d["out"] = f"{base}/{args.prefix}.{args.method}.{args.kind}.{i}"
        path = os.path.join(args.outdir, f"{args.method}.{args.kind}.{i}.json")
        with open(path, "w") as f:
            json.dump(d, f, indent=2)

    if args.all:
        for i, (start, end) in enumerate(cs):
            write_one(i, start, end)
        print(len(cs))
    else:
        start, end = cs[args.index]
        write_one(args.index, start, end)
        print(os.path.join(args.outdir, f"{args.method}.{args.kind}.{args.index}.json"))


if __name__ == "__main__":
    main()