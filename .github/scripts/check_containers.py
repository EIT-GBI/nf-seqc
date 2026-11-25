#!/usr/bin/env python3
"""Fail if any process resolves to an unpinned container image.

Reads `nextflow inspect` JSON on stdin. An image is unpinned if it is tagged
`:latest` or carries no tag at all, either of which makes a run
non-reproducible: the same pipeline revision can silently get different
software on a later day.
"""
import json
import sys


def is_unpinned(image: str) -> str | None:
    if not image:
        return "no container declared"
    # the tag lives after the last ':' in the final path segment, so that a
    # registry host with a port (host:5000/img) is not mistaken for a tag
    last = image.rsplit("/", 1)[-1]
    if ":" not in last:
        return "no tag"
    if last.rsplit(":", 1)[1] == "latest":
        return ":latest"
    return None


def main() -> int:
    raw = sys.stdin.read()
    start = raw.find("{")
    if start < 0:
        print("::error::no JSON found in `nextflow inspect` output")
        return 1
    data, _ = json.JSONDecoder().raw_decode(raw[start:])

    bad = []
    for proc in sorted(data.get("processes", []), key=lambda p: p["name"]):
        why = is_unpinned(proc.get("container", ""))
        if why:
            bad.append((proc["name"], proc.get("container", ""), why))
        else:
            print(f"  ok  {proc['name']:26s} {proc['container']}")

    if bad:
        print()
        for name, image, why in bad:
            print(f"::error::{name} uses an unpinned image ({why}): {image or '<none>'}")
        print(f"\n{len(bad)} process(es) are not pinned to a released version.")
        return 1

    print(f"\nAll {len(data.get('processes', []))} processes pinned to a specific image version.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
