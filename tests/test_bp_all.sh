#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

default_prefix="$tmpdir/default"
all_prefix="$tmpdir/all"

"$repo_root/ped-sim" \
  -d "$repo_root/tests/bp_all.def" \
  -m "$repo_root/refined_mf.simmap" \
  -o "$default_prefix" \
  --pois \
  --bp

"$repo_root/ped-sim" \
  -d "$repo_root/tests/bp_all.def" \
  -m "$repo_root/refined_mf.simmap" \
  -o "$all_prefix" \
  --pois \
  --bp_all

default_bp="$default_prefix.bp"
all_bp="$all_prefix.bp"

default_lines="$(wc -l < "$default_bp")"
all_lines="$(wc -l < "$all_bp")"

if [[ "$default_lines" -ne 4 ]]; then
  echo "expected default BP output to have 4 lines, got $default_lines" >&2
  exit 1
fi

if [[ "$all_lines" -ne 8 ]]; then
  echo "expected --bp_all alone to produce 8 BP lines, got $all_lines" >&2
  exit 1
fi

if grep -q '^bpall1_g1-b1-' "$default_bp"; then
  echo "default BP output unexpectedly included generation 1 samples" >&2
  exit 1
fi

if ! grep -q '^bpall1_g1-b1-i1 ' "$all_bp"; then
  echo "expected --bp_all BP output to include the generation 1 founder" >&2
  exit 1
fi

if ! grep -q '^bpall1_g1-b1-s1 ' "$all_bp"; then
  echo "expected --bp_all BP output to include the generation 1 spouse" >&2
  exit 1
fi

if ! grep -q '^bpall1_g2-b1-i1 ' "$default_bp"; then
  echo "default BP output did not include printed generation 2 samples" >&2
  exit 1
fi
