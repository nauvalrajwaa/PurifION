#!/usr/bin/env python3
"""
check_joins.py
--------------
Parse a GFA assembly graph, extract terminal windows from every contig,
run all-vs-all pairwise alignment, and report candidate contig joins
that could form a circular path.

Usage:
    python check_joins.py <assembly.gfa> [options]

Options:
    --window    Terminal window size in bp (default: 500)
    --min-id    Minimum alignment identity 0-100 (default: 75)
    --min-ovlp  Minimum overlap length in bp (default: 100)
    --out       Output TSV file (default: candidate_joins.tsv)
"""
import argparse
import sys
import itertools
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional

# Biopython pairwise aligner
from Bio import Align


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class Segment:
    name: str
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)

    def head(self, n: int) -> str:
        """First n bp (5'-end of contig)."""
        return self.seq[:n]

    def tail(self, n: int) -> str:
        """Last n bp (3'-end of contig)."""
        return self.seq[-n:]

    def revcomp(self, seq: str) -> str:
        comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        return seq.translate(comp)[::-1]


@dataclass
class CandidateJoin:
    seg_a: str
    end_a: str          # 'head' or 'tail'
    seg_b: str
    end_b: str          # 'head' or 'tail'
    orientation: str    # '+/+', '+/-', etc.
    identity: float
    overlap_len: int
    combined_len: int   # seg_a.length + seg_b.length - overlap_len
    note: str = ""


# ---------------------------------------------------------------------------
# GFA parser
# ---------------------------------------------------------------------------

def parse_gfa(path: str) -> Dict[str, Segment]:
    segments = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("S\t"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            name, seq = parts[1], parts[2]
            if seq == "*":
                continue
            segments[name] = Segment(name=name, seq=seq)
    return segments


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def make_aligner(mode: str = "overlap") -> Align.PairwiseAligner:
    """Configure a semi-global / end-gap-free aligner for overlap detection."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"           # local catches partial overlaps too
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    return aligner


def align_pair(
    aligner: Align.PairwiseAligner,
    query: str,
    target: str,
) -> Tuple[float, int]:
    """
    Return (identity, aligned_length) for the best local alignment.
    identity = matches / aligned_length.
    Returns (0.0, 0) if no alignment found.
    """
    if not query or not target:
        return 0.0, 0
    alignments = aligner.align(query, target)
    try:
        best = next(iter(alignments))
    except StopIteration:
        return 0.0, 0

    score = best.score
    # aligned length = length of aligned region in query
    qstart, qend = best.aligned[0][0][0], best.aligned[0][-1][-1]
    aligned_len = qend - qstart
    if aligned_len == 0:
        return 0.0, 0
    # identity = score / aligned_len  (since match=1, mismatch=-1 → approx)
    identity = score / aligned_len * 100
    return max(0.0, identity), aligned_len


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def extract_terminals(
    segments: Dict[str, Segment],
    window: int,
) -> Dict[str, Dict[str, str]]:
    """
    For each segment return {'head': seq, 'tail': seq, 'head_rc': ..., 'tail_rc': ...}.
    Segments shorter than window are used in full.
    """
    terminals = {}
    for name, seg in segments.items():
        w = min(window, seg.length)
        h = seg.seq[:w]
        t = seg.seq[-w:]
        terminals[name] = {
            "head": h,
            "tail": t,
            "head_rc": revcomp(h),
            "tail_rc": revcomp(t),
        }
    return terminals


def find_candidate_joins(
    segments: Dict[str, Segment],
    terminals: Dict[str, Dict[str, str]],
    aligner: Align.PairwiseAligner,
    min_identity: float,
    min_overlap: int,
) -> List[CandidateJoin]:
    """
    For every ordered pair (A, B) where A != B, test four join orientations:
      1. tail(A) vs head(B)  → A → B  (linear extension, circular-compatible)
      2. tail(A) vs head_rc(B) → A → B_rev
      3. tail_rc(A) vs head(B) → A_rev → B  (equivalent to head(A_rc) → B)
      4. tail(A) vs tail_rc(B)  (both point inward: A→ ←B)

    We focus on orientation (1) and (2) as the primary circular-path signals.
    All four are reported so the user can reason about graph topology.
    """
    candidates = []
    names = list(segments.keys())
    seen = set()  # deduplicate symmetric pairs

    orientations = [
        # (query_key, target_key, orientation_label, note)
        ("tail",    "head",    "+/+", "A_tail → B_head (A→B linear)"),
        ("tail",    "head_rc", "+/-", "A_tail → B_head_rc (A→B_rev)"),
        ("tail_rc", "head",    "-/+", "A_tail_rc → B_head (A_rev→B)"),
        ("tail_rc", "head_rc", "-/-", "A_tail_rc → B_head_rc (both rev)"),
    ]

    for a_name, b_name in itertools.permutations(names, 2):
        pair_key = tuple(sorted([a_name, b_name]))
        seg_a = segments[a_name]
        seg_b = segments[b_name]
        t_a = terminals[a_name]
        t_b = terminals[b_name]

        for qk, tk, orient_label, note in orientations:
            # Avoid duplicate (+/+ and -/- are symmetric across pairs)
            dedup_key = (pair_key, orient_label)
            if dedup_key in seen:
                continue
            seen.add(dedup_key)

            query = t_a[qk]
            target = t_b[tk]
            identity, ovlp_len = align_pair(aligner, query, target)

            if identity >= min_identity and ovlp_len >= min_overlap:
                combined = seg_a.length + seg_b.length - ovlp_len
                candidates.append(CandidateJoin(
                    seg_a=a_name,
                    end_a="tail" if "tail" in qk else "head",
                    seg_b=b_name,
                    end_b="head" if "head" in tk else "tail",
                    orientation=orient_label,
                    identity=round(identity, 2),
                    overlap_len=ovlp_len,
                    combined_len=combined,
                    note=note,
                ))

    return candidates


def infer_circular_paths(
    candidates: List[CandidateJoin],
    segments: Dict[str, Segment],
    expected_size: int = 494_000,
    size_tolerance: float = 0.20,
) -> List[List[str]]:
    """
    Simple greedy DFS to find chains of joins that could close a circle.
    Accepts paths whose total length is within ±size_tolerance of expected_size.
    """
    from collections import defaultdict

    # Build adjacency: tail of A → head of B (only +/+ for simplicity)
    adj = defaultdict(list)
    for c in candidates:
        if c.orientation == "+/+" and c.end_a == "tail" and c.end_b == "head":
            adj[c.seg_a].append((c.seg_b, c.overlap_len, c.identity))

    lo = expected_size * (1 - size_tolerance)
    hi = expected_size * (1 + size_tolerance)

    circular_paths = []

    def dfs(path: List[str], length: int, visited: set):
        node = path[-1]
        for nxt, ovlp, _ in adj.get(node, []):
            if nxt == path[0]:
                # potential circle closure
                total = length - ovlp
                if lo <= total <= hi:
                    circular_paths.append(list(path))
                continue
            if nxt in visited:
                continue
            new_len = length + segments[nxt].length - ovlp
            if new_len > hi:
                continue
            visited.add(nxt)
            path.append(nxt)
            dfs(path, new_len, visited)
            path.pop()
            visited.remove(nxt)

    for start in segments:
        dfs([start], segments[start].length, {start})

    return circular_paths


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def write_tsv(candidates: List[CandidateJoin], path: str):
    header = [
        "seg_a", "end_a", "seg_b", "end_b",
        "orientation", "identity_pct", "overlap_bp",
        "combined_len_bp", "note",
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for c in sorted(candidates, key=lambda x: -x.identity):
            fh.write("\t".join(str(getattr(c, f)) for f in [
                "seg_a", "end_a", "seg_b", "end_b",
                "orientation", "identity", "overlap_len",
                "combined_len", "note",
            ]) + "\n")


def print_summary(
    candidates: List[CandidateJoin],
    circular_paths: List[List[str]],
    segments: Dict[str, Segment],
):
    print(f"\n{'='*60}")
    print(f"  TERMINAL OVERLAP ANALYSIS — CANDIDATE JOINS")
    print(f"{'='*60}")
    print(f"  Segments analysed : {len(segments)}")
    print(f"  Candidate joins   : {len(candidates)}")
    print()

    if not candidates:
        print("  ❌  No candidate joins found above thresholds.")
        print("      → Contigs do not share terminal sequence at this identity.")
        print("      → This supports the conclusion that data depth/read length")
        print("        is the binding limit, not just Flye parameters.")
        return

    # Top joins by identity
    top = sorted(candidates, key=lambda x: -x.identity)[:20]
    print(f"  Top {min(20, len(top))} candidate joins (by identity):")
    print(f"  {'Seg A':<12} {'End':<6} {'Seg B':<12} {'End':<6} "
          f"{'Orient':<6} {'ID%':>6} {'Ovlp bp':>8} {'Note'}")
    print(f"  {'-'*90}")
    for c in top:
        seg_a_len = segments[c.seg_a].length
        seg_b_len = segments[c.seg_b].length
        print(f"  {c.seg_a:<12} {c.end_a:<6} {c.seg_b:<12} {c.end_b:<6} "
              f"{c.orientation:<6} {c.identity:>6.1f} {c.overlap_len:>8}   "
              f"{c.note}  [lens: {seg_a_len}/{seg_b_len}]")

    print()
    if circular_paths:
        print(f"  ✅  Potential circular paths ({len(circular_paths)} found):")
        for i, path in enumerate(circular_paths[:5], 1):
            total = sum(segments[n].length for n in path)
            print(f"  Path {i}: {' → '.join(path)} → (circle)  ~{total:,} bp")
    else:
        print("  ⚠️   No circular paths found via greedy DFS.")
        print("      Candidate joins exist but do not chain into a circle")
        print("      of the expected genome size (~494 kb).")

    print(f"\n{'='*60}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Check terminal contig overlaps in a GFA for circular join candidates."
    )
    parser.add_argument("gfa", help="Input GFA file (assembly graph)")
    parser.add_argument("--window", type=int, default=500,
                        help="Terminal window size in bp (default: 500)")
    parser.add_argument("--min-id", type=float, default=75.0,
                        help="Minimum alignment identity %% (default: 75)")
    parser.add_argument("--min-ovlp", type=int, default=100,
                        help="Minimum overlap length in bp (default: 100)")
    parser.add_argument("--expected-size", type=int, default=494_000,
                        help="Expected genome size for circular path DFS (default: 494000)")
    parser.add_argument("--size-tolerance", type=float, default=0.20,
                        help="±fraction tolerance for circular path size (default: 0.20)")
    parser.add_argument("--out", default="candidate_joins.tsv",
                        help="Output TSV file (default: candidate_joins.tsv)")
    args = parser.parse_args()

    print(f"[1/5] Parsing GFA: {args.gfa}")
    segments = parse_gfa(args.gfa)
    print(f"      → {len(segments)} segments loaded")
    total_bp = sum(s.length for s in segments.values())
    print(f"      → Total assembly: {total_bp:,} bp")

    print(f"[2/5] Extracting terminal windows ({args.window} bp each end)")
    terminals = extract_terminals(segments, args.window)

    print(f"[3/5] Building aligner (local, match=1 mismatch=-1 gap=-2/-0.5)")
    aligner = make_aligner()

    n_pairs = len(segments) * (len(segments) - 1)
    print(f"[4/5] Running all-vs-all terminal alignments ({n_pairs} ordered pairs × 4 orientations)…")
    candidates = find_candidate_joins(
        segments, terminals, aligner,
        min_identity=args.min_id,
        min_overlap=args.min_ovlp,
    )

    print(f"[5/5] Inferring circular paths (expected size: {args.expected_size:,} ±{int(args.size_tolerance*100)}%)")
    circular_paths = infer_circular_paths(
        candidates, segments,
        expected_size=args.expected_size,
        size_tolerance=args.size_tolerance,
    )

    write_tsv(candidates, args.out)
    print(f"      → Results written to: {args.out}")

    print_summary(candidates, circular_paths, segments)


if __name__ == "__main__":
    main()
