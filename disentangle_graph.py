#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set, Optional

from Bio import Align


@dataclass(frozen=True)
class Segment:
    name: str
    seq: str

    @property
    def length(self) -> int:
        return len(self.seq)


@dataclass(frozen=True)
class Edge:
    src_seg: str
    src_orient: str
    dst_seg: str
    dst_orient: str
    identity: float
    overlap: int
    src_len: int
    dst_len: int

    @property
    def score(self) -> float:
        return self.identity * self.overlap


@dataclass(frozen=True)
class PathResult:
    nodes: Tuple[Tuple[str, str], ...]
    total_len: int
    mean_identity: float
    min_identity: float
    total_overlap: int
    score: float
    kind: str


@dataclass(frozen=True)
class SelfTerminalHit:
    seg: str
    orient: str
    identity: float
    overlap: int
    length: int

    @property
    def score(self) -> float:
        return self.identity * self.overlap


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]


def parse_gfa(path: str) -> Dict[str, Segment]:
    out: Dict[str, Segment] = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("S\t"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[1]
            seq = parts[2]
            if seq == "*":
                continue
            out[name] = Segment(name=name, seq=seq)
    return out


def make_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    return aligner


def best_identity_and_len(aligner: Align.PairwiseAligner, query: str, target: str) -> Tuple[float, int]:
    if not query or not target:
        return 0.0, 0
    try:
        aln = next(iter(aligner.align(query, target)))
    except StopIteration:
        return 0.0, 0
    qstart, qend = aln.aligned[0][0][0], aln.aligned[0][-1][-1]
    aligned_len = qend - qstart
    if aligned_len <= 0:
        return 0.0, 0
    score = float(getattr(aln, "score"))
    identity = max(0.0, (score / aligned_len) * 100.0)
    return identity, aligned_len


def oriented_seq(seg: Segment, orient: str) -> str:
    return seg.seq if orient == "+" else revcomp(seg.seq)


def scan_self_terminals(
    segments: Dict[str, Segment],
    window: int,
) -> List[SelfTerminalHit]:
    aligner = make_aligner()
    hits: List[SelfTerminalHit] = []
    for seg in segments.values():
        w = min(window, seg.length)
        for orient in ["+", "-"]:
            seq = oriented_seq(seg, orient)
            head = seq[:w]
            tail = seq[-w:]
            identity, ovlp = best_identity_and_len(aligner, tail, head)
            hits.append(
                SelfTerminalHit(
                    seg=seg.name,
                    orient=orient,
                    identity=round(identity, 2),
                    overlap=ovlp,
                    length=seg.length,
                )
            )
    hits.sort(key=lambda x: (x.score, x.identity, x.overlap), reverse=True)
    return hits


def build_edges(
    segments: Dict[str, Segment],
    window: int,
    min_id: float,
    min_ovlp: int,
    top_out: int,
) -> List[Edge]:
    aligner = make_aligner()
    names = list(segments.keys())
    edges: List[Edge] = []

    terminals: Dict[Tuple[str, str], Tuple[str, str]] = {}
    for name, seg in segments.items():
        w = min(window, seg.length)
        for orient in ["+", "-"]:
            seq = oriented_seq(seg, orient)
            terminals[(name, orient)] = (seq[:w], seq[-w:])

    for src in names:
        for dst in names:
            if src == dst:
                continue
            for src_orient in ["+", "-"]:
                for dst_orient in ["+", "-"]:
                    src_head, src_tail = terminals[(src, src_orient)]
                    dst_head, _ = terminals[(dst, dst_orient)]
                    identity, ovlp = best_identity_and_len(aligner, src_tail, dst_head)
                    if identity >= min_id and ovlp >= min_ovlp:
                        edges.append(
                            Edge(
                                src_seg=src,
                                src_orient=src_orient,
                                dst_seg=dst,
                                dst_orient=dst_orient,
                                identity=round(identity, 2),
                                overlap=ovlp,
                                src_len=segments[src].length,
                                dst_len=segments[dst].length,
                            )
                        )

    grouped: Dict[Tuple[str, str], List[Edge]] = {}
    for e in edges:
        grouped.setdefault((e.src_seg, e.src_orient), []).append(e)

    reduced: List[Edge] = []
    for key in grouped:
        best = sorted(grouped[key], key=lambda x: (x.score, x.identity, x.overlap), reverse=True)[:top_out]
        reduced.extend(best)
    return reduced


def build_adjacency(edges: List[Edge]) -> Dict[Tuple[str, str], List[Edge]]:
    adj: Dict[Tuple[str, str], List[Edge]] = {}
    for e in edges:
        adj.setdefault((e.src_seg, e.src_orient), []).append(e)
    for k in adj:
        adj[k].sort(key=lambda x: (x.score, x.identity, x.overlap), reverse=True)
    return adj


def cycle_key(nodes: Tuple[Tuple[str, str], ...]) -> Tuple[Tuple[str, str], ...]:
    body = list(nodes[:-1])
    if not body:
        return nodes
    rots: List[Tuple[Tuple[str, str], ...]] = []
    n = len(body)
    for i in range(n):
        r = tuple(body[i:] + body[:i])
        rots.append(r)
    rev = list(reversed(body))
    for i in range(n):
        r = tuple(rev[i:] + rev[:i])
        rots.append(r)
    best = min(rots)
    return best + (best[0],)


def search_paths(
    segments: Dict[str, Segment],
    edges: List[Edge],
    expected_size: int,
    size_tolerance: float,
    max_depth: int,
    max_cycles: int,
    max_linear: int,
) -> Tuple[List[PathResult], List[PathResult]]:
    lo = int(expected_size * (1.0 - size_tolerance))
    hi = int(expected_size * (1.0 + size_tolerance))
    adj = build_adjacency(edges)

    cycles: List[PathResult] = []
    linear: List[PathResult] = []
    seen_cycles: Set[Tuple[Tuple[str, str], ...]] = set()

    for start in adj.keys():
        start_len = segments[start[0]].length
        stack: List[Tuple[Tuple[str, str], List[Edge], int, int, Set[str]]] = [
            (start, [], start_len, 0, {start[0]})
        ]

        while stack:
            node, used_edges, total_len, depth, used_segments = stack.pop()

            if lo <= total_len <= hi and used_edges:
                ids = [e.identity for e in used_edges]
                ov = sum(e.overlap for e in used_edges)
                linear.append(
                    PathResult(
                        nodes=tuple([start] + [(e.dst_seg, e.dst_orient) for e in used_edges]),
                        total_len=total_len,
                        mean_identity=round(sum(ids) / len(ids), 2),
                        min_identity=min(ids),
                        total_overlap=ov,
                        score=round(sum(e.score for e in used_edges), 2),
                        kind="linear",
                    )
                )

            if depth >= max_depth:
                continue

            for e in adj.get(node, []):
                nxt = (e.dst_seg, e.dst_orient)
                nxt_len = segments[e.dst_seg].length

                if nxt == start and depth >= 1:
                    cycle_len = total_len - e.overlap
                    nodes = tuple([start] + [(x.dst_seg, x.dst_orient) for x in used_edges] + [start])
                    key = cycle_key(nodes)
                    if key in seen_cycles:
                        continue
                    seen_cycles.add(key)

                    ids = [x.identity for x in used_edges + [e]]
                    ov = sum(x.overlap for x in used_edges + [e])
                    if lo <= cycle_len <= hi:
                        cycles.append(
                            PathResult(
                                nodes=nodes,
                                total_len=cycle_len,
                                mean_identity=round(sum(ids) / len(ids), 2),
                                min_identity=min(ids),
                                total_overlap=ov,
                                score=round(sum(x.score for x in used_edges + [e]), 2),
                                kind="cycle",
                            )
                        )
                    continue

                if e.dst_seg in used_segments:
                    continue

                new_total = total_len + nxt_len - e.overlap
                if new_total > hi:
                    continue

                stack.append(
                    (
                        nxt,
                        used_edges + [e],
                        new_total,
                        depth + 1,
                        used_segments | {e.dst_seg},
                    )
                )

    linear_sorted = sorted(
        linear,
        key=lambda p: (abs(p.total_len - expected_size), -p.mean_identity, -p.score),
    )[:max_linear]
    cycles_sorted = sorted(
        cycles,
        key=lambda p: (abs(p.total_len - expected_size), -p.mean_identity, -p.score),
    )[:max_cycles]
    return cycles_sorted, linear_sorted


def search_cycles_any_size(
    segments: Dict[str, Segment],
    edges: List[Edge],
    max_depth: int,
    max_cycles: int,
) -> List[PathResult]:
    adj = build_adjacency(edges)
    cycles: List[PathResult] = []
    seen_cycles: Set[Tuple[Tuple[str, str], ...]] = set()

    for start in adj.keys():
        start_len = segments[start[0]].length
        stack: List[Tuple[Tuple[str, str], List[Edge], int, int, Set[str]]] = [
            (start, [], start_len, 0, {start[0]})
        ]
        while stack:
            node, used_edges, total_len, depth, used_segments = stack.pop()
            if depth >= max_depth:
                continue
            for e in adj.get(node, []):
                nxt = (e.dst_seg, e.dst_orient)
                if nxt == start and depth >= 1:
                    cycle_len = total_len - e.overlap
                    nodes = tuple([start] + [(x.dst_seg, x.dst_orient) for x in used_edges] + [start])
                    key = cycle_key(nodes)
                    if key in seen_cycles:
                        continue
                    seen_cycles.add(key)
                    ids = [x.identity for x in used_edges + [e]]
                    ov = sum(x.overlap for x in used_edges + [e])
                    cycles.append(
                        PathResult(
                            nodes=nodes,
                            total_len=cycle_len,
                            mean_identity=round(sum(ids) / len(ids), 2),
                            min_identity=min(ids),
                            total_overlap=ov,
                            score=round(sum(x.score for x in used_edges + [e]), 2),
                            kind="cycle_any",
                        )
                    )
                    continue

                if e.dst_seg in used_segments:
                    continue

                nxt_len = segments[e.dst_seg].length
                stack.append(
                    (
                        nxt,
                        used_edges + [e],
                        total_len + nxt_len - e.overlap,
                        depth + 1,
                        used_segments | {e.dst_seg},
                    )
                )

    return sorted(cycles, key=lambda p: (-p.mean_identity, -p.score, p.total_len))[:max_cycles]


def path_to_str(nodes: Tuple[Tuple[str, str], ...]) -> str:
    return " -> ".join([f"{s}{o}" for s, o in nodes])


def build_edge_index(edges: List[Edge]) -> Dict[Tuple[Tuple[str, str], Tuple[str, str]], Edge]:
    idx: Dict[Tuple[Tuple[str, str], Tuple[str, str]], Edge] = {}
    for e in sorted(edges, key=lambda x: (x.score, x.identity, x.overlap), reverse=True):
        key = ((e.src_seg, e.src_orient), (e.dst_seg, e.dst_orient))
        if key not in idx:
            idx[key] = e
    return idx


def stitch_path_sequence(
    segments: Dict[str, Segment],
    edge_index: Dict[Tuple[Tuple[str, str], Tuple[str, str]], Edge],
    path: PathResult,
) -> Optional[str]:
    if not path.nodes:
        return None
    is_cycle = path.kind.startswith("cycle") and len(path.nodes) >= 2 and path.nodes[0] == path.nodes[-1]
    walk_nodes = path.nodes[:-1] if is_cycle else path.nodes
    if not walk_nodes:
        return None
    seq = oriented_seq(segments[walk_nodes[0][0]], walk_nodes[0][1])
    for i in range(1, len(walk_nodes)):
        prev = walk_nodes[i - 1]
        cur = walk_nodes[i]
        e = edge_index.get((prev, cur))
        if e is None:
            return None
        cur_seq = oriented_seq(segments[cur[0]], cur[1])
        seq += cur_seq[e.overlap:]
    if is_cycle and len(walk_nodes) >= 2:
        close_edge = edge_index.get((walk_nodes[-1], walk_nodes[0]))
        if close_edge is None:
            return None
        if close_edge.overlap > 0 and close_edge.overlap < len(seq):
            seq = seq[:-close_edge.overlap]
    return seq


def write_edges(path: str, edges: List[Edge]) -> None:
    with open(path, "w") as fh:
        fh.write("src_seg\tsrc_orient\tdst_seg\tdst_orient\tidentity_pct\toverlap_bp\tscore\n")
        for e in sorted(edges, key=lambda x: (x.score, x.identity, x.overlap), reverse=True):
            fh.write(
                f"{e.src_seg}\t{e.src_orient}\t{e.dst_seg}\t{e.dst_orient}\t{e.identity:.2f}\t{e.overlap}\t{e.score:.2f}\n"
            )


def write_paths(path: str, paths: List[PathResult]) -> None:
    with open(path, "w") as fh:
        fh.write("kind\ttotal_len\tmean_identity\tmin_identity\ttotal_overlap\tscore\tnodes\n")
        for p in paths:
            fh.write(
                f"{p.kind}\t{p.total_len}\t{p.mean_identity:.2f}\t{p.min_identity:.2f}\t{p.total_overlap}\t{p.score:.2f}\t{path_to_str(p.nodes)}\n"
            )


def write_self_hits(path: str, hits: List[SelfTerminalHit], min_id: float, min_ovlp: int) -> None:
    with open(path, "w") as fh:
        fh.write("seg\torient\tlength\tidentity_pct\toverlap_bp\tscore\tcircular_possible\n")
        for h in hits:
            possible = "yes" if (h.identity >= min_id and h.overlap >= min_ovlp) else "no"
            fh.write(
                f"{h.seg}\t{h.orient}\t{h.length}\t{h.identity:.2f}\t{h.overlap}\t{h.score:.2f}\t{possible}\n"
            )


def write_fasta(path: str, records: List[Tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def summarize_contig_support(
    segments: Dict[str, Segment],
    self_hits: List[SelfTerminalHit],
    edges: List[Edge],
    cycles_any: List[PathResult],
    min_id: float,
    min_ovlp: int,
) -> List[Tuple[str, str, str, int, int, float, int, str]]:
    best_self: Dict[str, SelfTerminalHit] = {}
    for h in self_hits:
        prev = best_self.get(h.seg)
        if prev is None or h.score > prev.score:
            best_self[h.seg] = h

    edge_count: Dict[str, int] = {k: 0 for k in segments}
    for e in edges:
        edge_count[e.src_seg] = edge_count.get(e.src_seg, 0) + 1
        edge_count[e.dst_seg] = edge_count.get(e.dst_seg, 0) + 1

    cycle_members: Dict[str, int] = {k: 0 for k in segments}
    for c in cycles_any:
        members = set(seg for seg, _ in c.nodes[:-1])
        for m in members:
            cycle_members[m] = cycle_members.get(m, 0) + 1

    rows: List[Tuple[str, str, str, int, int, float, int, str]] = []
    for seg in sorted(segments.keys()):
        h = best_self.get(seg)
        if h is None:
            rows.append((seg, "na", "no", segments[seg].length, 0, 0.0, edge_count.get(seg, 0), "none"))
            continue
        self_possible = "yes" if (h.identity >= min_id and h.overlap >= min_ovlp) else "no"
        cyc = "yes" if cycle_members.get(seg, 0) > 0 else "no"
        rows.append((seg, cyc, self_possible, h.length, h.overlap, h.identity, edge_count.get(seg, 0), f"{cycle_members.get(seg, 0)}"))
    return rows


def write_report(
    report_path: str,
    gfa: str,
    segments: Dict[str, Segment],
    edges: List[Edge],
    self_hits: List[SelfTerminalHit],
    cycles: List[PathResult],
    linear: List[PathResult],
    cycles_any: List[PathResult],
    contig_rows: List[Tuple[str, str, str, int, int, float, int, str]],
    min_id: float,
    min_ovlp: int,
) -> None:
    with open(report_path, "w") as fh:
        fh.write("Disentangle/Circularization Report\n")
        fh.write("================================\n\n")
        fh.write(f"Input GFA: {gfa}\n")
        fh.write(f"Segments: {len(segments)}\n")
        fh.write(f"Terminal overlap edges: {len(edges)}\n")
        fh.write(f"Expected-size cycles: {len(cycles)}\n")
        fh.write(f"Expected-size linear paths: {len(linear)}\n")
        fh.write(f"Any-size cycles: {len(cycles_any)}\n")
        fh.write(f"Self-terminal circular threshold: identity>={min_id}, overlap>={min_ovlp}\n\n")

        fh.write("Top Terminal Edges\n")
        fh.write("------------------\n")
        for e in sorted(edges, key=lambda x: (x.score, x.identity, x.overlap), reverse=True)[:15]:
            fh.write(
                f"{e.src_seg}{e.src_orient} -> {e.dst_seg}{e.dst_orient} | "
                f"id={e.identity:.2f}% ovlp={e.overlap} score={e.score:.2f}\n"
            )
        fh.write("\n")

        fh.write("Top Self-Terminal Hits\n")
        fh.write("----------------------\n")
        for h in self_hits[:15]:
            possible = "yes" if (h.identity >= min_id and h.overlap >= min_ovlp) else "no"
            fh.write(
                f"{h.seg}{h.orient} | len={h.length} id={h.identity:.2f}% ovlp={h.overlap} "
                f"score={h.score:.2f} circular_possible={possible}\n"
            )
        fh.write("\n")

        fh.write("Top Any-Size Cycles\n")
        fh.write("-------------------\n")
        for p in cycles_any[:20]:
            fh.write(
                f"len={p.total_len} mean_id={p.mean_identity:.2f}% min_id={p.min_identity:.2f}% "
                f"score={p.score:.2f} nodes={path_to_str(p.nodes)}\n"
            )
        fh.write("\n")

        fh.write("Per-Contig Circularization Possibility\n")
        fh.write("--------------------------------------\n")
        fh.write("contig\tcycle_possible\tself_terminal_possible\tlen\tself_ovlp\tself_id\tedge_degree\tcycle_count\n")
        for row in contig_rows:
            fh.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]:.2f}\t{row[6]}\t{row[7]}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description="Terminal-overlap graph disentangling and circularization search")
    parser.add_argument("gfa", help="Input GFA")
    parser.add_argument("--window", type=int, default=500)
    parser.add_argument("--min-id", type=float, default=75.0)
    parser.add_argument("--min-ovlp", type=int, default=100)
    parser.add_argument("--top-out", type=int, default=5)
    parser.add_argument("--expected-size", type=int, default=494000)
    parser.add_argument("--size-tolerance", type=float, default=0.20)
    parser.add_argument("--max-depth", type=int, default=8)
    parser.add_argument("--max-cycles", type=int, default=30)
    parser.add_argument("--max-linear", type=int, default=50)
    parser.add_argument("--max-cycles-any", type=int, default=200)
    parser.add_argument("--fasta-top", type=int, default=20)
    parser.add_argument("--out-prefix", default="disentangle")
    args = parser.parse_args()

    segments = parse_gfa(args.gfa)
    if not segments:
        raise SystemExit("No segment records found in GFA")

    self_hits = scan_self_terminals(segments=segments, window=args.window)
    edges = build_edges(
        segments=segments,
        window=args.window,
        min_id=args.min_id,
        min_ovlp=args.min_ovlp,
        top_out=args.top_out,
    )
    cycles, linear = search_paths(
        segments=segments,
        edges=edges,
        expected_size=args.expected_size,
        size_tolerance=args.size_tolerance,
        max_depth=args.max_depth,
        max_cycles=args.max_cycles,
        max_linear=args.max_linear,
    )
    cycles_any = search_cycles_any_size(
        segments=segments,
        edges=edges,
        max_depth=args.max_depth,
        max_cycles=args.max_cycles_any,
    )

    edge_path = f"{args.out_prefix}.terminal_overlaps.tsv"
    path_path = f"{args.out_prefix}.paths.tsv"
    cycle_any_path = f"{args.out_prefix}.cycles_any.tsv"
    self_path = f"{args.out_prefix}.self_terminal.tsv"
    fasta_path = f"{args.out_prefix}.path_sequences.fasta"
    report_path = f"{args.out_prefix}.report.txt"

    write_edges(edge_path, edges)
    write_paths(path_path, cycles + linear)
    write_paths(cycle_any_path, cycles_any)
    write_self_hits(self_path, self_hits, min_id=args.min_id, min_ovlp=args.min_ovlp)

    edge_index = build_edge_index(edges)
    top_paths = cycles_any[: args.fasta_top]
    if len(top_paths) < args.fasta_top:
        need = args.fasta_top - len(top_paths)
        top_paths = top_paths + (cycles + linear)[:need]
    records: List[Tuple[str, str]] = []
    for i, p in enumerate(top_paths, start=1):
        seq = stitch_path_sequence(segments=segments, edge_index=edge_index, path=p)
        if seq is None or not seq:
            continue
        # Name after the merged contigs (unique segment names in traversal order)
        merged_segs = []
        seen_seg_names: set = set()
        for seg_name, _ in p.nodes:
            if seg_name not in seen_seg_names:
                merged_segs.append(seg_name)
                seen_seg_names.add(seg_name)
        contig_name = "__".join(merged_segs)
        name = f"{contig_name}|{p.kind}|len={p.total_len}|mean_id={p.mean_identity:.2f}"
        records.append((name, seq))
    write_fasta(fasta_path, records)

    contig_rows = summarize_contig_support(
        segments=segments,
        self_hits=self_hits,
        edges=edges,
        cycles_any=cycles_any,
        min_id=args.min_id,
        min_ovlp=args.min_ovlp,
    )
    write_report(
        report_path=report_path,
        gfa=args.gfa,
        segments=segments,
        edges=edges,
        self_hits=self_hits,
        cycles=cycles,
        linear=linear,
        cycles_any=cycles_any,
        contig_rows=contig_rows,
        min_id=args.min_id,
        min_ovlp=args.min_ovlp,
    )

    print(f"Segments: {len(segments)}")
    print(f"Edges passing thresholds: {len(edges)}")
    print(f"Expected-size cycles: {len(cycles)}")
    print(f"Expected-size linear paths: {len(linear)}")
    print(f"Any-size cycles: {len(cycles_any)}")
    print(f"Per-contig self-terminal entries: {len(self_hits)}")
    print(f"Wrote: {edge_path}")
    print(f"Wrote: {path_path}")
    print(f"Wrote: {cycle_any_path}")
    print(f"Wrote: {self_path}")
    print(f"Wrote: {fasta_path}")
    print(f"Wrote: {report_path}")


if __name__ == "__main__":
    main()
