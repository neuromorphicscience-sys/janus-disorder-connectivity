"""Compare every feedback edge with two independent, unpruned full-graph oracles.

Run from the repository root:
    python scripts/validate_full_graph_classifier.py

The fixed suite includes every loopless directed graph on four ordered vertices
and 26 complete Janus model graphs, rather than induced subsets. No graph is
selected after inspecting its labels. These are algorithm checks, not additions
to the frozen physical ensembles used in the figures.
"""
from __future__ import annotations

import argparse
from collections import deque
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from janus_connectivity.feedback_order import classify_graph, make_filtered_graph
from janus_connectivity.graph_build import build_directed_graph
from janus_connectivity.model import generate_orientations, generate_points


def state_space_oracle(graph, y):
    """Exhaust all (vertex, upward-cost<=1) return states, excluding test edge.

    This uses the full adjacency. It has no geometric bounds, endpoint buckets,
    production reachability helpers, or early exit on finding a return.
    """
    src = np.repeat(np.arange(len(y)), np.diff(graph.indptr))
    upward = y[graph.indices] > y[src]
    labels = np.zeros(graph.nnz, dtype=np.uint8)
    for tested in np.flatnonzero(upward):
        start, goal = int(graph.indices[tested]), int(src[tested])
        seen = np.zeros((2, len(y)), dtype=bool)
        seen[0, start] = True
        queue = deque([(start, 0)])
        while queue:
            vertex, cost = queue.popleft()
            for edge in range(graph.indptr[vertex], graph.indptr[vertex + 1]):
                if edge == tested:
                    continue
                target = int(graph.indices[edge])
                new_cost = cost + int(upward[edge])
                if new_cost <= 1 and not seen[new_cost, target]:
                    seen[new_cost, target] = True
                    queue.append((target, new_cost))
        labels[tested] = 1 if seen[0, goal] else 2 if seen[1, goal] else 3
    return labels


def all_partner_oracle(graph, y):
    """Compute unpruned D closure and enumerate every upward partner.

    The dense closure is deliberately confined to small validation graphs.
    Identity entries include zero-length downward paths. No spatial bound or
    production DAG-search helper is used.
    """
    src = np.repeat(np.arange(len(y)), np.diff(graph.indptr))
    upward = y[graph.indices] > y[src]
    up = np.flatnonzero(upward)
    u, v = src[up], graph.indices[up]
    reach = np.eye(len(y), dtype=bool)
    for vertex in np.argsort(y):
        for edge in range(graph.indptr[vertex], graph.indptr[vertex + 1]):
            target = graph.indices[edge]
            if y[target] < y[vertex]:
                reach[vertex] |= reach[target]
    labels = np.zeros(graph.nnz, dtype=np.uint8)
    for index, edge in enumerate(up):
        partners = reach[v[index], u] & reach[v, u[index]]
        partners[index] = False
        labels[edge] = (1 if reach[v[index], u[index]] else
                        2 if partners.any() else 3)
    return labels, reach


def check_graph(graph, y):
    graph = graph.tocsr()
    graph.sort_indices()
    result = classify_graph(graph, y)
    actual = np.zeros(graph.nnz, dtype=np.uint8)
    actual[result["up_ix"]] = result["labels"]
    expected = state_space_oracle(graph, y)
    enumerated, reach = all_partner_oracle(graph, y)
    if not (np.array_equal(actual, expected) and np.array_equal(actual, enumerated)):
        bad = np.flatnonzero((actual != expected) | (actual != enumerated))
        raise AssertionError(f"Edge-label disagreement at CSR indices {bad.tolist()}")
    u, v = result["u"], result["v"]
    for i, label in enumerate(result["labels"]):
        j = int(result["witness"][i])
        if label == 1:
            assert j == i and reach[v[i], u[i]]
        elif label == 2:
            assert j != i and 0 <= j < len(u)
            assert reach[v[i], u[j]] and reach[v[j], u[i]]
    # Check the exact retained edges and complete SCC partitions at both orders.
    for order in (1, 2):
        keep = (expected == 0) | (expected <= order)
        reference = graph.copy()
        reference.data = keep.astype(np.uint8)
        reference.eliminate_zeros()
        filtered = make_filtered_graph(result, include_order2=(order == 2))
        assert (filtered != reference).nnz == 0
        _, a = connected_components(filtered, directed=True, connection="strong")
        _, b = connected_components(reference, directed=True, connection="strong")
        assert np.array_equal(a[:, None] == a[None, :], b[:, None] == b[None, :])
    counts = np.bincount(actual, minlength=4)
    return dict(N=len(y), edges=graph.nnz, downward=int(counts[0]),
                order1=int(counts[1]), order2=int(counts[2]),
                gt2_or_infinite=int(counts[3]), mismatches=0), actual, expected


def four_vertex_graphs():
    pairs = [(u, v) for u in range(4) for v in range(4) if u != v]
    for mask in range(1 << len(pairs)):
        edges = [pair for bit, pair in enumerate(pairs) if mask & (1 << bit)]
        src, dst = zip(*edges) if edges else ([], [])
        graph = csr_matrix((np.ones(len(edges), np.uint8), (src, dst)), shape=(4, 4))
        yield mask, graph, np.arange(4, dtype=float)


def model_cases():
    # Fixed before any oracle comparison; include every case, even without cycles.
    for L in (2, 4, 6):
        for W in (0.0, 0.805, 0.9, 1.6):
            for replicate in range(2):
                yield dict(L=L, W=W, N=24 * L * L, q=7, R=1,
                           geometry_seed=2026101000 + replicate,
                           orientation_seed=2026101100 + replicate)
    for W in (0.805, 0.9):
        yield dict(L=8, W=W, N=1536, q=7, R=1,
                   geometry_seed=2026101000, orientation_seed=2026101100)


def recheck_archived_subgraphs(output):
    """Rebuild all 80 archived induced graphs, including their isolated nodes."""
    folder = ROOT / "data/final_bridge"
    with (folder / "CLASSIFIER_VALIDATION_COVERAGE.csv").open() as handle:
        coverage = list(csv.DictReader(handle))
    records = []
    for case in coverage:
        name, n = case["name"], int(case["N"])
        with (folder / "classifier_edges" / f"{name}.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        y = np.full(n, np.nan)
        for row in rows:
            for endpoint in ("source", "target"):
                node, height = int(row[endpoint]), float(row[endpoint + "_y"])
                assert np.isnan(y[node]) or y[node] == height
                y[node] = height
        # Isolated vertices have no role in any witness. Give unrecorded heights
        # distinct values outside the observed positive coordinate interval.
        missing = np.flatnonzero(np.isnan(y))
        y[missing] = -1. - np.arange(len(missing))
        src = [int(row["source"]) for row in rows]
        dst = [int(row["target"]) for row in rows]
        graph = csr_matrix((np.ones(len(rows), np.uint8), (src, dst)), shape=(n, n))
        record, actual, _ = check_graph(graph, y)
        old = {(int(row["source"]), int(row["target"])): int(row["production_label"]) for row in rows}
        ordered_src = np.repeat(np.arange(n), np.diff(graph.indptr))
        expected = np.array([old[int(u), int(v)] for u, v in zip(ordered_src, graph.indices)])
        assert np.array_equal(actual, expected), name
        _, translated, _ = check_graph(graph, y - 512.)
        assert np.array_equal(translated, expected), name
        records.append(dict(case=name, **record, archived_labels_unchanged=True,
                            negative_translation_unchanged=True))
    assert len(records) == 80
    with (output / "ARCHIVED_SUBGRAPH_RECHECK.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    return dict(graphs=80, edges=sum(row["edges"] for row in records),
                order1=sum(row["order1"] for row in records),
                order2=sum(row["order2"] for row in records),
                archived_labels_unchanged=True, negative_translation_unchanged=True)


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    cases = list(model_cases())
    (output / "FIXED_CASES.json").write_text(json.dumps(cases, indent=2) + "\n")
    rows = []
    small_totals = dict(graphs=0, edges=0, downward=0, order1=0, order2=0,
                        gt2_or_infinite=0, mismatches=0)
    for mask, graph, y in four_vertex_graphs():
        record, _, _ = check_graph(graph, y)
        rows.append(dict(case=f"four_vertex_{mask:04d}", family="all_four_vertex", **record))
        small_totals["graphs"] += 1
        for key in small_totals:
            if key != "graphs":
                small_totals[key] += record[key]
    print(f"All {small_totals['graphs']} four-vertex graphs agree.", flush=True)
    model_rows = []
    with (output / "MODEL_EDGE_COMPARISONS.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["case", "source", "target", "source_y", "target_y",
                         "production_label", "state_space_label"])
        for case in cases:
            points = generate_points(case["N"], case["L"], case["geometry_seed"])
            theta = generate_orientations(case["N"], case["W"], case["orientation_seed"])
            graph, _, _ = build_directed_graph(points, theta, q=case["q"], R=case["R"])
            name = f"L{case['L']}_W{case['W']}_seed{case['geometry_seed']}"
            record, actual, expected = check_graph(graph, points[:, 1])
            record = dict(case=name, family="complete_model", **record,
                          **{key: value for key, value in case.items() if key != "N"},
                          graph_sha256=hashlib.sha256(graph.indptr.tobytes() + graph.indices.tobytes()).hexdigest())
            model_rows.append(record)
            src = np.repeat(np.arange(case["N"]), np.diff(graph.indptr))
            for edge, (u, v) in enumerate(zip(src, graph.indices)):
                writer.writerow([name, u, v, points[u, 1], points[v, 1], actual[edge], expected[edge]])
            print(f"{name}: N={case['N']}, edges={graph.nnz}, order2={record['order2']}; agree", flush=True)
    for name, records in [("FOUR_VERTEX_CASES.csv", rows), ("MODEL_CASES.csv", model_rows)]:
        with (output / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
    totals = {key: sum(row[key] for row in model_rows) for key in
              ("edges", "downward", "order1", "order2", "gt2_or_infinite", "mismatches")}
    totals.update(graphs=len(model_rows), min_N=min(row["N"] for row in model_rows),
                  max_N=max(row["N"] for row in model_rows))
    bound_graph = csr_matrix((np.ones(3, np.uint8), ([0, 1, 2], [1, 2, 0])), shape=(3, 3))
    _, bound_labels, _ = check_graph(bound_graph, np.array([-1.5, -1., .2]))
    assert bound_labels.tolist() == [2, 2, 0]
    archived = recheck_archived_subgraphs(output)
    summary = dict(all_comparisons_passed=True, oracle_count=2,
                   production_classifier_sha256=hashlib.sha256((ROOT / "src/janus_connectivity/feedback_order.py").read_bytes()).hexdigest(),
                   exhaustive_four_vertex=small_totals, complete_model=totals,
                   archived_subgraphs=archived, negative_coordinate_ceiling_regression=True,
                   witnesses_checked=True, filtered_edges_and_scc_partitions_checked=True,
                   scope="Algorithm validation only; frozen manuscript ensembles unchanged")
    (output / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/full_graph_classifier")
    run(parser.parse_args().output)
