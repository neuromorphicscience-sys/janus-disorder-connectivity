# Completeness of the order-two classifier

Let `D` contain strictly downward edges and `F` the upward edges, using the
same progress coordinate `y` as the production classifier. Every `D` path
strictly decreases `y`; reachability includes a zero-length path. Edges with
equal endpoint heights are rejected by the strict sign decomposition.

For `i=(u_i,v_i)`, order one means that `v_i` reaches `u_i` in `D`.
Such a path never goes below `y(u_i)`, so the lower slab in `exact_f1_mask`
does not exclude any return. For an edge without an order-one return, an
order-two witness exists exactly when some distinct `j=(u_j,v_j)` satisfies

    v_i ->*_D u_j   and   v_j ->*_D u_i.

The two paths and two upward edges concatenate to a closed return. Conversely,
a return with only one additional upward edge splits into these two downward
paths. Paths may have length zero, and the partner may itself have order one.

Let `h_j = y(v_j)-y(u_j)` and let `H` be the maximum upward-edge height gain in
the graph. Monotonicity gives `y(u_j)<=y(v_i)` and `y(v_j)>=y(u_i)`. Since
`h_j<=H`, every valid partner satisfies

    y(u_i)-H <= y(u_j) <= y(v_i),
    y(u_i)   <= y(v_j) <= y(v_i)+H.

Therefore the reverse search from `u_i` may stop above `y(v_i)+H`: the whole
reversed downward path to any valid `v_j` remains below this ceiling. The
implementation first expands the computed maximum height gain with
`nextafter(H,+inf)` and then expands the resulting sum with another
`nextafter(...,+inf)`. This covers both subtraction and addition rounding.
Endpoint buckets include every upward edge at each visited target, without
filtering by the partner's own order.

Let `C_i` be the distinct upward edges with targets found by that search. If
`C_i` is empty, no witness exists. Otherwise define
`m_i = min_{j in C_i} y(u_j)`. A downward path from `v_i` to any valid `u_j`
never goes below `y(u_j)>=m_i`; thus the forward search may discard vertices
below `m_i`. Inclusive comparisons retain equality at the floor/ceiling.
Finding an eligible source is exactly the stated two-path condition. These
bounds lose no valid partner; they affect search work, not classification.

## Direct reproduction

From the repository root after installing `requirements.txt`:

```sh
python scripts/validate_full_graph_classifier.py
python -m pytest tests/test_feedback_order.py tests/test_full_graph_classifier.py
```

The script compares every edge with two independent references:

1. Full-adjacency traversal of every reachable state `(vertex, upward cost)`
   with cost zero or one, excluding the tested upward edge. No geometric bound,
   endpoint bucket, production search helper, or early successful-return exit
   is used. A reached source at cost zero/one gives order one/two.
2. Complete unpruned downward reachability, followed by enumeration of every
   distinct upward partner for every tested edge. This validation-only dense
   closure is never used by the production algorithm.

The fixed suite is declared before any comparisons:

- All 4096 loopless directed graphs on four vertices with strictly ordered
  heights. Relabeling covers any strict height ordering of a four-vertex graph.
- Twenty-four complete Janus graphs at `L/R=2,4,6`, `sigma R^2=24`, `q=7`,
  `W=0,0.805,0.9,1.6`, and two fixed seed pairs.
- Two complete graphs at `L/R=8`, `W=0.805,0.9`, and the first seed pair.

The complete model graphs have `N=96,384,864,1536`. Each is constructed directly
from its parameters, rather than extracted as an induced subset. Geometry seeds
are 2026101000 and 2026101001; orientation seeds are 2026101100 and 2026101101.
The suite includes every declared case, without choosing cases by their labels
or agreement. It does not augment the frozen physical ensembles in the figures.

Results: all 4096 topologies agree on 24,576 directed-edge occurrences, including
1984 order-two labels. All 26 complete model graphs agree on 96,768 directed
edges, including 3278 order-one and 4395 order-two labels. Returned witnesses are checked against the unpruned downward closure. The
retained edge sets and complete SCC partitions of `G1` and `G2` agree with the
state-space reference. These are deterministic correctness checks, not statistical
confidence estimates or a substitute for the proof above.

Outputs in `outputs/full_graph_classifier` include all parameters, seeds,
per-case counts, CSR graph hashes, model-edge comparisons and the production
classifier's source hash. Reference outputs are archived in
`data/classifier_full_graph`. Regression tests separately exercise ceiling and
floor equality, shared endpoints and an order-one partner. Existing toy and
induced-subgraph validation remains available.

## Floating-point ceiling regression

Rounding only the final sum outward is insufficient when negative coordinates
cause cancellation. A three-node cycle with heights `[-1.5,-1.0,0.2]` and edges
`0->1, 1->2, 2->0` has two order-two upward edges. However, the computed gain
`fl(0.2-(-1.0))` is slightly below the exact gain. Even
`nextafter(fl(-1.0+fl(1.2)),+inf)` remains below the stored value `0.2`, so a
ceiling computed this way misses the witness for `0->1`.

The production ceiling now rounds the gain up before adding it, then rounds the
sum up. For correctly rounded finite arithmetic, the first step bounds the exact
height gain of every upward edge, and the second bounds the exact endpoint sum.
This is a conservative expansion of the search only: the unpruned two-path test
still determines every accepted witness. Regression tests include this example
and positive scaling/coordinate translations. The command also rechecks all 80
archived induced graphs, reconstructing their incident-vertex heights from the
edge archive, and translates their coordinates by -512 without changing any
edge direction. The original archived labels remain unchanged. Heights of
isolated vertices were not recorded in the edge CSV; distinct unused values
are assigned, which cannot affect any witness.
