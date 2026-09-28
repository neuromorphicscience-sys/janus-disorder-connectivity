"""Exact order<=2 feedback classification without materializing H or closure."""
from __future__ import annotations

import numpy as np
from numba import njit
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

ORDER1 = np.uint8(1)
ORDER2 = np.uint8(2)
GT2_OR_INF = np.uint8(3)


@njit(cache=True)
def filtered_row_pointer(indptr, keep):
    out = np.empty(len(indptr), np.int64)
    out[0] = 0
    for row in range(len(indptr) - 1):
        count = 0
        for k in range(indptr[row], indptr[row + 1]):
            if keep[k]:
                count += 1
        out[row + 1] = out[row] + count
    return out


@njit(cache=True)
def exact_f1_mask(d_indptr, d_indices, y, u, v):
    """Exact D-only return searches; the target-y slab is logically exact."""
    n = len(y)
    closed = np.zeros(len(u), np.bool_)
    marks = np.zeros(n, np.int32)
    queue = np.empty(n, np.int32)
    for i in range(len(u)):
        stamp = i + 1
        lower = y[u[i]]
        head = 0
        tail = 1
        queue[0] = v[i]
        marks[v[i]] = stamp
        while head < tail:
            x = queue[head]
            head += 1
            if x == u[i]:
                closed[i] = True
                break
            for k in range(d_indptr[x], d_indptr[x + 1]):
                z = d_indices[k]
                if y[z] >= lower and marks[z] != stamp:
                    marks[z] = stamp
                    queue[tail] = z
                    tail += 1
    return closed


@njit(cache=True)
def classify_order2_matrix_free(d_ptr, d_idx, r_ptr, r_idx, y, u, v,
                                f1, source_ptr, source_edges,
                                target_ptr, target_edges, max_up_dy):
    """Classify each F edge via two exact, bounded-slab D searches.

    For i=(u_i,v_i), a partner j exists iff u_j is D-reachable from v_i and
    v_j can D-reach u_i. We reverse-search from u_i to mark eligible target
    events, then forward-search from v_i for an eligible source. The y bounds
    follow from strict D descent and the graph's maximum positive dy, so they
    remove no valid reciprocal pair. Only O(N+M+E_D) arrays are held.
    """
    n = len(y)
    m = len(u)
    labels = np.full(m, GT2_OR_INF, np.uint8)
    witness = np.full(m, -1, np.int64)
    for i in range(m):
        if f1[i]:
            labels[i] = ORDER1
            witness[i] = i
    marks_r = np.zeros(n, np.int32)
    marks_f = np.zeros(n, np.int32)
    eligible = np.zeros(m, np.int32)
    queue = np.empty(n, np.int32)
    total_reverse_visits = 0
    total_forward_visits = 0
    max_reverse_visits = 0
    max_forward_visits = 0

    for i in range(m):
        if f1[i]:
            continue
        stamp = i + 1
        # Any partner source u_j satisfies y(u_j)>=y(u_i)-max_up_dy.
        # Its target v_j satisfies y(v_j)<=y(v_i)+max_up_dy.
        high_edge = y[v[i]] + max_up_dy
        high_edge = np.nextafter(high_edge, np.inf)

        # D-reverse search enumerates every v_j that reaches u_i.
        head = 0
        tail = 1
        queue[0] = u[i]
        marks_r[u[i]] = stamp
        min_candidate_source = np.inf
        while head < tail:
            x = queue[head]
            head += 1
            for p in range(target_ptr[x], target_ptr[x + 1]):
                j = target_edges[p]
                if j != i:
                    eligible[j] = stamp
                    if y[u[j]] < min_candidate_source:
                        min_candidate_source = y[u[j]]
            for k in range(r_ptr[x], r_ptr[x + 1]):
                z = r_idx[k]
                if y[z] <= high_edge and marks_r[z] != stamp:
                    marks_r[z] = stamp
                    queue[tail] = z
                    tail += 1
        total_reverse_visits += tail
        if tail > max_reverse_visits:
            max_reverse_visits = tail

        if min_candidate_source == np.inf:
            continue
        # D-forward search stops at its first exact reciprocal witness.
        head = 0
        tail = 1
        queue[0] = v[i]
        marks_f[v[i]] = stamp
        found = -1
        while head < tail and found < 0:
            x = queue[head]
            head += 1
            for p in range(source_ptr[x], source_ptr[x + 1]):
                j = source_edges[p]
                if j != i and eligible[j] == stamp:
                    found = j
                    break
            if found >= 0:
                break
            for k in range(d_ptr[x], d_ptr[x + 1]):
                z = d_idx[k]
                if y[z] >= min_candidate_source and marks_f[z] != stamp:
                    marks_f[z] = stamp
                    queue[tail] = z
                    tail += 1
        total_forward_visits += head
        if head > max_forward_visits:
            max_forward_visits = head
        if found >= 0:
            labels[i] = ORDER2
            witness[i] = found

    return (labels, witness, total_reverse_visits, total_forward_visits,
            max_reverse_visits, max_forward_visits)


def _event_buckets(n: int, endpoint: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ptr = np.zeros(n + 1, dtype=np.int64)
    np.add.at(ptr, endpoint.astype(np.int64) + 1, 1)
    np.cumsum(ptr, out=ptr)
    ids = np.empty(len(endpoint), dtype=np.int32)
    cursor = ptr[:-1].copy()
    for i, x in enumerate(endpoint):
        p = cursor[x]
        ids[p] = i
        cursor[x] += 1
    return ptr, ids


def decompose_graph(G: csr_matrix, y: np.ndarray) -> dict:
    G = G.tocsr(copy=False)
    G.sort_indices()
    y = np.ascontiguousarray(y, dtype=np.float64)
    n = int(G.shape[0])
    src = np.repeat(np.arange(n, dtype=np.int32), np.diff(G.indptr))
    dst = np.asarray(G.indices, dtype=np.int32)
    dy = y[dst] - y[src]
    if np.any(dy == 0.0):
        raise ValueError(f"strict sign decomposition failed: dy=0 count={np.count_nonzero(dy == 0.0)}")
    down = np.ascontiguousarray(dy < 0.0)
    up_ix = np.flatnonzero(dy > 0.0)
    u = np.ascontiguousarray(src[up_ix], dtype=np.int32)
    v = np.ascontiguousarray(dst[up_ix], dtype=np.int32)
    d_ptr = filtered_row_pointer(G.indptr, down)
    D = csr_matrix((G.data[down], G.indices[down], d_ptr), shape=G.shape, copy=False)
    D.sort_indices()
    R = D.transpose().tocsr()
    R.sort_indices()
    source_ptr, source_edges = _event_buckets(n, u)
    target_ptr, target_edges = _event_buckets(n, v)
    return dict(G=G, y=y, src=src, dst=dst, dy=dy, down=down, up_ix=up_ix,
                u=u, v=v,
                d_ptr=np.ascontiguousarray(D.indptr, dtype=np.int64),
                d_idx=np.ascontiguousarray(D.indices, dtype=np.int32),
                r_ptr=np.ascontiguousarray(R.indptr, dtype=np.int64),
                r_idx=np.ascontiguousarray(R.indices, dtype=np.int32),
                max_up_dy=float(dy[up_ix].max()) if len(up_ix) else 0.0,
                source_ptr=source_ptr, source_edges=source_edges,
                target_ptr=target_ptr, target_edges=target_edges)


def classify_graph(G: csr_matrix, y: np.ndarray) -> dict:
    d = decompose_graph(G, y)
    f1 = exact_f1_mask(d["d_ptr"], d["d_idx"], d["y"], d["u"], d["v"])
    result = classify_order2_matrix_free(
        d["d_ptr"], d["d_idx"], d["r_ptr"], d["r_idx"], d["y"],
        d["u"], d["v"], f1, d["source_ptr"], d["source_edges"],
        d["target_ptr"], d["target_edges"], d["max_up_dy"])
    labels, witness, rev_visits, fwd_visits, rev_max, fwd_max = result
    d.update(f1=f1, labels=labels, witness=witness,
             reverse_reachability_visits=int(rev_visits),
             forward_reachability_visits=int(fwd_visits),
             max_reverse_reachability_visits=int(rev_max),
             max_forward_reachability_visits=int(fwd_max))
    return d


def make_filtered_graph(d: dict, include_order2: bool) -> csr_matrix:
    keep = d["down"].copy()
    labels = d["labels"]
    if include_order2:
        keep[d["up_ix"]] = labels <= ORDER2
    else:
        keep[d["up_ix"]] = labels == ORDER1
    ptr = filtered_row_pointer(d["G"].indptr, keep)
    out = csr_matrix((d["G"].data[keep], d["G"].indices[keep], ptr),
                     shape=d["G"].shape, copy=False)
    return out


def scc_summary(G: csr_matrix) -> tuple[int, int, int, np.ndarray]:
    ncomp, labels = connected_components(G, directed=True, connection="strong")
    sizes = np.bincount(labels, minlength=ncomp)
    order = np.sort(sizes)[::-1]
    largest = int(order[0])
    second = int(order[1]) if len(order) > 1 else 0
    return int(ncomp), largest, second, labels
