# Model and definitions

## Spatial graph

N = round(sigma L^2) independent points are sampled uniformly in an open square of side L. Source i has quenched direction e_i = (sin theta_i, -cos theta_i), with theta_i = W z_i modulo 2 pi and z_i standard normal. Candidates lie at Euclidean distance smaller than 2R. The source retains the q candidates with largest projection e_i dot (r_j-r_i), or every candidate if fewer than q are present.

S is the fraction of vertices in the largest strongly connected component of the full graph G.

## Feedback-order filtration

D0 contains edges with strictly decreasing progress coordinate y and is acyclic. F contains the remaining upward edges. For f=(u,v) in F, c(f) is the minimum number of feedback edges in a directed cycle containing f. G_k contains D0 and feedback edges with c(f) <= k.

The order-one test asks whether u is reachable from v using D0 only. The order-two test seeks a distinct feedback edge g=(s,t) such that s is D0-reachable from v and u is D0-reachable from t. The public classifier performs exact paired searches on D0 and its reverse with vertex marks and endpoint buckets. Production storage is O(N+M+E_D); dense reachability is used only by tiny test oracles.

## Internal recurrence

For the strict core I, G_k[I] is induced after feedback labels are assigned in the full frozen graph. M_k^int is the number of vertices in nontrivial SCCs of G_k[I], and C_k^int is its largest nontrivial SCC.

    f_k,rec^int = M_k^int / |I|
    eta_k^int = |C_k^int| / M_k^int

Thus |C_k^int|/|I| = f_k,rec^int eta_k^int graph by graph.

## Surface observables

Downstream depth is measured from the boundary selected by the mean direction. d90 is the depth containing 90 percent of a component. Occupancy profiles divide component members in a strip by all vertices in the same strip.
