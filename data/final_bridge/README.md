# Final matched bridge source data

All analyses use existing graph identities only. The paired manifest contains 48 original/reassigned pairs sharing 16 geometry/orientation streams across three W values. Every CSR hash has been verified.

S, S1 and S2 are full-domain largest-SCC fractions, including singleton components. Low-order feedback labels are assigned in the complete graph before strict 4R core restriction. f_int_rec is the mass in nontrivial core SCCs divided by retained nodes; eta_int is the largest nontrivial mass divided by recurrent mass. S_int is the largest nontrivial core SCC fraction. If recurrent mass is zero, eta_int is undefined (blank/NaN); the separate zero_convention column records the software convention only. Paired concentration differences/ratios are undefined in that case. No pseudocount is used.

Orientation correlation uses the same 512 prespecified source anchors (seed 2026100801), all targets within 8R, and four radial bins. C_conn subtracts the measured global polarization squared. The exact same-multiset permutation expectation is (N P^2-1)/(N-1), also supplied. Pair averages are anchor sampled. No correlation length is fitted.

Indegree histograms contain all nodes. Reciprocal-edge fraction is the fraction of directed edges having a reverse edge. WCC is the largest weakly connected fraction.

Classifier coverage retains 40 historical induced subsets and 40 deterministic feedback-seeded local subsets from the same five parents. The historical rule uses two uniform and six weak-neighborhood BFS subsets per parent. All edges, vertex indices, production labels and independent oracle labels are supplied.

Finite-size uncertainty uses 4,000 joint shared-stream bootstrap draws per size (seed [2026100802,L]); percentile intervals are conditional on uniquely bracketed crossings. Successful fractions are explicit. No critical exponent, thermodynamic threshold, smoothing or collapse is fitted.

The clean profile retains raw per-graph/bin source-node and upward-edge counts on all 24 archived clean graphs. It averages over the full lateral interval, with ensemble means/SEM matching a mean prediction. Main Fig.4(c) uses all 192 deletion graphs at 4R; each pair uses the same retained node count.

Reproduction from the repository root:

    python scripts/replay_submission_bridge.py --workers 3
    python scripts/audit_submission_classifier.py
    python scripts/replay_clean_profile.py
    python scripts/bootstrap_submission_crossovers.py
    python scripts/reproduce_bridge.py
    python scripts/reproduce_fig4.py

These replay commands rebuild only the supplied frozen identities and verify hashes. They are optional for figure reproduction, which reads the supplied CSVs directly. Derived outputs go under outputs/final_bridge. No new seed production is exposed by these scripts.
