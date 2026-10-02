"""Draw Figure 1(a) entirely from explicit vector geometry."""
from collections import Counter

import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch


def draw_schematic(ax, definition, colors):
    width, height = definition["canvas"]
    ax.set(xlim=(0, width), ylim=(height, 0), aspect="equal")
    ax.axis("off")
    neutral = colors["neutral"]

    def arrow(p, q, color=neutral, linewidth=.7, head=4.1, style="-|>", zorder=3):
        artist = FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=head,
                                 linewidth=linewidth, color=color,
                                 shrinkA=0, shrinkB=0, zorder=zorder)
        ax.add_patch(artist)
        return artist

    def node(p, radius=11.5, color=neutral, linewidth=.8, filled=False):
        ax.add_patch(Circle(p, radius, facecolor=color if filled else "white",
                            edgecolor=color, linewidth=linewidth, zorder=5))

    def link(p, q, radius=13.5, both=False, linewidth=.7, head=4.1):
        p, q = np.asarray(p, float), np.asarray(q, float)
        direction = (q - p) / np.linalg.norm(q - p)
        arrow(p + radius * direction, q - radius * direction,
              linewidth=linewidth, head=head, style="<|-|>" if both else "-|>")

    local = definition["local"]
    source = np.asarray(local["source"], float)
    candidates = np.asarray(local["candidates"], float)
    direction = np.asarray(local["preferred_direction"], float)
    projections = (candidates - source) @ direction
    selected = set(np.argsort(projections)[-local["q"]:].tolist())
    ax.add_patch(Circle(source, local["radius"], fill=False, edgecolor="#AEB4BA",
                        linewidth=.85, linestyle=(0, (4, 4)), zorder=1))
    for i, target in enumerate(candidates):
        if i in selected:
            link(source, target, radius=18, linewidth=1.15, head=6.5)
        else:
            ax.plot([source[0], target[0]], [source[1], target[1]],
                    color="#BCC1C6", linewidth=.65, linestyle=(0, (3, 3)), zorder=1)
        node(target, radius=12, color="#555A5E", linewidth=.85)
    node(source, radius=16, filled=True)
    arrow(*local["direction_arrow"], color=colors["direction"], linewidth=2, head=9)
    ax.text(211, 168, r"$\mathbf{e}_i$", ha="center", va="center", fontsize=6.5,
            color=colors["direction"])

    for branch in definition["aligned_branches"]:
        points = np.asarray(branch["nodes"], float)
        for u, v in branch["edges"]:
            link(points[u], points[v])
        for point in points:
            node(point)
        arrow(*branch["direction_arrow"], color=colors["direction"], linewidth=1, head=5.8)

    restored = definition["restored"]
    points = np.asarray(restored["nodes"], float)
    edges = {tuple(edge) for edge in restored["edges"]}
    # Reciprocal arrows share one path to keep both arrowheads readable.
    drawn = set()
    for u, v in sorted(edges):
        if (u, v) in drawn:
            continue
        reciprocal = (v, u) in edges
        link(points[u], points[v], both=reciprocal, radius=15, linewidth=.6, head=4.2)
        drawn.add((u, v))
        if reciprocal:
            drawn.add((v, u))
    for point in points:
        node(point, radius=12.5, color=colors["full"], linewidth=1.1)
    center = points.mean(axis=0)
    for u in restored["external_direction_nodes"]:
        point = points[u]
        outgoing = np.array([points[v] - point for s, v in sorted(edges) if s == u])
        preferred = (outgoing / np.linalg.norm(outgoing, axis=1)[:, None]).sum(axis=0)
        preferred /= np.linalg.norm(preferred)
        radial = point - center
        radial /= np.linalg.norm(radial)
        arrow_center = point + 65 * radial
        padding = 18 + np.abs(25 * preferred)
        arrow_center = np.clip(arrow_center, padding, np.array([width, height]) - padding)
        arrow(arrow_center - 25 * preferred, arrow_center + 25 * preferred,
              color=colors["direction"], linewidth=1, head=5.8)

    for p, q in definition["connectors"]:
        arrow(p, q, color=colors["reference"], linewidth=.9, head=6)
    labels = definition["labels"]
    ax.text(*labels["panel"], "(a)", ha="left", va="top", fontsize=8.5, fontweight="bold")
    for name, text, fontsize in [
        ("same_rule", "same\ntop-$q$ rule", 7),
        ("clean", r"$W=0$", 7.5), ("disorder", r"$W\uparrow$", 7.5),
        ("q", r"$q=3$", 7), ("local", "local rule", 7.3),
        ("fragmented", "fragmented", 7.3), ("restored", "restored GSCC", 7.3),
    ]:
        ax.text(*labels[name], text, ha="center", va="center", fontsize=fontsize,
                linespacing=1.12, zorder=6)

    # Record the logical properties of these fixed illustrative objects.
    degree = Counter(u for u, _ in edges)
    reachable = []
    for source in range(len(points)):
        seen, todo = {source}, [source]
        while todo:
            u = todo.pop()
            for s, v in edges:
                if s == u and v not in seen:
                    seen.add(v)
                    todo.append(v)
        reachable.append(len(seen))
    return {"local_selected_indices": sorted(selected), "local_projection_scores": projections.tolist(),
            "all_candidates_inside_disk": bool((np.linalg.norm(candidates - np.asarray(local["source"]), axis=1) <= local["radius"]).all()),
            "aligned_edges_strictly_downstream": all(branch["nodes"][v][1] > branch["nodes"][u][1]
                                                     for branch in definition["aligned_branches"] for u, v in branch["edges"]),
            "restored_nodes": len(points), "restored_edges": len(edges),
            "restored_max_outdegree": max(degree.values()),
            "restored_all_nodes_mutually_reachable": all(n == len(points) for n in reachable)}
