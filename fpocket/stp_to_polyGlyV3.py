"""
fpocket_to_polyglycine_v3.py

Key insight: fpocket alpha-spheres are already a connected spatial graph.
Large gaps in v1/v2 were caused by over-aggressive DBSCAN clustering
(eps=1.5) destroying that connectivity by collapsing legitimate neighbors.

Strategy:
  1. Use eps=0.3 to remove only true duplicates (identical positions <0.5 Å)
  2. Build an explicit adjacency graph where spheres within CONNECT_DIST are
     neighbors — this recovers the Voronoi connectivity
  3. Find the longest path through that graph (approximated via DFS from the
     most-exposed endpoint) — no interpolation needed for well-connected pockets
  4. Only interpolate residual gaps that remain after graph traversal (rare)
  5. PCA orientation + both directions as before
  6. Optional protein clash check carried over from v2
"""

import sys
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from scipy.spatial import KDTree
from collections import defaultdict

# ── Constants ─────────────────────────────────────────────────────────────────
DEDUP_EPS    = 1    # only merge near-identical spheres (<0.3 Å apart)
CONNECT_DIST = 4.5    # max Å between sphere centers to be graph-adjacent
                      # (alpha-sphere typical neighbor distance is 2–4 Å)
CA_DIST      = 3.8    # ideal Cα–Cα distance for gap-filling
CA_VDW       = 1.7
CLASH_PAD    = 0.5
MAX_NUDGE    = 50
NUDGE_STEP   = 0.2
VDW = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'P':1.80,'H':1.20,'default':1.70}

# ── Parse PQR ─────────────────────────────────────────────────────────────────
def parse_pqr(path):
    coords, radii = [], []
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            p = line.split()
            coords.append([float(p[5]), float(p[6]), float(p[7])])
            radii.append(float(p[9]))
    return np.array(coords), np.array(radii)

# ── Parse protein PDB ─────────────────────────────────────────────────────────
def parse_protein(path):
    coords, radii = [], []
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            elem = line[76:78].strip().upper()
            if elem in ('H','') and line[13:14] == 'H':
                continue
            try:
                x,y,z = float(line[30:38]),float(line[38:46]),float(line[46:54])
            except ValueError:
                continue
            e = (elem[0] if elem else 'C')
            coords.append([x,y,z])
            radii.append(VDW.get(e, VDW['default']))
    return np.array(coords), np.array(radii)

# ── Clash checker ─────────────────────────────────────────────────────────────
class ClashChecker:
    def __init__(self, prot_coords, prot_radii):
        self.tree   = KDTree(prot_coords)
        self.radii  = prot_radii
        self.coords = prot_coords
    def clashes(self, pt):
        max_r = self.radii.max() + CA_VDW + CLASH_PAD
        for i in self.tree.query_ball_point(pt, max_r):
            if np.linalg.norm(pt - self.coords[i]) < self.radii[i]+CA_VDW+CLASH_PAD:
                return True
        return False
    def nudge(self, pt):
        p = pt.copy()
        for _ in range(MAX_NUDGE):
            if not self.clashes(p):
                return p, True
            _, idx = self.tree.query(p)
            d = p - self.coords[idx]
            n = np.linalg.norm(d)
            p += NUDGE_STEP * (d/n if n>1e-6 else np.array([1.,0.,0.]))
        return p, False

# ── Step 1: Deduplicate only ──────────────────────────────────────────────────
def dedup(coords, radii, eps=DEDUP_EPS):
    db     = DBSCAN(eps=eps, min_samples=1).fit(coords)
    labels = db.labels_
    cen, rad = [], []
    for lbl in sorted(set(labels)):
        if lbl == -1: continue
        mask = labels == lbl
        cen.append(coords[mask].mean(axis=0))
        rad.append(radii[mask].mean())
    return np.array(cen), np.array(rad)

# ── Step 2: Build adjacency graph ─────────────────────────────────────────────
def build_graph(centroids, connect_dist=CONNECT_DIST):
    """
    Each centroid is a node. Edges connect pairs within connect_dist Å.
    Returns adjacency list.
    """
    tree  = KDTree(centroids)
    graph = defaultdict(list)
    for i, pt in enumerate(centroids):
        neighbors = tree.query_ball_point(pt, connect_dist)
        for j in neighbors:
            if j != i:
                dist = np.linalg.norm(centroids[i] - centroids[j])
                graph[i].append((j, dist))
    # Sort each adjacency list by distance (prefer closer neighbors)
    for i in graph:
        graph[i].sort(key=lambda x: x[1])
    return graph

# ── Step 3: DFS longest path from start node ─────────────────────────────────
def dfs_longest_path(graph, n_nodes, start):
    """
    Greedy DFS: at each step visit the nearest unvisited neighbor.
    This is O(N) and sufficient for small pocket graphs (~20-70 nodes).
    Returns ordered list of node indices.
    """
    visited = [False] * n_nodes
    path    = []

    def dfs(node):
        visited[node] = True
        path.append(node)
        for (nbr, _) in graph[node]:
            if not visited[nbr]:
                dfs(nbr)

    sys.setrecursionlimit(10000)
    dfs(start)

    # If graph has disconnected components, visit remaining nodes
    # by jumping to the nearest unvisited node from current tail
    centroids_ref = None   # set externally below
    return path

def greedy_path(centroids, graph, start):
    """
    Greedy nearest-unvisited traversal that PREFERS graph edges.
    Falls back to any unvisited node if current node has no unvisited neighbors.
    """
    n       = len(centroids)
    visited = [False] * n
    path    = [start]
    visited[start] = True

    for _ in range(n - 1):
        cur  = path[-1]
        # First: try graph neighbors in order of distance
        moved = False
        for (nbr, _) in graph[cur]:
            if not visited[nbr]:
                path.append(nbr)
                visited[nbr] = True
                moved = True
                break
        if not moved:
            # No graph neighbors left: jump to globally nearest unvisited
            dists = np.linalg.norm(centroids - centroids[cur], axis=1)
            dists[[i for i,v in enumerate(visited) if v]] = np.inf
            nxt = int(np.argmin(dists))
            path.append(nxt)
            visited[nxt] = True

    return path

# ── Step 4: PCA orient ────────────────────────────────────────────────────────
def pca_orient(pts, path):
    proj = PCA(n_components=1).fit_transform(pts[path])[:,0]
    return path if proj[-1] >= proj[0] else list(reversed(path))

# ── Step 5: Minimal gap-filling interpolation ─────────────────────────────────
def interpolate(pts, path, checker=None):
    result = []
    n_clash = n_fail = 0

    def add(p):
        nonlocal n_clash, n_fail
        if checker and checker.clashes(p):
            n_clash += 1
            p, ok = checker.nudge(p)
            if not ok: n_fail += 1
        result.append(p)

    add(pts[path[0]])
    for i in range(1, len(path)):
        a, b = result[-1], pts[path[i]]
        d    = np.linalg.norm(b - a)
        if d > CA_DIST:
            n_ins = int(np.ceil(d / CA_DIST)) - 1
            for k in range(1, n_ins + 1):
                t = k / (n_ins + 1)
                add(a + t*(b-a))
        add(b)

    return np.array(result), n_clash, n_fail

# ── Write PDB ─────────────────────────────────────────────────────────────────
def write_pdb(ca, fpath, chain):
    with open(fpath, 'w') as f:
        for i,(x,y,z) in enumerate(ca, 1):
            f.write(f"ATOM  {i:5d}  CA  GLY {chain}{i:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        f.write("TER\nEND\n")

def report(ca, label, checker=None):
    d = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    n_b = (d > 4.2).sum()
    n_c = sum(1 for p in ca if checker and checker.clashes(p))
    print(f"  [{label}] {len(ca)} residues | "
          f"max gap={d.max():.2f} Å  mean={d.mean():.2f} Å | "
          f"chain-breaks: {n_b} | Cα clashes: {n_c}")

# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    pqr_file     = sys.argv[1] if len(sys.argv) > 1 else "pocket348_vert.pqr"
    protein_file = sys.argv[2] if len(sys.argv) > 2 else None

    print(f"\nLoading {pqr_file}")
    raw_coords, raw_radii = parse_pqr(pqr_file)
    print(f"  Raw alpha-spheres: {len(raw_coords)}")

    # Step 1: conservative dedup only
    cen, rad = dedup(raw_coords, raw_radii, eps=DEDUP_EPS)
    print(f"  After dedup (eps={DEDUP_EPS} Å): {len(cen)} nodes  "
          f"(removed {len(raw_coords)-len(cen)} true duplicates)")

    # Step 2: build connectivity graph
    graph = build_graph(cen, connect_dist=CONNECT_DIST)
    n_edges = sum(len(v) for v in graph.values()) // 2
    print(f"  Graph: {len(cen)} nodes, {n_edges} edges "
          f"(connect_dist={CONNECT_DIST} Å)")

    # Check connectivity
    isolated = [i for i in range(len(cen)) if not graph[i]]
    if isolated:
        print(f"  Warning: {len(isolated)} isolated nodes "
              f"(no neighbors within {CONNECT_DIST} Å)")

    # Optional protein clash checker
    checker = None
    if protein_file:
        print(f"\nLoading protein: {protein_file}")
        pc, pr = parse_protein(protein_file)
        checker = ClashChecker(pc, pr)
        print(f"  {len(pc)} heavy atoms loaded")

    # Step 3: path — start from most-exposed centroid (largest alpha-sphere radius)
    start = int(np.argmax(rad))
    path_raw = greedy_path(cen, graph, start)

    # Step 4: PCA orient
    path_fwd = pca_orient(cen, path_raw)
    path_rev = list(reversed(path_fwd))

    # Step 5: gap-fill only residual large gaps
    print()
    ca_fwd, nc_f, nf_f = interpolate(cen, path_fwd, checker)
    ca_rev, nc_r, nf_r = interpolate(cen, path_rev, checker)

    if nc_f: print(f"  FWD: {nc_f} interpolated points clashed, {nf_f} unresolved")
    if nc_r: print(f"  REV: {nc_r} interpolated points clashed, {nf_r} unresolved")

    report(ca_fwd, "FWD", checker)
    report(ca_rev, "REV", checker)

    write_pdb(ca_fwd, "polyglycine_fwd.pdb", "A")
    write_pdb(ca_rev, "polyglycine_rev.pdb", "B")

    with open("polyglycine_both.pdb", "w") as f:
        for i,(x,y,z) in enumerate(ca_fwd, 1):
            f.write(f"ATOM  {i:5d}  CA  GLY A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        f.write("TER\n")
        for i,(x,y,z) in enumerate(ca_rev, 1):
            f.write(f"ATOM  {i:5d}  CA  GLY B{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
        f.write("TER\nEND\n")

    print("\n  → polyglycine_fwd.pdb  polyglycine_rev.pdb  polyglycine_both.pdb")
    print("Done.\n")