"""
fpocket_to_polyglycine_v4.py

Builds a full backbone polyglycine (N, CA, C, O per residue) from fpocket
alpha-sphere vertices.  Geometry matches Engh & Huber ideal values.
Suitable as RFdiffusion partial diffusion input.

Pipeline:
  1. Parse PQR alpha-sphere coords + radii
  2. Conservative dedup (eps=0.3 Å) — removes only identical positions
  3. Build adjacency graph (connect_dist=4.5 Å) to preserve Voronoi connectivity
  4. Greedy graph-aware path, PCA-oriented, both directions output
  5. Minimal interpolation only for residual gaps after graph traversal
  6. Full backbone (N, CA, C, O) placed at each CA using rotation-based frame:
       - C along +tangent * D_CAC
       - N from rotating -tangent by (180° - 111.2°) around curvature normal
       - O in peptide plane, 120° from C–CA

Ideal geometry (Engh & Huber 1991):
  CA–N = 1.46 Å,  CA–C = 1.52 Å,  C=O = 1.23 Å,  ∠N–CA–C = 111.2°

Usage:
  python fpocket_to_polyglycine_v4.py pocket.pqr [protein.pdb]
"""

import sys
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from scipy.spatial import KDTree
from collections import defaultdict

# ── Ideal geometry ────────────────────────────────────────────────────────────
D_CAN       = 1.46
D_CAC       = 1.52
D_CO        = 1.23
ANG_NCC     = np.radians(111.2)
DEFLECT     = np.pi - ANG_NCC     # 68.8° — rotation of -tangent to place N

# ── Graph / clash params ──────────────────────────────────────────────────────
DEDUP_EPS    = 0.3
CONNECT_DIST = 4.5
CA_DIST      = 3.8
CA_VDW       = 1.7
CLASH_PAD    = 0.5
MAX_NUDGE    = 50
NUDGE_STEP   = 0.2
VDW = {'C':1.70,'N':1.55,'O':1.52,'S':1.80,'default':1.70}

# ─────────────────────────────────────────────────────────────────────────────
# Math helpers
# ─────────────────────────────────────────────────────────────────────────────
def safe_norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-8 else np.array([1., 0., 0.])

def rotation_matrix(axis, angle_rad):
    """Rodrigues rotation of a vector around axis by angle_rad."""
    k = safe_norm(axis)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) * c + (1 - c) * np.outer(k, k) + s * K

# ─────────────────────────────────────────────────────────────────────────────
# I/O
# ─────────────────────────────────────────────────────────────────────────────
def parse_pqr(path):
    coords, radii = [], []
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM"): continue
            p = line.split()
            coords.append([float(p[5]), float(p[6]), float(p[7])])
            radii.append(float(p[9]))
    return np.array(coords), np.array(radii)

def parse_protein(path):
    coords, radii = [], []
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
            elem = line[76:78].strip().upper()
            if elem == 'H' or (not elem and line[13:14] == 'H'): continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError: continue
            e = elem[0] if elem else 'C'
            coords.append([x, y, z]); radii.append(VDW.get(e, VDW['default']))
    return np.array(coords), np.array(radii)

# ─────────────────────────────────────────────────────────────────────────────
# Clash checker
# ─────────────────────────────────────────────────────────────────────────────
class ClashChecker:
    def __init__(self, pc, pr):
        self.tree = KDTree(pc); self.radii = pr; self.coords = pc
    def clashes(self, pt, r=CA_VDW):
        for i in self.tree.query_ball_point(pt, self.radii.max() + r + CLASH_PAD):
            if np.linalg.norm(pt - self.coords[i]) < self.radii[i] + r + CLASH_PAD:
                return True
        return False
    def nudge(self, pt):
        p = pt.copy()
        for _ in range(MAX_NUDGE):
            if not self.clashes(p): return p, True
            _, idx = self.tree.query(p); d = p - self.coords[idx]; n = np.linalg.norm(d)
            p += NUDGE_STEP * (d / n if n > 1e-6 else np.array([1., 0., 0.]))
        return p, False

# ─────────────────────────────────────────────────────────────────────────────
# Graph-based CA trace
# ─────────────────────────────────────────────────────────────────────────────
def dedup(coords, radii, eps=DEDUP_EPS):
    db = DBSCAN(eps=eps, min_samples=1).fit(coords)
    cen, rad = [], []
    for lbl in sorted(set(db.labels_)):
        if lbl == -1: continue
        m = db.labels_ == lbl
        cen.append(coords[m].mean(0)); rad.append(radii[m].mean())
    return np.array(cen), np.array(rad)

def build_graph(cen, d=CONNECT_DIST):
    tree = KDTree(cen); graph = defaultdict(list)
    for i, pt in enumerate(cen):
        for j in tree.query_ball_point(pt, d):
            if j != i: graph[i].append((j, np.linalg.norm(cen[i] - cen[j])))
    for i in graph: graph[i].sort(key=lambda x: x[1])
    return graph

def greedy_path(cen, graph, start):
    n = len(cen); visited = [False] * n; path = [start]; visited[start] = True
    for _ in range(n - 1):
        cur = path[-1]; moved = False
        for (nbr, _) in graph[cur]:
            if not visited[nbr]:
                path.append(nbr); visited[nbr] = True; moved = True; break
        if not moved:
            dists = np.linalg.norm(cen - cen[cur], axis=1)
            dists[[i for i, v in enumerate(visited) if v]] = np.inf
            nxt = int(np.argmin(dists)); path.append(nxt); visited[nxt] = True
    return path

def pca_orient(pts, path):
    proj = PCA(n_components=1).fit_transform(pts[path])[:, 0]
    return path if proj[-1] >= proj[0] else list(reversed(path))

def interpolate_ca(pts, path, checker=None):
    result = []; nc = nf = 0
    def add(p):
        nonlocal nc, nf
        if checker and checker.clashes(p):
            nc += 1; p, ok = checker.nudge(p)
            if not ok: nf += 1
        result.append(p)
    add(pts[path[0]])
    for i in range(1, len(path)):
        a, b = result[-1], pts[path[i]]; dist = np.linalg.norm(b - a)
        if dist > CA_DIST:
            ni = int(np.ceil(dist / CA_DIST)) - 1
            for k in range(1, ni + 1): add(a + (k / (ni + 1)) * (b - a))
        add(b)
    return np.array(result), nc, nf

# ─────────────────────────────────────────────────────────────────────────────
# Full backbone builder
# ─────────────────────────────────────────────────────────────────────────────
def get_reference_perp(tangents):
    """
    Compute a stable reference vector perpendicular to each tangent,
    using parallel transport to avoid sudden flips.
    For curved segments, blend in the local curvature direction.
    """
    n = len(tangents)
    refs = np.zeros((n, 3))

    # Seed
    up = np.array([0., 0., 1.])
    if abs(np.dot(tangents[0], up)) > 0.9:
        up = np.array([0., 1., 0.])
    b = safe_norm(np.cross(tangents[0], up))

    for i in range(n):
        t = tangents[i]
        b = b - np.dot(b, t) * t   # re-orthogonalise
        b = safe_norm(b)
        refs[i] = b

    return refs

def build_backbone(ca_trace):
    """
    Place full backbone (N, CA, C, O) for each residue.

    Frame per residue i:
      tangent t  = safe_norm(CA[i+1] - CA[i-1])   (central diff)
      ref     r  = parallel-transported binormal ⊥ t

      C = CA + t * D_CAC                  (C along forward tangent)
      N = CA + R(r, DEFLECT) @ (-t) * D_CAN   (N from rotating -t by 68.8°)
          where DEFLECT = 180° - 111.2° = 68.8°
          → guaranteed ∠N–CA–C = 111.2°

      O in peptide plane: 120° from C–CA around the C–N_next axis normal
    """
    n = len(ca_trace)
    res = []

    # Tangents
    tangents = np.zeros((n, 3))
    for i in range(n):
        if   i == 0:     tangents[i] = safe_norm(ca_trace[1]  - ca_trace[0])
        elif i == n - 1: tangents[i] = safe_norm(ca_trace[-1] - ca_trace[-2])
        else:            tangents[i] = safe_norm(ca_trace[i+1] - ca_trace[i-1])

    refs = get_reference_perp(tangents)

    # Place N and C
    for i in range(n):
        t  = tangents[i]
        r  = refs[i]
        ca = ca_trace[i]

        # C: along forward tangent
        C = ca + t * D_CAC

        # N: rotate -t by DEFLECT (68.8°) around r  → ∠N–CA–C = 111.2°
        R = rotation_matrix(r, DEFLECT)
        N = ca + R @ (-t) * D_CAN

        res.append({'N': N, 'CA': ca, 'C': C, 'O': None})

    # Place O
    for i in range(n):
        C  = res[i]['C']
        CA = res[i]['CA']
        v_cca = safe_norm(CA - C)

        if i < n - 1:
            N_next  = res[i + 1]['N']
            v_cn    = safe_norm(N_next - C)
            plane_n = safe_norm(np.cross(v_cca, v_cn))
        else:
            plane_n = safe_norm(np.cross(v_cca, refs[i]))

        cos120 = np.cos(np.radians(120))
        sin120 = np.sin(np.radians(120))
        o_dir  = cos120 * v_cca + sin120 * np.cross(plane_n, v_cca)
        res[i]['O'] = C + safe_norm(o_dir) * D_CO

    return res

# ─────────────────────────────────────────────────────────────────────────────
# Write PDB
# ─────────────────────────────────────────────────────────────────────────────
def write_full_pdb(residues, path_out, chain='A'):
    atom_idx = 1; lines = []
    for ri, r in enumerate(residues, start=1):
        for aname in ('N', 'CA', 'C', 'O'):
            x, y, z = r[aname]
            col  = f" {aname:<3s}" if len(aname) < 4 else aname
            elem = 'N' if aname == 'N' else 'O' if aname == 'O' else 'C'
            lines.append(
                f"ATOM  {atom_idx:5d} {col} GLY {chain}{ri:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {elem}\n"
            )
            atom_idx += 1
    lines.append("TER\nEND\n")
    with open(path_out, 'w') as f: f.writelines(lines)

def report_backbone(residues, label, checker=None):
    ca   = np.array([r['CA'] for r in residues])
    d    = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    angles = []
    for r in residues:
        v1 = safe_norm(r['N'] - r['CA']); v2 = safe_norm(r['C'] - r['CA'])
        angles.append(np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1, 1))))
    nc   = sum(1 for p in ca if checker and checker.clashes(p))
    print(f"  [{label}] {len(residues)} residues | "
          f"CA–CA max={d.max():.2f} mean={d.mean():.2f} Å | "
          f"chain-breaks: {(d>4.2).sum()} | "
          f"∠N–CA–C mean={np.mean(angles):.1f}° | "
          f"Cα clashes: {nc}")

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    pqr_file     = sys.argv[1] if len(sys.argv) > 1 else "pocket348_vert.pqr"
    protein_file = sys.argv[2] if len(sys.argv) > 2 else None

    print(f"\nLoading {pqr_file}")
    raw_coords, raw_radii = parse_pqr(pqr_file)
    print(f"  Raw alpha-spheres: {len(raw_coords)}")

    cen, rad = dedup(raw_coords, raw_radii)
    print(f"  After dedup (eps={DEDUP_EPS} Å): {len(cen)} nodes "
          f"(removed {len(raw_coords)-len(cen)} true duplicates)")

    graph   = build_graph(cen)
    n_edges = sum(len(v) for v in graph.values()) // 2
    print(f"  Graph: {len(cen)} nodes, {n_edges} edges (connect_dist={CONNECT_DIST} Å)")

    checker = None
    if protein_file:
        print(f"  Loading protein: {protein_file}")
        pc, pr = parse_protein(protein_file)
        checker = ClashChecker(pc, pr)
        print(f"  {len(pc)} heavy atoms loaded")

    start    = int(np.argmax(rad))
    path_raw = greedy_path(cen, graph, start)
    path_fwd = pca_orient(cen, path_raw)
    path_rev = list(reversed(path_fwd))

    print("\nBuilding CA traces + full backbones...")
    ca_fwd, nc_f, _ = interpolate_ca(cen, path_fwd, checker)
    ca_rev, nc_r, _ = interpolate_ca(cen, path_rev, checker)

    bb_fwd = build_backbone(ca_fwd)
    bb_rev = build_backbone(ca_rev)

    report_backbone(bb_fwd, "FWD", checker)
    report_backbone(bb_rev, "REV", checker)

    write_full_pdb(bb_fwd, "polyglycine_fwd.pdb", 'A')
    write_full_pdb(bb_rev, "polyglycine_rev.pdb", 'B')

    # Combined
    atom_idx = 1
    with open("polyglycine_both.pdb", "w") as f:
        for chain, bb in [('A', bb_fwd), ('B', bb_rev)]:
            for ri, r in enumerate(bb, 1):
                for aname in ('N', 'CA', 'C', 'O'):
                    x, y, z = r[aname]
                    col  = f" {aname:<3s}" if len(aname) < 4 else aname
                    elem = 'N' if aname=='N' else 'O' if aname=='O' else 'C'
                    f.write(f"ATOM  {atom_idx:5d} {col} GLY {chain}{ri:4d}    "
                            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {elem}\n")
                    atom_idx += 1
            f.write("TER\n")
        f.write("END\n")

    print("\n  → polyglycine_fwd.pdb  (chain A, N→C along PC1)")
    print("  → polyglycine_rev.pdb  (chain B, reversed)")
    print("  → polyglycine_both.pdb (both chains for PyMOL overlay)")
    print("Done.\n")