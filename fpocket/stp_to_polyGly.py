import numpy as np
from sklearn.cluster import DBSCAN
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue, Atom

# 1. Parse PQR coordinates
coords = []
with open("pocket348_vert.pqr") as f:
    for line in f:
        if line.startswith("ATOM"):
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            coords.append([x, y, z])

coords = np.array(coords)

# 2. Cluster to remove redundancy
db = DBSCAN(eps=1.5, min_samples=1).fit(coords)
centroids = np.array([coords[db.labels_ == i].mean(axis=0) 
                      for i in set(db.labels_) if i != -1])

# 3. Nearest-neighbor path through centroids
def nn_path(pts):
    visited, path = [False]*len(pts), []
    idx = 0
    for _ in range(len(pts)):
        path.append(idx)
        visited[idx] = True
        dists = np.linalg.norm(pts - pts[idx], axis=1)
        dists[visited] = np.inf
        idx = np.argmin(dists) if not all(visited) else idx
    return path

order = nn_path(centroids)
ca_coords = centroids[order]

# 4. Write polyglycine PDB (Cα only for simplicity)
with open("polyglycine_pocket348.pdb", "w") as out:
    for i, (x, y, z) in enumerate(ca_coords):
        out.write(f"ATOM  {i+1:5d}  CA  GLY A{i+1:4d}    "
                  f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
    out.write("END\n")