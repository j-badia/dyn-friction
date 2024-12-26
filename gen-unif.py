import numpy as np
import sys

if len(sys.argv) != 6:
    print("Usage: python gen-unif.py dist half-count big-pos big-vel fname")
    exit()

d, N, r, v, fname = sys.argv[1:]
d = float(d)
N = int(N)
r = float(r)
v = float(v)

coord = np.linspace(-(N-0.5)*d, (N-0.5)*d, 2*N)
x, y, z = np.meshgrid(coord, coord, coord)
small_pos = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T # cartesian product

with open(fname, "w") as f:
    f.write(f"{r} 0 0 ")
    for pos in small_pos[:-1]:
        f.write(f"{pos[0]} {pos[1]} {pos[2]} ")
    f.write(f"{small_pos[-1][0]} {small_pos[-1][1]} {small_pos[-1][2]}\n")
    f.write(f"0 {v} 0 ")
    f.write((len(small_pos)-1)*"0 0 0 ")
    f.write("0 0 0")