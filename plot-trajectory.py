import numpy as np
from matplotlib import pyplot as plt
import sys

fname = sys.argv[1]
x, y, z = np.loadtxt(fname, usecols=[1, 2, 3], unpack=True)

fig, ax = plt.subplots()
ax.plot(x, y, color="black", linestyle="-", linewidth=1)
ax.set_aspect(1)
plt.show()