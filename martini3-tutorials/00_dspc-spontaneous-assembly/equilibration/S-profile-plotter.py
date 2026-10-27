import numpy as np
import matplotlib.pyplot as plt

bonds = ["NC3-PO4","PO4-GL1","GL1-GL2","GL1-C1A","GL2-C1B",
         "C1A-C2A","C2A-C3A","C3A-C4A","C1B-C2B","C2B-C3B","C3B-C4B"]

# per-frame values (skip the header, separator and summary lines)
rows = [l.split() for l in open("order.dat") if l.split() and l.split()[0].isdigit()]
S = np.array(rows, float)[:, 1:]
mean = S.mean(0)
sem = S.std(0, ddof=1) / np.sqrt(len(S))   # assumes roughly independent frames

x = np.arange(len(bonds))
plt.errorbar(x, mean, yerr=sem, fmt="o-", capsize=3)
plt.axhline(0, color="gray", lw=0.5)
plt.xticks(x, bonds, rotation=45, ha="right")
plt.ylabel("S")
plt.title("DSPC bond order parameters")
plt.tight_layout()
plt.savefig("S-profile.png", dpi=300)
plt.show()