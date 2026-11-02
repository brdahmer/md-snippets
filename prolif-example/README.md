# ProLIF interaction fingerprints over an MD trajectory

`mTOR_RSV36.ipynb` computes protein–ligand interaction fingerprints with [ProLIF](https://github.com/chemosim-lab/ProLIF) over a 300 ns molecular dynamics trajectory of mTOR in complex with resveratrol (residue name `RSV`). `mTOR_RSV36.html` is a static export of the same notebook with all outputs.

The notebook is kept with its outputs so it can be read directly on GitHub.

## Running

The input files (`mtor-rsv36.gro` and `md-rsv36-fitted.xtc`) are not included because of their size. To reuse the notebook, replace them with your own topology/structure and a fitted trajectory, and change the ligand residue name in the selections.

## Reference

Bouysset, C. & Fiorucci, S. ProLIF: a library to encode molecular interactions as fingerprints. *J. Cheminform.* **13**, 72 (2021). [doi:10.1186/s13321-021-00548-6](https://doi.org/10.1186/s13321-021-00548-6)
