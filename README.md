Tools for computing alpha beta gamma euler angles to describe the orientation of RNA motifs.



## How to install 

```python
python setup.py install
```



## How to compute the euler angles

Should generate the module `abg` and the script `calc_abg`

If you only want to compute the alpha beta gamma angles use the script

```bash
# using 3 base pairs above and below the motif of interest
calc_abg -pdb examples/ens-1.pdb -helix_1_res "3,4,5,25,26,27" -helix_2_res "10,11,12,21,22,23"

# using 2
calc_abg -pdb examples/ens-1.pdb -helix_1_res "4,5,25,26" -helix_2_res "10,11,22,23"

# using 1
calc_abg -pdb examples/ens-1.pdb -helix_1_res "5,25" -helix_2_res "10,23"

```

Output looks like

```bash
68.421   20.051  -19.221    0.442    0.657
```

In the format: alpha beta gamma rmsd1 rmsd2 

where rmsd1 and rmsd 2 is how well your bottom and top stem can be overlayed with an idealized helix.

## How to call the package in another script

```python
import abg.compute

pdb_path = "examples/ens-1.pdb"
helix_1_resi = [3, 4, 5, 25, 26, 27]
helix_2_resi = [10, 11, 12, 21, 22, 23]

abg_computer = abg.compute.ABGComputer()
r = abg_computer.compute(pdb_path, helix_1_resi, helix_2_resi)

print(r.a, r.b, r.g, r.rmsd1, r.rmsd2)

```









