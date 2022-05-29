Tools for computing alpha beta gamma euler angles to describe the orientation of RNA motifs.



## How to install 

```python
python setup.py install
```



## How to compute the euler angles

Should generate the module `abg` and the script `calc_abg`

If you only want to compute the alpha beta gamma angles use the script. Use 3 basepairs whenever possible.

```bash
# using 3 base pairs above and below the motif of interest
calc_abg -helix_1_res "3,4,5,25,26,27" -helix_2_res "10,11,12,21,22,23" examples/ens-1.pdb

# using 2
calc_abg -helix_1_res "4,5,25,26" -helix_2_res "10,11,22,23" examples/ens-1.pdb

# using 1
calc_abg -helix_1_res "5,25" -helix_2_res "10,23" examples/ens-1.pdb
```

Output looks like

```bash
2022-05-28 19:13:18,043 - abg.main - INFO - residues in helix 1: [3, 4, 5, 25, 26, 27]
2022-05-28 19:13:18,043 - abg.main - INFO - residues in helix 2: [10, 11, 12, 21, 22, 23]
2022-05-28 19:13:18,043 - abg.main - INFO - examples/ens-1.pdb is determined to be in PDB format
   a     b      g    rmsd1    rmsd2      x     y      z
----  ----  -----  -------  -------  -----  ----  -----
68.4  20.1  -19.2     0.44     0.66  -5.65  8.92  -11.5
```

where rmsd1 and rmsd 2 is how well your bottom and top stem can be overlayed with an idealized helix.


You can also supply a csv file and calculate abg on every pdb in the csv file. Ensure your pdb path is supplied the 'pdb' column

```shell
calc_abg -helix_1_res "3,4,5,25,26,27" -helix_2_res "10,11,12,21,22,23" examples/test.csv 
```

## How to call the package in another script

```python
import abg.compute

pdb_path = "examples/ens-1.pdb"
helix_1_resi = [3, 4, 5, 25, 26, 27]
helix_2_resi = [10, 11, 12, 21, 22, 23]

abg_computer = abg.compute.ABGComputer()
r = abg_computer.compute(pdb_path, helix_1_resi, helix_2_resi)

print(r.a, r.b, r.g, r.rmsd1, r.rmsd2, r.x, r.y, r.z)

```









