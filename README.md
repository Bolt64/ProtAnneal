# ProtAnneal

## About ProtAnneal
ProtAnneal tries to construct a protein with a given backbone shape. It
does that by starting out with a skeleton protein, and goes on mutating the
amino acid and tries the reduce the Zdope score, which is an indicator of the
energy of the protein. This optimized protein is might be the one we're looking
for.

## Dependencies
1. naccess - Can be downloaded from [this
   page](http://www.bioinf.manchester.ac.uk/naccess/).
2. Modeller - Can be downloaded from [this page](https://salilab.org/modeller/).
3. Dunbrack Rotamer library - Can be downloaded from [this
   page](http://dunbrack.fccc.edu/bbdep2010/DownloadPageExample.php).
4. reduce - Can be downloaded from [this
   page](http://kinemage.biochem.duke.edu/software/).

Note that the above dependencies are licensed software, and a license has to be
obtained before downloading them.

## Installation and usage
After installing the dependencies, download ProtAnneal and edit the `config.py`
file in the `source` directory. Modify the location constants appropriately.

Before performing simulated annealing on a pdb file, two things need to be kept
in mind. First, ProtAnneal can't handle pdb files with more than one subunit,
so you need to get only one subunit out of a pdb file. Next, before starting
annealing, all the buried residues need to be turned into alanine.

To do the first, use the `get_subunit.py` script to get out the subunit.
```
python get_subunit.py <subunit> <pdb_file> > <output_file>
```

Next use the `create_primer.py` to change all the inner residues to alanine.
```
python create_primer.py <pdb_file> > <output_file>
```

Finally, use the `main.py` to do the annealing.
```
python main.py <pdb_file> <number_of_iterations>
```
