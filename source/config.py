"""
Contains the constants needed for proper operation
of the program.
"""

# Location constants
residue_path = "/home/bolt/protanneal/database/residues/{0}_naiveH.pdb"
naccess_location = "/home/bolt/protanneal/third-party/binaries/naccess2.1.1/naccess"
reduce_location = "/home/bolt/protanneal/third-party/binaries/reduce"
Zdope_command = "/usr/local/bin/Zdope"
dunbrack_file = "/home/bolt/protanneal/other/dunbrack2.csv"
output_directory = "../workspace/temp_pdb"
result_directory = "../workspace/result_pdb"


# Other constants
allowed_residues = ["VAL", "LEU", "ILE", "MET", "ALA"]
ATOMS_TO_LOOK = (" C",)# " H", " S")
THRESHOLD = 10 # Angstroms
TO_CONSIDER = (' CD1',)# " CD2")
RESIDUE = "LEU"
