import os

def get_lib_path():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-1])
    return base_dir

AFORM_HELIX_PDB_PATH = get_lib_path() + "/resources/iAformRNA.pdb"