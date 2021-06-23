import sys
import numpy as np
import math
import argparse

from abg import pdb_parser, matvec, settings


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", help="target_pdb", required=True)
    parser.add_argument("-helix_1_res", help="residue nums in helix 1", required=True)
    parser.add_argument("-helix_2_res", help="residue nums in helix 2", required=True)
    args = parser.parse_args()
    return args


def get_helix_coord_matrix(mol, target_res_nums):
    atns = ("O3'", "P", "O5'", "C5'", "C4'", "C3'", "C2'", "C1'", "O4'", "O2'")
    atom_coords = {}
    count = 1
    for res_num in target_res_nums:
        mres1 = pdb_parser.getres(mol, res_num)
        for atn in atns:
            mat1 = mres1.getat(atn)
            # add warning
            if mat1 is None:
                continue
            atom_coords[str(count) + "," + atn] = mat1.r
        count += 1

    return atom_coords


def get_only_coords_shared(coords_1, coords_2):
    coords_1_reduced = []
    coords_2_reduced = []

    for key in coords_1.keys():
        coords_1_reduced.append(coords_1[key])
        coords_2_reduced.append(coords_2[key])

    return coords_1_reduced, coords_2_reduced


def compute_abg_from_rotation_matrix(rot):
    a = np.arctan2(rot[2, 1], rot[2, 0])
    g = np.arctan2(rot[1, 2], -rot[0, 2])
    b = np.arctan2(rot[1, 2] / np.sin(g), rot[2, 2])
    # solve the degeneracy issue
    A = np.array([a, a - math.pi, a + math.pi, a + math.pi, a - math.pi])
    B = np.array([b, -b, -b, -b, -b])
    G = np.array([g, g + math.pi, g - math.pi, g + math.pi, g - math.pi])
    M = A ** 2 + B ** 2 + G ** 2
    idx = np.argmin(M)
    a = -A[idx]
    b = B[idx]
    g = -G[idx]
    return (a, b, g)


class ABGComputer(object):
    # return object
    class ABGResults(object):
        def __init__(self, a, b, g, rmsd1, rmsd2, x, y, z):
            self.a, self.b, self.g = a, b, g
            self.rmsd1, self.rmsd2 = rmsd1, rmsd2
            self.x, self.y, self.z = x, y, z

    def __init__(self):
        self.ref_mol = pdb_parser.Mol(settings.AFORM_HELIX_PDB_PATH)
        self.__reset_ref_resi()

    def compute(self, pdb_path, target_resi_1, target_resi_2):
        """
        :param pdb_path: path to pdb with target RNA of interest
        :param target_resi_1: a list of residue numbers that contains the lower helical stem
        :param target_resi_2: a list of residue numbers that contain the upper helical stem
        :return: ABGResults object, containing a,b,g angles and rmsds of the idealized helix compared to both stems

        Example RNA

        G1 - C8
        A2 - U9
        C3 - G10
        A4   |
        A5 - U11
        G6 - C12
        U7 - A14

        target_resi_1 would be [1, 2, 3, 8, 9, 10]
        target_resi_2 would be [5, 6, 7, 11, 12, 14]

        """
        target_resi_1.sort()
        target_resi_2.sort()
        self.__reset_ref_resi()
        target_mol = pdb_parser.Mol(pdb_path)

        # only 2, 4 or 6 residues is a valid arg
        if not (
            len(target_resi_1) == 2
            or len(target_resi_1) == 4
            or len(target_resi_1) == 6
        ):
            raise ValueError(
                "2,4,6 residues are the only valid number of residues allowed"
            )
        if not (
            len(target_resi_2) == 2
            or len(target_resi_2) == 4
            or len(target_resi_2) == 6
        ):
            raise ValueError(
                "2,4,6 residues are the only valid number of residues allowed"
            )

        # check if resi are valid
        for resi in target_resi_1 + target_resi_2:
            if pdb_parser.getres(target_mol, resi) is None:
                raise ValueError(str(resi) + " not contained in pdb: " + pdb_path)

        # set number of residues in the reference helix to those in the target helix regions
        self.ref_resi_1 = self.__set_ref_resi_1(self.ref_resi_1, len(target_resi_1))
        self.ref_resi_2 = self.__set_ref_resi_2(self.ref_resi_2, len(target_resi_2))

        M1, ref1 = self.__get_helix_1_alignment(target_mol, target_resi_1)
        rot, rmsd1 = matvec.lsqfit(M1, ref1)
        M2, ref2 = self.__get_helix_2_alignment(target_mol, target_resi_2, rot, M1)
        rot, rmsd2 = matvec.lsqfit(M2, ref2)

        a, b, g = compute_abg_from_rotation_matrix(rot)
        if len(target_resi_1) == 6:
            last_bp_res_1 = target_resi_1[2:-2]
            last_bp_res_2 = target_resi_2[2:-2]
        elif len(target_resi_1) == 4:
            last_bp_res_1 = target_resi_1[1:-1]
            last_bp_res_2 = target_resi_2[1:-1]
        else:
            last_bp_res_1 = target_resi_1
            last_bp_res_2 = target_resi_2
        first_bp_coords = get_helix_coord_matrix(target_mol, last_bp_res_1)
        avg_bp_coords_1 = np.array([0.0, 0.0, 0.0])
        avg_bp_coords_1 += np.array(first_bp_coords["1,C1'"])
        avg_bp_coords_1 += np.array(first_bp_coords["2,C1'"])
        avg_bp_coords_1 /= 2
        second_bp_coords = get_helix_coord_matrix(target_mol, last_bp_res_2)
        avg_bp_coords_2 = np.array([0.0, 0.0, 0.0])
        avg_bp_coords_2 += np.array(second_bp_coords["1,C1'"])
        avg_bp_coords_2 += np.array(second_bp_coords["2,C1'"])
        avg_bp_coords_2 /= 2
        diff = avg_bp_coords_2 - avg_bp_coords_1
        x, y, z = diff

        return ABGComputer.ABGResults(a, b, g, rmsd1, rmsd2, x, y, z)

    # private methods
    def __reset_ref_resi(self):
        self.ref_resi_1 = [9, 10, 11, 34, 35, 36]
        self.ref_resi_2 = [12, 13, 14, 31, 32, 33]

    def __set_ref_resi_1(self, ref_resi, target_length):
        if target_length == 6:
            return ref_resi
        elif target_length == 4:
            return ref_resi[1:-1]
        elif target_length == 2:
            return ref_resi[2:-2]

    def __set_ref_resi_2(self, ref_resi, target_length):
        if target_length == 6:
            return ref_resi
        elif target_length == 4:
            return ref_resi[:2] + ref_resi[-2:]
        elif target_length == 2:
            return [ref_resi[0], ref_resi[-1]]

    def __get_helix_1_alignment(self, target_mol, target_res_1):
        coords_1 = get_helix_coord_matrix(target_mol, target_res_1)
        coords_2 = get_helix_coord_matrix(self.ref_mol, self.ref_resi_1)
        coords_1_reduced, coords_2_reduced = get_only_coords_shared(coords_1, coords_2)

        M1 = np.array(coords_1_reduced)
        M1 -= np.mean(M1, axis=0)
        ref1 = np.array(coords_2_reduced)
        ref1 -= np.mean(ref1, axis=0)

        return M1, ref1

    def __get_helix_2_alignment(self, target_mol, target_res_2, rot, M1):
        coords_1 = get_helix_coord_matrix(target_mol, target_res_2)
        coords_2 = get_helix_coord_matrix(self.ref_mol, self.ref_resi_2)
        coords_1_reduced, coords_2_reduced = get_only_coords_shared(coords_1, coords_2)

        M2 = np.array(coords_1_reduced)
        M2 -= np.mean(M1, axis=0)
        M2 = np.dot(M2, rot)
        M2 -= np.mean(M2, axis=0)
        ref2 = np.array(coords_2_reduced)
        ref2 -= np.mean(ref2, axis=0)
        return M2, ref2


def main():
    args = parse_args()
    target_res_1 = [int(x) for x in args.helix_1_res.split(",")]
    target_res_2 = [int(x) for x in args.helix_2_res.split(",")]
    abg_computer = ABGComputer()
    r = abg_computer.compute(args.pdb, target_res_1, target_res_2)
    print(
        "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f"
        % (
            r.a * 180.0 / math.pi,
            r.b * 180.0 / math.pi,
            r.g * 180.0 / math.pi,
            r.rmsd1,
            r.rmsd2,
            r.x,
            r.y,
            r.z
        )
    )


if __name__ == "__main__":
    main()
