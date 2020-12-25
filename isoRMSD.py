#!/usr/bin/env python
"""
calculates RMSD differences between 2 conformation with different atom names.

@author: JC <yangjincai@nibs.ac.cn>
"""
# %%
import os
from os.path import split
import sys
import math

# rdkit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import AlignMol


def GetBestRMSD(probe, ref, refConfId=-1, probeConfId=-1, maps=None):
    """Returns the optimal RMS for aligning two molecules, taking
    symmetry into account. As a side-effect, the probe molecule is
    left in the aligned state.

    Arguments:
      - ref: the reference molecule
      - probe: the molecule to be aligned to the reference
      - refConfId: (optional) reference conformation to use
      - probeConfId: (optional) probe conformation to use
      - maps: (optional) a list of lists of (probeAtomId,refAtomId)
        tuples with the atom-atom mappings of the two molecules.
        If not provided, these will be generated using a substructure
        search.

    Note:
    This function will attempt to align all permutations of matching atom
    orders in both molecules, for some molecules it will lead to 'combinatorial
    explosion' especially if hydrogens are present.
    Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing the
    atom order.

    """
    # When mapping the coordinate of probe will changed!!!
    ref.pos = orginXYZ(ref)
    probe.pos = orginXYZ(probe)

    if not maps:
        matches = ref.GetSubstructMatches(probe, uniquify=False)
        if not matches:
            # raise ValueError(
            #     "mol %s does not match mol %s"
            #     % (ref.GetProp("_Name"), probe.GetProp("_Name"))
            # )
            pass
        if len(matches) > 1e6:
            warnings.warn(
                "{} matches detected for molecule {}, this may lead to a performance slowdown.".format(
                    len(matches), probe.GetProp("_Name")
                )
            )
        maps = [list(enumerate(match)) for match in matches]
    bestRMS = 1000.0
    bestRMSD = 1000.0
    bestMap = None
    for amap in maps:
        rms = AlignMol(probe, ref, probeConfId, refConfId, atomMap=amap)
        rmsd = RMSD(probe, ref, amap)
        if rmsd < bestRMSD:
            bestRMSD = rmsd
        if rms < bestRMS:
            bestRMS = rms
            bestMap = amap

    # finally repeate the best alignment :
    if bestMap != amap:
        AlignMol(probe, ref, probeConfId, refConfId, atomMap=bestMap)

    return bestRMS, bestRMSD


# Map is probe -> ref
# [(1:3),(2:5),...,(10,1)]
def RMSD(probe, ref, amap):
    rmsd = 0.0
    # print(amap)
    atomNum = ref.GetNumAtoms() + 0.0
    for (pi, ri) in amap:
        posp = probe.pos[pi]
        posf = ref.pos[ri]
        rmsd += dist_2(posp, posf)
    rmsd = math.sqrt(rmsd / atomNum)
    return rmsd


def dist_2(atoma_xyz, atomb_xyz):
    dis2 = 0.0
    for i, j in zip(atoma_xyz, atomb_xyz):
        dis2 += (i - j) ** 2
    return dis2


def orginXYZ(mol):
    mol_pos = {}
    for i in range(0, mol.GetNumAtoms()):
        pos = mol.GetConformer().GetAtomPosition(i)
        mol_pos[i] = pos
    return mol_pos


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    from oddt.toolkits import rdk as toolkit

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="reference mol in .pdb, .mol2 or .sdf format",
    )
    parser.add_argument(
        "-p", "--probe", required=True, help="probe mols in .pdb, .mol2, or .sdf"
    )
    parser.add_argument(
        "-o",
        "--output_csv",
        default="rmsd.csv",
        help="output rmsd in a csv file, default: rmsd.csv",
    )
    args = parser.parse_args()

    ref_fmt = args.reference.split(".")[-1]
    ref_oddt = next(toolkit.readfile(ref_fmt, args.reference))
    ref_rdk = Chem.RemoveHs(ref_oddt.Mol)

    probe_fmt = args.probe.split(".")[-1]
    probe_oddt_supp = toolkit.readfile(probe_fmt, args.probe)

    for probe_oddt in probe_oddt_supp:
        if probe_oddt is None:
            continue
        probe_rdk = Chem.RemoveHs(probe_oddt.Mol)
        ref = AllChem.AssignBondOrdersFromTemplate(probe_rdk, ref_rdk)
        rms, rmsd = GetBestRMSD(probe_rdk, ref)
        print("\nBest_RMSD: %.3f\nBest_Not_Fit_RMSD: %.3f\n" % (rms, rmsd))
