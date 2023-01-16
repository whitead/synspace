import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.rdFMCS as FMCS
import rdkit.rdBase as rdBase
import rdkit.Chem.rdmolops as rdmolops
import functools
import numpy as np
from rdkit.DataStructs.cDataStructs import BulkTanimotoSimilarity
from .reos import REOS


def flatten(s):
    if len(s) == 0:
        return s
    return functools.reduce(lambda a, b: a + b, s)


def moldiff(template, query):
    """Compare the two rdkit molecules.
    :param template: template molecule
    :param query: query molecule
    :return: list of modified atoms in query, list of modified bonds in query
    """
    r = FMCS.FindMCS([template, query], timeout=5)
    if r is None:
        return [], []
    substructure = Chem.MolFromSmarts(r.smartsString)
    raw_match = query.GetSubstructMatches(substructure)
    if len(raw_match) == 0:
        return [], []
    template_match = template.GetSubstructMatches(substructure)
    # flatten it
    match = list(raw_match[0])
    template_match = list(template_match[0])

    # need to invert match to get diffs
    inv_match = [i for i in range(query.GetNumAtoms()) if i not in match]

    # get bonds
    bond_match = []
    for b in query.GetBonds():
        if b.GetBeginAtomIdx() in inv_match or b.GetEndAtomIdx() in inv_match:
            bond_match.append(b.GetIdx())

    # now get bonding changes from deletion
    def neigh_hash(a):
        return "".join(sorted([n.GetSymbol() for n in a.GetNeighbors()]))

    for ti, qi in zip(template_match, match):
        if neigh_hash(template.GetAtomWithIdx(ti)) != neigh_hash(
            query.GetAtomWithIdx(qi)
        ):
            inv_match.append(qi)

    return inv_match, bond_match


def atom_match(template, query):
    r = FMCS.FindMCS([template, query], timeout=5)
    if r is None:
        return ()
    substructure = Chem.MolFromSmarts(r.smartsString)
    raw_match = query.GetSubstructMatches(substructure)
    if len(raw_match) == 0:
        return ()
    template_match = template.GetSubstructMatches(substructure)
    # flatten it
    match = list(raw_match[0])
    template_match = list(template_match[0])

    # need to invert match to get diffs
    template_idx = [i for i in range(template.GetNumAtoms()) if i not in match]
    return template_idx


def match_cluster(mols, props, cmax=3):
    clusters = dict()
    rmol = []
    rprop = []
    for m, p in zip(mols, props):
        mi = p["match"]
        if mi not in clusters:
            clusters[mi] = 1
            rmol.append(m)
            rprop.append(p)
        elif clusters[mi] < cmax:
            clusters[mi] += 1
            rmol.append(m)
            rprop.append(p)
    return rmol, rprop


reos = REOS()


def reos_filter(mols, props):
    result = [
        (x, p) for x, p in zip(mols, props) if reos.process_mol(x) == ("ok", "ok")
    ]
    return [r[0] for r in result], [r[1] for r in result]


def extract(products):
    for p in products:
        lp = list(p)
        lp.sort(key=lambda m: m.GetNumAtoms())
        success = True
        with rdBase.BlockLogs():
            for l in lp:
                try:
                    Chem.SanitizeMol(l)
                except ValueError:
                    success = False
                if "*" in Chem.MolToSmiles(l):
                    success = False
                # strong condition
                if rdmolops.GetFormalCharge(l) != 0:
                    success = False
                if not success:
                    break
        if success:
            yield lp[-1], ".".join([Chem.MolToSmiles(l) for l in lp])


def remove_dups(mols, props):
    if len(mols) == 0:
        return mols, props
    smis = set()
    out = []
    out_props = []
    for m, p in zip(mols, props):
        s = Chem.MolToSmiles(m, canonical=True)
        if s in smis:
            continue
        out.append(m)
        out_props.append(p)
        smis.add(s)
    return out, out_props


def get_fp(mol):
    return AllChem.GetMorganFingerprint(mol, 2, useFeatures=True)
    # return AllChem.RDKFingerprint(mol)


def sort_mols(mols, props, base, threshold=0):
    if len(mols) == 0:
        return mols, props
    fps = [get_fp(m) for m in mols]
    base_fp = get_fp(base)
    M = BulkTanimotoSimilarity(base_fp, fps)
    order = np.argsort(M)
    order = [i for i in order[::-1] if M[i] > threshold or abs(1 - M[i]) < 0.01]
    for p in props:
        if "similarity" in p:
            del p["similarity"]
    props = [dict(**props[i], similarity=np.round(M[i], 2)) for i in order]
    return [mols[i] for i in order], props


def find_prop(smi, mols, props):
    # canonicalize smiles
    smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)
    for m, p in zip(mols, props):
        if Chem.MolToSmiles(m, canonical=True) == smi:
            return p
    return None
