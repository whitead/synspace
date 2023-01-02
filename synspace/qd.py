from .synspace import forward, retro, mannifold_retro
from .utils import get_fp, remove_dups, extract_one, flatten
import os
import vdict
import numpy as np
from .data import get_reactions, get_blocks
from rdkit import Chem
import random


def embed_mols(mols, extra_x, extra_d, dims=8):
    from sklearn.decomposition import PCA
    from rdkit.DataStructs.cDataStructs import BulkTanimotoSimilarity

    fps = [get_fp(m) for m in mols]
    M = np.array([BulkTanimotoSimilarity(f, fps) for f in fps])
    dist_mat = 1 - M
    dist_mat = np.concatenate([dist_mat, np.array(extra_x)], axis=1)
    pca = PCA(n_components=dims)
    proj_dmat = pca.fit_transform(dist_mat)
    # ensure each column is bounded from 0 to 1
    proj_dmat = (proj_dmat - proj_dmat.min(axis=0)) / (
        proj_dmat.max(axis=0) - proj_dmat.min(axis=0)
    )
    result = vdict.vdict()
    for i, m in enumerate(mols):
        result[proj_dmat[i]] = m, extra_d[i]
    return result


def reverse_blocks(blocks):
    out = {}
    for name, r_pos in blocks.items():
        for i, r_pos_i in enumerate(r_pos):
            for block in r_pos_i:
                for r in block:
                    out.setdefault(r, []).append((name, i))
    return out


def qd_chemical_space(
    mol,
    steps=(1, 1),
    threshold=0.2,
    blocks=None,
    rxns=None,
    use_mannifold=None,
    strict=None,
    nblocks=25,
    num_samples=250,
    _pbar=None,
    embed_dim=8,
):
    if type(mol) == str:
        mol = Chem.MolFromSmiles(mol)
    if type(steps) == int:
        steps = (0, steps)
    if blocks is None:
        blocks = get_blocks()
    if rxns is None:
        rxns = get_reactions()
    if use_mannifold is None:
        use_mannifold = os.environ.get("POSTERA_API_KEY") is not None
    if use_mannifold:
        if _pbar:
            _pbar.set_description("⚗️Synspace Retrosynthesis (Mannifold)⚗️")
        mols, props = mannifold_retro(mol)
    else:
        if _pbar:
            _pbar.set_description("⚗️Synspace Retrosynthesis...⚗️")
        mols, props = retro(mol, rxns=rxns, strict=False if strict is None else strict)
        for _ in range(steps[0] - 1):
            to_add = []
            for m, p in zip(mols, props):
                ms, ps = retro(
                    m,
                    rxns=rxns,
                    strict=False if strict is None else strict,
                    start_props=p,
                )
                to_add.append((ms, ps))
            for m, p in to_add:
                mols.extend(m)
                props.extend(p)
                if _pbar:
                    _pbar.update(len(mols))
        mols, props = remove_dups(mols, props)
    eretro = embed_mols(mols, dims=embed_dim)
    rblocks = reverse_blocks(blocks)
    # need to get the extra_x, extra_d
    eblocks = embed_mols(list(rblocks.keys()), dims=embed_dim)

    def _simulate(x):
        x0 = x[:embed_dim]
        m = eretro[x0]
        result = 0, ()
        for i in range(steps[1]):
            x1 = x[embed_dim * (i + 1) : embed_dim * (i + 2)]
            m1 = eblocks[x1]
            # TODO: Cannot work because random choice is not deterministic
            name, pos = random.choice(rblocks[m1])
            rxn = rxns[name][0]
            reactants = [None for _ in range(rxn.GetNumReactantTemplates())]
            reactants[pos] = m1
            j = reactants.index(None)
            reactants[j] = m
            match = flatten(mol.GetSubstructMatches(t))
            p = rxn.RunReactants(reactants)
            if len(p) == 1 or (len(p) > 0 and not strict):
                m, _ = extract_one(p)

    return mols, props
