import bz2
import pickle
import os

_BLOCKS = None
_REACTIONS = None

def get_reactions():
    global _REACTIONS
    if _REACTIONS is None:
        from importlib_resources import files
        import syngen.rxn_data
        import json
        reactions_path = files(syngen.rxn_data).joinpath("rxns.json")
        with bz2.open(reactions_path, "rb") as f:
            _REACTIONS = pickle.load(f)
    return _REACTIONS

def get_blocks():
    global _BLOCKS
    if _BLOCKS is None:
        from importlib_resources import files
        import syngen.rxn_data
        blocks_path = files(syngen.rxn_data).joinpath("blocks.pk.bz2")
        with bz2.open(blocks_path, "rb") as f:
            _BLOCKS = pickle.load(f)
    return _BLOCKS

def make_blocks(reactions, smi_path):
    import tqdm
    from rdkit import Chem
    results = dict()
    matchers = dict()
    for n, (fr, _) in reactions.items():
        matchers[n] = []
        results[n] = []
        for i in range(fr.GetNumReactantTemplates()):
            t = fr.GetReactantTemplate(i)
            matchers[n].append(t)
            results[n].append([])
    with Chem.SmilesMolSupplier(
        smi_path, titleLine=False, delimiter="\t"
    ) as mol_blocks:
        for i, x in enumerate(tqdm.tqdm(mol_blocks)):
            if x is None:
                continue
            for n, ms in matchers.items():
                for i, m in enumerate(ms):
                    if x.HasSubstructMatch(m):
                        results[n][i].append(x)
    return results

