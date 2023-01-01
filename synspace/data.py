import bz2
import pickle
import os
from rdkit import rdBase
import rdkit.Chem as Chem
import rdkit.Chem.rdChemReactions as rdChemReactions

_BLOCKS = None
_REACTIONS = None


def reverse(rxn):
    with rdBase.BlockLogs():
        rxn2 = rdChemReactions.ChemicalReaction()
        for i in range(rxn.GetNumReactantTemplates()):
            rxn2.AddProductTemplate(rxn.GetReactantTemplate(i))
        for i in range(rxn.GetNumProductTemplates()):
            rxn2.AddReactantTemplate(rxn.GetProductTemplate(i))
        rxn2.Initialize()
    return rxn2


def get_reactions():
    global _REACTIONS
    if _REACTIONS is None:
        from importlib_resources import files
        import synspace.rxn_data
        import json

        reactions_path = files(synspace.rxn_data).joinpath("rxns.json")
        with open(reactions_path, "r") as f:
            data = json.load(f)
        _REACTIONS = dict()
        for n, s in data.items():
            r = rdChemReactions.ReactionFromSmarts(s)
            rr = reverse(r)
            _REACTIONS[n] = (r, rr)
    return _REACTIONS


def get_blocks():
    global _BLOCKS
    if _BLOCKS is None:
        from importlib_resources import files
        import synspace.rxn_data

        blocks_path = files(synspace.rxn_data).joinpath("blocks.pk.bz2")
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
