import requests
import random
import json
import os
from .data import get_reactions, get_blocks
from rdkit import Chem
from .utils import remove_dups, sort_mols, extract, flatten, atom_match

def chemical_space(mol, steps=(1, 1), threshold=0.2, blocks=None, rxns=None,  use_mannifold=None,  strict=False, samples=10):
    if type(mol) == str:
        mol = Chem.MolFromSmiles(mol)
    if type(steps) == int:
        steps = (0, steps)
    if use_mannifold is None:
        use_mannifold = os.environ.get("POSTERA_API_KEY") is not None
    if use_mannifold:
        mols, props = mannifold_retro(mol)
    else:
        mols, props = retro(mol, rxns=rxns, strict=strict)
        for _ in range(steps[0] - 1):
            to_add = []
            for m,p in zip(mols, props):
                ms, ps = retro(m, rxns=rxns, strict=strict, start_props=p)
                to_add.append((ms, ps))
            for m,p in to_add:
                mols.extend(m)
                props.extend(p)
        mols,props = remove_dups(mols, props)
    for _ in range(steps[1]):
        to_add = []
        for m,p in zip(mols, props):
            ms, ps = forward(m, blocks=blocks, samples=samples, rxns=rxns,
                            threshold=threshold, strict=strict, start_props=p)
            to_add.append((ms, ps))
        for m,p in to_add:
            mols.extend(m)
            props.extend(p)
    print(mol, mols[0])
    return sort_mols(*remove_dups(mols, props), mol, threshold=threshold)

def mannifold_retro(query_mol):
    # try to get the API Key
    api_key = os.environ.get("POSTERA_API_KEY")
    if api_key is None:
        raise RuntimeError("Please set the POSTERA_API_KEY environment variable.")
    smi = Chem.MolToSmiles(query_mol)
    response = requests.post(
        "https://api.postera.ai/api/v1/retrosynthesis/",
        headers={
            "X-API-KEY": api_key,
        },
        json={
            "smiles": smi,
            "maxSearchDepth": 4,
            # 'catalogs': ['generic']
        },
    )
    response = json.loads(response.text)
    mols, props = [], []
    for route in response["routes"]:
        for mol in route["molecules"]:
            # if not mol['isBuildingBlock']:
            s = mol["smiles"]
            mols.append(Chem.MolFromSmiles(s))
            props.append({"rxn-name": "mannifold", "rxn": "", "match": ()})
    mols, props = remove_dups(*sort_mols(mols, props, mol, threshold=0.2))
    for m, p in zip(mols, props):
        match = atom_match(query_mol, m)
        p["match"] = tuple(match)
    return mols, props


def retro(mol, threshold=0.5, strict=True, start_props=None, rxns=None):
    if rxns is None:
        rxns = get_reactions()
    out = []
    props = []
    smi = Chem.MolToSmiles(mol)
    for n, (_, br) in rxns.items():
        t = br.GetReactantTemplate(0)
        # collapse it because I can't think of a better idea
        match = flatten(mol.GetSubstructMatches(t))
        p = br.RunReactants((mol,))
        if len(p) == 1 or (len(p) > 0 and not strict):
            for m, s in extract(p):
                out.append(m)
                props.append(
                    {
                        "rxn-name": "retro-" + n.replace("_", " "),
                        "rxn": f"{s}>>{smi}",
                        "match": match,
                    }
                )
                if start_props:
                    props[-1]["rxn-name"] = (
                        start_props["rxn-name"] + "," + props[-1]["rxn-name"]
                    )
                    props[-1]["rxn"] = start_props["rxn"] + ":" + props[-1]["rxn"]
                    props[-1]["match"] = start_props["match"] + props[-1]["match"]
    return remove_dups([mol] + out, [{"rxn-name": "", "rxn": "", "match": ()}] + props)


def forward(mol, blocks=None, rxns=None, samples=1000, threshold=0.5, strict=True, start_props=None):
    if blocks is None:
        blocks = get_blocks()
    if rxns is None:
        rxns = get_reactions()
    out = []
    props = []
    smi = Chem.MolToSmiles(mol)
    for n, (fr, _) in rxns.items():
        for i in range(fr.GetNumReactantTemplates()):
            t = fr.GetReactantTemplate(i)
            # collapse it because I can't think of a better idea
            match = flatten(mol.GetSubstructMatches(t))
            if match:
                break
        if not match:
            continue
        for _ in range(samples):
            reactants = []
            success = True
            for j in range(fr.GetNumReactantTemplates()):
                if i != j:
                    if len(blocks[n][j]) == 0:
                        success = False
                        break
                    selection = blocks[n][j]
                    b = random.choice(selection)
                    reactants.append(b)
                else:
                    reactants.append(mol)
            if not success:
                continue
            p = fr.RunReactants(reactants)
            if len(p) == 1 or (len(p) > 0 and not strict):
                rsmi = ".".join([Chem.MolToSmiles(r) for r in reactants])
                for m, s in extract(p):
                    out.append(m)
                    props.append(
                        {
                            "rxn-name": n.replace("_", " "),
                            "rxn": f"{s}>>{smi}",
                            "match": match,
                        }
                    )
                    if start_props:
                        props[-1]["rxn-name"] = (
                            start_props["rxn-name"] + "," + props[-1]["rxn-name"]
                        )
                        props[-1]["rxn"] = start_props["rxn"] + ":" + props[-1]["rxn"]
                        props[-1]["match"] = start_props["match"] + props[-1]["match"]
    return [mol] + out, [{"rxn-name": "", "rxn": "", "match": ()}] + props
