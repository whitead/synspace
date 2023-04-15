import requests
import numpy as np
import json
import os
from .data import get_reactions, get_blocks
from rdkit import Chem
from .utils import remove_dups, sort_mols, extract, flatten, atom_match, reos_filter


def chemical_space(
    mol,
    steps=(1, 1),
    threshold=0.2,
    blocks=None,
    rxns=None,
    use_mannifold=None,
    strict=None,
    nblocks=25,
    num_samples=250,
    filter=True,
    _pbar=None,
):
    if type(mol) == str:
        mol = Chem.MolFromSmiles(mol)
    if type(steps) == int:
        steps = (0, steps)
    mols, props = None, None
    if use_mannifold is None:
        use_mannifold = os.environ.get("POSTERA_API_KEY") is not None
    if use_mannifold:
        if _pbar:
            _pbar.set_description("⚗️Synspace Retrosynthesis (Mannifold)⚗️")
        mols, props = mannifold_retro(mol)
    else:
        if _pbar and steps[0] > 0:
            _pbar.set_description("⚗️Synspace Retrosynthesis...⚗️")
        if steps[0] > 0:
            mols, props = retro(mol, rxns=rxns, strict=False if strict is None else strict)
        for _ in range(steps[0]-1):
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
    if _pbar:
        _pbar.set_description("⚗️Forward synthesis...⚗️")
    if mols is None:
        mols, props = [mol], [{"rxn-name": "", "rxn": "", "match": ()}]
    for i in range(steps[1]):
        to_add = []
        for m, p in zip(mols, props):
            ms, ps = forward(
                m,
                blocks=blocks,
                samples=nblocks,
                rxns=rxns,
                threshold=threshold,
                strict=True if strict is None else strict,
                start_props=p,
            )
            to_add.append((ms, ps))
        for m, p in to_add:
            mols.extend(m)
            props.extend(p)
        # poor estimate that we can allow 1/steps *  threshold for each step
        mols, props = sort_mols(
            *remove_dups(mols, props),
            mol,
            threshold=threshold / steps[1] if i < steps[1] - 1 else threshold,
        )
        if filter:
            mols, props = reos_filter(mols, props)
        mols, props = mols[:num_samples], props[:num_samples]
        if _pbar:
            _pbar.update(len(mols))
    return mols, props


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
            # sometimes happens?
            s = s.replace("~", "")
            mols.append(Chem.MolFromSmiles(s))
            props.append({"rxn-name": "mannifold", "rxn": "", "match": ()})
    if len(mols) < 2:
        raise RuntimeError("No retrosynthetic routes found.")
    mols, props = remove_dups(mols, props)
    for m, p in zip(mols, props):
        match = atom_match(query_mol, m)
        p["match"] = tuple(match)

    return mols, props


def merge_props(p0, p1):
    if p0["rxn"] != "":
        p1["rxn"] = p0["rxn"] + ":" + p1["rxn"]
        p1["rxn-name"] = p0["rxn-name"] + "," + p1["rxn-name"]
    p1["match"] = p0["match"] + p1["match"]


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
                    merge_props(start_props, props[-1])
    return remove_dups([mol] + out, [{"rxn-name": "", "rxn": "", "match": ()}] + props)


def prob_transform(x):
    x = x - np.min(x)
    x = x / np.max(x)
    p = np.exp(-x)
    return p / np.sum(p)


def forward(
    mol,
    blocks=None,
    rxns=None,
    samples=1000,
    threshold=0.5,
    strict=True,
    start_props=None,
):
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
                    # with decreasing probability, since ordered by weight
                    if len(blocks[n][j]) == 1:
                        b = selection[0]
                    else:
                        b = np.random.choice(
                            selection, p=prob_transform(np.arange(len(selection)))
                        )
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
                            "rxn": f"{rsmi}>>{s}",
                            "match": match,
                        }
                    )
                    if start_props:
                        merge_props(start_props, props[-1])
    return [mol] + out, [{"rxn-name": "", "rxn": "", "match": ()}] + props
