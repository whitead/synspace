import synspace
import rdkit.Chem as Chem


def test_chemical_space():
    smi = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
    # smi = 'CC1([C@@H]2[C@H]1[C@H](N(C2)C(=O)[C@H](C(C)(C)C)NC(=O)C(F)(F)F)C(=O)N[C@@H](C[C@@H]3CCNC3=O)C#N)C'
    # smi = 'O=C(C)Oc1ccccc1C(=O)O'
    mols, props = synspace.chemical_space(smi, use_mannifold=False)
    assert Chem.MolToSmiles(mols[0], canonical=True) == Chem.MolToSmiles(
        Chem.MolFromSmiles(smi), canonical=True
    )
    assert len(mols) == len(props)
    assert len(mols) > 100

    smi = "CC(C)c4nc(CN(C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](Cc1ccccc1)C[C@H](O)[C@H](Cc2ccccc2)NC(=O)OCc3cncs3)cs4"
    mols, props = synspace.chemical_space(
        smi, use_mannifold=False, steps=(2, 1), nblocks=25, threshold=0.5
    )
    assert Chem.MolToSmiles(mols[0], canonical=True) == Chem.MolToSmiles(
        Chem.MolFromSmiles(smi), canonical=True
    )
    assert len(mols) > 100


def test_mannifold():
    smi = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
    mols, props = synspace.chemical_space(smi, use_mannifold=True)
    assert Chem.MolToSmiles(mols[0], canonical=True) == Chem.MolToSmiles(
        Chem.MolFromSmiles(smi), canonical=True
    )
    assert len(mols) == len(props)
    assert len(mols) > 100


def test_find_prop():
    smi = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
    # smi = 'CC1([C@@H]2[C@H]1[C@H](N(C2)C(=O)[C@H](C(C)(C)C)NC(=O)C(F)(F)F)C(=O)N[C@@H](C[C@@H]3CCNC3=O)C#N)C'
    # smi = 'O=C(C)Oc1ccccc1C(=O)O'
    mols, props = synspace.chemical_space(smi, use_mannifold=False)
    assert synspace.find_prop(smi, mols, props) == props[0]


def test_draw_rxn():
    smi = "O=C(C)Oc1ccccc1C(=O)O"
    mols, props = synspace.chemical_space(smi, use_mannifold=False)
    p = props[-1]
    synspace.draw_rxn(p["rxn"])
    synspace.draw_rxn(p["rxn"], p["rxn-name"])
    synspace.draw_rxn(p)


def test_draw_rxn():
    smi = "O=C(C)Oc1ccccc1C(=O)O"
    mols, props = synspace.chemical_space(smi, use_mannifold=True)
    p = props[-1]
    synspace.draw_rxn(p["rxn"])
    synspace.draw_rxn(p["rxn"], p["rxn-name"])
    synspace.draw_rxn(p)
