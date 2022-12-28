import syngen
import rdkit.Chem as Chem

def test_chemical_space():
    smi = 'Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C'
    # smi = 'CC1([C@@H]2[C@H]1[C@H](N(C2)C(=O)[C@H](C(C)(C)C)NC(=O)C(F)(F)F)C(=O)N[C@@H](C[C@@H]3CCNC3=O)C#N)C'
    #smi = 'O=C(C)Oc1ccccc1C(=O)O'
    #smi = 'CC(C)c4nc(CN(C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](Cc1ccccc1)C[C@H](O)[C@H](Cc2ccccc2)NC(=O)OCc3cncs3)cs4'
    mols, props = syngen.chemical_space(smi)
    assert Chem.MolToSmiles(mols[0]) == smi
    assert len(mols) == len(props)
    assert len(mols) > 100