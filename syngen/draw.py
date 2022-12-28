import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.Draw as Draw
import rdkit.Chem.rdChemReactions as rdChemReactions
import matplotlib as mpl
from .utils import moldiff
import skunk


delete_color = mpl.colors.to_rgb("#F06060")
modify_color = mpl.colors.to_rgb("#1BBC9B")


def draw(m, dos, ref, size=(-1, -1)):
    AllChem.Compute2DCoords(m)
    aidx, bidx = moldiff(ref, m)
    AllChem.GenerateDepictionMatching2DStructure(ref, m, acceptFailure=True)
    kwargs = dict(
        highlightAtoms=aidx,
        highlightBonds=bidx,
        highlightAtomColors={k: modify_color for k in aidx}
        if len(bidx) > 0
        else {k: delete_color for k in aidx},
        highlightBondColors={k: modify_color for k in bidx}
        if len(bidx) > 0
        else {k: delete_color for k in bidx},
    )
    d = Draw.rdMolDraw2D.MolDraw2DSVG(*size)
    d.SetDrawOptions(dos)
    d.DrawMolecule(m, **kwargs)
    d.FinishDrawing()
    return d.GetDrawingText()


def draw_grid(mols, props=None):
    dos = Draw.MolDrawOptions()
    dos.useBWAtomPalette()
    # dos.minFontSize = 12
    if props is None:
        props = [{}] * len(mols)
    svgs = []
    for m, p in zip(mols, props):
        si = []
        labels = []
        si.append(draw(m, dos, mols[0]))
        if "similarity" in p:
            labels.append(str(p["similarity"]))
        if "rxn" in p and p["rxn"]:
            si.append(draw_rxn(p["rxn"]))
            labels.append(p["rxn-name"])
            svgs.append(skunk.layout_svgs(si, labels=labels))
        else:
            svgs.append(si[0])
    return skunk.layout_svgs(svgs)


def draw_rxn(s):
    rxn = rdChemReactions.ReactionFromSmarts(s, useSmiles=True)
    return Draw.ReactionToImage(rxn, useSVG=True)
