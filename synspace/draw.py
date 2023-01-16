import rdkit.Chem.Draw as Draw
import rdkit.Chem.rdChemReactions as rdChemReactions
import skunk


def draw_rxn(s, names=None):
    if type(s) is dict:
        names = s["rxn-name"]
        s = s["rxn"]
    rxns = s.split(":")
    svgs = []
    for rxn in rxns:
        rxn = rdChemReactions.ReactionFromSmarts(rxn, useSmiles=True)
        svgs.append(Draw.ReactionToImage(rxn, useSVG=True, subImgSize=(200, 200)))
    labels = None
    if names is not None:
        labels = names.split(",")
    return skunk.layout_svgs(svgs, labels=labels)
