import sys
import bz2
import pickle
from synspace.data import get_reactions, make_blocks


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: python build_blocks.py <smi_path> <additional smi paths>")
        sys.exit(1)
    rxns = get_reactions()
    all_blocks = dict()
    for path in sys.argv[1:]:
        blocks = make_blocks(rxns, path)
        # merge them into all_blocks
        for name, rxn in blocks.items():
            if name not in all_blocks:
                all_blocks[name] = [[] for _ in rxn]
            for i, rxn_block in enumerate(rxn):
                all_blocks[name][i].extend(rxn_block)
    # now sort them by molecular weight
    for name, rxn in all_blocks.items():
        for i, rxn_block in enumerate(rxn):
            rxn_block.sort(key=lambda x: x.GetNumHeavyAtoms())
    # write them to bz2 file
    with bz2.open("blocks.pk.bz2", "wb") as f:
        pickle.dump(all_blocks, f)
    print("done")
