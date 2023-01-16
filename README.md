# synspace

**This is early-stage code that is in progress. It is in flux**

This package generates a local chemical space around a given molecule using retro and forward synthesis rules. The reactions used are the 50 robust medchem reactions proposed by [Hartenfeller et al.](https://pubs.acs.org/doi/10.1021/ci200379p). The retrosynthesis is done either via [PostEra Mannifold](https://postera.ai/) if you have an API key, or by reversing the 50 robust reactions. The purchasable building blocks come from the [Purchasable Mcule supplier building block catalogs](https://mcule.com/database/). All of these things can be customized though. 

<img src="https://user-images.githubusercontent.com/908389/212747862-577182f2-a880-46ec-a0e1-33fe7f5d9987.png" alt="chemical space showing a molecule space that is predicted not to cross the blood brain barrier along with three synthetically feasible modifications" width=500>


## Installation

```sh
pip install synspace
```

## Usage

Generate local chemical space given a SMILES string
```py
mols, props = synspace.chemical_space('CCC=O')
```
`props` contains information like the synthesis route for the molecules. Note that all synthesis routes are relative to the given molecule (it is assumed to be synthetically feasible). 


## Citation

TODO

Also, this idea is similar to [Dolfus et al.](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00246)

## NOTICE

This product includes software developed by Pat Walters (MIT Licensed)
https://github.com/PatWalters/useful_rdkit_utils
Copyright (c) 2022 Pat Walters
