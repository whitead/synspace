# MIT License

# Copyright (c) 2021 PatWalters

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import sys
from pathlib import Path
import importlib.resources as pkg_resources

import pandas as pd
import os
from rdkit import Chem


class REOS:
    """REOS - Rapid Elimination Of Swill\n
    Walters, Ajay, Murcko, "Recognizing molecules with druglike properties"\n
    Curr. Opin. Chem. Bio., 3 (1999), 384-387\n
    https://doi.org/10.1016/S1367-5931(99)80058-1
    """

    def __init__(self, active_rules=None):
        if active_rules is None:
            active_rules = ["Glaxo"]

        # Use the bundled alert_collection.csv file
        try:
            # For Python 3.9+
            with pkg_resources.as_file(
                pkg_resources.files("synspace.reos.data").joinpath(
                    "alert_collection.csv"
                )
            ) as path:
                self.rule_path = path
                self.rule_df = pd.read_csv(self.rule_path)
        except (ImportError, AttributeError):
            # Fallback for older Python versions
            package_dir = os.path.dirname(os.path.abspath(__file__))
            self.rule_path = os.path.join(package_dir, "data", "alert_collection.csv")
            self.rule_df = pd.read_csv(self.rule_path)

        self.read_rules(active_rules)

    def parse_smarts(self):
        """Parse the SMARTS strings in the rules file to molecule objects and check for validity

        :return: True if all SMARTS are parsed, False otherwise
        """
        smarts_mol_list = []
        smarts_are_ok = True
        for idx, smarts in enumerate(self.rule_df.smarts, 1):
            mol = Chem.MolFromSmarts(smarts)
            if mol is None:
                smarts_are_ok = False
                print(f"Error processing SMARTS on line {idx}", file=sys.stderr)
            smarts_mol_list.append(mol)
        self.rule_df["pat"] = smarts_mol_list
        return smarts_are_ok

    def read_rules(self, active_rules=None):
        """Read a rules file

        :param rules_file: name of the rules file
        :param active_rules: list of active rule sets, all rule sets are used if
            this is None
        :return: None
        """
        if self.parse_smarts():
            self.active_rule_df = self.rule_df.query("rule_set_name in @active_rules")
            if len(self.active_rule_df) == 0:
                available_rules = sorted(list(self.rule_df["rule_set_name"].unique()))
                raise ValueError(
                    f"Supplied rules: {active_rules} not available. Please select from {available_rules}"
                )

        else:
            print(
                "Error reading rules, please fix the SMARTS errors reported above",
                file=sys.stderr,
            )
            sys.exit(1)
        if active_rules is not None:
            self.active_rule_df = self.rule_df.query("rule_set_name in @active_rules")
        else:
            self.active_rule_df = self.rule_df

    def set_active_rule_sets(self, active_rules=None):
        """Set the active rule set(s)

        :param active_rules: list of active rule sets
        :return: None
        """
        assert active_rules
        self.active_rule_df = self.rule_df.query("rule_set_name in @active_rules")

    def get_available_rule_sets(self):
        """Get the available rule sets in rule_df

        :return: a list of available rule sets
        """
        return self.rule_df.rule_set_name.unique()

    def get_active_rule_sets(self):
        """Get the active rule sets in active_rule_df

        :return: a list of active rule sets
        """
        return self.active_rule_df.rule_set_name.unique()

    def get_rule_file_location(self):
        """Get the path to the rules file as a Path

        :return: Path for rules file
        """
        return self.rule_path

    def process_mol(self, mol):
        """Match a molecule against the active rule set

        :param mol: input RDKit molecule
        :return: the first rule matched or "ok" if no rules are matched
        """
        cols = ["description", "rule_set_name", "pat", "max"]
        for desc, rule_set_name, pat, max_val in self.active_rule_df[cols].values:
            if len(mol.GetSubstructMatches(pat)) > max_val:
                return rule_set_name, desc
        return "ok", "ok"

    def process_smiles(self, smiles):
        """Convert SMILES to an RDKit molecule and call process_mol

        :param smiles: input SMILES
        :return: process_mol result or None if the SMILES can't be parsed
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error parsing SMILES {smiles}")
            return None
        return self.process_mol(mol)
