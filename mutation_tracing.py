"""
Classes and functions for determining pairs of homoplasies on a phylogenetic tree

Authors:
  Anna G. Green
"""

import pandas as pd
import numpy as np
import sys
from genotype_matrix import *
from collections import defaultdict
import glob
import re
from copy import deepcopy
import random
from evcouplings.align import Alignment

def name_branches(t1, lineage):
    """
    Assigns uniquen ames to each branch in forumlation: {lineage}-{anc_node}-{des_node}

    Parameters
    ----------
    t1: ete3.Tree
        the phylogeny of strains
    lineage: int or str
        name to prepend to branch names

    Returns
    -------
    list of str

    """
    branches = []

    for node in t1.traverse("levelorder"):
        if len(node.children) > 1:

            for idx, child_node in enumerate(node.children):
                name_anc = node.name.split("_")[0]
                name_des = child_node.name.split("_")[0]

                branch_name = f"{lineage}-{name_anc}-{name_des}"
                branches.append(branch_name)

    return branches


def get_ancestral_sequence(ancestral_seq_file, all_positions):
    """
    Fetches the ancestral sequence of H37Rv for all positions in input array
    all_positions MUST be in H37Rv numbering (same as ancestral seq file)

    Parameters
    ----------
    ancestral_seq_file: str
        path to file to input
    all_positions: list of int
        list of positions to grab from input

    Returns
    np.array of str
    """
    # subtract 1 or get wrecked. H37Rv positions are 1-indexed
    all_positions = [int(x) - 1 for x in all_positions]

    ancestral_seq = Alignment.from_file(open(ancestral_seq_file, "r"))

    ancestral_seq_subset = np.array([ancestral_seq.matrix[0, x] for x in all_positions])

    return ancestral_seq_subset


class HomoplasyAnalysis:

    def __init__(self, mutation_event_table, phylo, lineage, ancestral_seq_file, phylo_aln):

        self.mutation_table = mutation_event_table
        self.positions = sorted(list(self.mutation_table.Position.unique()))
        self.pos_to_index = {pos: idx for idx, pos in enumerate(self.positions)}

        self.tree = phylo
        self.lineage_label = lineage

        self.branches = name_branches(self.tree, self.lineage_label)
        self.idx_to_branch = {idx: x for idx, x in enumerate(self.branches)}

        self.aln = phylo_aln

        self.ancestral_seq = get_ancestral_sequence(ancestral_seq_file, self.positions)

        self.mutation_event_matrix = None
        self.reversion_event_matrix = None
        self.background_matrix = None

        indices_to_keep = []
        for idx, row in self.mutation_table.iterrows():
            if row.Ancestor_Call == self.aln.matrix[self.aln.id_to_index[row.Ancestor_Node],
                                                    self.pos_to_index[row.Position]]:
                indices_to_keep.append(idx)
        print(f"had to eliminate {len(self.mutation_table) - len(indices_to_keep)} rows due to outside error")
        self.mutation_table = self.mutation_table.loc[indices_to_keep, :]

    def filter_positions(self, new_positions):
        """
        Parameters
        ----------
        new_positions: list of int
        """

        if not (set(new_positions).issubset(set(self.positions))):
            raise ValueError("positions provided are not a subset of positions in the model")

        new_positions = sorted(new_positions)

        new_pos_to_idx = {pos: idx for idx, pos in enumerate(new_positions)}
        indices_to_keep = [self.pos_to_index[x] for x in new_positions]

        self.positions = new_positions
        self.pos_to_index = new_pos_to_idx

        self.ancestral_seq = self.ancestral_seq[indices_to_keep]

        self.mutation_table = self.mutation_table.query("Position in @new_positions")

        print(self.aln.matrix.shape)
        self.aln = self.aln.select(columns=indices_to_keep)
        print(self.aln.matrix.shape)

        if self.mutation_event_matrix is not None:
            self.mutation_event_matrix = self.mutation_event_matrix[:, indices_to_keep]

        if self.reversion_event_matrix is not None:
            self.reversion_event_matrix = self.reversion_event_matrix[:, indices_to_keep]

        if self.background_matrix is not None:
            self.background_matrix = self.background_matrix[:, indices_to_keep]

    def create_mutation_matrix(self, mode="mutation"):
        """
        Creates a numpy matrix where each branch gets a row, and each position gets a column.
        Contains a 1 if a mutation (or reversion, in reversion mode) happened on that branch

        """

        # initialize mutation event matrix
        mutation_event_matrix = np.zeros((len(self.branches), len(self.positions)))

        for _, row in self.mutation_table.iterrows():

            # look up row in mutation event matrix based on branch name
            branch_name = f"{self.lineage_label}-{row.Ancestor_Node}-{row.Derived_Node}"
            branch_index = self.branches.index(branch_name)

            # look up column in mutation eevn matrix based on position
            position_index = self.pos_to_index[row.Position]

            # function can either find all mutations AWAY FROM reference
            # or all mutations TO reference (reversions)
            if mode == "mutation":
                if row.Ancestor_Call != self.ancestral_seq[position_index]:
                    continue
            elif mode == "reversion":
                if row.Derived_Call != self.ancestral_seq[position_index]:
                    continue
            else:
                print(f"select either mutation or reversion for mode, you selected {mode}")
                return None

            mutation_event_matrix[branch_index, position_index] = 1

        return mutation_event_matrix

    def create_background_matrix(self):

        """
        Creates a numpy matrix where each branch gets a row, and each position gets a column.
        Contains a 0 if the ancestor node of that branch had the ancestral state, and a 1 otherwise

        """

        # create a dictionary of ancestral node to branches descending from that node
        node_to_branch_names = defaultdict(list)
        for x in self.branches:
            node_to_branch_names[x.split("-")[1]].append(x)

        # initialize mutation event matrix
        mutation_background_matrix = np.zeros((len(self.branches), len(self.positions)))

        # For each node in the tree
        for node in self.tree.traverse("levelorder"):
            # ignore if the node is a leaf
            if len(node.children) < 1:
                continue

            # get the name of the ancestral node and its sequence
            name_anc = node.name.split("_")[0]
            seq_anc = self.aln.matrix[self.aln.id_to_index[name_anc], :]

            # find positions where different from LUCA sequence
            differences = (self.ancestral_seq != seq_anc).astype(int)

            # for every branch descending from this node
            for b in node_to_branch_names[name_anc]:
                # save positions where there's a difference as 1, else 0
                mutation_background_matrix[self.branches.index(b)] = differences

        return mutation_background_matrix

    def dump_equals(self, filename):
        """
        Write a file of all position pairs that both mutated on the same branch

        Parameters
        ----------
        filename: str
            path to file to write

        """

        def _dump_equality(filehandle, positions, branch_index):

            position_names = [self.positions[x] for x in positions]

            b = self.idx_to_branch[branch_index]

            for idx, p1 in enumerate(position_names[:-1]):
                for p2 in position_names[idx + 1::]:
                    string = f"{p1},{p2},{b},EQUAL\n{p2},{p1},{b},EQUAL\n"
                    of.write(string)

        with open(filename, "w") as of:
            # iterate through rows
            for row_i in range(self.mutation_event_matrix.shape[0]):
                # find each column that has a 1
                hits = np.where(self.mutation_event_matrix[row_i, :] == 1)[0]
                # write the informaiton on those positions
                _dump_equality(of, hits, row_i)

    def dump_subset(self, filename):
        """
        Write a file of all position pairs where position j mutated after position i

        Parameters
        ----------
        m: numpy.ndarray
            2-d array with each row representing a branch and each column representing a position
        x: numpy.ndarray
            2-d background matrix (0 if branch started with ancestor, 1 else)
        all_positions: list of int
            positions to include (H37Rv numbering)
        idx_to_branch: dict of int:str
            indices of first dimension of m and their corresponding name
        filename: str
            path to file to write

        """

        def _dump_to_file(filehandle, positions_with_mutation,
                          positions_with_background, branch, CODE
                          ):

            positions_1 = [self.positions[x] for x in positions_with_mutation]
            positions_2 = [self.positions[x] for x in positions_with_background]
            branch_name = self.idx_to_branch[branch]

            for p1 in positions_1:
                for p2 in positions_2:
                    string = f"{p1},{p2},{branch_name},{CODE}\n"
                    of.write(string)

        with open(filename, "w") as of:
            # test for equal mutation events
            for row_i in range(self.mutation_event_matrix.shape[0]):
                # Get all positions that mutated on that branch
                positions_with_mutation = np.where(self.mutation_event_matrix[row_i, :] == 1)[0]

                # Get all positions that had already been mutated when we got to that branch
                positions_with_background = np.where(self.background_matrix[row_i, :] == 1)[0]

                _dump_to_file(
                    of, positions_with_mutation,
                    positions_with_background, row_i, "SUBSET"
                )


def run(t1, aln, mutation_table, lineage,
        prefix, n_events_required=2,
        ancestral_file="MTB_ancestor_reference.fasta", mapping_cutoff=0.9,
        agreement_cutoff=0.9, n_shuffling_permutations=100,
        quality_table_path='../co-homoplasy/data_tables/empirical_agreement_df.csv'
        ):
    """
    Runs analysis of homoplasies

    Parameters
    ----------
    t1: ete3.Tree
    aln: evcouplings.align.Alignment
    mutation_table: pd.DataFrame
    lineage: int
    prefix: str
    """

    ha = HomoplasyAnalysis(mutation_table, t1, lineage, ancestral_file, aln)

    # Filter out unwanted positions
    empirical_agreement_df = pd.read_csv(quality_table_path, index_col=0)
    above_filters = empirical_agreement_df.query("agreement > @agreement_cutoff and mappability > @mapping_cutoff")

    positions_to_keep = list(above_filters.pos)
    positions_to_keep = list(set(positions_to_keep).intersection(set(ha.positions)))
    print(f"keeping {len(positions_to_keep)} positions based on mapping filter, of {len(ha.positions)}")
    ha.filter_positions(positions_to_keep)

    # create the three matrices
    ha.mutation_event_matrix = ha.create_mutation_matrix(mode="mutation")
    ha.reversion_matrix = ha.create_mutation_matrix(mode="reversion")
    ha.background_matrix = ha.create_background_matrix()

    # Create a df with information on branches
    n_muts = ha.mutation_event_matrix.sum(axis=1)
    n_revs = ha.reversion_matrix.sum(axis=1)
    n_total = n_muts + n_revs
    branch_table = pd.DataFrame({
        "branch_name": ha.branches,
        "n_muts_from_numpy": n_muts,
        "n_reversions_from_numpy": n_revs,
        "total_mutations_from_numpy": n_total
    })
    branch_table.to_csv(f"{prefix}_branch_info_table.csv")

    # Create a df with informaiton on each position
    n_muts = ha.mutation_event_matrix.sum(axis=0)
    n_revs = ha.mutation_event_matrix.sum(axis=0)
    mutation_info = pd.DataFrame({
        "i": [idx for idx, i in enumerate(ha.positions)],
        "genomic_position": ha.positions,
        "n_mutations_from_anc": n_muts,
        "n_reversion_to_anc": n_revs,
        "total_n_mutations": n_muts + n_revs
    })
    mutation_info.to_csv(f"{prefix}_mutation_count_table.csv")

    # Filter for positions with more than n_events_required hits
    kept_position_indices = n_muts >= n_events_required
    kept_positions = np.array(ha.positions)[kept_position_indices]
    ha.filter_positions(kept_positions)

    print(f"{len(kept_position_indices)} pre-filtering, {len(ha.positions)} post filtering")

    for i in range(0, n_shuffling_permutations):
        print("shuffling iteration", i)

        new_m = deepcopy(ha.mutation_event_matrix)
        for col in range(new_m.shape[1]):
            new_m[:, col] = random.sample(list(ha.mutation_event_matrix[:, col]), len(ha.mutation_event_matrix[:, col]))

        ha.mutation_event_matrix = new_m
        # Save the mutations that happen on the same branch
        ha.dump_equals(f"{prefix}_shuffle{i}_table_equals.csv")

        # Save the mutations where one subsets another
        ha.dump_subset(f"{prefix}_shuffle{i}_table_subset.csv")
