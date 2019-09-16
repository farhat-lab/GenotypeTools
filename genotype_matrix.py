"""
Tools for creating a matrix of strains by genomes. 
The GenotypeMatrix class is inspired by and modified from EVCouplings/align Alignment class, 
which is (c) 2017 EVcouplings development team and freely available under the MIT license

Authors:
Anna G. Green
"""

import pandas as pd
import numpy as np

ALPHABET_NOGAP = "ACGT"
GAP_CHAR = "-"
ALPHABET = GAP_CHAR + ALPHABET_NOGAP

def map_matrix(matrix, map_):
    """
    Map elements in a numpy array using alphabet
    Parameters
    ----------
    matrix : np.array
        Matrix that should be remapped
    map_ : defaultdict
        Map that will be applied to matrix elements
    Returns
    -------
    np.array
        Remapped matrix
    """
    return np.vectorize(map_.__getitem__)(matrix)


def calculate_frequencies(matrix, num_symbols):
    """
    Calculate single-site frequencies of symbols in alignment
    Parameters
    ----------
    matrix : np.array
        N x L matrix containing N sequences of length L.
        Matrix must be mapped to range(0, num_symbols) using
        map_matrix function
    num_symbols : int
        Number of different symbols contained in alignment
    Returns
    -------
    np.array
        Matrix of size L x num_symbols containing relative
        column frequencies of all characters
    """
    N, L = matrix.shape
    fi = np.zeros((L, num_symbols))
    for s in range(N):
        for i in range(L):
            fi[i, matrix[s, i]] += 1
    return np.divide(fi, N)

class GenotypeMatrix:
    """
    Takes inspiration from the evcouplings python package Alignment class
    """
    def __init__(self, strain_list, position_list, genotype_data=None, matrix_data=None):

        self.positions = position_list
        self.position_to_index = {p:idx for idx, p in enumerate(self.positions)}
        self.strains = strain_list
        self.strain_to_index = {p:idx for idx, p in enumerate(self.strains)}

        # If user provided genotype data, do some sanity checks
        self.genotype_data = genotype_data
        # Make sure that every position and every genotype is represented in the genotype_data
        if not self.genotype_data is None:
            assert set(self.positions).issubset(genotype_data.POS)
            assert set(self.strains).issubset(genotype_data.STRAIN)

        # if user provided matrix data, set it up
        if not matrix_data is None:
            assert matrix_data.shape[0] == len(self.strains)
            assert matrix_data.shape[1] == len(self.positions)
            self.matrix = matrix_data
        else:
            self.matrix = np.zeros((len(self.strains), len(self.positions)), dtype=str)

        self.alphabet_map = {
            c: i for i, c in enumerate(ALPHABET)
        }
        self.num_symbols = len(self.alphabet_map)

        self.matrix_mapped = None
        self._frequencies = None
        self._major_alleles = None

    @classmethod
    def from_df(cls, file):
        """
        Initializes a GenotypeMatrix using an input pd.DataFrame, where rows = strains and columns=positions, and
        the data will become the contents of self.matrix
        """
        df = pd.read_csv(file, index_col=0, header=0)

        return cls(
            list(df.index), list(df.columns), genotype_data=None, matrix_data=df.values
        )

    def make_matrix(self):
        '''
        Populates the matrix of strains by positions with genotype data from
        the input genotype data file via the following heuristic:

        0. Initialize every position in every genome as the reference allele
        1. For each position in each strain, determine if there is any information
            in the genotype data input. If not, leave that position as Reference. If yes,
            proceed to step 2.
        2. If the alternate allele is an insertion or a non-standard character, insert a gap
            at that position.
        3. If the alternate allele did not pass quality threshold (determined by QUAL flag being 0),
            insert a gap at that position
        4. If alternate allele is a SNP and is not low quality, replace reference with alternate allale

        Returns
        -------
        np.ndarray:
            an N_strains x N_positions array of strings
        '''

        for position in self.positions:
            subset = self.genotype_data.query("POS==@position")
            allele = subset.loc[subset.index[0], "REF"]
            self.matrix[: , self.position_to_index[position]] = allele

        # iterate through each strain
        strain_group = self.genotype_data.query("STRAIN in @self.strains").groupby("STRAIN")


        for strain, subset in strain_group:
            print(strain)

            # iterate through each position
            subset = subset.query("POS in @self.positions")
            position_groups = subset.groupby("POS")
            for position, subsubset in position_groups:

                assert(len(subsubset) < 2)

                # if we didn't find any information for this position, leave as reference allele
                if len(subsubset) == 0:
                    continue

                subsubset = subsubset.iloc[0,:]

                # If we did find information but its an indel, treat as missing
                if len(subsubset.ALT) > 1 or  not subsubset.ALT in ALPHABET:
                    self.matrix[self.strain_to_index[strain],self.position_to_index[position]] = "-"

                # If we found information but it didn't pass quality filter, treat as missing
                elif subsubset.QUAL==0:
                    self.matrix[self.strain_to_index[strain],self.position_to_index[position]] = "-"

                # else, change to alterante allele
                else:
                    self.matrix[self.strain_to_index[strain],self.position_to_index[position]] = subsubset.ALT


    def __ensure_mapped_matrix(self):
        """
        Ensure self.matrix_mapped exists
        """
        if self.matrix_mapped is None:
            self.matrix_mapped = map_matrix(
                self.matrix, self.alphabet_map
            )

    @property
    def major_alleles(self):
        """
        Returns/calculates identity of major allele for each position in alignment.
        Also sets self._major_alleles member variable for later reuse.

        Returns
        -------
        np.array
            Reference to self._major_alleles
        """
        if self._major_alleles is None:

            self._major_alleles = np.argmax(
                self.frequencies, axis=1
            )

        return self._major_alleles

    @property
    def frequencies(self):
        """
        Returns/calculates single-site frequencies of symbols in alignment.
        Also sets self._frequencies member variable for later reuse.

        Returns
        -------
        np.array
            Reference to self._frequencies
        """
        if self._frequencies is None:
            self.__ensure_mapped_matrix()

            self._frequencies = calculate_frequencies(
                self.matrix_mapped, self.num_symbols
            )

        return self._frequencies

    def write_fasta(self, fileobj, strain_list=None):
        """
        Creates a fasta-formatted file of pseudo-sequences for strains and positions in self.matrix

        Parameters
        ----------
        fileobj: filehandle
            fileobject to which to write
        strain_list: list of str, optional (default=None)
            List of strain names, in order, which will be written to file. If None,
            will write all strains in self.strains

        """
        if strain_list is None:
            strain_list = self.strains

        for strain in strain_list:
            seq = "".join(self.matrix[self.strain_to_index[strain],:])
            fileobj.write(f"{strain}\n>{seq}\n")

    def write_GEMMA(self, fileobj, major_allele_string="\"A\"", minor_allele_string="TRUE",
        major_allele_code='1', minor_allele_code='0',gap_code="NA"):
        """
        Creates a file of loci in correct input format for GEMMA

        User warning: Changing the optional variables to values other than their defaults may result in errors when
        using the program GEMMA. Verify that all inputs are still valid before proceeding.

        Gemma input format:
        position_name,major_allele_string,minor_allele_string,{allele code for strain 1, ..., allele code for strain N}

        Parameters
        ----------
        fileobj: filehandle
            file to which to write
        major_allele_string: str, optional (default "A")
            string used to represent the major allele
        minor_allele_string: str, optional (default "T")
            string used to represent the minor allele
        major_allele_code: str, optional (default "1")
            string used to code for the major allele in the input
        minor_allele_code: str, optional (default "0")
            string used to code for the major allele in the input
        gap_code: str, optional (default "NaN")
            string used to code for the gap/missing data in the input

        """
        self.__ensure_mapped_matrix()

        for position in self.positions:
            # slice corresponding to all strains
            vec = self.matrix_mapped[:, self.position_to_index[position]]

            is_major_allele = np.equal(vec, self.major_alleles[self.position_to_index[position]])
            is_gap = np.equal(vec, 0)

            string = np.array([minor_allele_code] * len(vec), dtype=object)
            string[is_major_allele] = major_allele_code
            string[is_gap] = gap_code

            fileobj.write(f'{position},{major_allele_string},{minor_allele_string},{",".join(string)}\n')

    def to_df(self):
        """
        Creates a pd.DataFrame, with rows corresponding to strains, columns corresponding to positions,
        and entries corresponding to the alleles in self.matrix
        Returns
        -------
        pd.DataFrame
        """
        df = pd.DataFrame(self.matrix, index=self.strains, columns=self.positions)
        return df
