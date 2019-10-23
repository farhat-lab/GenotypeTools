"""
Tools for reading and writing inputs and outputs of GEMMA

Authors:
Anna G. Green
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def prepare_phenotpye_file(prefix, phenotype_file, strain_file, outcfg):
    """
    Prepare input phenotype file for GEMMA

    Parameters
    ----------
    prefix: str
        prefix to prepend to output file
    phenotype_file: str
        file containing two columns: the first called 'strain'
        with strains to use, and the second containing corresponding phenotypes
    strain_file: str
        file containing names of strains to use.
    outcfg: dict of str:str

    Returns
    -------
    outcfg

    """

    # write the input phenotype file
    strain_table = pd.read_csv(strain_file)
    strain_list = list(strain_table.strain)

    phenotypes = pd.read_csv(phenotype_file, index_col='strain', header=0)

    phenotypes = phenotypes.loc[strain_list,:]

    gemma_input_phenotype_file = f"{prefix}.phenotypes"

    phenotypes.iloc[:,-1].to_csv(gemma_input_phenotype_file, header=False, index=None, na_rep="NA")
    outcfg["gemma_input_phenotype_file"] = gemma_input_phenotype_file
    return outcfg

def prepare_loci_file(prefix, genotype_file,
    position_file, outcfg_prefix, gap_code, outcfg):
    """
    Prepare input loci file for GEMMA
    """

    verify_resources(
        "Input file does not exist",
        genotype_file, position_file
    )

    # Read in list of positions to consider
    position_data = pd.read_csv(position_file)
    positions = sorted(list(position_data.POS.unique()))
    outcfg["num_positions_" + outcfg_prefix] = len(positions)

    # read in genotype data
    genotype_data = pd.read_csv(genotype_file)

    strains = genotype_data.STRAIN.unique()
    print("analyzing N strains:", len(strains))

    #initialize genotype GenotypeMatrix
    gm = GenotypeMatrix(strains, positions, genotype_data)

    #make real mutations
    print("filling in the genoype matrix")
    gm.make_matrix()

    print('preparing gemma input files')
    # write the input loci file
    gemma_input_loci_file = f"{prefix}_{outcfg_prefix}.loci"
    of = open(gemma_input_loci_file, "w")
    gm.write_GEMMA(of, gap_code=gap_code)
    outcfg["gemma_loci_file_" + outcfg_prefix] = gemma_input_loci_file

    # write the major allele file
    major_allele_file = f"{prefix}_{outcfg_prefix}_major_alleles.csv"
    maj_alleles = pd.DataFrame(gm.major_alleles, index=gm.positions, columns=["major_alleles"])
    maj_alleles.loc[:,"major_alleles"] = [ALPHABET[x] for x in maj_alleles.major_alleles]
    maj_alleles.to_csv(major_allele_file)
    outcfg["major_allele_file_" + outcfg_prefix] = major_allele_file

    # write the strain file
    strain_file= f"{prefix}_strains.csv"
    strains = pd.DataFrame(gm.strains, columns=["strain"])
    strains.to_csv(strain_file, index=None)
    outcfg["strain_file"]=strain_file

    return outcfg

def execute(incfg,  **kwargs):
    """
    Command line wrapper for executing GEMMA
    """

    outcfg = incfg

    prefix = str(kwargs["prefix"])
    gemma_executable = kwargs["gemma"]
    relatedness_loci_file = incfg["gemma_loci_file_relatedness"]
    trait_loci_file = incfg["gemma_loci_file_trait"]
    pheno_file = incfg["gemma_input_phenotype_file"]
    outfile = prefix + "_gemma_out"
    maf_GRM = incfg["maf_GRM"]
    maf_LMM = incfg["maf_LMM"]

    gemma_output = outfile

    #calculate relatedness matrix
    string = f"{gemma_executable} -g {relatedness_loci_file} -p {pheno_file} -gk -maf {maf_GRM} -o {outfile}"
    print(string)
    os.system(string)

    # run the actual gwas
    string = f"{gemma_executable} -g {trait_loci_file} -p {pheno_file} -k output/{outfile}.cXX.txt -maf {maf_LMM} -lmm -notsnp -o {outfile}"
    print(string)
    os.system(string)
    outcfg["gemma_output"] = gemma_output
    return outcfg

def read_GEMMA_assoc_table(assoc_file, snp_annotation_file):
    """
    assoc_file: str
        path to input file
    snp_annotation_file: str
        path to file of annotations for SNPs
    """

    d = pd.read_csv(assoc_file, sep="\t", dtype=object)
    snp_annotation = pd.read_csv(snp_annotation_file, index_col=0, dtype=object)

    data = d.merge(snp_annotation, left_on="rs", right_on="POS")
    data = data[[
        "POS","AA_pos", "name",  'SNP_type', "type","beta", 'p_wald',
        "accession", 'description',  'percent_hi_quality',
        'promoter_coding_gene', 'reference_allele', "se", "logl_H1", "l_remle"
    ]]
    data["p_wald"] = data["p_wald"].astype(float)
    return data.sort_values("p_wald")

def qq_plot(df, p_val_column, outfile="QQ.pdf", ylim=None):

    """
    makes a QQ plot of the data
    """

    df = df.sort_values(p_val_column, ascending=False)

    length = len(df)

    increment = 1./length

    expected_p_values = np.arange(1,0,-increment)

    fig = plt.figure(figsize = (5,5))
    ax = fig.gca()

    p_val_plot_expected = -np.log10(expected_p_values)
    p_val_plot_real = -np.log10(df[p_val_column])

    ax.scatter(p_val_plot_expected, p_val_plot_real, s=8, color="black")

    ax.plot(p_val_plot_expected, p_val_plot_expected, color="grey")
    sns.despine()
    ax.set_xlabel("-log10 expected P values")
    ax.set_ylabel("-log10 real P values")
    if ylim:
        ax.set_ylim(ylim)
    plt.savefig(outfile)

    return expected_p_values, df
