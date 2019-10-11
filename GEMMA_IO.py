"""
Tools for reading and writing inputs and outputs of GEMMA

Authors:
Anna G. Green
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

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
