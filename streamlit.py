!pip install streamlit biopython requests numpy pandas

import streamlit as st
import pandas as pd
import numpy as np
import biopython.seqio as seqio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import requests


def parse_vcf(vcf_file):
    """
    Parses a VCF file and extracts the relevant information.

    Args:
        vcf_file: A file-like object containing the VCF file.

    Returns:
        A Pandas DataFrame containing the variant information.
    """
    vcf_reader = pysam.VariantFile(vcf_file)
    variant_info = []
    
    for record in vcf_reader:
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt = ','.join(record.alts)
        qual = record.qual
        filter_ = ','.join(record.filter.keys()) if record.filter.keys() else "PASS"
        info = record.info

        variant_info.append([chrom, pos, ref, alt, qual, filter_, info])

    variant_df = pd.DataFrame(variant_info, columns=["Chromosome", "Position", "Reference", "Alternate", "Quality", "Filter", "Info"])

    return variant_df

def get_clinical_significance(variant_df):
    """
    Retrieves clinical significance information from public databases.

    Args:
        variant_df (Pandas DataFrame): A DataFrame containing the variant information.

    Returns:
        A Pandas DataFrame containing the clinical significance information.
    """
    clinical_info = []
    for _, row in variant_df.iterrows():
        chrom = row["Chromosome"]
        pos = row["Position"]
        ref = row["Reference"]
        alt = row["Alternate"]

        # Construct the Ensembl VEP API URL for variant annotation
        vep_url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}:{ref}/{alt}?"
        vep_response = requests.get(vep_url, headers={"Content-Type": "application/json"})

        if vep_response.status_code == 200:
            vep_data = vep_response.json()
            if vep_data:
                gene = vep_data[0].get("gene_symbol", "N/A")
                clin_sig = vep_data[0].get("most_severe_consequence", "N/A")
                clinical_info.append([gene, chrom, pos, clin_sig])
            else:
                clinical_info.append(["N/A", chrom, pos, "No data"])
        else:
            clinical_info.append(["N/A", chrom, pos, "Error fetching data"])

    clinical_df = pd.DataFrame(clinical_info, columns=["Gene", "Chromosome", "Position", "Clinical Significance"])
    return clinical_df

def main():
    # Set the page title
    st.title("Genomic Variant Analyzer")

    # Allow the user to upload a VCF file
    uploaded_file = st.file_uploader("Upload a VCF file", type="vcf")

    if uploaded_file:
        # Parse the VCF file
        variant_df = parse_vcf(uploaded_file)

        # Retrieve the clinical significance information
        clinical_df = get_clinical_significance(variant_df)

        # Display the variant table
        st.subheader("Variant Table")
        st.write(variant_df)

        # Display the clinical significance table
        st.subheader("Clinical Significance")
        st.write(clinical_df)

if __name__ == "__main__":
    main()

