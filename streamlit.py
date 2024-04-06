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
        vcf_file (str): The path to the VCF file.

    Returns:
        A Pandas DataFrame containing the variant information.
    """
    # Parse the VCF file using Biopython
    handle = open(vcf_file)
    records = list(seqio.parse(handle, "vcf"))
    handle.close()

    # Extract the relevant information and create a Pandas DataFrame
    variant_info = []
    for record in records:
        chrom = record.id.split(":")[0]
        pos = record.id.split(":")[1]
        ref = record.seq.ungap("").upper()
        alt = str(record.annotations["ALT"]).upper()
        qual = record.annotations["QUAL"]
        filter_ = record.annotations["FILTER"]
        info = record.annotations["INFO"]

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
    # Initialize an empty DataFrame to store the clinical significance information
    clinical_df = pd.DataFrame(columns=["Gene", "Chromosome", "Start", "End", "ClinVar", "OMIM"])

    # Loop through each variant and retrieve the clinical significance information
    for index, row in variant_df.iterrows():
        # Extract the variant information
        chrom = row["Chromosome"]
        pos = row["Position"]
        ref = row["Reference"]
        alt = row["Alternate"]

        # Construct the ClinVar API URL
        clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/api/v2/variation/?chrom={chrom}&start={pos-1}&end={pos}&reference_allele={ref}&alternate_allele={alt}"

        # Retrieve the ClinVar data
        try:
            clinvar_response = requests.get(clinvar_url)
            clinvar_data = clinvar_response.json()

            # Extract the relevant information
            clinvar_gene = clinvar_data["results"][0]["gene_symbol"]
            clinvar_chrom = clinvar_data["results"][0]["chrom"]
            clinvar_start = clinvar_data["results"][0]["start"]
            clinvar_end = clinvar_data["results"][0]["end"]
            clinvar_clinical_significance = clinvar_data["results"][0]["clinical_significance"]

            # Construct the OMIM API URL
            omim_url = f"https://api.omim.org/api/entry/?search={clinvar_gene}"

            # Retrieve the OMIM data
            omim_response = requests.get(omim_url)
            omim_data = omim_response.json()

            # Extract the relevant information
            omim_phenotype = omim_data["entry_list"][0]["phenotype_description"]

            # Add the clinical significance information to the DataFrame
            clinical_df = clinical_df.append({"Gene": clinvar_gene, "Chromosome": clinvar_chrom, "Start": clinvar_start, "End": clinvar_end, "ClinVar": clinvar_clinical_significance, "OMIM": omim_phenotype}, ignore_index=True)

        except:
            continue

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
