# Calculating heteroplasmic burden scores with mitoHPC results
This repository includes R codes to extract and compute number of heteroplasmic variants and mitochondrial local constraint score sum (MSS) with mitoHPC results.


**Computing heteroplasmic burden scores contains following steps:**
1. Extract SNVs and create three matrixs for allele, frequency, and read depth <-- "1_extract_vcf.R"
2. Calculate the maximum AAF for mtDNA loci <-- "2_aaf.R"
3. Sum all heteroplasmic variants and apply mitochondrial local constraint score (MLC score) <-- "3_het_score.R"
