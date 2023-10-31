# IDUL
IDUL (iterative dispersion update to fit linear mixed model) is designed for multi-omics analysis where each SNPs are tested for association with many phenotypes. IDUL has both theoretical and practical advantages over the Newton-Raphson method. 

## Current version 
Version 0.51 was compiled on 7 Aug 2023. Linux exectuable can be downloaded from: http://www.haplotype.org.  
You may also choose to compile from the source code in src/, but you want to install eigen and boost packages first. A preprint documenting IDUL can be found here: https://biorxiv.org/cgi/content/short/2023.10.25.563975v1. 

## Usage exmpales  
1) This example takes vcf file as genotypes input.  The output contains likelihood ratio test p-values, effect size, likelihood, etc.
    
       ./idul -i input.vcf.gz -p phenotypes.gz -k kinship.txt -o pref 

3) This example takes bimbam mean genotype as genotype input. Since -b is invoked, pref.pval.gz only contains a p-value for each SNP each phenotype.
   
       ./idul -g input.bimbam.mgt.gz -p phenotypes.gz -k kinship.gz -o pref -b 

## Input and options  
    y = W a + x b + u + e 
This is the standard linear mixed model. The following options specify each components (except e). 

    -p phenotype_file (y). 
    -c covariate_file (W). 
    -i vcf_file  (x).  
    -g bimbam_mean_genotypes (x). 
    -k kinship_file (u).  
  
Phenotypes file is expected to have multiple columns (a single column is okay too) with one column for each phenotype. Phenotypes should have no missing values. Covariates file is a square, with the same number of rows as phenotypes, and however many columns needed. Missing values are allowed for vcf files (replaced with the mean genotypes), but not allowed for bimbam mean genotypes. Kinship is a square and symmetric kinsihp matrix, such as one estimated by [kindred](https://github.com/haplotype/kindred). Rows of bimbam mean genotype looks like below. 

    rs123 A G 0 1 2 2 0 1.1 0 
    rs456 C T 1 0 1 0.05 0 2 0 
  
### The options related to input files: 

    -e eigenvector_file.   
    -v eigenvalues_file.  
    -f minor_allele_freq. 
  
In eigenvector file, each column is an eigenvector. Eigenvalues file contains one row with n numbers, its jth number should correspond to j-th column of the eigenvectors. Two options -e and -v come together, and mutually exclusive with -k. 
SNPs whose minor allele frequencies are below the threshold will be removed. 

## Output and options

     -o output_pref. 
     -b 
    
There are three output files: pref.log (a txt document contains log), pref.snpinfo.txt.gz, and pref.pval.gz file, where pref is specified by -o. 
There are four columns in pref.snpinfo.txt.gz, SNP ID, A-allele, B-allele, and MAF (minor allele frequency). 
By default pref.pval.gz contains p-values, beta, sigma, l0 (loglikelihood of null), l1 (loglikelihood of alternative), h (eta/(eta+1)), and niter (number of iterations used) in optimization. 
If -b was invoked, pref.pval.gz only contains pvalues.  


## Other options
       -t  number_of_threads.  
       -w  compute wald test p-values instead of likelihood ratio test p-values, which is default. 

