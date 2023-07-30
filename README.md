# IDUL
IDUL (iterative dispersion update to fit linear mixed model) is designed for multi-omics analysis where each SNPs are tested for association with many phenotypes. IDUL has both theoretical and practical advantages over the Newton-Raphson method. Below is the abstract of a manuscript that document IDUL with a title "Iterative dispersion update achieves an exact optimum for the linear mixed model in genetic association studies".   

In genetic association studies, the linear mixed model (LMM) has become a standard practice to account for population stratification and relatedness in the sample, which if unaccounted for tend to produce false positives. Despite recent progress, fitting LMMs to achieve an exact optimum remains computationally demanding, particularly in multi-omics studies where tens of thousands of phenotypes are tested for association with millions of genetic markers.  
The state-of-the-art method GEMMA used Newton-Raphson algorithm for optimization to achieve ``effectively exact" optimum.   
Here we present IDUL, an approach that uses iterative dispersion updates to fit LMMs. The dispersion update requires no evaluation of derivatives of the likelihood function, rather, it fits two weighted least square regressions in each iteration.  Applied to Sabre protein assay from the Framingham Heart Study,  IDUL converged in a few iterations and provided consistent estimates between different initial values, showed superiority over the Newton-Raphson method. We analyzed IDUL through a theoretical lens to show that it is as efficient as the Newton-Raphson method near the optimum, and with a sufficiently large sample size, IDUL update  always increases the likelihood, even outside the neighborhood of the optimum, which ensures its consistency. Most significantly, with sufficiently large sample, IDUL converges to the global optimum with a probability of one, achieving the exact optimum.  Software implementing IDUL can be downloaded at http://www.haplotype.org. 

## Current version 
Version 0.51 was compiled on 28 July 2023. Exectuables for Linux and Mac can be found here: http://www.haplotype.org/software.html.  
You may also choose to compile from the source code in src/. 

## Input and options  
    y = W a + x b + Z u + e 
  -p phenotype_file (y). Phenotypes file is expected to have multiple columns (a single column is okay too) with one column for each phenotype. Phenotypes should have no missing values. 
  
  -c covariate_file (W). Covariates file is a square, with the same number of rows as phenotypes, and however many columns needed. 
  
  -i vcf_file  (x).  Missing values are allowed for vcf files (replaced with the mean genotypes). 
  
  -g bimbam_mean_genotypes (x). Missing value is not allowed. Rows of bimbam mean genotype looks like below. 

    rs123 A G 0 1 2 2 0 1.1 0 
    rs456 C T 1 0 1 0.05 0 2 0 
  
  -k kinship_file (u).   Kinship is a square and symmetric kinsihp matrix, such as one estimated by [kindred](https://github.com/haplotype/kindred)  


The options related to input files: 

  -e eigenvector_file.   Each column is an eigenvector.  
  
  -v eigenvalues_file.   One row with n numbers, assumed jth number of the eigenvalue correpsonding to j-th column of the eigenvectors. These two options come together, and mutually exclusive with -k. 

  -f minor_allele_freq.  SNPs whose maf below the threshold will be removed. 

## Output and options
There are three output files: pref.log (a txt document contains log), pref.snpinfo.txt.gz, and pref.pval.gz file, where pref is specified by -o. 
There are four columns in pref.snpinfo.txt.gz, SNP ID, A-allele, B-allele, and MAF (minor allele frequency). 
By default pref.pval.gz only contains p-values, one column per phenotypes and one SNP per row. 
If -b is invoked, then additional informations will be writen, include SNP effect size (beta), sigma of beta, eta/(1+eta), and number of IDUL iterations to obtain eta estimates. 

## Usage exmpales  
1) This example takes vcf file as genotypes input. The output pref.pval.gz only contains a p-value for each SNP each phenotype.
   
       ./idul -i input.vcf.gz -p phenotypes.gz -k kinship.txt -o pref 

2) This example takes bimbam mean genotype as genotype input. Since -b is invoked, the pref.pval.gz contains more columns in addition to p-values, such as beta and sigma, etc.
   
       ./idul -g input.bimbam.mgt.gz -p phenotypes.gz -k kinship.gz -o pref -b 
