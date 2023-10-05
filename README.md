# EPIC-unmix

Deconvolution of bulk expression profiles with a two-step Bayesian method.

Note that our package is built upon [bMIND](https://github.com/randel/MIND).
All functions in the `bMIND` package are included with slightly modifications, which allows users to directly use the `bMIND` function.
We highly recommend users not loading both `bMIND` package and our `EPICunmix` package at the same time to avoid conflicts.

## Introduction of EPIC-unmix

EPIC-unmix is a two-step Bayesian method designed for cell-type-specific (CTS) analysis built upon the bMIND method. 
Similarly, it can provide sample-level CTS expression from the bulk RNA-seq data, 
with the aid of a reference panel derived from single-cell/single-nuclei RNA-seq data.
The second layer of Bayesian method, which distinguishes EPIC-unmix from bMIND, 
could adjust for the heterogeneity between reference and target samples, leading to more enhanced and robust performance.

## Instructions about running EPIC-unmix

### Installation

1. First, install the `devtools` package if you haven't already. 

		install.packages("devtools")

2. Load the devtools package and install through Github.

		library(devtools)
		install_github("quansun98/EPICunmix")

### Checklist before running the package

Data needed:

* bulk RNA-seq data of the study samples;
* single-cell or single-nuclei RNA-seq data as the reference;
* estimated cell type proportions of bulk RNA-seq data, 
which could be obtained using [MuSiC](https://xuranw.github.io/MuSiC/articles/MuSiC.html),
or estimated through the [`est_frac()`](https://randel.github.io/MIND/reference/est_frac.html) 
function which is also included in this package.

### QC of bulk data

We recommend performing QC of bulk data before inferring cell-type-specific (CTS) profiles.
For example, only keep the genes that are expressed in more than 80% of the samples.

		library(EPICunmix)
		bulk = read.table("example/by_id.csv",sep=',',header=T, row.name=1,check.names=F)
		n = dim(bulk)[2]
		n_express = apply(bulk, 1, function(x) sum(x > 0))
		bulk = bulk[(n_express > n*0.8),]

### Read in cell type fraction

In this example, we provided cell type fraction, a sample by cell type matrix.

		frac <- read.table("example/fracs.csv", sep = ',', header=T, check.names=F) 

### Read in sc/snRNA-seq reference

		sc_ref <- read.table("example/mat_train.csv",sep = ',',header=T,check.names=F)
		sc_meta <- read.table("example/meta_train.csv",sep = ',',header=T,check.names=F)
		
		# match cells
		index <- intersect(colnames(sc_ref),sc_meta$sample_name)
		sc_ref <- sc_ref[,index]
		sc_meta <- sc_meta[which(sc_meta$sample_name %in% index),]
		
		# match genes
		genes <- rownames(bulk)
		sc_ref <- sc_ref[genes,]
		
		# make sure the genes between bulk and reference are exactly the same
		bulk = bulk[order(rownames(bulk)),]
		sc_ref = sc_ref[order(rownames(sc_ref)),]
		colnames(sc_meta) <- c('sample_name','sample','cell_type')

### First-step Bayesian inference 

This step is the same as running the `bMIND` function.

		# get prior from reference
		prior = get_prior(sc_ref, meta_sc = sc_meta)

		# match bulk from prior because some genes would be dropped in this process
		bulk_sub <- bulk[rownames(bulk) %in% rownames(prior$profile),]
		bulk_sub <- bulk_sub[,order(colnames(bulk_sub))]
		frac <- frac[,order(colnames(frac))]
		celltype <- colnames(frac)
		colnames(frac) = paste0('c', 1:ncol(frac))
		dimnames(prior$cov)[[2]] = colnames(frac)
		dimnames(prior$cov)[[3]] = colnames(frac)
		colnames(prior$profile) = colnames(frac)

		# bMIND inference
		posterior = bMIND(bulk_sub, frac = frac, profile = prior$profile, covariance = prior$cov)
		## Note that if you have covariates that needed to be corrected, use the function below
		## posterior = bmind_de(bulk_sub, frac = frac, profile = prior$profile, covariance = prior$cov, covariate = covariate, noRE = F)
		
Note that the original `bmind_de` function in the `bMIND` package do not output the CTS profiles. 
We modified it to generate output that could be used in our second-step Bayesian inference.

### Second-step Bayesian inference

The second-step Bayesian inference takes the output from `bMIND()` or `bmind_de()`.

		epic_unmix = run_epic_unmix(bulk, frac, posterior, out_prefix = "example/epic_unmix")



