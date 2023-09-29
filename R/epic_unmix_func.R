#' EPIC-unmix for sample-level cell-type-specific gene expression analysis
#'
#' A two step Bayesian method. In the first step, a normal Bayesian framework was adopted, similar to the bMIND algorithm. In the second step, we update the prior using results from the first step.
#'
#' @param bulk bulk gene expression (gene x sample). We recommend log2-transformed data for better performance. If the max(bulk) > 50, bulk will be transformed to log2(count per million + 1).
#' @param frac sample-specific cell type fraction (sample x cell type). This could be estimated by external software, e.g., MuSiC, bMIND, etc.
#' @param input_cts an R object, standard output from `bMIND` or `bMIND_de` in this package. It will be used to infer the second layer of prior.
#' @param outf logical, write output (CTS profile) to files or not, one file for each cell type. Default TRUE
#' @param out_prefix file output name prefix. Default value "epic_unmix"
#' @param out_ct selected output cell types. Default all the cell types
#' @param nstop number of iteration before early stop. In practice, we do not recommend iterating too many times, since it will significantly reduce the number of genes 
#' (it will drop genes during each iteration). In our evaluations, 1 time iteration would already give satisfying performance.
#' @param delta stoping criteria. The maximum mean difference between two iterations. Default 0.1. We do not recommend setting this to be very stringent, 
#'which will significantly increase the number of iterations needed and will drop genes.
#' @param nu0 hyper-parameter for the prior covariance matrix initially derived from the reference. The larger the nu0, the higher the certainty about the information in covariance, and the more 
#' informative is the distribution. The default is 50.
#' @param nu1 hyper-parameter for the prior covariance matrix derived from the last iteration. The larger the nu1, the more informative after each iteration. Default is 50.
#' @param seed random seed used. Default 1.
#' @param ncore number of cores to run in parallel. The default is all available cores.
#'
#' @return Same as bMIND output, a list containing of the following elements:
#' \item{A}{the deconvolved cell-type-specific gene expression (gene x cell type x sample). Note that if outf == TRUE, it will be written to files.}
#' \item{SE}{the standard error of cell-type-specific gene expression (gene x cell type x sample).}
#' \item{Sigma_c}{the covariance matrix for the deconvolved cell-type-specific expression (gene x cell type x cell type).}
#' \item{mu}{the estimated profile matrix (gene x cell type).}
#' \item{frac}{the estimated cell type fractions (sample x cell type) if fractions are not provided.}
#'
#' @examples
#'
#' file_dir <- 'example/'
#' bulk <- read.table(paste0(file_dir,'by_id.csv'),sep=',',header=T, row.name=1,check.names=F)
#' frac <- read.table(paste0(file_dir,'fracs.csv'),sep = ',',header=T,check.names=F)
#' cts = readRDS("bmind.rds")
#' #note that the example already matched genes, otherwise genes must be exactly matched
#' epic_unmix = run_epic_unmix(bulk, frac, cts, out_prefix = paste0(file_dir,'epic_unmix'))
#'
#' @export run_epic_unmix
#'
#' @importFrom matrixcalc is.positive.definite
#' 

run_epic_unmix = function(bulk, frac, input_cts, outf = TRUE, out_prefix = NULL, out_ct = NULL, nstop = 1, 
			  delta = 0.1, nu0 = 50, nu1 = 50, seed = 1, ncore = NULL){

	current_cts = input_cts
	for(iter in 1:nstop){

		n_gene = dim(current_cts$A)[1]
		n_ct = dim(current_cts$A)[2]
		n_sample = dim(current_cts$A)[3]

		bp = apply(current_cts$A, c(1,2), mean) 
		bp_cov = array(NA, dim = c(n_gene, n_ct, n_ct))
		rownames(bp_cov) = rownames(current_cts$A)
		colnames(bp_cov) = colnames(current_cts$A)
		for(i in 1:n_gene) {
    			tmp.cov = cov(t(current_cts$A[i,,] + 0.1), use = 'pairwise')
    			if(matrixcalc::is.positive.definite(tmp.cov)){
       				bp_cov[i,,] = tmp.cov
    			}
		}

		na_index = is.na(bp_cov[,1,1])
		bp_cov = bp_cov[!na_index,,]
		bp = bp[!na_index,]
		new_bulk_sub = bulk[rownames(bulk) %in% rownames(bp_cov),]

		if((dim(new_bulk_sub)[1] == dim(bp)[1]) && (dim(bp)[1] == dim(bp_cov)[1])){
			out_cts = bmind_de(new_bulk_sub, frac = frac, profile = bp, covariance = bp_cov, noRE = F, nu = nu0 + nu1*iter, seed = seed, ncore = ncore)
		}else{
			stop("bulk data dimension and priors are not the same, please check")
		}

		origA = current_cts$A[rownames(current_cts$A) %in% rownames(out_cts$A),,]

		e=c()
		for(j in 1:n_ct){
			e = c(e, sum(rowMeans((origA[,j,]-out_cts$A[,j,])^2)/nrow(out_cts$A)))
		}

		message(paste(iter,"iteration: maximum CTS difference between the last two iterations is", max(e)))

		if(max(e) < delta){
			break
		}

		current_cts = out_cts
		bulk = new_bulk_sub[rownames(new_bulk_sub) %in% rownames(current_cts$A),]

	}

	message(paste("Analysis finished!"))

	if(outf == T){
		out_prefix = ifelse(is.null(out_prefix), "epic_unmix", out_prefix)
		if(is.null(out_ct)) {out_ct = colnames(out_cts$A)}
		index = match(out_ct, colnames(out_cts$A))
		for(j in index){
			write.table(out_cts$A[,j,], paste0(out_prefix,"_",colnames(out_cts$A)[j],"_CTS.txt"), col.names = T, row.names = T, quote = F, sep = "\t")
		}
	}

	return(out_cts)

}


