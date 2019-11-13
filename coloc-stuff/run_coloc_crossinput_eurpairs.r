#!/usr/bin/env Rscript

# ypar
# reference for main analysis: https://github.com/chr1swallace/coloc


args = commandArgs(trailingOnly=TRUE)
print(args)

library(dplyr)
library(data.table)
library(coloc)

Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}

logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

logdiff <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}

approx.bf.p <- function(p,f,type, N, s, suffix=NULL) {
  if(type=="quant") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ret <- data.frame(V,z,r,lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)
}

approx.bf.estimates <- function (z, V, type, suffix=NULL, sdY=1) {
  sd.prior <- if (type == "quant") { 0.15*sdY } else { 0.2 }
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

combine.abf <- function(l1, l2, p1, p2, p12) {
  lsum <- l1 + l2
  lH0.abf <- 0
  lH1.abf <- log(p1) + logsum(l1)
  lH2.abf <- log(p2) + logsum(l2)
  lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
  lH4.abf <- log(p12) + logsum(lsum)

  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  print(signif(pp.abf,3))
  print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  return(pp.abf)
}

sdY.est <- function(vbeta, maf, n) {
    warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
    oneover <- 1/vbeta
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    return(sqrt(cf))
}

process.dataset <- function(d, suffix) {
  #message('Processing dataset')

  nd <- names(d)
  if (! 'type' %in% nd)
    stop("dataset ",suffix,": ",'The variable type must be set, otherwise the Bayes factors cannot be computed')

  if(!(d$type %in% c("quant","cc")))
      stop("dataset ",suffix,": ","type must be quant or cc")

  if(d$type=="cc") {
      if(! "s" %in% nd)
          stop("dataset ",suffix,": ","please give s, proportion of samples who are cases")
      if("pvalues" %in% nd && !( "MAF" %in% nd))
          stop("dataset ",suffix,": ","please give MAF if using p values")
      if(d$s<=0 || d$s>=1)
          stop("dataset ",suffix,": ","s must be between 0 and 1")
  }

  if(d$type=="quant") {
      if(!("sdY" %in% nd || ("MAF" %in% nd && "N" %in% nd )))
          stop("dataset ",suffix,": ","must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
  }

  if("beta" %in% nd && "varbeta" %in% nd) {  ## use beta/varbeta.  sdY should be estimated by now for quant
    if(length(d$beta) != length(d$varbeta))
      stop("dataset ",suffix,": ","Length of the beta vectors and variance vectors must match")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$beta))
    if(length(d$snp) != length(d$beta))
      stop("dataset ",suffix,": ","Length of snp names and beta vectors must match")

    if(d$type=="quant" && !('sdY' %in% nd))
          d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
   df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df$snp <- as.character(d$snp)
    return(df)
  }

  if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) { ## no beta/varbeta: use p value / MAF approximation
    if (length(d$pvalues) != length(d$MAF))
      stop('Length of the P-value vectors and MAF vector must match')
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$pvalues))
    df <- data.frame(pvalues = d$pvalues,
                     MAF = d$MAF,
                     snp=as.character(d$snp))
    colnames(df)[-3] <- paste(colnames(df)[-3], suffix, sep=".")
    df <- subset(df, df$MAF>0 & df$pvalues>0) # all p values and MAF > 0
    abf <- approx.bf.p(p=df$pvalues, f=df$MAF, type=d$type, N=d$N, s=d$s, suffix=suffix)
    df <- cbind(df, abf)
    return(df)
  }

  stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}

coloc.abf <- function(dataset1, dataset2, MAF=NULL,
                      p1=1e-4, p2=1e-4, p12=1e-5) {

  if(!is.list(dataset1) || !is.list(dataset2))
    stop("dataset1 and dataset2 must be lists.")
  if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF

  df1 <- process.dataset(d=dataset1, suffix="df1")
  df2 <- process.dataset(d=dataset2, suffix="df2")
  # the na.omit is the fix here:
  merged.df <- na.omit(merge(df1,df2))

   if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)

  pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)

  output<-list(summary=results, results=merged.df)
  return(output)
}

############################## above here are general coloc package functions


# example file
# due to constraints in the computing environment,
# individual genes were tested first before parallelization

# head -2 input_per_gene/Whole_Blood/Coronary_Artery_Disease/ENSG00000000457.13_for_coloc.txt
# gene_id	eqtl_variant_id	tss_distance	eqtl_ma_samples	eqtl_ma_count	eqtl_maf	eqtl_pval_nominal	eqtl_slope	eqtl_slope_se	tissue	eqtl_sample_size	gwas_variant_id	chromosome	position	eqtl_effect_allele	eqtl_non_effect_allele	current_build	frequency	gwas_sample_size	gwas_zscore	gwas_pvalue	gwas_effect_size	gwas_standard_error	gwas_imputation_status	gwas_n_cases	trait
# ENSG00000000457.13	chr1_168894411_A_T_b38	-999856	30	30	0.0223881	0.4447680000000001	0.0509097	0.0665773	Whole_Blood	670	rs114383479	chr1	168894411	T	A	hg38	0.023088023088023088	184305	1.3586797714233398	0.17424808484624	NA	NA	imputed	60801	Coronary_Artery_Disease

# this script assumes input files contains per-gene or per-ldblock variants only

# inputfile = 'ENSG00000000457.13_for_coloc.txt'
inputfile = args[1]
tissueid = args[2]
traitid = args[3]

workdir = '/project/chrbrolab/gtex/coloc_v8/'
# inputdir = paste(workdir, 'input_per_gene/', sep='')
outputdir = paste(workdir, 'output_per_gene_eur/', tissueid, '/', traitid, '/', sep='')
# file names are changed to reflect new pheno labels
priorsfile = 'ypark_enloc_enrichment_and_prior_estimation_results.txt'

geneid = unlist(strsplit(inputfile, '_for_'))[1]

print(c(tissueid, traitid, geneid))

# before organizing
# inputdir = paste(workdir, 'input_per_gene/', sep='')
# after organizing
inputdir = paste(workdir, 'input_per_gene_eur/', tissueid, '/', traitid, '/', sep='')
outfile = paste(outputdir, geneid, '_coloc_output_file.txt', sep='')
output = c()

inputdata = tbl_df(fread(paste(inputdir, inputfile, sep='')))
inputdata = inputdata[inputdata$frequency > 0,]
inputdata$gwas_varbeta = 1

inputdata = inputdata %>% filter(complete.cases(gwas_zscore)) %>% mutate(max_gwas_zscore = max(abs(gwas_zscore)))
inputdata = inputdata %>% arrange(desc(abs(gwas_zscore)))

priors = tbl_df(fread(paste(workdir, priorsfile, sep='')))
priorsused = priors %>% filter(trait == traitid & tissue == tissueid)
p1used = priorsused$p1_trait_prior
p2used = priorsused$p2_tissue_prior
p12used = priorsused$p12

# eqtlinput = list(beta = inputdata$eqtl_slope, varbeta = (inputdata$eqtl_slope_se)**2, N = inputdata$eqtl_sample_size, type = 'quant', snp = inputdata$eqtl_variant_id, MAF = inputdata$eqtl_maf)
eqtlinput = list(pvalues = inputdata$eqtl_pval_nominal, N = inputdata$eqtl_sample_size, type = 'quant', snp = inputdata$eqtl_variant_id, MAF = inputdata$frequency)


caseproportion = inputdata$gwas_n_cases/inputdata$gwas_sample_size
datatype = ifelse(is.na(caseproportion[1]), 'quant', 'cc')

if(datatype == 'quant') {
  gwasinput = list(beta = inputdata$gwas_zscore, varbeta = inputdata$gwas_varbeta, N = inputdata$gwas_sample_size, type = datatype, snp = inputdata$eqtl_variant_id, MAF = inputdata$frequency)
} else {
  gwasinput = list(beta = inputdata$gwas_zscore, varbeta = inputdata$gwas_varbeta, N = inputdata$gwas_sample_size, type = datatype, snp = inputdata$eqtl_variant_id, MAF = inputdata$frequency, s = caseproportion)
}

# if(datatype == 'quant') {
#   gwasinput = list(pvalues = inputdata$gwas_pvalue, N = inputdata$gwas_sample_size, type = datatype, snp = inputdata$eqtl_variant_id, MAF = inputdata$frequency)
# } else {
#   gwasinput = list(pvalues = inputdata$gwas_pvalue, N = inputdata$gwas_sample_size, type = datatype, snp = inputdata$eqtl_variant_id, MAF = inputdata$frequency, s = caseproportion)
# }

# priors were stored in a table like this:
# tissue	trait	p1_trait_prior	p2_tissue_prior	p12	pm1	pm2	enloc_enrich_intercept_effect(inteffect)	inteffect_95ci_lower	inteffect_95ci_upper	enloc_enrich_log_odds_ratio_effect(logor)	logor_95ci_lower	logor_95ci_upper
# Adipose_Subcutaneous	Adiponectin	7.31680E-06	3.32503E-03	3.05385E-08	7.34734E-06	3.32506E-03	-11.822	-11.849	-11.795	0.224	-2.384	2.832
# Adipose_Subcutaneous	Alzheimers_Disease	1.32258E-05	3.32490E-03	1.60605E-07	1.33864E-05	3.32506E-03	-11.230	-11.409	-11.052	1.292	-2.956	5.540

runcoloc = coloc.abf(gwasinput, eqtlinput, p1=p1used, p2=p2used, p12=p12used)
output = rbind(output, c(unlist(inputdata[1,]), runcoloc$summary))
output = as.data.table(output)

output = mutate(output, PP.H1.abf = ifelse(PP.H0.abf == 1, 0, PP.H1.abf))
output = mutate(output, PP.H2.abf = ifelse(PP.H0.abf == 1, 0, PP.H2.abf))
output = mutate(output, PP.H3.abf = ifelse(PP.H0.abf == 1, 0, PP.H3.abf))
output = mutate(output, PP.H4.abf = ifelse(PP.H0.abf == 1, 0, PP.H4.abf))

outputdata = tbl_df(as.data.frame(output)) %>% mutate(PP.H4.abf = as.numeric(as.character(PP.H4.abf))) %>% arrange(desc(PP.H4.abf))
write.table(outputdata, file=outfile, quote=F, col.names=T, row.names=F, sep='\t')




