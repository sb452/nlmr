#' Instrumental variable analysis using fractional polynomials based on summary data
#'
#' @description frac_poly_summ_mr performs instumental variable analysis by fitting fractional polynomial models to localised average causal effects using meta-regression.
#'
#' Please note that if you provide synthetic data that are too regular (eg the by associations are all zero or the xmean values are exactly 1, 2, 3, ...), the function may error as several fractional polynomials provide the same fit.
#'
#' @param by vector of gene-outcome associations.
#' @param bx vector of gene-exposure associations.
#' @param byse vector of standard errors of gene-outcome associations.
#' @param bxse vector of standard errors of gene-exposure associations.
#' @param xmean average value of the exposure in each stratum (or whatever summary of the exposure level in the stratum is desired).
#' @param method meta-regression method parsed to the rma package. The default is fixed-effects ("FE").
#' @param d fractional polynomial degree. The default is degree 1. The other options are: 1, 2, or "both".
#' @param pd p-value cut-off for choosing the best-fitting fractional polynomial of degree 2 over the best-fitting fractional polynomial degree 1. This option is only used if d="both". The default is 0.05.
#' @param ci the type of 95\% confidence interval. There are three options: (i) using the model standard errors ("model_se"), (ii) using bootstrap standard errors ("bootstrap_se"), (iii) using bootstrap percentile confidence intervals ("bootstrap_per"). The default is the model standard errors.
#' @param nboot the number of bootstrap replications (if required). The default is 100 replications.
#' @return model the model specifications. The first column is the number of quantiles (q); the second column is the position used to relate x to the LACE in each quantiles (xpos); the third column is the type of confidence interval constructed (ci); the fourth column is the number of bootstrap replications performed (nboot).
#' @return powers the powers of the chosen polynomial.
#' @return coefficients the regression estimates. The first column is the regression coefficients (beta); the second column is the standard errors of regression coefficients (se); the third column is the lower confidence interval (lci); the fourth column is the upper confidence interval (uci); the fifth column is the p-value (pval).
#' @return lace the localised average causal effect estimate in each quantile. The first column is the regression coefficients (beta); the second column is the standard errors of regression coefficients (se); the third column is the lower confidence interval (lci); the fourth column is the upper confidence interval (uci); the fifth column is the p-value (pval).
#' @return xcoef the association between the instrument and the exposure in each quantile. The first column is the regression coefficients (beta); the second column is the standard errors of regression coefficients (se).
#' @return p_tests the p-value of the non-linearity tests. The first column is the p-value of the test between the fractional polynomial degrees (fp_d1_d2);
#'         the second column is the p-value from the fractional polynomial non-linearity test (fp);
#'         the third column is the p-value from the quadratic test (quad); the fourth column is the p-value from the Cochran Q test (Q).
#' @return p_heterogeneity the p-value of heterogeneity. The first column is the p-value of the Cochran Q heterogeneity test (Q); the second column is the p-value from the trend test (trend).
#' @author Stephen Burgess <sb452@medschl.cam.ac.uk>, leading heavily on James R Staley <js16174@bristol.ac.uk>
#' @export
frac_poly_summ_mr <- function(by, bx, byse, bxse, xmean, method="FE", d=1, pd=0.05, ci="model_se", nboot=100){
  
  ##### Error messages #####
  if(!(d==1 | d==2 | d=="both")) stop('the degree has to be equal to 1, 2 or "both"')
  if(!(ci=="model_se" | ci=="bootstrap_se" | ci=="bootstrap_per")) stop('the confidence intervals must be one of "model_se", "bootstrap_se" and "bootstrap_per"')
  
      frac_coef    = by
      frac_se      = byse
      xcoef_sub    = bx
      xcoef_sub_se = bxse
      xcoef = sum(bx*bxse^-2)/sum(bxse^-2)
      q = length(by)
  
  ##### Test of IV-exposure assumption #####
  p_het <- 1- pchisq(rma(xcoef_sub, vi=(xcoef_sub_se)^2)$QE, df=(q-1))
  p_het_trend <- rma.uni(xcoef_sub ~ xmean, vi=xcoef_sub_se^2, method=method)$pval[2]
  
  ##### Best-fitting fractional polynomial of degree 1 #####
  powers<-c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
  p<-NULL
  ML<-NULL
  j<-1
  for(p1 in powers){
    if(p1==-1){x1 <- xmean^p1}else{x1 <- (p1+1)*xmean^p1}
    p[j]<-p1
    ML[j] <- summary(rma(frac_coef/xcoef ~ -1 + x1, vi=(frac_se/xcoef)^2, method=method))$fit.stats[1,1]
    j<-j+1
  }
  fits <- data.frame(p, ML)
  fits$max <- 0
  fits$max[fits$ML==max(fits$ML)] <- 1
  p_ML <- fits$p[fits$max==1]
  
  ##### Best-fitting fractional polynomial of degree 2 #####
  if(d==1 | d==2 | d=="both"){
    powers1 <- c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
    powers2 <- c(0, -3, -2, -1.5, -1, -0.5, 1, 2)
    p1<-NULL
    p2 <-NULL
    ML <- NULL
    j <- 1
    for(p11 in powers1){
      if(p11==-1){x1 <- xmean^p11}else{x1 <- (p11+1)*xmean^p11}
      for(p21 in powers2){
        if(p11==p21){if(p21==-1){x2 <- 2*(xmean^p21)*log(xmean)}else{x2 <- ((p21+1)*(xmean^p21)*log(xmean) + xmean^p21)}}
        else{if(p21==-1){x2 <- xmean^p21}else{x2 <- (p21+1)*xmean^p21}}
        p1[j]<-p11
        p2[j]<-p21
        cc <- try(rma(frac_coef/xcoef ~ -1 + x1 + x2, vi=(frac_se/xcoef)^2, method=method), silent=TRUE)
        if(is(cc, "try-error")==T){ML[j] <- NA}
        if(is(cc, "try-error")==F){ML[j] <- summary(rma(frac_coef/xcoef ~ -1 + x1 + x2, vi=(frac_se/xcoef)^2, method=method))$fit.stats[1,1]}
        j<-j+1
      }
      powers2<-powers2[-1]
    }
    fits <- data.frame(p1, p2, ML)
    fits$max <- 0
    fits$max[fits$ML==max(fits$ML, na.rm=T)] <- 1
    p1_ML <- fits$p1[fits$max==1]
    p2_ML <- fits$p2[fits$max==1]
  }
  
  ##### Best-fitting fractional polynomial of either degree 1 or degree 2 #####
  p_d1_d2 <- NA
  if(d==1 | d==2 | d=="both"){
    if(p_ML==-1){x1<-xmean^p_ML}else{x1 <- (p_ML+1)*xmean^p_ML}
    best_fracp_d1 <- rma(frac_coef/xcoef ~ -1 + x1, vi=(frac_se/xcoef)^2, method=method)
    dev_best_fracp_d1 <- best_fracp_d1$fit.stats[2,1]
    if(p1_ML==-1){x1 <- xmean^p1_ML}else{x1 <- (p1_ML+1)*xmean^p1_ML}
    if(p1_ML==p2_ML){if(p2_ML==-1){x2 <- 2*(xmean^p2_ML)*log(xmean)}else{x2 <- ((p2_ML+1)*(xmean^p2_ML)*log(xmean) + xmean^p2_ML)}}
    else{if(p2_ML==-1){x2 <- xmean^p2_ML}else{x2 <- (p2_ML+1)*xmean^p2_ML}}
    best_fracp_d2 <- rma(frac_coef/xcoef ~ -1 + x1 + x2, vi=(frac_se/xcoef)^2, method=method)
    dev_best_fracp_d2 <- best_fracp_d2$fit.stats[2,1]
    p_d1_d2 <- 1 - pchisq((dev_best_fracp_d1 - dev_best_fracp_d2), df=2)
    if(p_d1_d2>=pd){d1 <- 1}else{d1 <- 2}
    if(d=="both"){d <- d1}
  }
  
  ##### Model #####
  if(d==1){if(p_ML==-1){x1<-xmean^p_ML}else{x1 <- (p_ML+1)*xmean^p_ML}; model <- rma(frac_coef/xcoef ~ -1 + x1, vi=(frac_se/xcoef)^2, method=method)}
  if(d==2){if(p1_ML==-1){x1<-xmean^p1_ML}else{x1 <- (p1_ML+1)*xmean^p1_ML}; if(p1_ML==p2_ML){if(p2_ML==-1){x2 <- 2*(xmean^p2_ML)*log(xmean)}else{x2 <- ((p2_ML+1)*(xmean^p2_ML)*log(xmean) + xmean^p2_ML)}}else{if(p2_ML==-1){x2 <- xmean^p2_ML}else{x2 <- (p2_ML+1)*xmean^p2_ML}}; model <- rma(frac_coef/xcoef ~ -1 + x1 + x2, vi=(frac_se/xcoef)^2, method=method)}
  
  ##### Bootstrap #####
  if(ci=="bootstrap_per" | ci=="bootstrap_se"){
    if(d==1){
      frac_coef_boot <- NULL
      for(i in 1:nboot){
            frac_coef1 <- by + rnorm(q, 0, byse)
            frac_se1 <- byse
         xmean1 <- xmean
        if(p_ML==-1){x111<-xmean1^p_ML}else{x111 <- (p_ML+1)*xmean1^p_ML}
        mod <- rma.uni(frac_coef1/xcoef ~ -1 + x111, vi=(frac_se1/xcoef)^2, method=method)
        frac_coef_boot[i] <- mod$b[1]
      }
    }
    if(d==2){
      frac_coef_boot <- matrix(, nrow = nboot, ncol = 2)
      for(i in 1:nboot){
            frac_coef1 <- by + rnorm(q, 0, byse)
            frac_se1 <- byse
         xmean1 <- xmean
        if(p1_ML==-1){x111<-xmean1^p1_ML}else{x111 <- (p1_ML+1)*xmean1^p1_ML}
        if(p1_ML==p2_ML){if(p2_ML==-1){x211 <- 2*(xmean1^p2_ML)*log(xmean1)}else{x211 <- ((p2_ML+1)*(xmean1^p2_ML)*log(xmean1) + xmean1^p2_ML)}}
        else{if(p2_ML==-1){x211 <- xmean1^p2_ML}else{x211 <- (p2_ML+1)*xmean1^p2_ML}}
        mod <- rma.uni(frac_coef1/xcoef ~ -1 + x111 + x211, vi=(frac_se1/xcoef)^2, method=method)
        frac_coef_boot[i,1] <- mod$b[1]
        frac_coef_boot[i,2] <- mod$b[2]
      }
    }
  }
  
  ##### Fractional polynomial degree 1 test against linearity #####
  if(p_ML==-1){x1<-xmean^p_ML}else{x1 <- (p_ML+1)*xmean^p_ML}
  linear <- rma(frac_coef/xcoef ~ 1, vi=(frac_se/xcoef)^2, method=method)
  dev_linear <- linear$fit.stats[2,1]
  best_fracp_d1 <- rma(frac_coef/xcoef ~ -1 + x1, vi=(frac_se/xcoef)^2, method=method)
  dev_best_fracp_d1 <- best_fracp_d1$fit.stats[2,1]
  p_fp <- 1 - pchisq((dev_linear - dev_best_fracp_d1), df=1)
  
  ##### Other tests #####
  p_quadratic <- rma(frac_coef/xcoef ~ xmean, vi=(frac_se/xcoef)^2, method=method)$pval[2]
  p_Q <- 1 - pchisq(rma(frac_coef/xcoef, vi=(frac_se/xcoef)^2)$QE, df=(q-1))
  
  ##### Results #####
  beta <- as.numeric(model$b)
  if(ci=="model_se"){if(d==1){powers <- p_ML + 1}; if(d==2){powers <- c(p1_ML, p2_ML); powers <- powers + 1}; cov <- model$vb; se <- model$se; lci <- beta - 1.96*se; uci <- beta + 1.96*se; pval <- 2*pnorm(-abs(beta/se))}
  if(ci=="bootstrap_se"){if(d==1){powers <- p_ML + 1; cov <- var(frac_coef_boot); se <- sqrt(cov)}; if(d==2){powers <- c(p1_ML, p2_ML); powers <- powers + 1; cov <- cov(frac_coef_boot); se <- sqrt(diag(cov))}; lci <- beta - 1.96*se; uci <- beta + 1.96*se; pval <- 2*pnorm(-abs(beta/se))}
  if(ci=="bootstrap_per"){if(d==1){powers <- p_ML + 1; se <- NA; lci <- quantile(frac_coef_boot, probs=0.025); uci <- quantile(frac_coef_boot, probs=0.975); pval <- NA}; if(d==2){powers <- c(p1_ML, p2_ML); powers <- powers + 1; se <- rep(NA, 2); lci <- NULL; uci <- NULL; pval <- NULL; lci[1] <- quantile(frac_coef_boot[,1], probs=0.025); lci[2] <- quantile(frac_coef_boot[,2], probs=0.025); uci[1] <- quantile(frac_coef_boot[,1], probs=0.975); uci[2] <- quantile(frac_coef_boot[,2], probs=0.975); pval <- rep(NA,2)}}
  lci <- as.numeric(lci); uci <- as.numeric(uci)
  if(ci=="model_se"){nboot<-NA}
  
  ##### Return #####
  model <- as.matrix(data.frame(q=q, xpos="user", ci_type=ci, nboot=nboot))
  coefficients <- as.matrix(data.frame(beta=beta, se=se, lci=lci, uci=uci, pval=pval))
  rownames(coefficients) <- powers
  lace <- as.matrix(data.frame(beta=(frac_coef/xcoef), se=(abs(frac_se/xcoef)), lci=(frac_coef/xcoef - 1.96*(abs(frac_se/xcoef))), uci=(frac_coef/xcoef + 1.96*(abs(frac_se/xcoef))), pval=(2*pnorm(-abs(frac_coef/frac_se)))))
  rownames(lace) <- 1:nrow(lace)
  xcoef_quant <- as.matrix(data.frame(beta=xcoef_sub, se=xcoef_sub_se))
  rownames(xcoef_quant) <- 1:nrow(xcoef_quant)
  p_tests <- as.matrix(data.frame(fp_d1_d2=p_d1_d2, fp=p_fp, quad=p_quadratic, Q=p_Q))
  p_heterogeneity <- as.matrix(data.frame(Q=p_het, trend=p_het_trend))

  results <- list(n=NA, model=model, powers=powers, coefficients=coefficients, lace=lace, xcoef=xcoef_quant, p_tests=p_tests, p_heterogeneity=p_heterogeneity)
  class(results) <- "frac_poly_mr"
  return(results)
}
