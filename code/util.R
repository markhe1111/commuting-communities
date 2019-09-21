
bh_reject_p <- function (pvals, alpha,p) {
  pvals_adj <- length(pvals) * pvals / rank(pvals)^(1+p)
  if (sum(pvals_adj <= alpha) > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}

bh_reject <- function (pvals, alpha) {
  pvals_adj <- length(pvals) * pvals / rank(pvals)
  
  sum.pval = sum(pvals_adj <= alpha) 
  if(is.na(sum.pval)) sum.pval <- 0
  
  if (sum.pval > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}

#bh_reject <- function (pvals, alpha) {
#  pvals_adj <- length(pvals) * pvals / rank(pvals)
#  if (sum(pvals_adj <= alpha) > 0) {
#    thres <- max(pvals[pvals_adj <= alpha])
#    return(which(pvals <= thres))
#  } else {
#    return(integer(0))
#  }
#}
bh_rejectR_old <- function (pvals, alpha, conserv = T) {
  
  m <- length(pvals)
  
  if (!conserv) {
    pvals_adj <- m * pvals / rank(pvals, ties.method = "first")
  } else {
    mults <- sum(1 / c(1:m))
    pvals_adj <- mults * m * pvals / rank(pvals, ties.method = "first")
  }
  
  if (sum(pvals_adj <= alpha) > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
  
} # orders p values, including maximum, and finds threshold, and returns pvals below thresh


bh_reject_sq <- function (pvals, alpha) {
  
  pvals_adj <- (length(pvals)) * (pvals)^2 / rank(pvals)
  
  sum.pval = sum(pvals_adj <= alpha) 
  if(is.na(sum.pval)) sum.pval <- 0
  
  if (sum.pval > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}


bh_reject_sqrt <- function (pvals, alpha) {
  
  pvals_adj <- (length(pvals)) * sqrt(pvals) / rank(pvals)
  
  sum.pval = sum(pvals_adj <= alpha) 
  if(is.na(sum.pval)) sum.pval <- 0
  
  if (sum.pval > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}




bh_rejectR <- function (pvals, alpha,  conserv = F) {
  pvals_adj <- length(pvals) * pvals / rank(pvals)
  
  sum.pval = sum(pvals_adj <= alpha) 
  if(is.na(sum.pval)) sum.pval <- 0
  
  if (sum.pval > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}

bhy <- function(pvals, alpha = 0.05){
    
    # Sorted p-vals
    sp = sort(pvals)
    
    # Save original order of p-vals
    ord = order(pvals)
    
    # Find bhy cutoff
    nums = 1:length(pvals)
    cms = cumsum(1/nums)
    
    # Find which p-vals are less than bh cutoff
    under = sp < (nums/(length(pvals)*cms)*alpha)
    
    # Return indices of significant p-vals
    if(sum(under) == 0){
      return(c())
    }else{
      cutoff = max(which(under))
      return(ord[1:cutoff])
    }
  }



symdiff <- function (s1, s2) {
  return(union(setdiff(s1, s2), setdiff(s2, s1)))
}
jaccard <- function (s1, s2) {
  return(length(symdiff(s1, s2)) / length(union(s1, s2)))
}


make.log.dedup<-function(init.yr){
  #z.thr =   init.yr$Z.THR[[1]]
  z.thr =   init.yr$Z.THR
  
  pop = z.thr$weight
  z.thr$weight = log(pop)  
  
  z.dedup = z.thr[!duplicated(z.thr),]
  z.dedup
} # this might be raw pop

make.log10<-function(init.yr){
  #z.thr =   init.yr$Z.THR[[1]]
  z.thr =   init.yr$Z.THR
  
  pop = z.thr$weight
  z.thr$weight = log10(pop)  
  z.thr
} # this might be raw pop
make.log10.dedup<-function(init.yr){
  #z.thr =   init.yr$Z.THR[[1]]
  z.thr =   init.yr$Z.THR
  
  pop = z.thr$weight
  z.thr$weight = log10(pop)  
  
  z.dedup = z.thr[!duplicated(z.thr),]
  z.dedup
} # this might be raw pop
sqrt.dedup<-function(init.yr){
  #z.thr =   init.yr$Z.THR[[1]]
  z.thr =   init.yr$Z.THR
  
  pop = z.thr$weight
  z.thr$weight = sqrt(pop)  
  
  z.dedup = z.thr[!duplicated(z.thr),]
  z.dedup
} # this might be raw pop
power.dedup<-function(init.yr, pow){
  #z.thr =   init.yr$Z.THR[[1]]
  z.thr =   init.yr$Z.THR
  
  pop = z.thr$weight
  z.thr$weight = (pop)^pow  
  
  z.dedup = z.thr[!duplicated(z.thr),]
  z.dedup
} # this might be raw pop



make.ln<-function(init.yr){
  #z.thr =   init.yr$Z.THR[[1]]
  z.thr =   init.yr$Z.THR
  
  pop = z.thr$weight
  z.thr$weight = log(pop)  
  z.thr
} # this might be raw pop

nonlog<-function(init.yr){
  #z.thr =   init.yr$Z.THR[[1]]
  z.thr =   init.yr$Z.THR
  z.thr$weight = (z.thr$weight)
  z.thr
} # this might be raw pop

