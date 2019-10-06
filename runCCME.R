# 2009 crucial year: California, Arizona, Florida unstable:  Massachu, Oreg, NC More Stable
rm(list = ls())

x <- c( "rgdal", "rgeos", "maptools", 'Matrix', 'parallel'); 
lapply(x, library, character.only = TRUE) # load

source('code/util.R')
EDC = read.csv('data/C2010_withDC.csv')


make.Adj.mat<-function(E, n=n){
  
  A <-matrix(nc=n, nr = n, data = 0)
  
  for(k in 1:nrow(E)){
    a =  (E[k,])
    A[a[1], a[2]] = a[3]
    A[a[2], a[1]] = a[3]
  }
  return(A)
  
}
register_EDC2<-function(EDC, key){
  library(dplyr)
  
  # get region split data
  
  run_regionsplit_100<-function( ){
    
    setwd('/Users/markhe 1/Dropbox/regiondemarcation/Regionsplitdata/')
    setwd('threshold_100/')
    
    csv_files = list.files()
    #n = nrow(key)
    
    regional_result = list()
    keys = list()
    
    Es = list()
    
    for(i in seq_along(csv_files)){
      
      csv = read.csv(csv_files[i])
      csv$from <- NULL
      csv$to <- NULL
      
      #  A = make.ad
      Es[[i]] = unique (csv )
    }
    EE = unique( do.call(rbind, Es) )
    EE
  }
  
  
  library(plyr)
  
  EE = run_regionsplit_100()
  make.key<-function(EE){
    
    from = select( EE , from_name, from_ctyname )
    to = select( EE , to_name, to_ctyname )
    
    from.key = unique(from)
    to.key = unique(to)
    
    colnames(from.key) <-colnames(to.key) <- c('fips','name')
    key = unique(rbind(from.key, to.key))
    key = data.frame(node=seq(nrow(key)), name=key$name, fips= key$fips)
    sorted = key[ order(key$fips),]
    sorted$node = seq(nrow(key))
    sorted
  }
  key <<- make.key(EE)
  
  
  
  J = EDC[EDC$commuters>99,]
  
  J$fromfips = J$node1
  J$tofips = J$node2
  
  
  howmany  = unique(c(J$node1, J$node2))
  
  unprocessed = setdiff(  unique(c(J$node1, J$node2)), key$fips )
  
  
  J[J$node1 %in% unprocessed,] <- NA # remove unprocessed for now!
  J[J$node2 %in% unprocessed,] <-NA
  
  
  Z = na.omit(J)
  
  
  key1 = key
  
  key2 = key
  
  
  key1$fromfips = key$fips
  key2$tofips = key$fips
  
  newN01 = join(Z, key1)
  newN11 = select(newN01,   f1=node1, f2=node2, commuters, node1 = node, tofips)
  
  newN2 = join(newN11, key2)
  edgelist = select(newN2, node1 , node2=node, commuters)
  
  edgelist 
  
} #uses join
DC.2010 = register_EDC2(EDC, key)
n = nrow(key)
init.yr  = DC.2010
E = ( ( DC.2010 ))

init = E[E$node1==E$node2 & E$commuters>20000, ] #subsets nodes that have large self loops, to initialize

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


#getting monads 

get_island_sample  <- function(E, sig = .05, key_j ){
  
  n = nrow(key_j)
  
  A = make.Adj.mat(as.matrix(E), n = nrow(key_j))
  A0 = A
  diag(A0)<- 0
  
  d <- colSums(A0 > 0)
  S <- colSums(A)
  
  sT = sum(S)
  dT = sum(d)
  SL = diag(A)
  
  
  rho =SL/S 
  p = sum(SL)/ sum(A)
  
  
  # get good a
  sample_var = 1/(n-1) *sum((na.omit(rho) - p)^2)
  
  a = -p*(p^2 - p - sample_var) /sample_var
  b = a*(1-p)/p
  
  
  J <<- na.omit( data.frame(key, S, SL, rho) )
  
  
  sample_hub<-function(J){
    
    node = J$node
    random_rho = rbeta(a,b,n = nrow(J))
    
    F_sl = J$SL -  (random_rho* J$S)
    
    qtl = quantile(F_sl , 1-sig)
    
    
    emp =  J$SL -  (p*J$S) 
    
    which_over = J$SL -  (p*J$S) > qtl
    
    node[which_over]
  }
  
  
  hub_samples = lapply(1:10000, function(x) sample_hub(J))
  
  hub_sam = table( unlist(hub_samples))
  good = cbind( key[ as.numeric( names(hub_sam)), ] , hub_sam)
  
  hubs = good[,c(-1, -4)]
  good
  #hubs_gooder = good[good$Freq>5000,]
  
  #hubs_gooder
  
}
monads = get_island_sample( E , sig = .05, key_j=key)
monads_05 = monads[monads$Freq==10000,]

key_j = key

Comm_detect   <- function ( E,  alpha  , init_nodes, key_j ) {
  
  each.B<-function(adjMat, B){
    if(length(B)==1){ 
      nodesToB <-  which( ( adjMat[B, ] ) > 0 ) 
    }else{
      nodesToB <- which(colSums(adjMat[B, ]) > 0)
    }
    nodesToB
  }
  
  
  n = nrow(key_j)
  
  A = make.Adj.mat(as.matrix(E), n = nrow(key_j))
  A0 = A
  diag(A0)<- 0
  
  d <- colSums(A0 > 0)
  S <- colSums(A)
  
  sT = sum(S)
  dT = sum(d)
  SL = diag(A)
  rho =SL/S 
  p = sum(SL)/ sT

  
  nn = length(na.omit(rho))
  sample_var = 1/(nn-1) *sum( (na.omit(rho) - p)^2)
  a = -p*(p^2 - p - sample_var) /sample_var
  
  
  var_at_a_for_beta <-function(a,p){
    b = a*(1-p)/p
    numerator = a*b
    denominator = (a+b)^2*(a+b+1)
    numerator/denominator
  }
  
  get_K_nonSL_B <-function(A, S, d, p ){
    
    zeta <-function(u, v){
      a = (  S[u]*S[v] / sT )
      b =  ( d[u]*d[v] / dT )
      if(b==0) b = 1/dT
      a/b
    }
    
    numerator_matrix<- matrix(nrow = n, ncol=n, data = 0)
    denominator_matrix<- matrix(nrow = n, ncol=n, data = 0)
    
    for(u in 1:n){
      for(v in 1:n){
        if(v !=u){
          zeta_uv = zeta(u,v)
          numerator_matrix [u,v] = (A[u,v] - (1-p)*zeta(u,v))^2
          denominator_matrix [u,v] = (1-p)^2*zeta(u,v)^2
        }
      }
    }
    
    Kappa_notself = sum(numerator_matrix)  / sum( denominator_matrix)
    return(Kappa_notself) 
  }
  
  get_sum_of_zeta <-function(S,d,sT,dT,A ){
    
    zeta <-function(u, v){
      a = (  S[u]*S[v] / sT )
      b =  ( d[u]*d[v] / dT )
      if(b==0) b = 1/dT
      a/b
    }
    
    sumsq_E_matrix<- matrix(nrow = n, ncol=n, data = 0)
    zetasq_matrix<- matrix(nrow = n, ncol=n, data = 0)
    
    for(u in 1:n){
      for(v in 1:n){
        if(v !=u){
          zeta_uv = zeta(u,v)
          sumsq_E_matrix [u,v] = (A[u,v] - (1-p)*zeta(u,v))^2
          zetasq_matrix [u,v] = zeta(u,v)^2
        }
      }
    }
    sum_zeta_sq = sum(zetasq_matrix)
    sum_sq_E_sq = sum(sumsq_E_matrix)
    
    list(sum_zeta_sq=sum_zeta_sq,
         sum_sq_E_sq=sum_sq_E_sq)
  }
  
  sums = get_sum_of_zeta(S,d,sT,dT,A)
  sum_sq_E_sq =  sums$sum_sq_E_sq
  sum_zeta_sq =  sums$sum_zeta_sq
  
  get_K_nonSL <-function(  p, variance_rho ,   sum_sq_E_sq, sum_zeta_sq ){
    
    numerator = sum_sq_E_sq - variance_rho * sum_zeta_sq
    denominator = (variance_rho+(1-p)^2 )  * sum_zeta_sq
    
    Kappa_notself = numerator/denominator
    return(Kappa_notself) 
  }
  get_K_SL<-function(SL, p, S, variance_rho ){
    
    sample_var = sum( (SL - p*S)^2 )
    numerator = sample_var - sum( S^2*variance_rho)
    denom_u =    S^2*(variance_rho + p^2 )
    denominator = sum(   na.omit(denom_u))
    Kappa = numerator/denominator 
    
    return(Kappa)
    
  }

  variance_rho =sample_var 
  
  
  K_nSL = get_K_nonSL( p, variance_rho , sum_sq_E_sq, sum_zeta_sq)
  K_SL = get_K_SL ( SL,p,S,variance_rho)
  
  
  init_Set =  lapply(init_nodes, function(B) each.B(A, B))

  
  pvalFun<- function (B,   A ) {
    
    if(length(B)==1){ 
      nodesToB <-  which( ( A [B, ] ) > 0 ) 
    }else{
      nodesToB <-  which(colSums(A [B, ]) > 0)
    }
    
    stats <- rowSums(  A [nodesToB, B, drop = FALSE])
    
    
    meanfun <-function( B, nodesToB , p, S ,sT){
      
      meansB = list()
      
      for( u in nodesToB ){
        
        meansToB = list()
        
        for(v in B){
          if(v!=u) meansToB[[v]] =  S[v]  
        }
        
        other_sum  = (1-p) *S[u]* sum(unlist(meansToB))/sT
        # self = S[u]*p
        meansB[[u]] = other_sum # + self
      }
      
      unlist(meansB)
      
    }
    
    varsfun <-function( B,nodesToB, 
                        p, S, d, sT, dT,
                        K_SL, K_nSL,
                        variance_rho){
      
      varsB = list()
      
      for( u in nodesToB ){
        
        varsToB = list()
        
        for(v in B){
          if(v!=u) {
            
            d_uv = (d[u]*d[v])/dT 
            
            r_uv = (S[u]*S[v] /sT )^2 / d_uv 
            
            varsToB[[v]] =  r_uv  *  (  (variance_rho+ (1-p)^2) *K_nSL + variance_rho + 1 - d_uv )   
          }
        }
        
        other_sum  = (1-p)^2 * sum(  unlist(varsToB) )
        # self = S[u]^2*( K_SL*(variance_rho +  p) + variance_rho)
        varsB[[u]] = other_sum   #+ self
      }
      unlist(varsB)
    }
    
    
    means = meanfun(B, nodesToB, p,S,sT)
    vars = varsfun(B,nodesToB, 
                   p, S, d, sT, dT,   K_SL, K_nSL,  variance_rho)
    if(length(stats)>0){ 
      pvalsToB <- pnorm(stats, means, sqrt(vars), lower.tail = FALSE)
    }else{
      pvalsToB = 1
    }
    
    pvals <- rep(1, n)
    pvals[ nodesToB ] <- pvalsToB
    pvals
  }
  
  extract <- function (i,init_Set, A) {
    
    #    doTheNode <- !inNodes[i] %in% unionC
    B0 <- init_Set[[i]]
    #unionInitial <<- union(unionInitial, B0)
    B_old <- B0
    B_new <- 1:n
    chain <- list(B_old)
    consec_jaccards <- NULL
    found_cycle <- found_break <- NULL
    mean_jaccards <- NULL ## add all these to update_info
    itCount <- 0
    cycledSets <- NULL
    
    cat("Beginning Updates\n", i, init_Set[[i]])
    
    while ( length(B_new) > 1) {
      
      itCount <- itCount + 1
      
      pvals <- pvalFun(B_old, A = A) # ############ pvalfun here
      
      B_new <- bh_reject(pvals, alpha)
      
      consec_jaccard <- jaccard(B_new, B_old)
      jaccards <- unlist(lapply(chain, function (B) jaccard(B, B_new)))
      found_cycle <- c(found_cycle, FALSE)
      found_break <- c(found_break, FALSE)
      consec_jaccards <- c(consec_jaccards, consec_jaccard)
      mean_jaccards <- c(mean_jaccards, mean(jaccards))
      cat(paste0("Update ", itCount,   " is size ", length(B_new), ", ","jaccard to last is ", round(consec_jaccard, 3), ", ","mean jaccard along chain is ", round(mean(jaccards), 3), "\n", sep=""))
      
      # Checking for cycles (4.4.1 in paper)
      if ( jaccard(B_new, B_old ) > 0) { # Otherwise loop will end naturally
        jaccards <- unlist(lapply(chain, function (B) jaccard(B, B_new)))
        if (sum(jaccards == 0) > 0) { # Cycle has been found
          found_cycle[itCount] <- TRUE
          Start <- max(which(jaccards == 0))
          cycle_chain <- chain[Start:length(chain)]
          
          # Checking for cycle break (4.4.1a)
          seq_pair_jaccards <- rep(0, length(cycle_chain) - 1)
          for (j in seq_along(seq_pair_jaccards)) {
            seq_pair_jaccards[j] <- jaccard(chain[[Start + j - 1]], chain[[Start + j]])
          }
          if (sum(seq_pair_jaccards == 1) > 0) {# then break needed
            cat(" ---- Break found\n")
            found_break[itCount] <- TRUE
            B_new <- NULL
            break
          }
          
          # Create conglomerate set (and check, 4.4.1b)
          B_J <- unique(unlist(cycle_chain))
          B_new <- B_J
          B_J_check <- unlist(lapply(chain, function (B) jaccard(B_J, B)))
          if (sum(B_J_check == 0) > 0) 
            break
        } # From checking jaccards to cycle_chain
      } else { # From checking B_new to B_old; if B_new = B_old: 
        break
      }
      B_old <- B_new
      chain <- c(chain, list(B_new))
      
      if(length(chain) >500 ) break
    } # From Updates
    
    unionC <<- union(B_new, unionC)
    cat(paste0(length(unionC), " vertices in communities.\n"))
    update_info <- list("mean_jaccards" = mean_jaccards,
                        "consec_jaccards" = consec_jaccards,
                        "found_cycle" = found_cycle,
                        "found_break" = found_break
    )
    return(list("comm" = B_new, "update_info" = update_info))
  } # no self, bh_y instead of bh_reject
  
  # Extractions - 
  unionC <- unionInitial <- integer(0)
  extractRes <- lapply( 1: length(init_Set) , function(x) extract(x, init_Set = init_Set, A  )) #### IMPORTANT PART
  
  
  comms <- lapply(extractRes, function (L) L$comm)
  update_info <- lapply(extractRes, function (L) L$update_info)
  
  #  Clean-up and return --
  
  nonNullIndxs <- which( unlist(lapply(comms, length)) > 0)     # pt of contention: 0 or 1?
  comms0 <- comms[nonNullIndxs]


  returnList <- list( "communities_before_OLfilt" = comms0,
                     "nonNullIndxs" = nonNullIndxs )
  return(returnList)
} # previous global baseline model as described in manuscript

filter_OL_posthoc  <- function (clus, E, tau, n, bysize=F, key_j) {
  
  
  commus = Filter(function(x) length(x)>1,    clus) 

  A = make.Adj.mat(as.matrix(E),  n = n)
  K <- length(commus)
  
  A0 = A
  diag(A0) = 0
  
  comm_str =  sapply(commus,  function(J)  {
    sum(  A0[J,J]  /2 + diag(A[J,J]))
    })
  
  AA = A>0
  comm_connections =  sapply(commus,  function(J)  {
    sum(  AA[J,J]  /2 + diag(AA[J,J]))
  })
  
  if(bysize==F){
    scores = comm_str/ comm_connections
  }else{
    scores = 1/comm_connections
  }
#  scores = comm_str/ sapply(commus, length)
  
  
  jaccard_mat0 <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      jaccard_mat0[i, j] <- length( intersect(commus[[i]], commus[[j]])) / length(commus[[i]])
    }
  }
  
  jaccard_mat <- jaccard_mat0
  jaccard_mat[is.nan(jaccard_mat)] <- 0
  diag(jaccard_mat) <- 0 # delete this bc no self loops??
  max_jacc <- max(jaccard_mat)
  deleted_comms <- integer(0)
  
  while (max_jacc > tau) {
    inds <- which(jaccard_mat == max_jacc, arr.ind = TRUE)[1, ]
    # keep comm with larger score
    delete_comm <- inds[which.min(c(scores[inds[1]], scores[inds[2]]))]
    jaccard_mat[delete_comm, ] <- 0
    jaccard_mat[, delete_comm] <- 0
    deleted_comms <- c(deleted_comms, delete_comm)
    max_jacc <- max(jaccard_mat)
  }
  
  kept_comms <- setdiff(1:K, deleted_comms)
  return(list("final_comms" = commus[kept_comms], "kept_comms" = kept_comms))
  
}



Comm_detec = Comm_detect(  E,  alpha=.01,  init_nodes=init$node1  , key_j = key  ) # this one 
 
  
