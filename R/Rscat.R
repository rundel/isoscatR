MCMC_Chains <- function(geneofile, locfile, boundfile, options=list(), 
                        nChains=2, nIter=1000, nThin=50, nBurn=500, adjbound=TRUE, 
                        cv_indivs=c(), cv_locs=c()) {
    
    res = foreach(i = 1:nChains) %do% {
        
        options["FILESUFFIX"] = paste("_",i,sep="")
        
        z= MCMC_Main(geneofile, locfile, boundfile, options, 
                     nIter, nThin, nBurn, adjbound, 
                     cv_indivs, cv_locs, i)
        
        return(z)
    }
    
    
    
    l1 = length(res)
    l2 = length(res[[1]])
    
    out = list()
    for(i in 1:l2) {
        
        mcl = list()
        for(j in 1:l1) {
            mcl[[j]] = res[[j]][[i]]
        }
        out[[i]] = mcmc.list(mcl)
    }
     
    return( out )
}


MCMC_Main <- function(geneofile, locfile, boundfile, options=list(), 
                      nIter=1000, nThin=50, nBurn=500, adjbound=TRUE, 
                      cv_indivs=c(), cv_locs=c(), chain=1) {
    
    run_opts = scat_options
    for(n in names(options)) {
        run_opts[n] = options[n]
    }
    
    locdata = read_location_file(locfile)
    geneodata = read_geneotype_file(geneofile)
    bound = read_boundary_file(boundfile)
    
    locs = as.matrix(locdata[,5:6])
    #locnames = locdata[,1]
    locindex = locdata[,2]
    
    geneo = geneodata[[1]]
    indivlocs = geneodata[[2]]
    nAllele = geneodata[[3]]
    #indivID = geneodata[[4]]
    
    stopifnot(all(is.numeric(bound)))
    stopifnot(all(is.numeric(locs)))
    stopifnot(all(is.numeric(nAllele)))
    stopifnot(all(is.numeric(geneo)))
    #stopifnot(all(is.character(locnames)))
    #stopifnot(all(is.character(indivID)))
    
    if (adjbound) {
        bound = adjust_boundary(bound, locs)
    }
    
    nIndivs = nrow(geneo)
    nLocs = nrow(locs)
    if (length(cv_locs) != 0) {
        cv_locs = unique(cv_locs)
        stopifnot(all(cv_locs %in% 1:nLocs))
        
        cv_indivs = c(cv_indivs, which(indivlocs %in% cv_locs) )
    }
    
    cv_geneo = numeric()
    if (length(cv_indivs) != 0) {
        cv_indivs = sort(unique(cv_indivs))
        stopifnot(all(cv_indivs %in% 1:nIndivs))
        
        rows = sort(c(2*cv_indivs-1,2*cv_indivs))
        cv_geneo = geneo[rows,]
        geneo = geneo[-rows,]
        indivlocs = indivlocs[-cv_indivs]
        
        locs = locs[which( (1:nLocs %in% unique(indivlocs)) ),]
        
        i=1
        newid = rep(NA,length(indivlocs))
        for(l in sort(unique(indivlocs))) {
            newid[indivlocs == l] = i
            i=i+1
        }
        indivlocs = newid
        
        run_opts["CROSSVALIDATE"] = TRUE
        run_opts["LOCATE"] = FALSE
    }
    
    if (run_opts["CROSSVALIDATE"] == TRUE | run_opts["LOCATE"] == TRUE) {
        # if we are locating we need to create a directory for the tmp data
        if (!file.exists(run_opts[["TMPDIR"]]))
            dir.create(run_opts[["TMPDIR"]], recursive=TRUE)
    }

    
    if (run_opts$PERMUTE) {
        n = nrow(locs)
        perm = sample(1:n,n)
        
        locs = locs[perm,]
        #locnames = locnames[perm]
        locindex = locindex[perm]
        
        newid = rep(NA,length(indivlocs))
        for(i in 1:n) {
            newid[ indivlocs==perm[i] ] = i 
        }
        indivlocs = newid
    }
    
    z= .Call("mcmc_main",
              chain,
              bound,
              locs,
              geneo,
              indivlocs-1,
              nAllele,
              nIter,
              nThin,
              nBurn,
              cv_indivs,
              cv_geneo,
              run_opts,
              PACKAGE = "Rscat" )
  
    mcmc_list = lapply(z,function(z) {
            
                    d = z[["values"]]
                    colnames(d) = z[["names"]]
            
                    mcmc(data = d, start = (nBurn+1), end = (nBurn+nIter), thin=1)
                })
            
        
    return( mcmc_list )
    
}