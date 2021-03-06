MCMC_Chains <- function(genofile, locfile, boundfile, options=list(),
                        nChains=2, nIter=1000, nThin=50, nBurn=500, adjbound=TRUE,
                        locate_indivs=c(), cv_locs=c()) {

    #res = foreach(i = 1:nChains) %do% {
    #
    #    options["FILESUFFIX"] = paste("_",i,sep="")
    #
    #    z= MCMC_Main(genofile, locfile, boundfile, options,
    #                 nIter, nThin, nBurn, adjbound,
    #                 locate_indivs, cv_locs, i)
    #
    #    return(z)
    #}

    res = list()
    for(i in 1:nChains) {

        options["FILESUFFIX"] = paste("_",i,sep="")

        z= MCMC_Main(genofile, locfile, boundfile, options,
                     nIter, nThin, nBurn, adjbound,
                     locate_indivs, cv_locs, i)

        res[[i]] = z
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


MCMC_Main <- function(genofile, locfile, boundfile, options=list(),
                      nIter=1000, nThin=50, nBurn=500, adjbound=TRUE,
                      locate_indivs=c(), cv_locs=c(), chain=1) {

    run_opts = scat_options
    for(n in names(options)) {
        run_opts[n] = options[n]
    }

    locdata = read_location_file(locfile)
    genodata = read_genotype_file(genofile)
    bound = read_boundary_file(boundfile)

    locs = as.matrix(locdata[,5:6])

    geno = genodata[[1]]
    indivlocs = genodata[[2]]
    nAllele = genodata[[3]]
    #indivID = genodata[[4]]

    stopifnot(all(is.numeric(bound)))
    stopifnot(all(is.numeric(locs)))
    stopifnot(all(is.numeric(nAllele)))
    stopifnot(all(is.numeric(geno)))
    #stopifnot(all(is.character(locnames)))
    #stopifnot(all(is.character(indivID)))

    if (adjbound) {
        bound = adjust_boundary(bound, locs)
    }

    nIndivs = nrow(geno)/2
    nLocs = nrow(locs)
    if (length(cv_locs) != 0) {
        cv_locs = unique(cv_locs)
        stopifnot(all(cv_locs %in% 1:nLocs))

        locate_indivs = c(locate_indivs, which(indivlocs %in% cv_locs) )
    }

    locate_geno = numeric()
    if (length(locate_indivs) != 0) {
        locate_indivs = sort(unique(locate_indivs))
        stopifnot(all(locate_indivs %in% 1:nIndivs))

        rows = sort(c(2*locate_indivs-1,2*locate_indivs))
        locate_geno = geno[rows,]
        geno = geno[-rows,]
        indivlocs = indivlocs[-locate_indivs]

        # if we remove all indivs from a loc drop that loc
        locs = locs[which( (1:nLocs %in% unique(indivlocs)) ),]

        i=1
        newid = rep(NA,length(indivlocs))
        for(l in sort(unique(indivlocs))) {
            newid[indivlocs == l] = i
            i=i+1
        }
        indivlocs = newid
    }

    if (!file.exists(run_opts[["TMPDIR"]]))
        dir.create(run_opts[["TMPDIR"]], recursive=TRUE)

    if (run_opts[["LOCATE"]]) {
        if (length(locate_indivs) == 0) {
            locate_indivs = 1:nIndivs
            locate_geno = geno
        }
    }


    if (run_opts[["PERMUTE"]]) {
        n = nrow(locs)
        perm = sample(1:n,n)

        locs = locs[perm,]

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
              geno,
              indivlocs-1, # change to 0 indexing
              nAllele,
              nIter,
              nThin,
              nBurn,
              locate_indivs,
              locate_geno,
              run_opts,
              PACKAGE = "scatR" )

    mcmc_list = lapply(z,function(z) {

                    d = z[["values"]]
                    colnames(d) = z[["names"]]

                    mcmc(data = d, start = (nBurn+1), end = (nBurn+nIter), thin=1)
                })


    return( mcmc_list )

}
