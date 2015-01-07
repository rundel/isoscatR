emp_variog = function(genofile, locfile,transform) {

    if (missing(transform))
        transform = function(x) x

    locdata = read_location_file(locfile)
    genodata = read_genotype_file(genofile, nrow(locdata))

    locs = as.matrix(locdata[,c(5,6)])
    locnames = locdata[,1]

    geno = as.matrix( genodata[[1]] )
    nAllele = genodata[[2]]
    indivID = genodata[[3]]


    res = list()

    for(i in 1:length(nAllele)) {
        res[[i]] = matrix(0, nrow=length(locnames), ncol=nAllele[i])

        for(j in 1:length(locnames)) {
            row_sub = geno[,1]==j
            for(k in which(row_sub)) {
                if (geno[k,i+1] != -999)
                    res[[i]][j,geno[k,i+1]+1] = res[[i]][j,geno[k,i+1]+1]+1
            }
            res[[i]][j,] = transform( res[[i]][j,] / (sum(geno[row_sub,i+1] != -999)) )
        }
    }

    dist = distance_mat(locs)
    ndists = choose(nrow(dist),2)

    variog = list(dist = rep(NA,ndists), data=matrix(NA,length(nAllele),ndists))

    l=1
    for(i in 2:length(locnames)) {
        for(j in 1:(i-1)) {

            variog[["dist"]][l] = dist[i,j]
            for(k in 1:length(nAllele)) {
                variog[["data"]][k,l] = (1/nAllele[k]) * sum( (res[[k]][i,]-res[[k]][j,])^2 )
            }
            l=l+1
        }
    }

    return(variog)
}


bin_variog = function(variog,breaks,step=500,add=FALSE,...) {

    if (missing(breaks)) {
        breaks = seq(min(variog[["dist"]])+step,max(variog[["dist"]])+step,step)
    }

    data = variog[["data"]]
    dist = variog[["dist"]]

    res = rep(NA,length(breaks))
    for(i in 1:length(breaks)) {
        sub = dist < breaks[i]
        if (sum(sub) != 0)
            res[i] = mean(data[,sub])
        else
            res[i] = NA

        data = data[,!sub,drop=FALSE]
        dist = dist[!sub]
    }

    if (add)
        points(breaks,res,...)
    else
        plot(breaks,res,ylim=c(0,max(res*(1.1),na.rm=TRUE)),...)

    print(res)
}


calc_cov_matern = function(dist, sigma2, phi, nu, tau2,uselog=FALSE) {

    stopifnot(is.numeric(dist))
    stopifnot(is.numeric(nu))
    stopifnot(is.numeric(phi))
    stopifnot(is.numeric(tau2))
    stopifnot(is.numeric(sigma2))
    stopifnot(is.logical(uselog))

    distmat = !is.null(dim(dist))

    res = .Call("R_cov_matern", sigma2, nu, phi, tau2, dist, distmat, uselog, PACKAGE = "scatR");

    return(res)
}

calc_cov_powered_exp = function(dist,sigma2,phi,kappa,tau2) {
    return ( sigma2 * exp(-(dist/phi)^kappa) + tau2*(dist==0) )
}

plot_cov = function(params, d, cov_func, CI=95, variog = FALSE) {

    if (inherits(params,"mcmc.list")) {
        nchains = length(params)
    } else {
        nchains = 1
        params = list(params)
    }


    ymax = round(2*sum(summary(params)$quantiles[,3]),1)
    ymax = 4
    plot(0,0,type='n',xlab="Distance",ylab="cov",xlim=range(d),ylim=c(0,ymax))#,main=paste("covariance",suff,sep=""))

    up = rep(0,length(d))
    low = rep(0,length(d))
    up_vario = rep(0,length(d))
    low_vario = rep(0,length(d))

    for(i in 1:nchains) {

        ma = apply(params[[i]],2,median)
        med_cov = cov_func(d,ma[1],ma[2],ma[3],ma[4] )

        ma = apply(params[[i]],2,mean)
        mean_cov = cov_func(d,ma[1],ma[2],ma[3],ma[4])

        pos_up = nrow(params[[i]])*((100-CI)/100/2)
        pos_low = nrow(params[[i]])-pos_up

        for(j in 1:length(d)) {
            cov =  sort( cov_func(d[j],params[[i]][,1], params[[i]][,2], params[[i]][,3], params[[i]][,4]) )

            up[j] = cov[pos_up]
            low[j] = cov[pos_low]

            up_vario[j] = up[1]-cov[pos_low]
            low_vario[j] = low[1]-cov[pos_up]
        }

        if (variog) {
            med_vario = med_cov[1]-med_cov

            lines(d,med_vario,col=i)

            lines(d,up_vario,col=i,lty=2)
            lines(d,low_vario,col=i,lty=2)
        } else {
            lines(d,med_cov,col=i)

            lines(d,up,col=i,lty=2)
            lines(d,low,col=i,lty=2)
        }
    }
}
