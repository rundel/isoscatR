plot_allele_freq = function(allele_file, sep="", loci = NULL, return_af = FALSE, cols = rainbow(length(locs))) {
    
    data = read.table(allele_file,sep=sep,stringsAsFactors=FALSE)
    
    locs = unique(data[,2])
    
    n_loci = ncol(data)-2
    if (is.null(loci)) {
        loci = 1:n_loci
    }
    
    stopifnot(all(loci %in% 1:n_loci))
    
    res = calc_allele_freq(allele_file, sep)
    
    nc = round(sqrt(n_loci),0)
    nr = round(n_loci/nc,0)
    
    par(mfrow=c(nr,nc))
    for(l in loci) {
        n = nrow(res[[l]])
        
        
        plot(c(1,n),c(0,1),type='n',xlab="",ylab="Allele Freq",main=paste("Locus",l))
        for(loc in locs) {
            lines(1:n, res[[l]][,loc], col = cols[which(loc == locs)])
        }
        legend("topright", paste("Loc", locs),col=cols,lty=1,ncol=4,lwd=1,cex=1/2)
    }

    if (return_af)
        return(res)
}


generate_reports = function(m, dir="./", prefix = "", pts = 1000, ac.lag = 500, trace=TRUE, autocorr=TRUE, cov=TRUE, cov_func=calc_cov_powered_exp) {
    
    if (inherits(m,"mcmc.list"))
        m = list(m)
    
    stopifnot(all(sapply(m,function(x) inherits(x,"mcmc.list"))))
    
    thin = niter(m[[1]])/pts  
    w=lapply(m, function(x) window(x,thin=thin))

    if (trace) {
        file = paste(prefix,"mcmc_chains.pdf",sep="")
        pdf( file.path(dir,file) )
        for(i in 1:length(m)) {
            plot_mcmc(w[[i]],)
        }
        dev.off()
    }

    if (autocorr) {
        file = paste(prefix,"mcmc_autocorr.pdf",sep="")
        pdf( file.path(dir,file))
        for(i in 1:length(m)) {
            for(j in 1:length(w[[i]])) 
                autocorr.plot(w[[i]][[j]],ac.lag)
        }
        dev.off()
    }

    if (cov) {
        file = paste(prefix,"mcmc_model_cov.pdf",sep="")
        pdf( file.path(dir,file))
        plot_cov(w[[1]],seq(0,2000,5), cov_func)
        plot_cov(w[[1]],seq(0,2000,5), cov_func,variog=TRUE)
        dev.off()
    }
}


plot_mcmc = function(mcmc) {
    
    par(mfrow=c(4,2))
    for(i in 1:nvar(mcmc)) {
        plot_mcmc_var(mcmc,i)
    }

}

plot_mcmc_var = function(mcmc, var) {
    plot(mcmc[,var],density=FALSE,auto.layout=FALSE)
    title(varnames(mcmc)[var])
    plot_density(mcmc[,var])
}

plot_locs = function(boundary_file,location_file,plot_bound = FALSE, plot_map = TRUE, map_color="black", adjust=FALSE, ...) {

    poly = as.matrix(read.table(boundary_file,sep=" "))[,2:1]
    
    loc_data = read.table(location_file,sep="\t")
    locs = as.matrix(loc_data[,4:3])
    
    if (adjust)
        poly = adjust_boundary(poly, locs)

    if (plot_bound) type = "l"
    else type = "n"

    plot(poly, asp=1,type=type,main="",xlab="",ylab="",axes=FALSE)
    
    if (plot_map) {
        if (require(maptools)) {
            data(wrld_simpl)
            plot(wrld_simpl,add=T,lwd=0.5,col=map_color)
        } else {
            warning("maptools library is not available, unable to plot map.")
        }
    }

    loc_names = paste("(",loc_data[,2],") ",loc_data[,1],sep="")
    
    text(locs,labels=loc_names, ...)            
}



plot_density = function (x, show.obs = TRUE, show.med = TRUE, bwf, main = "", col, ylim, xlim, alpha=0.5, ...) {
    if (alpha < 1) {
        alpha = round(255*alpha,0)
    }
    
    if (missing(col))
        col = 1:length(x)
        
    stopifnot(length(col) == length(x))
    
    if (missing(bwf)) {
        bwf = function(x) {
            x = x[!is.na(as.vector(x))]
            return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
        }
    }
    
    y=as.matrix(x)
    if (max(abs(y - floor(y))) == 0) {
        hist(y, prob = TRUE, main = main, ...)
        return()
    }
    
    if (missing(ylim)) ylim = c(0, 0)
    if (missing(xlim)) xlim = c(0, 0)
    
    dens_list = list()
    for (i in 1:length(x)) {
        y = x[[i]]

        bw = bwf(y)
        width = 4 * bw

        scale = "open"
        if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
            if (min(y) >= 0 && min(y) < 2 * bw) {
                scale = "proportion"
                y = c(y, -y, 2 - y)
            }
        } else if (min(y) >= 0 && min(y) < 2 * bw) {
            scale = "positive"
            y = c(y, -y)
        } else {
            scale = "open"
        }
        if (width == 0) {
            dens = NA
        } else {
            dens = density(y, width = width)
        
            if (scale == "proportion") {
                dens$y = 3 * dens$y[dens$x >= 0 & dens$x <= 1]
                dens$x = dens$x[dens$x >= 0 & dens$x <= 1]
            } else if (scale == "positive") {
                dens$y = 2 * dens$y[dens$x >= 0]
                dens$x = dens$x[dens$x >= 0]
            }
         
            ylim = c(0, max(ylim[2],dens$y))
            xlim = c(min(xlim[1],dens$x), max(xlim[2],dens$x) )
        }
        dens_list[[i]] = dens
    }
    
    bws = paste(sapply(dens_list, function(x) { ifelse(class(x)=="density", formatC(x$bw,digits=2), "") }),collapse=", ")
    
    plot(0,0,type='n',ylab="",xlab = paste("N = ", niter(x), "  Bandwidths = (", bws, ")",sep=""),xlim=xlim,ylim=ylim,...)
    for (i in 1:length(x)) {
        
        #curcol = col2rgb(col[i],TRUE)
        #curcol[4] = alpha
        if (class(dens_list[[i]]) == "density") {
            lines(dens_list[[i]], col = col[i])
        
            if (show.obs) 
                lines(x[[i]], rep(ylim[2]/100, niter(x)), type = "h",col=col[i])
                
            if (show.med)
                abline(v=median(x[[i]]),col=col[i])
        }
    }
}