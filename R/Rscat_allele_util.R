find_ind_file = function(root, ind, chain="") {
    d = dir(root, "Loc",full.names=TRUE)
    if (length(d) != 0) {
        file = dir(d, pattern = paste("CV_Ind",ind,"_",chain,sep=""),full.names=TRUE)
    } else {
        d = dir(root,paste("Ind",ind,sep=""),full.names=TRUE)
        file = dir(d, pattern = paste("CV_Ind",ind,"_",chain,sep=""), full.names=TRUE)
    }
    
    return( file )
}

calc_step = function(locs) {

    if (min(locs) > -2*pi & max(locs) < 2*pi)
        locs=locs*180/pi

    step = min(dist(locs))

    if (step < 1)
        step = 1/round(1/step,0)

    return(step)
}

create_raster = function(locs,step) {

    require(raster)

    if (min(locs) > -2*pi & max(locs) < 2*pi)
        locs=locs*180/pi

    xr = range(locs[,1]) + c(-1,1)*step/2
    yr = range(locs[,2]) + c(-1,1)*step/2

    nr = (yr[2]-yr[1])/step
    nc = (xr[2]-xr[1])/step

    r = raster(nrow=nr,ncol=nc,xmn=xr[1], xmx=xr[2], ymn=yr[1], ymx=yr[2])

    return( r )
}

allele_raster_from_file = function(root, ind = NULL, chain = 1) {
    
    if (is.null(ind))
        ind=".*"
    
    file = find_ind_file(root, ind, chain)
    
    if (length(file) == 0) {
        stop("Error unable to locate allele file in: ", paste(head(d),collapse=" "), " ...")
    }    
    if (length(file) != 1) {
        if (ind == ".*")
            file = file[1]
        else
            stop("Error multiple allele files found: ",file)
    }
    
    file_dir = dirname(file)
    pred_file = dir(file_dir, pattern=paste("pred_coords[0-9]+_",chain,".mat",sep=""), full.names=TRUE)
    stopifnot(length(pred_file)==1)
    
    pred_locs = matrix(scan(pred_file),ncol=2,byrow=TRUE)*180/pi
    
    step = calc_step(pred_locs)
    r = create_raster(pred_locs,step)
    r[]=NA
    
    cells = cellFromXY(r,pred_locs)
    stopifnot(all(!is.na(cells)))

    z = read_allele_file(file, nr=nrow(pred_locs))
    r[cells] = exp(apply(z,1,median))
    
    return(r)
}

read_allele_file = function(file, nr, nc, byrow=FALSE) {
    
    stopifnot(!missing(file))
    stopifnot(is.character(file))
    stopifnot(length(file)==1)
    
    stopifnot(file.exists(file))
        
    vec = .Call("read_allele_file", file, PACKAGE = "Rscat" )

    if (missing(nr) & missing(nc)) {
        return(vec)
    } else if (missing(nr)) {
        nr = length(vec) / nc
    } else if (missing(nc)) {
        nc = length(vec) / nr
    }
    
    stopifnot(nc == floor(nc))
    stopifnot(nr == floor(nr))
    
    res = matrix(vec, nrow=nr, ncol=nc, byrow=byrow)
    
    return(res)
}


calc_model_allele_freq = function(result_dir, nChains, nIter, loc_file, sep="\t", probs = 0.5) {
    
    stopifnot(file.exists(result_dir))
    
    al_files = dir(result_dir, pattern = "Al[0-9]+-[0-9]+_[0-9]+\\.gz")
    stopifnot(length(al_files)!=0)

    pred_file = file.path(result_dir,"pred_coords.mat")
    stopifnot(file.exists(pred_file))
    
    pred_locs = matrix(scan(pred_file),ncol=2,byrow=TRUE)*180/pi
    step = calc_step(pred_locs)
    r = create_raster(pred_locs,step)
    r[]=NA
    
    cells = cellFromXY(r,pred_locs)
    stopifnot(all(!is.na(cells)))
    
    
    locs = read_location_file(loc_file,sep)[,4:3]
    samp_cells = cellFromXY(r,locs)
    
    res = list()
    for(c in 1:nChains) {
        res[[c]] = list()
        
        l=1
        repeat {
            sub = grepl(paste("Al",l,"-[0-9]+_",c,"\\.gz",sep=""),al_files)
            
            if (sum(sub)==0)
                break
            
            res[[c]][[l]] = array(NA,c(sum(sub),nrow(locs),length(probs)))
            for(a in 1:sum(sub)) {
                file = paste("Al",l,"-",a,"_",c,".gz",sep="")
                print(file)
                stopifnot(file.exists(file.path(result_dir,file)))
                
                z = read_allele_file(file.path(result_dir,file), nr=nrow(pred_locs),nc=nIter)
                
                quants = t( matrix(apply(z,1,function(x) quantile(x,probs)),nrow=length(probs)) )
                
                loc_sub = sapply(samp_cells, function(x) which(x==cells))
                
                res[[c]][[l]][a,,] = quants[loc_sub,]
            }
            
            l=l+1
        }
    }

    return(res)
}

calc_allele_freq = function(allele_file, sep="") {
    
    data = read.table(allele_file,sep=sep,stringsAsFactors=FALSE)
    
    n_loci = ncol(data)-2
    loci = 1:n_loci
    
    locs = unique(data[,2])
    
    alleles = data[,loci+2]
        
    res = list()
    for(l in loci) {
        vals = alleles[,l]
        u = sort(unique(vals[vals!=-999]))

        allele_names =  as.character(u)
        
        res[[l]] = matrix(NA,nrow=length(u),ncol=length(locs))
        rownames(res[[l]]) = u
        colnames(res[[l]]) = paste("Loc",locs)
        
        for(loc in locs) {
            sub = data[,2] == loc
            res[[l]][,loc] = sapply(u, function(x) sum(vals[sub]==x)/sum(vals[sub]!=-999))
        }
    }

    return(res)
}