adjust_boundary = function(bound, locs, perc=0.2) {
    
    require(rgeos)
    
    xr = range(locs[,1])
    yr = range(locs[,2])
    
    bbx = xr + c(-1,1)*(xr[2]-xr[1])*perc
    bby = yr + c(-1,1)*(yr[2]-yr[1])*perc
    
    coords = cbind(x = c(bbx[1],bbx[2],bbx[2],bbx[1],bbx[1]),
                   y = c(bby[1],bby[1],bby[2],bby[2],bby[1]))

    bbox = SpatialPolygons(list(Polygons(list(Polygon(coords)),ID="box")))
    boundsp = SpatialPolygons(list(Polygons(list(Polygon(bound)),ID="box")))
    
    return( gIntersection(bbox,boundsp)@polygons[[1]]@Polygons[[1]]@coords )
     
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
    
    if (min(locs) > -2*pi & max(locs) < 2*pi)
        locs=locs*180/pi
    
    xr = range(locs[,1]) + c(-1,1)*step/2
    yr = range(locs[,2]) + c(-1,1)*step/2
    
    nr = (yr[2]-yr[1])/step
    nc = (xr[2]-xr[1])/step
    
    r = raster(nrow=nr,ncol=nc,xmn=xr[1], xmx=xr[2], ymn=yr[1], ymx=yr[2])
        
    return( r )
}

read_allelefile = function(file, nr, nc, byrow=FALSE) {
    
    stopifnot(!missing(file))
    stopifnot(is.character(file))
    stopifnot(length(file)==1)
    
    stopifnot(file.exists(file))
        
    vec = .Call("read_allelefile", file, PACKAGE = "Rscat" )

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

prec_sum = function(..., na.rm=TRUE) {
    
    vec = unlist(list(...))

    if (length(vec)==0)
        return(0)

    if (na.rm==FALSE & any(is.na(vec)))
        return(NA)
    
    if (na.rm==TRUE)
        vec = vec[!is.na(vec)]
        
    return (.Call("prec_sum",vec,PACKAGE="Rscat"))

}

calc_distance = function(x,y) {
    
    if(missing(y)) {
        stopifnot(ncol(x) == 2)
        y=x[,2]
        x=x[,1]
    }

    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    
    res = .Call("R_calc_distance", x, y, PACKAGE = "Rscat" )
    
    return(res)
}

calc_distance_to_point = function(px,py,x,y) {
    
    if(missing(y)) {
        stopifnot(ncol(x) == 2)
        y=x[,2]
        x=x[,1]
    }

    stopifnot(is.numeric(px))
    stopifnot(is.numeric(py))
    stopifnot(length(px)==1)
    stopifnot(length(py)==1)

    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    
    res = .Call("R_calc_distance_to_point", px, py, x, y, PACKAGE = "Rscat" )
    
    return(res)
}