raster_fast_aggregate = function(rast, fact=2) {
    if (length(fact)==1) {
        fact = as.integer(round(fact))
        if (fact < 2) { stop('fact should be > 1') }
        rstep = cstep = fact
    } else if (length(fact)==2) {
        rstep = as.integer(round(fact[[1]]))
        cstep = as.integer(round(fact[[2]]))
        if (rstep < 2) { stop('fact[[1]] should be > 1') }
        if (cstep < 2) { stop('fact[[2]] should be > 1') }
    } else {
        stop('length(fact) should be 1 or 2')
    }

    nrow = dim(rast)[2]
    ncol = dim(rast)[1]
    nlayer = nlayers(rast)

    z = .Call("R_fast_aggregate", rast[], rstep, cstep, nrow, ncol, nlayer, PACKAGE = "scatR" )

    out = aggregate(brick(rast, values=FALSE), fact)
    out[] = z

    return(out)
}

create_allele_data = function(allele_file, location_file)
{
    all = read.table(allele_file)[,1:2]
    all = all[seq(1,nrow(all),2),]

    loc = read.table(location_file)

    return( merge(all,loc,by.x=2,by.y=2)[,-1] )
}

create_allelel_isotope_data = function(allele_file, location_file, isotope_file) {
    iso = read.table(isotope_file,sep=",")#[,-2]

    all = create_allele_data(allele_file, location_file)

    return( merge(all,iso,by.x=1,by.y=1) )
}


raster_area = function(rast, R = 6378.388) {

    cell_area = function(lat1, lat2, lon_delta) {
        2*pi * R^2 * abs( sin(lat1*pi/180)-sin(lat2*pi/180) ) * abs( lon_delta ) / 360
    }

    lon_res = xres(rast)
    lat_res = yres(rast)

    areas = sapply(1:nrow(rast), function(i){
        lats = yFromRow(rast,i) + c(-0.5,0.5)*lat_res
        return( sum(!is.na(rast[i,])) * cell_area(lats[1],lats[2], lon_res) )
    })

    return(sum(areas))
}



normalize_raster = function(rast) {
    rast[] = rast[] / sum(rast[], na.rm=TRUE)
    return(rast)
}



adjust_boundary = function(bound, locs, perc=0.2, poly=FALSE) {

    require(rgeos)

    xr = range(locs[,1])
    yr = range(locs[,2])

    bbx = xr + c(-1,1)*(xr[2]-xr[1])*perc
    bby = yr + c(-1,1)*(yr[2]-yr[1])*perc

    coords = cbind(x = c(bbx[1],bbx[2],bbx[2],bbx[1],bbx[1]),
                   y = c(bby[1],bby[1],bby[2],bby[2],bby[1]))

    bbox = SpatialPolygons(list(Polygons(list(Polygon(coords)),ID="box")))
    boundsp = SpatialPolygons(list(Polygons(list(Polygon(bound)),ID="box")))

    adj = gIntersection(bbox,boundsp)

    if (poly)
        return( adj )
    else
        return( adj@polygons[[1]]@Polygons[[1]]@coords )

}


calc_distance = function(x,y) {

    if(missing(y)) {
        stopifnot(ncol(x) == 2)
        y=x[,2]
        x=x[,1]
    }

    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))

    res = .Call("R_calc_distance", x, y, PACKAGE = "scatR" )

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

    res = .Call("R_great_circle_dist_point", px, py, x, y, PACKAGE = "scatR" )

    return(res)
}