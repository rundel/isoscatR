mask_precip_map = function(map.mask) {
    data(map.precip.dat,package="isoscape")
    
    map.precip.rast = raster(map.precip.dat)
    map.precip.rast[map.precip.rast[]==-999] = NA

    stopifnot( !is.null(map.mask) )
    
    if(is.character(map.mask)) { 
        map.mask = readShapeSpatial(mask)
    }

    if (class(map.mask) %in% c("Extent","matrix","SpatialLines","SpatialPoints","SpatialPolygons","SpatialPolygonsDataFrame")) {
        map.mask = rasterize(map.mask, map.precip.rast)
    }

    stopifnot(inherits(map.mask,"RasterLayer"))

    return( crop(mask(map.precip.rast, map.mask), map.mask) )    
}

mask_precip_with_boundary = function(boundary_file) {
    bound = read_boundary_file(boundary_file)
    r = range(bound)
    if (!(r[1] < -2*pi || r[2] > 2*pi))
        bound = bound*180/pi
    
    bound.poly = SpatialPolygons(list(Polygons(list(Polygon(bound)),ID=1)))

    return( mask_precip_map(bound.poly) )
}

plot_feather_model = function(data, poly_order=1, cv_ind=NULL) {
    
    if (!is.null(cv_ind))
        data = data[-cv_ind,]
    
    lm.res = lm(d2 ~ poly(map.d2,poly_order),data=data)
    
    plot(data$map.d2,data$d2,type='n',xlab="Precipitation Map d2", ylab="Feather d2")
    text(data$map.d2,data$d2,data$loc,cex=0.5)
    
    x=data$d2
    if(poly_order == 1)
        abline(lm.res)
}

calc_feather_map = function(ind, feather.data, poly_order=1, plot=FALSE, pog = map.precip.rast) {
    
    tdat = feather.data[-ind,]

    lm.res = lm(d2 ~ poly(map.d2,poly_order),data=tdat)
    
    new = data.frame(map.d2 = pog[!is.na(pog[])])
    lm.est.mean = rep(NA,length(pog[]))
    lm.est.mean[!is.na(pog[])] = predict(lm.res,newdata=new)
    lm.est.sigma = sd(lm.res$residuals)
    
    
    pog[] = dnorm(feather.data$d2[ind], lm.est.mean, lm.est.sigma)
    
    return(pog)
}
    




