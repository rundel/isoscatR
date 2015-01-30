calc_feather_brick = function(ind_d, obs_d, locs, brick, progressbar = FALSE) {

    map_d = brick[cellFromXY(brick,locs)]

    data = brick[]

    if(progressbar) pb = txtProgressBar(style=3,max=ncol(map_d))

    for(c in 1:ncol(map_d)) {
        l = lm(obs_d ~ map_d[,c])

        coefs = l$coefficients
        subset = !is.na(data[,c])

        m = rep(NA,nrow(data))
        m[subset] = coefs[1]+coefs[2]*data[subset,c]
        s = sd(l$residuals)


        data[,c] = dnorm(ind_d, m, s)

        if(progressbar) setTxtProgressBar(pb, c)
    }

    if(progressbar) close(pb)

    b=brick
    values(b) = data

    return(b)
}


calc_feather_map = function(d2, data, poly_order=1, pog = map.precip.rast) {


    l = lm(d2 ~ poly(map.d2,poly_order),data=data)

    lm.est.mean = rep(NA,length(pog[]))
    lm.est.mean[!is.na(pog[])] = predict(l,newdata=data.frame(map.d2 = pog[!is.na(pog[])]))

    lm.est.sigma = sd(l$residuals)


    pog[] = dnorm(d2, lm.est.mean, lm.est.sigma)

    return(pog)
}

mask_precip_map = function(map.mask) {
    data(map.precip.dat,package="scatR")

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


