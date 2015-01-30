###
### Preprocessing of data files to use in isoscape model
###
### Note some data versions may be out of date (newer versions available with different naming schemes)
### and different column names
###

library(raster)
library(lubridate)
library(maptools)


####
#### ETOPO Data Prep
####
#### ETOPO Data - http://www.ngdc.noaa.gov/mgg/global/global.html
#### gshhs Data - http://www.soest.hawaii.edu/pwessel/gshhg/
####

r = raster(file.path(data_dir,"etopo/ETOPO1_Ice_c_geotiff.tif"))

r = merge(r,shift(r,x=360))

res = c("f", "h", "i", "l", "c")[3]
s = Rgshhs(file.path(data_dir, paste("GSHHS21/gshhs_",res,".b",sep="")), level=1)[[3]]
proj4string(s) = CRS("")

rs = rasterize(s,r,mask=TRUE)

rs_a = aggregate(rs,fact=30)
cells = which( !is.na(rs_a[]) )

x = xFromCell(rs_a,cells)
y = yFromCell(rs_a,cells)

e = extent( floor(min(x)), ceiling(max(x)), floor(min(y)), ceiling(max(x)) )

rs_a = crop(rs_a,e)
gc()

writeRaster(crop(rs,e),file.path(data_dir,"etopo/etopo1_ice_land.grd"),overwrite=TRUE)
gc()

r1=crop(rs_a,extent(floor(min(x)),180,-90,90))
r2=crop(rs_a,extent(180,ceiling(max(x)),-90,90))
r2=shift(r2,x=-360)
r3=mosaic(r1,r2,fun=mean)

writeRaster(r3,file.path(data_dir,"etopo/etopo1_ice_land_0.5.grd"),overwrite=TRUE)


####
#### GNIP Data Prep
####
#### GNIP data - http://www-naweb.iaea.org/napc/ih/IHS_resources_isohis.html
####

gnip = read.table(file.path(data_dir,"GNIP_monthly.txt"),sep="\t",comment.char="",header=TRUE, stringsAsFactors=FALSE)
gnip = gnip[!is.na(gnip$H2) & gnip$month %in% months,]


####
#### CRU Data Prep
####
#### CRU TS Data (Wiser) - http://catalogue.ceda.ac.uk/uuid/3f8944800cc48e1cbc29a5ee12d8542d
####



types = c("pre","tmn","tmp","tmx")
cru_years = 1901:2009

for(type in types)
{
    in_file  = file.path(data_dir,"cru",paste0("cru_ts_3_10.1901.2009.",type,".dat.gz"))
    out_file = file.path(data_dir,"cru",paste0("cru_ts_3_10.1901.2009.",type,".Rdata"))

    x = scan(in_file,na.strings="-999")
    gc()

    x=array(x, c(720,360, 12, max(cru_years)-min(cru_years)+1))

    r = list()
    for(year in cru_years) {
        r[[year]] = brick(x[,360:1,,year-min(cru_years)+1],-180,180,-90,90,transpose=TRUE)
    }

    save(r,file=out_file)

    rm(x)
    rm(r)
    gc()
}


pred_rasts = list()

cru_scaling = function(type) {
    scaling = list(cld = 10,
                   dtr = 10,
                   frs = 100,
                   pre = 10,
                   tmp = 10,
                   tmn = 10,
                   tmx = 10,
                   vap = 10,
                   wet = 100,
                   pet = 10)

    return(scaling[[type]])
}


cru_units = function(type) {
    units = list(cld = "percentage",
                 dtr = "deg C",
                 frs = "days",
                 pre = "mm",
                 tmp = "deg C",
                 tmn = "deg C",
                 tmx = "deg C",
                 vap = "hecta-Pascals",
                 wet = "days",
                 pet = "mm" )

    return(units[[type]])
}

# for(type in c("pre","tmn","tmp","tmx")) {
#     file_name = paste0("cru_ts_3_10.1901.2009.",type,".Rdata")
#     load(file.path(cru_dir,"cru",file_name))

#     gnip$loc_index = cellFromXY(r[[min(years)]],gnip[,c("long","lat")])

#     v = apply(gnip[,c("years","month","loc_index")], 1,
#               function(x) {
#                   years = eval(parse(text=gsub("-", ":", x[1])))
#                   vals = sapply(years, function(y,index,month) r[[ y ]][ index ][ month-4 ], index=as.numeric(x[3]), month=as.numeric(x[2]))
#                   return( mean(vals,na.rm=TRUE))
#               })
#     gnip[[paste("cru_",type,sep="")]] = v / cru_scaling(type)

#     pred_rasts[[type]] = lapply(years, function(i) r[[i]] / cru_scaling(type))
# }

# r = raster(file.path(data_dir,"rasters/ETOPO1_Ice_0.5deg.grd")
# i = cellFromXY(r,gnip[,c("long","lat")])
# gnip$etopo_alt = r[i]
# pred_rasts[["etopo"]] = r

# r = raster(file.path(data_dir,"rasters/etopo1_ice_land_0.5.grd")
# i = cellFromXY(r,gnip[,c("long","lat")])
# gnip$etopo_alt_land = r[i]
# pred_rasts[["etopo_land"]] = r


# ####
# #### Collapse data across years and months
# ####

# model_params = c("H2", "alt", "cru_tmx", "cru_tmn", "precip", "cru_tmp", "lat", "long")
# gnip_mean = aggregate(gnip[,model_params],by=list(gnip[,"site"]),mean,na.rm=TRUE)
# model_sub = apply(gnip_mean[,model_params],1,function(x) !any(is.na(x)))

# gnip_mean = gnip_mean[model_sub,]

# save(gnip, gnip_mean, pred_rasts, file="data_monthly.Rdata")
