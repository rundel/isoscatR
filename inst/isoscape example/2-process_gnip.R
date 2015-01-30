library(maptools)
library(lubridate)
library(raster)
data(wrld_simpl)

exclude = c("DIEGO GARCIA ISLAND (INDIAN O.)",
            "DESTRUCTION ISLAND (PACIFIC O.)",
            "WEATHERSHIP E (ATLANTIC O.)",
            "BERMUDA ISLAND (ATLANTIC O.)",
            "SAN JUAN (PUERTO RICO)",
            "WEATHERSHIP V (PACIFIC O.)",
            "MIDWAY IS. (PACIFIC O.)",
            "TAGUAC GUAM IS. (PACIFIC O.)",
            "WAKE ISLAND (PACIFIC O.)",
            "JOHNSTON ISLAND (PACIFIC O.)",
            "HAWAII (HILO PACIFIC O.)",
            "TRUK (EASTERN CAROLINE IS. PACIFIC O.)",
            "PONAPE CAROLINE ISLAND (PACIFIC O.)",
            "MAJURO (MARSHALL ISLAND PACIFIC O.)",
            "YAP (WESTERN CAROLINE IS. PACIFIC O.)",
            "CANTON ISLAND (PACIFIC O.)")

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


gnip = read.table("Data/GNIP_NA.txt",sep=";",comment.char="",header=TRUE, stringsAsFactors=FALSE)
gnip = gnip[!(gnip[,1] %in% exclude | is.na(gnip$H2)),]

gnip$sample_date = ymd(strtrim(gnip$sample_date,10))

gnip = gnip[,c("sample_site", "sample_date", "H2", "longitude", "latitude", "altitude", "precipitation", "airtemp")]
names(gnip) = c("site", "date", "H2", "long", "lat", "alt", "precip", "temp")

gnip$year = year(gnip$date)
gnip$month = month(gnip$date)



# Select only data for months we care about
gnip = gnip[gnip$month %in% months, ]

locs = cbind(gnip$long, gnip$lat)


pred_rasts = list()

for(type in c("pre","tmn","tmp","tmx")) {
    file_name = paste("cru_ts_3_10.1901.2009.",type,".Rdata",sep="")
    load(file.path(data_dir,"cru",file_name))

    gnip$loc_index = cellFromXY(r[[1901]],gnip[,c("long","lat")])

    v = apply(gnip[,c("year","month","loc_index")], 1, function(x) { r[[ x[1] ]][ x[3] ][ x[2] ] })
    gnip[[paste("cru_",type,sep="")]] = v / cru_scaling(type)

    pred_rasts[[type]] = r[[pred_year]] / cru_scaling(type)
}
gnip$loc_index = NULL


r = raster("Data/rasters/ETOPO1_Ice_0.5deg.grd")
i = cellFromXY(r,gnip[,c("long","lat")])
gnip$etopo_alt = r[i]
pred_rasts[["etopo"]] = r


# Replace missing precip / temp with cru data
gnip$precip[is.na(gnip$precip)] = gnip$cru_pre[is.na(gnip$precip)]
gnip$temp[is.na(gnip$temp)] = gnip$cru_tmp[is.na(gnip$temp)]

# Fix -2.2 deg C in OTTAWA in august
gnip$temp[487] = gnip$cru_tmp[487]


sub = !(names(gnip) %in% c("site","date","year","month","loc_index"))
gnip_mean = data.frame( sapply(gnip[,sub],function(x) tapply(x,gnip[,1],mean,na.rm=TRUE)) )


pred_data = data.frame(alt     = pred_rasts[["etopo"]][],
                       cru_tmx = pred_rasts[["tmx"]][][,1],
                       cru_tmn = pred_rasts[["tmn"]][][,1],
                       temp    = pred_rasts[["tmp"]][][,1],
                       precip  = pred_rasts[["pre"]][][,1],
                       lat     = yFromCell(pred_rasts[["tmx"]],1:ncell(pred_rasts[["tmx"]])),
                       long    = xFromCell(pred_rasts[["tmx"]],1:ncell(pred_rasts[["tmx"]]))

save(gnip, gnip_mean, pred_data, pred_rasts, file = file.path(data_dir,"gnip_data.Rdata"))

