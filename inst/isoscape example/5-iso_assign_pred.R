library(raster)
library(scatR)

load(r,file="sp_lm_pred.Rdata")

results_dir = "Results/"

sp = "WIWA"

for(chain in 1:2) {

    out_dir = file.path(results_dir,"out")
    if(!file.exists(out_dir))
        dir.create(out_dir,recursive=TRUE)

    gen_r = create_raster( read.table(file.path(results_dir, paste0("pred_coords_",chain,".mat"))))

    iso_r = crop(r,gen_r)
    iso_r = aggregate(iso_r, round(res(gen_r)/res(iso_r)))

    if (!all(res(iso_r) == res(gen_r)))
        iso_r = raster::resample(iso_r,gen_r,"ngb")



    boundary_file = file.path(data_dir,"Data/birdboundary.txt")
    allele_file = file.path(data_dir,"Data",paste(sp,"_all.txt",sep=""))
    location_file = file.path(data_dir,"Data",paste(sp,"_loc.txt",sep=""))
    isotope_file = file.path(data_dir,"Data",paste(sp,"_isotope.csv",sep=""))

    all = read.table(allele_file)[,1:2]
    all = all[seq(1,nrow(all),2),]
    iso = read.table(isotope_file,sep=",")#[,-2]
    loc = read.table(location_file)

    locs = loc[,c(4,3)]

    data = merge(merge(all,loc,by.x=2,by.y=2)[,-1],iso,by.x=1,by.y=1)

    kfd = data[,-c(1,2)]
    names(kfd) = c("lat","long","d")
    kfd$loc = data[,2]

    for(i in 1:nrow(kfd)) {
        cat("Ind",i,"/",nrow(kfd),"\n")
        b = calc_feather_brick(kfd$d[i], kfd$d, kfd[,c("long","lat")], iso_r,progressbar=TRUE)
        writeRaster(b,filename=file.path(out_dir,paste("Ind",i,".grd",sep="")),overwrite=TRUE)

        #cat("CV - Ind",i,"/",nrow(kfd),"\n")
        b = calc_feather_brick(kfd$d[i], kfd$d[-i], kfd[-i,c("long","lat")], iso_r,progressbar=TRUE)
        writeRaster(b,filename=file.path(out_dir,paste("Ind",i,"_cv.grd",sep="")),overwrite=TRUE)

        #cat("CV Loc - Ind",i,"/",nrow(kfd),"\n")
        loc_sub = (kfd$loc[i] == kfd$loc)
        b = calc_feather_brick(kfd$d[i], kfd$d[!loc_sub], kfd[!loc_sub,c("long","lat")], iso_r,progressbar=TRUE)
        writeRaster(b,filename=file.path(out_dir,paste("Ind",i,"_cv_loc.grd",sep="")),overwrite=TRUE)
    }
}
