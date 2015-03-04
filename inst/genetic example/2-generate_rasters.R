library(scatR)
library(raster)
library(ROCR)
library(stringr)
library(sp)

sp = "WIWA"

root = "./"
output_dir = file.path(root,"Scratch")

boundary_file = file.path(root, "Data/birdboundary.txt")
allele_file   = file.path(root, "Data",paste(sp,"_all.txt",sep=""))
location_file = file.path(root, "Data",paste(sp,"_loc.txt",sep=""))

geno = read_genotype_file(allele_file)
al = geno[[1]]
indloc = geno[[2]]
ind_ids = geno[[4]]
nind = length(indloc)

loc = read.table(location_file)
locs = loc[,c(4,3)]

r = allele_raster_from_file(output_dir)

bound = read_boundary_file(boundary_file)*180/pi
bound = SpatialPolygons(list(Polygons(list(Polygon(bound,hole=FALSE)),ID="boundary")))
bound = rasterize(bound,r)

g = list()
for(c in 1:nChains)
    g[[c]] = list()

pb = txtProgressBar(max = nind, style = 3)

for(ind in 1:nind)
{
    for(c in 1:nChains)
    {
        r = mask(allele_raster_from_file(root, ind, c), bound)
        g[[c]][[ind]] = list(r, cellFromXY(r,locs[indloc[ind],]))

    }

    setTxtProgressBar(pb, ind)
}
close(pb)

save(g, file=file.path(output_dir,"genetic_rasts.Rdata"))
