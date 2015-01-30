library(spBayes)

load(file = file.path(data_dir,"gnip_data.Rdata"))

# This will be a very slow step, depends on n.samples
# Model based IsoMAP proposed model - probably worthwhile to revisit

n.samples = 5000
l = spLM(
    H2 ~ alt+cru_tmx+I(cru_tmx^2)+cru_tmn+I(cru_tmn^2)+precip+temp+I(temp^2)+lat+I(lat^2)+long+I(long^2),
    data = gnip_mean,
    coords = as.matrix(gnip_mean[,c("long","lat")]),
    starting = list( nu = 1,
                     phi = 1, 
                     sigma.sq = 0.08,
                     tau.sq = 0.02), 
    sp.tuning = list( nu = 0.2,
                      phi = 5,
                      sigma.sq = 0.05, 
                      tau.sq = 0.05), 
    priors = list( nu.Unif = c(0.5,1.5),
                   phi.Unif = c(0.01, 5*max.dist),
                   sigma.sq.IG = c(2, 0.08), 
                   tau.sq.IG = c(2, 0.5)), 
    cov.model = "matern",
    n.samples = n.samples,
    sub.samples = c(1000, n.samples, 10),
    verbose = TRUE,
    n.report = 100
)

save(l, file = file.path(data_dir,"splm_fit.Rdata"))