library(spBayes)

load(file = file.path(data_dir,"gnip_data.Rdata"))
load(file = file.path(data_dir,"splm_fit.Rdata"))


p = spPredict(l, pred.coords = pred_data[,c("long", "lat")], pred.covars=pred_data, start=10020, end=30000, thin=20)



r = brick( pred_rasts[[1]], nl = ncol(p$y.pred))
values(r) = p$y.pred

save(p,r,file="sp_lm_pred.Rdata")
