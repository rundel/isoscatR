root_dir = path.expand("Scratch/")

if(!file.exists(root_dir)) 
    dir.create(root_dir, recursive=TRUE, showWarnings=FALSE)

options = list( TMPDIR = root_dir,
                RETURNFIT = FALSE,
                LOCATE = TRUE,
                USEMATERN = FALSE,
                VERBOSE = FALSE,

                NULLPROB = 0.01,
                DELTA = 0.05,

                MAXCELL = 3000, #1 deg
                #MAXCELL = 2000,  #2 deg

                ALPHAMIN = c(0,1,0.01,0),
                ALPHAMAX = c(10,20000,2,10),

                ALPHASD = c(0.2,0.2, 0.2,0.05),
                THETASD = 0.1,
                ETASD = 1,
                XISD = 0.1,
                BETASD = 3,
                BETARANGE = c(0.01,25),

                PERMUTE = TRUE,

                OUTPUTALFREQ = FALSE,
                GZIPOUTPUT = FALSE
            )




nBurn = 10 #1000
nIter = 10 #1000
nThin = 1 #100
nChains=2

pts = 1000
ac.lag = 500


sp = "WIWA"

boundary_file = "Data/birdboundary.txt"
allele_file = paste("Data/",sp,"_all.txt",sep="")
location_file = paste("Data/",sp,"_loc.txt",sep="")

options["OUTPUTALFREQ"] = TRUE

m=MCMC_Chains(allele_file,location_file,boundary_file,
              nChains=nChains,nBurn=nBurn,nIter=nIter,nThin=nThin,
              options=options)

save(m,nBurn, nIter, nThin, nChains,
     file=file.path(root_dir, paste0(sp,"_mcmc.Rdata"))

