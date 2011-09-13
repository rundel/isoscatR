require(coda)
require(doMC)
registerDoMC()

scat_options = list( VERBOSE = FALSE,
                
                LOCATE = FALSE,
                MAXCELL = 2500,
                
                TMPDIR = "./",
                FILESUFFIX = "",
                
                ADAPT = FALSE,
                TUNEINTERVAL = 1000,
                CROSSVALIDATE = FALSE,
                
                RETURNFIT = FALSE,
                USEMATERN = FALSE,
                
                PSEUDOCOUNT = 0,
                
                FIXALPHA = c(FALSE,FALSE,FALSE,FALSE),
                ALPHA = c(1,1,1,0),
                ALPHAMIN = c(0.2,1,0.1,0),
                ALPHAMAX = c(100,100000,2,100),
                
                ALPHASD = c(0.3,0.3,0.25,0.3),
                
                ANGLE = 0,
                FIXANGLE = TRUE,
                ANGLESD = 0.1,
                
                RATIO = 1,
                FIXRATIO = TRUE,
                RATIOSD = 0.1,
                
                XIRANGE = c(-1000,1000),
                FIXXI = FALSE,
                XI = 1,
                XISD = 1,
                SIGMAXI = 10,
                
                FIXMU = FALSE,
                MU = 0,
                MUSD = 1.1,
                
                FIXETA = FALSE,
                ETA = 0,
                ETASD = 1.1,
                
                BETARANGE = c(0,25),
                FIXBETA = FALSE,
                BETA = 1,
                BETASD = 0.5,
                SIGMABETA = 5,
                
                XSD = 2.1,
                
                LOCALSD = 0.225,
                GLOBALSD = 0.225,
                
                NULLPROB = 0,
                DELTA = 0,
                
                PERMUTE = TRUE )