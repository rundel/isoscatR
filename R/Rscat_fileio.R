read_boundary_file = function(file,sep=" ") {
    data = as.matrix(read.table(file, sep=sep))
    
    stopifnot(all(is.numeric(data)))
    stopifnot(ncol(data)==2)
    
    stopifnot(all(data[,1] > -90 & data[,1] < 90))
    stopifnot(all(data[,2] > -180 & data[,2] < 180))
    
    l = nrow(data)
    
    if(!all(data[1,] == data[l,])) 
        data[l+1,] = data[1,]
    
    res = pi*data[,c(2,1)] / 180
    
    return(res)
}

read_location_file = function(file,sep="\t") {
    
    data = read.table(file,sep=sep,stringsAsFactors=FALSE)
    
    # Reorder the data if needed
    data = data[order(data[,2]),]
    
    stopifnot(ncol(data) == 4)
    stopifnot(all( 1:nrow(data) == data[,2] ))
    
    stopifnot(all(data[,3] > -90 & data[,3] < 90))
    stopifnot(all(data[,4] > -180 & data[,4] < 180))
    
    stopifnot( length(unique(data[,1])) == nrow(data)  )
    
    data[,5:6] = pi*data[,c(4,3)] / 180 # add new columns with pos in radians
    
    return(data)
}


read_geneotype_file = function(file, sep="\t") {
    
    data = read.table(file,sep=sep,stringsAsFactors=FALSE)
    odd_rows = seq(1,nrow(data),2) 
    
    stopifnot(nrow(data) %% 2 == 0)
    stopifnot(all(data[odd_rows,1] == data[-odd_rows,1]))
    stopifnot(all(data[odd_rows,2] == data[-odd_rows,2]))
    
    
    nAlleles = rep(NA,ncol(data)-2)
    
    allele_lookup = list()
    for (i in 3:ncol(data)) {
        unq = sort(unique( data[ data[,i] > 0, i] ))
        nAlleles[i-2] = length(unq)
        
        recoded = sapply(data[ data[,i] > 0, i], function(x) which(x==unq))-1 #so indexes start from 0
        data[ data[,i] > 0, i] = recoded
        
        allele_lookup[[i-2]] = unq
    }
    
    odd_rows = seq(1,nrow(data),2)
    indivID = data[odd_rows,1]
    indivlocs = data[odd_rows,2]
    
    alleles = as.matrix(data[,-(1:2)])
    
    return(list(alleles, indivlocs, nAlleles, indivID, allele_lookup))
}
