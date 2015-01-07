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


read_genotype_file = function(file, sep="\t", allele_lookup = NULL) {

    data = read.table(file,sep=sep,stringsAsFactors=FALSE)
    odd_rows = seq(1,nrow(data),2)

    stopifnot(nrow(data) %% 2 == 0)
    stopifnot(all(data[odd_rows,1] == data[-odd_rows,1]))
    stopifnot(all(data[odd_rows,2] == data[-odd_rows,2]))

    if (!is.null(allele_lookup)) {
        stopifnot(length(allele_lookup) == ncol(data)-2)
    } else {
        allele_lookup = list()
    }


    nAlleles = rep(NA,ncol(data)-2)

    alleles = as.matrix(data[,-(1:2)])
    for (i in 1:ncol(alleles)) {

        row_sub = alleles[,i] > 0

        if (length(allele_lookup) >= i) {
            unq = allele_lookup[[i]]
        } else {
            unq = sort(unique( alleles[row_sub, i] ))
        }


        nAlleles[i] = length(unique( alleles[row_sub, i] ))
        allele_lookup[[i]] = unq

        if (sum(row_sub) == 0) next

        alleles[row_sub, i] = sapply(alleles[row_sub, i], function(x) {
                                                            z=which(x==unq)-1 #so indexes start from 0
                                                            if(length(z)==0) z=-1
                                                            return(z)
                                                          })
    }

    indivID = data[odd_rows,1]
    indivlocs = data[odd_rows,2]

    return(list(alleles, indivlocs, nAlleles, indivID, allele_lookup))
}
