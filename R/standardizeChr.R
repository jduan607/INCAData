standardizeChr = function(chrom) {
    chrom = paste0('chr', gsub('^chr|chrom|chromosome', '', chrom, ignore.case=TRUE))
    return(chrom)
}
