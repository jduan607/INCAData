standardizeChr = function(chrom) {
    chrom = paste0('chr', gsub('^chr|chrom|chromosome', '', chrom, ignore.case=TRUE))
    return(chrom)
}


summarizeReadCounts = function(input, scale=1e+7, chrom=c(1:22,'X','Y'), output=NULL) {
    bam = scanBam(BamFile(input), param=ScanBamParam(what=c('rname','pos','qwidth')))[[1]]
    gr = GRanges(seqnames=as.character(bam$rname), ranges=IRanges(start=bam$pos, width=bam$qwidth))
    counts = coverage(gr) # the number of reads that cover each position
    names(counts) = standardizeChr(names(counts)) # standardize chromosome names

    if (!is.null(chrom)) {
        counts = counts[standardizeChr(chrom)] # specific chromosomes
    }

    if (scale > 0) {
        counts = counts / length(gr) * scale
    } else if (scale < 0) {
        print('Invalid scaling factor for normalization. Use default scale = 1,000,000.')
        scale = 1e+7
        counts = counts / length(gr) * scale
    }

    counts_dt = lapply(names(counts), function(x) data.table(Chr = x, Length = counts[[x]]@lengths,
                                                             Value = counts[[x]]@values))
    counts_dt = rbindlist(counts_dt)
    colnames(counts_dt)[3] = ifelse(scale!=0, 'NormRC', 'RC')

    if (!is.null(output)) {
        fwrite(counts_dt, output, sep='\t')
    }
    return(counts_dt)
}


summarizePeakSignals = function(input, suffix=NULL, chrom=c(1:22,'X','Y'), output=NULL) {
    peak = lapply(input, function(x) fread(x, select=c(1,2,3,7,8),
                                                col.names=c('Chr','Start','End','signalValue','p')))
    peak = lapply(peak, function(x) x[p.adjust(10^(-p), method='fdr') < 0.05, ]) # FDR control p-value
    peak = lapply(peak, function(x) x[, Chr := standardizeChr(Chr)]) # standardize chromosome names

    if (!is.null(chrom)) {
        peak = lapply(peak, function(x) x[Chr %in% standardizeChr(chrom),]) # specific chromosomes
    }

    if (length(input) > 1) {
        if (length(suffix) < length(input)) {
            print('Invalid suffix.')
            suffix = 1:length(input)
        }
        gr = do.call(c, lapply(peak, function(x) GRanges(seqnames=x[,Chr], ranges=IRanges(start=x[,Start], end=x[,End]))))
        gr = disjoin(gr, with.revmap=TRUE) # split to non-overlapping regions

        signals = do.call(rbind, lapply(1:length(peak), function(i) cbind(peak[[i]][,signalValue], i)))
        signals = t(sapply(mcols(gr)[,1], function(x) { s=rep(0,length(peak)); s[signals[x,2]] = signals[x,1]; s }))

        signals_dt = cbind(as.data.table(gr)[,c(1:3)], signals)
        colnames(signals_dt) = c('Chr','Start','End', paste('signalValue', suffix, sep='_'))
    } else {
        signals_dt = peak[[1]]
        colnames(signals_dt) = c('Chr','Start','End', paste(c('signalValue', suffix), collapse='_'))
    }

    if (!is.null(output)) {
        fwrite(signals_dt, output, sep='\t')
    }
    return(signals_dt)
}
