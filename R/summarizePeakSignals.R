summarizePeakSignals = function(input, suffix=NULL, chrom=c(1:22,'X','Y'), output=NULL) {
  peak = lapply(input, function(x) fread(x, select=c(1,2,3,7,8),
                                         col.names=c('Chr','Start','End','signalValue','p')))
  peak = lapply(peak, function(x) x[p.adjust(10^(-p), method='fdr') < 0.05,.(Chr,Start,End,signalValue)]) # FDR control p-value
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
