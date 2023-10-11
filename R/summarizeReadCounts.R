summarizeReadCounts = function(input, scale=1e+7, chrom=c(1:22,'X','Y'), output=NULL) {
  bam = Rsamtools::scanBam(Rsamtools::BamFile(input), param=Rsamtools::ScanBamParam(what=c('rname','pos','qwidth')))[[1]]
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
