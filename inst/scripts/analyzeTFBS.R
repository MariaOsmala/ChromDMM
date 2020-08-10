library(GenomicRanges)

parse.narrowPeak <- function (peaks.file, summit=T) {
  #load human chromosome lengths
  library(BSgenome.Hsapiens.UCSC.hg19)

  human.chromlens = seqlengths(Hsapiens)

  colnames <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                "signalValue", "pValue","qValue", "peak")

  peaks <- read.table(peaks.file)
  names(peaks) <- colnames[seq(ncol(peaks))]

  if ('peak' %in% names(peaks) && summit) {
    ranges.start <- peaks$chromStart + peaks$peak + 1
    ranges.end   <- ranges.start
  } else {
    ranges.start <- peaks$chromStart + 1
    ranges.end   <- peaks$chromEnd
  }

  all.GRanges <- GRanges(seqnames = peaks$chrom,
                         ranges = IRanges(start=ranges.start, end=ranges.end),
                         seqlengths = human.chromlens)

  all.GRanges <- all.GRanges[seqnames(all.GRanges) != "chrM"]

  if ('signalValue' %in% names(peaks))
    elementMetadata(all.GRanges) <- peaks$signalValue

  if (summit)
    all.GRanges <- resize(all.GRanges, 1, "center")
  return(all.GRanges)
}


test.significance <- function(gr,
                              labels,
                              tfbs.filename.list,
                              tfbs.details.tsv,
                              mappability,
                              gat.command='gat-run.py') {


  #parse narrowPeak files
  peaks <- lapply(tfbs.filename.list, parse.narrowPeak, summit=F)
  names(peaks) <- tfbs.filename.list

  #read peak detail files
  peakDetails <- read.table(tfbs.details.tsv)

  #store file names and antibody names as metadata columns
  peaks <- lapply(names(peaks),
                  function(pn){
                    mcols(peaks[[pn]])$antibody <- strsplit(toString(peakDetails[match(pn, peakDetails[,1]), 9]), '=')[[1]][2]
                    mcols(peaks[[pn]])$filename <- pn
                    peaks[[pn]]
                  })

  peaks <- do.call(c, peaks) #combine all TFBS peaks

  tfbs.file <- tempfile()

  #create BED files out of prepared GRanges object
  write.table(data.frame(seqnames(peaks),
                         start(peaks),
                         end(peaks),
                         mcols(peaks)$filename),
              file=tfbs.file,
              sep = '\t',
              quote = F,
              row.names = F,
              col.names = F)


  clust.file <- tempfile()

  #convert GRanges file to BED file
  write.table(data.frame(seqnames(gr),
                         start(resize(gr, 50, 'center')),
                         end(resize(gr, 50, 'center')),
                         paste0('cluster', labels)),
              file=clust.file,
              sep = '\t',
              quote = F,
              row.names = F,
              col.names = F)

  log.file.name <- tempfile()

  gat.args <- c('-s', tfbs.file,
                '-w', mappability.file,
                '-a', clust.file,
                '-L', log.file.name)
  gat.result <- system2(gat.command, gat.args, stdout=T)

  unlink(tfbs.file, clust.file, log.file.name)
  read.table(textConnection(paste0(gat.result, collapse = '\n')), header=T)
}

