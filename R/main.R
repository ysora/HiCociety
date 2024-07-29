options(scipen = 100)

#' Example download
#'
#' This function downloads an example Hi-C file from chromosome 19 of mice (mm10) Th1 cells.
#'
#' @description It downloads example Hi-C dataset.
#' @param destfile Path to save the example Hi-C file
#' @return Hi-C file download at destfile.
#' @examples
#' # Not run
#' # get_example_hic_file('./example.hic')
#' @export
get_example_hic_file <- function(destfile = "example.hic") {
  url <- "https://github.com/ysora/HiCociety/raw/main/example.hic"
  if (!file.exists(destfile)) {
    message("Downloading example.hic file...")
    utils::download.file(url, destfile, mode = "wb")
  } else {
    message("Example .hic file already exists.")
  }
  return(destfile)
}

#' Pakcage install and load
#'
#' Loading a package
#'
#' @author Sora Yoon, PhD
#' @description It loads a package.
#' @param package package name
#' @import BiocManager
#' @importFrom utils getFromNamespace install.packages
install_and_load <- function(package) {
  # Check if the package is installed
  if (!requireNamespace(package, quietly = TRUE)) {
    # If not installed, install it using BiocManager
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(package)
  }
  # Load the package
  library(package, character.only = TRUE)
}

#' Get contact Frequency
#'
#' Get Contact Frequncy from .hic file
#'
#' @author Sora Yoon, PhD
#' @description Retrieve contact frequency from .hic file using strawR package.
#' @param fname .hic data for any types of genome conformation capture data
#' @param chr Chromosome number of network extraction
#' @param resol DNA basepair Resolution. Default=10000
#' @import strawr
#' @examples
#' # Not run
#' # myhic=get_example_hic_file()
#' # getContactFrequency(myhic,19,5000)
#' @export
getContactFrequency <- function(fname, chr, resol){
  chr <- as.character(chr)
  cf <- straw("NONE", fname, chr, chr, "BP", resol)
  return(cf)
}

#' Contact probability
#'
#' Get Contact probablity
#'
#' @author Sora Yoon, PhD
#' @description It estimates contact probability based on the distance of a pair of a loci.
#' @param tab Output from getContactFrequency function.
#' @param farthest Maximum 1-D distance to search. Default=2Mb
#' @param resol Hi-C resolution for test. Default = 10000
#' @param prob Significance cutoff for negative binomial distribution. Default =0.975
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom Rcpp sourceCpp
#' @import fitdistrplus
#' @importFrom stats pnbinom
#' @examples
#' # Not run
#' # myhic=get_example_hic_file();
#' # mydf=getContactFrequency(myhic, 19, 5000);
#' # myprob=getContactProbability(mydf,farthest=2000000, resol=5000,prob=0.975);
#' @export
getContactProbability <- function(tab, farthest = 2000000, resol = 10000, prob) {
  bin1 = bin2 = cnt = c()
  AREA <- list()

  # Detect number of cores
  n_cores <- parallel::detectCores() - 1
  if (n_cores > 30) n_cores <- 30
  if (n_cores <= 0) n_cores <- 1

  # Create cluster and register parallel backend
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  # Define the function for unique pairs
  UniqueTwo <- function(vec1, vec2) {
    a <- paste(vec1, vec2, sep = "_")
    b <- unique(a)
    ab <- sapply(b, function(x) unlist(strsplit(x, split = "_")))
    res1 <- ab[1, ]
    res2 <- ab[2, ]
    return(list(uniVec1 = res1, uniVec2 = res2))
  }

  # Run the parallel processing
  A <- foreach::foreach(mylen = seq(0, farthest, resol), .combine = 'c', .init = list(list()),
                        .packages = c('fitdistrplus', 'stats')) %dopar% {
                          sameD <- which(tab$y - tab$x == mylen)
                          allCounts <- tab$counts[sameD]
                          fd <- fitdistrplus::fitdist(allCounts, distr = 'nbinom')
                          sz <- fd$estimate[1]
                          mu <- fd$estimate[2]
                          area <- stats::pnbinom(q = allCounts, size = sz, mu = mu, lower.tail = TRUE)
                          sigarea <- which(area > prob)
                          ABD <- UniqueTwo(allCounts, area)
                          AA <- list(list(contactfreq = ABD$uniVec1, prob = ABD$uniVec2))
                          list(Area = AA, bin1 = tab$x[sameD][sigarea], bin2 = tab$y[sameD][sigarea], cnt = area[sigarea])
                        }


  # Stop the cluster
  parallel::stopCluster(cl)
  sizeA = (length(A)-1)/4
  AREAindex = (0:(sizeA-1))*4+2
  AREA = A[AREAindex]
  names(AREA) <- seq(0,farthest,resol)

  bin1index = (0:(sizeA-1))*4+3
  bin1 = unlist(A[bin1index], use.names = F)
  bin2index = (0:(sizeA-1))*4+4
  bin2 = unlist(A[bin2index], use.names = F)
  cntindex = (0:(sizeA-1))*4+5
  cnt = unlist(A[cntindex],use.names=F)
  #names(AREA) = seq(0,farthest,resol)

  LB = length(bin1)

  return(list(AREA=AREA, original=data.frame(x=bin1,y=bin2,counts=cnt), len1=LB))
}

# #' Prune Network
# #'
# #' Prune Network
# #'
# pruneNet <- function(net){
#   halfOfAverageConnectivity <- length(E(net)) / length(V(net))
#   nodedegrees <- degree(net)
#   to_delete <- which(nodedegrees < halfOfAverageConnectivity)
#   if(length(to_delete)>0){
#     delname <- names(nodedegrees)[to_delete]
#     prune <- delete_vertices(net, delname)
#     return(prune)
#   }else{
#     return(net)
#   }
#
# }

#' HiC to network data format
#'
#' Convert HiC to network data format
#'
#' @author Sora Yoon, PhD
#' @description It converts Hi-C dataframe to network object.
#' @param ftab three-column data composed of locus1, locus2 and value
#' @examples
#' # Not run
#' # myhic=get_example_hic_file();
#' # ftab=getContactFrequency(paste0('./',myhic),19,5000);
#' # net = hic2network(ftab[100:200,]);
#' # plot(net)
#' @useDynLib HiCociety
#' @export
hic2network <- function(ftab){
  g = graph_from_data_frame(ftab, directed = FALSE, vertices = NULL)
  return(g)
}

#' Calculate Average Count within 5-pixel Padding
#'
#' This function calculates the average count within a 25kb padding
#' around each (x, y) coordinate pair.
#'
#' @param x Numeric vector of x-coordinates of contact frequency data frame.
#' @param y Numeric vector of y-coordinates of contact frequency data frame.
#' @param counts Numeric vector of contact frequency counts.
#' @param resol Integer specifying the HiC resolution.
#' @return A numeric vector of average counts.
#' @examples
#' # Nor run
#' # x <- c(1, 2, 3, 4, 5)
#' # y <- c(1, 2, 3, 4, 5)
#' # counts <- c(10, 20, 30, 40, 50)
#' # resol <- 10000
#' # calculate_avg_count(x, y, counts, resol)
#' @export
calculate_avg_count <- function(x, y, counts, resol) {
  .Call('_HiCociety_calculate_avg_count', x, y, counts, resol)
}



#' Retrieve chromosome names from .hic file
#'
#' To extract all chromosome names from .hic file
#'
#' @author Sora Yoon, PhD
#' @description It retrieves all chromosome names having longer than 2.5Mbp.
#' @param fname Path to .hic file
#' @examples
#' # Not run
#' # myhic=get_example_hic_file();
#' # get_all_chr_names(myhic)
#' @export
get_all_chr_names = function(fname)
{
  chromtab = strawr::readHicChroms(fname)
  filterChrom = chromtab$name[which(chromtab$length>5000*500)]
  filterChrom = sort(filterChrom)
  filterChrom = filterChrom[which(filterChrom!="ALL")]
  return(filterChrom)
}

#' Create module objects from the Hi-C data
#'
#' It generates a list of graph of significant interactions, module table and module elements.
#'
#' @author Sora Yoon, PhD
#' @description It generates a list of graph of significant interactions, module table and module elements.
#' @param fname  Path to .hic file
#' @param chr    chromosome numbersto run.
#' @param resol  Resolution of Hi-C data
#' @param nbprob  Negative binomial probability. Higher value gives smaller number of stronger interaction.
#' @param farthest  The maximum searching distance between two nodes
#' @param par.noise Parameter for noise removal. Default is 1, higher value gives more filtered interactions.
#' @param network.cluster.method Can select between 'louvain' as default and 'label.prop' which means the label propagation method.
#' @import igraph
#' @examples
#' # Not run
#' # myhic = get_example_hic_file();
#' # mycom = hic2community(myhic, "19", 5000, 0.975, 2000000, par.noise=1, 'louvain')
#' @export
hic2community <- function(fname, chr, resol, nbprob, farthest, par.noise = 1, network.cluster.method = 'louvain'){

  set.seed(1)
  chr = as.character(chr)
  totalNet = list()
  totalInfo  = data.frame(chr=character(), module_start = character(), module_end = character())
  totalModule = list()
  for(ch in chr)
  {
    cat("Chromosome: ",ch,"\n")
    tab <- getContactFrequency(fname, ch, resol)

    cat("--------Contact frequency table has been loaded (1/5)--------")
    fetab <- getContactProbability(tab, farthest, resol=resol, prob=nbprob)
    ftab_orig<-fetab$original
    cat("--------Contact probability estimation is completed (2/5)--------")
    ftab_orig_cnt= merge(tab, ftab_orig, by=c('x','y'))
    ftab_orig_cnt$counts = ftab_orig_cnt$counts.x
    ftab_orig_cnt$avg_count <- calculate_avg_count(ftab_orig_cnt$x, ftab_orig_cnt$y, ftab_orig_cnt$counts)

    min_c = min(tab$x)
    max_c = max(tab$y)
    width_hic = ((max_c - min_c)/resol) + 1
    height = farthest / resol
    area_hic = width_hic * height
    c_area = sum(tab$counts[which(abs(tab$x-tab$y) < farthest)])
    acf1 = c_area / area_hic
    df = ftab_orig_cnt[which(ftab_orig_cnt$avg_count> par.noise*acf1),]
    cat("--------Noise filtering is completed (3/5)--------")
    net <- hic2network(df)
    cat("--------Significant interaction network is generated (4/5)--------")
    if(network.cluster.method == 'louvain'){A = cluster_louvain(as.undirected(net))
    }else if(network.cluster.method== 'label_prop'){ A = cluster_label_prop(as.undirected(net))
    }else{stop("Choose network.cluster.method between 'louvain' and 'label_prop'.")}
    cat("--------Clustering data have been generated (5/5)--------")
    modu=sapply(1:max(A$membership), function(x) sort(as.numeric(A$names[which(A$membership==x)])))
    modu = modu[which(lapply(modu,length)>3)]
    modu2 <- lapply(modu, function(x){x = as.character(x); g1=subgraph(net,x);return(names(V(g1)))})
    EE = lapply(modu2, function(x){return(length(E(subgraph(net, vids = x))))})
    modu2 = modu2[order(unlist(EE),decreasing = T)]
    modu2 = lapply(modu2,function(x) sort(as.numeric(as.character(x))))
    totalNet = append(totalNet, list(net))
    totalInfo = rbind(totalInfo, data.frame(chr=ch, module_start = unlist(lapply(modu2, function(x) as.character(min(as.numeric(x))))),
                                            module_end = unlist(lapply(modu2, function(x) as.character(max(as.numeric(x)))))))
    totalModule = append(totalModule, modu2)
  }
  names(totalNet) = paste0("graph_",chr)

  # Gene name - add_Genes
  # Connectivity
  cat('Start to estimate connectivity, transitivity and centrality node.\n')
  pt = proc.time()
  CON = c()
  TRA = c()
  EV  = c()
  for(i in 1:nrow(totalInfo))
  {
    chrom = totalInfo$chr[i]
    idx = which(names(totalNet) == paste0('graph_',chrom))
    grp = totalNet[[idx]]
    ele = totalModule[[i]]
    subg = subgraph(grp, as.character(ele))
    conne = length(E(subg))
    trans = transitivity(subg, type='global')
    ev= eigen_centrality(subg)
    ev = names(ev$vector)[which.max(ev$vector)]
    CON = append(CON, conne)
    TRA = append(TRA, trans)
    EV  = append(EV , ev)
  }
  timetaken=proc.time() - pt
  cat('Elapsed time : ', timetaken[3], 's', sep="")
  totalInfo$connectivity = CON
  totalInfo$transitivity = TRA
  totalInfo$centrality_node=EV

  ord = order(totalInfo$connectivity, decreasing = T)

  totalInfo=totalInfo[ord,]
  totalModule=totalModule[ord]
  return(list(Graphs=totalNet, ModuleSummary = totalInfo, ModuleElements =totalModule))
}

#' All available Txdb
#'
#' Check all available Txdb package
#'
#' @author Sora Yoon, PhD
#' @description It finds all available Txdb packages used in add_Genes function.
#' @examples
#' get_txdb()
#' @importFrom AnnotationDbi mapIds
#' @export
get_txdb = function()
{
  available_packages <- BiocManager::available()
  txdb_packages <- available_packages[grep("TxDb", available_packages)]
  return(txdb_packages)
}

#' Add gene information
#'
#' Adding gene list to ModuleSummary data frame obtained from hic2community function.
#'
#' @author Sora Yoon, PhD
#' @description This function adds a column with a list of genes included in each locus to the ModuleSummary data frame of the hic2community function.
#' @param df The ModuleSummary data frame obtained by running hic2community function
#' @param speciesObj Any Txdb package name corresponding
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges findOverlaps pintersect makeGRangesFromDataFrame width
#' @examples
#' # Not run
#' # myhic = get_example_hic_file();
#' # mycom = hic2community(myhic, "19", 5000, 0.975, 2000000, par.noise=1, 'louvain')
#' # mycom$ModuleSummary = add_Genes(mycom$ModuleSummary, 'TxDb.Mmusculus.UCSC.mm10.knownGene')
#' @importFrom AnnotationDbi mapIds
#' @importFrom S4Vectors queryHits subjectHits
#' @export
add_Genes <- function(df, speciesObj) {
  # Prefix 'chr' to chromosome names if not present
  if (substr(df$chr[1], 1, 1) != "c") {
    df$chr <- paste0("chr", df$chr)
  }

  # Create a GRanges object
  gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = as.numeric(df$module_start), end = as.numeric(df$module_end))
  )

  # Load the TxDb package dynamically
  if (!requireNamespace(speciesObj, quietly = TRUE)) {
    stop(paste("Package", speciesObj, "not found"))
  }
  txdb <- getFromNamespace(speciesObj, asNamespace(speciesObj))

  # Define a mapping between TxDb and org packages
  org_mapping <- list(
    "TxDb.Hsapiens.UCSC.hg38.knownGene" = "org.Hs.eg.db",
    "TxDb.Mmusculus.UCSC.mm10.knownGene" = "org.Mm.eg.db",
    "TxDb.Rnorvegicus.UCSC.rn6.knownGene" = "org.Rn.eg.db",
    "TxDb.Dmelanogaster.UCSC.dm6.ensGene" = "org.Dm.eg.db",
    "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" = "org.Sc.sgd.db",
    "TxDb.Athaliana.BioMart.plantsmart28" = "org.At.tair.db",
    "TxDb.Btaurus.UCSC.bosTau9.refGene" = "org.Bt.eg.db",
    "TxDb.Cfamiliaris.UCSC.canFam3.refGene" = "org.Cf.eg.db",
    "TxDb.Celegans.UCSC.ce11.refGene" = "org.Ce.eg.db",
    "TxDb.Drerio.UCSC.danRer10.refGene" = "org.Dr.eg.db",
    "TxDb.EcoliK12.UCSC.ecoliK12.refGene" = "org.EcK12.eg.db",
    "TxDb.Ggallus.UCSC.galGal5.refGene" = "org.Gg.eg.db",
    "TxDb.Pfalcifarum.PlasmoDB.v46" = "org.Pf.plasmo.db",
    "TxDb.Sscrofa.UCSC.susScr11.refGene" = "org.Ss.eg.db",
    "TxDb.Xlaevis.UCSC.xenTro9.refGene" = "org.Xl.eg.db"
  )

  # Get the corresponding org package name
  orgPkg <- org_mapping[[speciesObj]]

  # Load the org package dynamically
  install_and_load(orgPkg)
  if (!is.null(orgPkg) && requireNamespace(orgPkg, quietly = TRUE)) {
    orgdb <- get(orgPkg)
  } else {
    warning(paste("Package", orgPkg, "not found. Gene symbols will not be annotated."))
    orgdb <- NULL
  }

  # Find genes and overlaps
  genes <- genes(txdb)
  overlaps <- findOverlaps(gr, genes)

  # Get gene IDs
  gene_ids <- genes$gene_id[subjectHits(overlaps)]

  if (!is.null(orgdb)) {
    # Get gene symbols if org package is loaded
    gene_symbols <- mapIds(orgdb, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  } else {
    # Use gene IDs as symbols if org package is not loaded
    gene_symbols <- gene_ids
  }

  # Add Genes column to the data frame
  df$Genes <- sapply(seq_len(nrow(df)), function(i) {
    overlapping_genes <- gene_symbols[queryHits(overlaps) == i]
    paste(overlapping_genes, collapse = ",")
  })

  return(df)
}


#' Rank difference between two conditions
#'
#' Rank difference between two conditions
#'
#' @author Sora Yoon, PhD
#' @description output table of rank difference of modules between cell types is generated.
#' @param wt hic2community result from condition 1
#' @param ko hic2community result from condition 2
#' @examples
#' # Not run
#' # RankDiff(mycom1, mycom2)
#' @export
RankDiff <- function(wt, ko)
{
  wt.dat = wt$ModuleSummary
  ko.dat = ko$ModuleSummary

  wt.range = GRanges(seqnames=wt.dat$chr, IRanges(start=as.numeric(wt.dat$module_start), end=as.numeric(wt.dat$module_end)), connectivity = wt.dat$connectivity, idx = 1:nrow(wt.dat), rank = order(wt.dat$connectivity, decreasing = T))
  ko.range = GRanges(seqnames=ko.dat$chr, IRanges(start=as.numeric(ko.dat$module_start), end=as.numeric(ko.dat$module_end)), connectivity = ko.dat$connectivity, idx = 1:nrow(ko.dat), rank = order(ko.dat$connectivity, decreasing = T))

  ov = findOverlaps(wt.range, ko.range)
  ov_ratio = pintersect(wt.range[queryHits(ov)], ko.range[subjectHits((ov))])
  width_ov = width(ov_ratio)
  width_wt = width(wt.range[queryHits(ov)])
  width_ko = width(ko.range[subjectHits(ov)])
  min_width = apply(cbind(width_wt,width_ko),1,min)
  ov_perc = width_ov/min_width
  ov_wt = width_ov/width_wt
  ov_ko = width_ov/width_ko
  overlap_wt = wt.range[queryHits(ov)]
  nooverlap_wt = setdiff(wt.range,overlap_wt)
  #nooverlap_wt = wt.range[which(wt.range %in% nooverlap_wt)]
  overlap_ko = ko.range[subjectHits(ov)]
  nooverlap_ko = setdiff(ko.range,overlap_ko)
  #nooverlap_ko = ko.range[which(ko.range %in% nooverlap_ko)]

  score_overlap = abs(log(overlap_wt$rank) - log(overlap_ko$rank))
  score_wt_only = abs(log(nooverlap_wt$rank)) # - log(length(wt.range)))
  score_ko_only = abs(log(nooverlap_ko$rank)) # - log(length(ko.range)))

  ovl = cbind(as.data.frame(overlap_wt), as.data.frame(overlap_ko), rank_diff=score_overlap, overlap_pixel = width_ov, overlap_perc_min = ov_perc, overlap_perc_cond1 = ov_wt, overlap_perc_cond2 = ov_ko)
  wtonly =cbind(as.data.frame(nooverlap_wt), data.frame(seqnames=NA,start=NA,end=NA, width=NA, strand=NA, connectivity = 3,idx = NA, rank = length(nooverlap_wt), rank_diff=score_wt_only),overlap_pixel = NA, overlap_perc_min = NA, overlap_perc_cond1 = NA, overlap_perc_cond2 = NA)
  koonly =cbind(data.frame(seqnames=NA,start=NA,end=NA, width=NA, strand=NA, connectivity = rep(3, length(nooverlap_ko)), idx = NA, rank = length(nooverlap_ko)), as.data.frame(nooverlap_ko), rank_diff=score_ko_only,overlap_pixel = NA, overlap_perc_min = NA, overlap_perc_cond1 = NA, overlap_perc_cond2 = NA)
  res = rbind(ovl, wtonly, koonly)
  res = res[order(res$rank_diff, decreasing=T),]
  return(res)
}

#' Numbers to color vector
#'
#' Numbers to color vector
#'
#' @author Sora Yoon, PhD
#' @description It maps numbers to colors.
#' @param vec Numeric vector to be converted to the color vector
#' @param col Colors to be mapped.
#' @param num Length of sequence where the min and max are those from the truncated / or the original vectors
#' @param range Truncated values' range
maptocolors <- function(vec,col,num=100,range=NULL)
{
  if (is.null(range) == TRUE)
  {
    breaks <- seq(min(vec), max(vec),length.out=num)
  }
  if (is.null(range) == FALSE)
  {
    vec[which(vec < range[1])] = range[1]
    vec[which(vec > range[2])] = range[2]
    breaks <- seq(range[1], range[2],length.out=num)
  }

  cols <- col(length(breaks) + 1)
  colvec = as.character(cut(vec, c(-Inf, breaks, Inf), labels=cols))
  return(colvec)
}

#' Visualization of module
#'
#' Visualization of module
#'
#' @author Sora Yoon, PhD
#' @description It draws a triangle heatmap and arcplot of a module
#' @param hicpath Path to the .hic file
#' @param HC.object The object name from hic2community result
#' @param moduleNum The row index of module to draw
#' @param resolution Resolution of HiC data
#' @param hic.norm Normalization method. If not, set 'NONE'
#' @param heatmap.color.range Min and max value of contact frequency, e.g., c(0,10)
#' @param heatmap.color Color for heatmap. For example, colorRampPalette(c("white","red))
#' @param arc.depth Height of arc plot
#' @param arc.color Arc color
#' @param nbnom.param Negative binomial probability cutoff. Higher cutoff gives less number of arcs.
#' @param txdb Character. One of Txdb list obtained from get_txdb().
#' @param highlight.centrality Boolean input to set if highlight eigenvector centrality node.
#' @param highlight.cent.col The color of arcs stemming from the centrality node.
#' @param highlight.node The coordiante of a node of which the user will highlight the arcs stemming from this node. Default=NULL
#' @param highlight.node.col The color of arcs stemming from the node which the user highlight.
#' @param show.sig.int Boolean. If TRUE, it marks significant contact on the triangle heatmap.
#' @param netinfo Boolean. If TRUE, it shows network information of the module as text in the plot.
#' @export
visualizeModule <- function(hicpath, HC.object, moduleNum, resolution, hic.norm, heatmap.color.range=NULL, heatmap.color = colorRampPalette(c('white','red')), arc.depth=10, arc.color = "gray80", nbnom.param=0.99, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,highlight.centrality, highlight.cent.col, highlight.node=NULL, highlight.node.col, show.sig.int=TRUE, netinfo){

  N = as.numeric(moduleNum) # 1
  r = resolution
  chr = HC.object$ModuleSummary$chr[N]
  chr = gsub("chr","",chr)
  chromstart = as.numeric(as.character(HC.object$ModuleSummary$module_start[N]))
  chromend = as.numeric(as.character(HC.object$ModuleSummary$module_end[N]))
  H = arc.depth
  max_y = (chromend - chromstart) / 2 / r + H + 20
  stepsize = r / 2
  hic = straw(norm = hic.norm, fname = hicpath, chr1loc = as.character(chr), chr2loc= as.character(chr), unit = "BP", binsize = r)
  myidx = which(hic$x >= chromstart & hic$y <= chromend)
  myhic = hic[myidx,]
  min_val = min(myhic$x)
  max_val = max(myhic$y)
  cols = rows = seq(min_val, max_val, r)

  graph.names = names(HC.object$Graphs)
  my.graph.name = paste('graph',chr,sep="_")
  graph.idx = which(graph.names == my.graph.name)
  tot.graph = HC.object$Graphs[[graph.idx]]
  my.compo = HC.object$ModuleElements[[N]]
  my.graph = subgraph(tot.graph, as.character(my.compo))
  my.graph.dat = as_data_frame(my.graph)
  filt.dat = my.graph.dat[which(my.graph.dat$counts.y > nbnom.param),]


  mymat = matrix(0, length(rows), length(cols))
  sigmat = matrix(FALSE, length(rows), length(cols))
  for(i in 1:nrow(myhic))
  {
    idx1 = which(cols == myhic$x[i])
    idx2 = which(rows == myhic$y[i])
    mymat[idx1,idx2] = mymat[idx2, idx1] = myhic$counts[i]
  }
  for(i in 1:nrow(filt.dat))
  {
    idx1 = which(cols == filt.dat$from[i])
    idx2 = which(rows == filt.dat$to[i])
    sigmat[idx1,idx2] = sigmat[idx2, idx1] = TRUE
  }
  hicm = mymat
  if(!is.null(heatmap.color.range)){
    max_z = heatmap.color.range[2]
    min_z = heatmap.color.range[1]
  }else{
    max_z = max(hicm)
    min_z = min(hicm)
  }
  hicmcol = matrix(maptocolors(hicm, heatmap.color, num=100, range = c(min_z, max_z)), nrow = nrow(hicm))

  # Empty plot
  plot(1,1,xlim=c(chromstart,chromend), ylim=c(0, max_y), type = 'n', xaxs = 'i', yaxs ='i', bty='n', xaxt='n', yaxt='n', xlab="",ylab="")
  # Fill plot
  for(rownum in 1:nrow(hicm)){
    y = H + 15 - 0.5
    x = chromstart + (rownum * 2 * stepsize) - (stepsize * 2)
    for(colnum in rownum:ncol(hicm))
    {
      x = x+stepsize
      y = y+.5
      xs = c(x-stepsize, x, x+stepsize, x, x-stepsize)
      ys = c(y, y+.5, y, y-.5, y)
      highlight.sig.int = sigmat[rownum,colnum]
      if(highlight.sig.int & show.sig.int){
        polygon(xs,ys, border='black', col = hicmcol[rownum,colnum])
      }else{
        polygon(xs,ys, border=NA, col = hicmcol[rownum,colnum])
      }
    }
  }
  # arcplot
  graph.names = names(HC.object$Graphs)
  my.graph.name = paste('graph',chr,sep="_")
  graph.idx = which(graph.names == my.graph.name)
  tot.graph = HC.object$Graphs[[graph.idx]]
  my.compo = HC.object$ModuleElements[[N]]
  my.graph = subgraph(tot.graph, as.character(my.compo))
  my.graph.dat = as_data_frame(my.graph)
  filt.dat = my.graph.dat[which(my.graph.dat$counts.y > nbnom.param),]
  cat("size of filtered data is:", nrow(filt.dat))
  start_points = list()
  end_points = list()
  for(i in 1:nrow(filt.dat))
  {
    start_points = append(start_points, list(c(as.numeric(as.character(my.graph.dat$from[i])), H+2)))
    end_points = append(end_points, list(c(as.numeric(as.character(my.graph.dat$to[i])), H+2)))
  }
  # Function to calculate the radius given the height H and chord length
  calculate_radius <- function(H, chord_length) {
    return((H^2 + (chord_length / 2)^2) / (2 * H))
  }
  # Calculate the radius for each arc and generate points
  arcs <- lapply(1:nrow(filt.dat), function(i) {
    start <- start_points[[i]]
    end <- end_points[[i]]
    # Calculate the chord length
    chord_length <- sqrt((end[1] - start[1])^2 + (end[2] - start[2])^2)
    # Calculate the radius
    radius <- calculate_radius(H, chord_length)
    # Calculate the midpoint of the chord
    mid_x <- (start[1] + end[1]) / 2
    mid_y <- (start[2] + end[2]) / 2
    # Calculate the center of the circle
    if (start[2] == end[2]) {
      center_x <- mid_x
      center_y <- mid_y + sqrt(radius^2 - (chord_length / 2)^2)
    } else if (start[1] == end[1]) {
      center_x <- mid_x + sqrt(radius^2 - (chord_length / 2)^2)
      center_y <- mid_y
    } else {
      slope <- (end[2] - start[2]) / (end[1] - start[1])
      perp_slope <- -1 / slope
      dx <- sqrt(radius^2 - (chord_length / 2)^2) / sqrt(1 + perp_slope^2)
      dy <- perp_slope * dx
      center_x <- mid_x + sign(H) * dx
      center_y <- mid_y + sign(H) * dy
    }
    # Calculate angles
    start_angle <- atan2(start[2] - center_y, start[1] - center_x)
    end_angle <- atan2(end[2] - center_y, end[1] - center_x)
    # Generate points for the arc
    theta <- seq(start_angle, end_angle, length.out = 100)
    x <- center_x + radius * cos(theta)
    y <- center_y + radius * sin(theta)
    list(x = x, y = y)
  })

  par(new=T)
  for (i in 1:nrow(filt.dat)) {
    lines(arcs[[i]]$x, arcs[[i]]$y, col = arc.color, lwd = 1)
  }
  # Optionally, add the start and end points
  for (i in 1:nrow(filt.dat)) {
    points(start_points[[i]][1], start_points[[i]][2], pch = 16, col = arc.color)
    points(end_points[[i]][1], end_points[[i]][2], pch = 16, col = arc.color)
  }

  if(highlight.centrality)
  {
    centra = HC.object$ModuleSummary$centrality_node[N]
    idx = which(filt.dat$from == centra | filt.dat$to == centra)
    for(i in idx){
      lines(arcs[[i]]$x, arcs[[i]]$y, col = highlight.cent.col, lwd = 1)
      points(start_points[[i]][1], start_points[[i]][2], pch = 16, col = highlight.cent.col)
      points(end_points[[i]][1], end_points[[i]][2], pch = 16, col = highlight.cent.col)
    }
  }

  if(!is.null(highlight.node))
  {
    hl = highlight.node
    idx = which(filt.dat$from == hl | filt.dat$to == hl)
    if(length(idx)==0){warning("No such nodes in the module interaction list.")
    }else{
      for(i in idx){
        lines(arcs[[i]]$x, arcs[[i]]$y, col = highlight.node.col, lwd = 1)
        points(start_points[[i]][1], start_points[[i]][2], pch = 16, col = highlight.node.col)
        points(end_points[[i]][1], end_points[[i]][2], pch = 16, col = highlight.node.col)
      }
    }
  }
  # add gene track
  if(substr(chr,1,1) !="c"){chr = paste0("chr",chr)}
  gr = GRanges(
    seqnames = chr,
    ranges = IRanges(start = as.numeric(chromstart), end = as.numeric(chromend))
  )
  if(is.null(txdb)){cat("Provide a proper txdb object."); exit(1)}
  org_mapping <- list(
    "TxDb.Hsapiens.UCSC.hg38.knownGene" = "org.Hs.eg.db",
    "TxDb.Mmusculus.UCSC.mm10.knownGene" = "org.Mm.eg.db",
    "TxDb.Rnorvegicus.UCSC.rn6.knownGene" = "org.Rn.eg.db",
    "TxDb.Dmelanogaster.UCSC.dm6.ensGene" = "org.Dm.eg.db",
    "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene" = "org.Sc.sgd.db",
    "TxDb.Athaliana.BioMart.plantsmart28" = "org.At.tair.db",
    "TxDb.Btaurus.UCSC.bosTau9.refGene" = "org.Bt.eg.db",
    "TxDb.Cfamiliaris.UCSC.canFam3.refGene" = "org.Cf.eg.db",
    "TxDb.Celegans.UCSC.ce11.refGene" = "org.Ce.eg.db",
    "TxDb.Drerio.UCSC.danRer10.refGene" = "org.Dr.eg.db",
    "TxDb.EcoliK12.UCSC.ecoliK12.refGene" = "org.EcK12.eg.db",
    "TxDb.Ggallus.UCSC.galGal5.refGene" = "org.Gg.eg.db",
    "TxDb.Pfalcifarum.PlasmoDB.v46" = "org.Pf.plasmo.db",
    "TxDb.Sscrofa.UCSC.susScr11.refGene" = "org.Ss.eg.db",
    "TxDb.Xlaevis.UCSC.xenTro9.refGene" = "org.Xl.eg.db"
  )

  # Get the corresponding org package name
  orgPkg <- org_mapping[[txdb]]

  install_and_load(orgPkg)
  if (!is.null(orgPkg) && requireNamespace(orgPkg, quietly = TRUE)) {
    orgdb <- get(orgPkg)
  } else {
    warning(paste("Package", orgPkg, "not found. Gene symbols will not be annotated."))
    orgdb <- NULL
  }
  txdb <- getFromNamespace(txdb, asNamespace(txdb))
  genes = genes(txdb)
  overlaps = findOverlaps(gr, genes)
  gene_ids = genes$gene_id[subjectHits(overlaps)]
  gene_rg = genes[subjectHits(overlaps)]
  gene_symbols = mapIds(orgdb, keys = gene_ids , column = "SYMBOL", keytype="ENTREZID", multiVals = "first")
  gene_rg$gene_symbol = gene_symbols
  abline(h=c(H + 2))
  coord_write = round(seq(chromstart/1000000, chromend/1000000, .1),1)
  coord_orig  = coord_write * 1000000
  segments(x0 = c(chromstart,coord_orig,chromend), y0 = H+1.5, y1 = H+2.5)
  text(x = coord_orig, y=H+4.5, labels = coord_write)
  text(x = chromstart+r, y = H+5.5, labels = chr)
  text(x = chromend-r, y = H+5.5, labels = "Mb")

  for(g in 1:length(gene_rg))
  {
    tra = gene_rg[g]
    st = start(tra)
    en = end(tra)
    gsy = tra$gene_symbol
    str = as.character(tra@strand@values)
    if(str == "+"){Arrows(x0 = st, x1 = en, y0 = H+8, y1 = H+8, arr.type = 'triangle', arr.lwd=3,lwd=6, col='purple')}
    if(str == "-"){Arrows(x0 = en, x1 = st, y0 = H+8, y1 = H+8, arr.type = 'triangle', arr.lwd=3,lwd=6, col='pink')}
    text(x = ((st + en)/2)-10000, y = H+11, labels = gsy)
  }
  if(netinfo)
  {
    x0 = chromstart+r*10
    y0 = max_y
    text(x0, y0-3, paste0("Rank:",N), pos=4)
    text(x0, y0-9, paste0("Connectivity:",HC.object$ModuleSummary$connectivity[N]), pos=4)
    text(x0, y0-15, paste0("Transitivity:",signif(HC.object$ModuleSummary$transitivity[N],3)), pos=4)
    text(x0, y0-21, paste0("Centrality node:",HC.object$ModuleSummary$centrality_node[N]), pos=4)
  }
}
