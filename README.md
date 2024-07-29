# HiCociety 
### (Ha-ee-so-sa-ee-uh-ti)
HiCociety is a compound word of **'Hi-C'** and **'Society'**. 
<br>
Through the network analysis of Hi-C or its variant data, it shows the chromatin elements interact each other and forms huge network.
<br>
This network splits into ***modules*** each of which is a cluster of network. 
This concept is similar as the topologically associating domain (TAD) but unlike TAD defined based on the isulation score change, it is defined based on the density of significant chromatin contact.
HiCociety, as a network-based analysis tool for 3D chromatin conformation capture data such as Hi-C, it finds ***modules*** and estimate their ***connectivity, transitivity and centrality node***.
In addition, it provides a function that **compares module connectivities*** from two Hi-C datasets.

## The Usage
<br>
# --------------------------------------
### 1. Finding modules
# --------------------------------------

myhic = get_example_hic_file();
mycom = hic2community(fname=myhic, chr="19", resol=5000, nbprob=0.975, farthest=2000000, par.noise = 1, network.cluster.method = 'louvain')
head(mycom$ModuleSummary)

where,
fname : Path to Hi-C data
chr : Chromosome numbers to test/ e.g., c(1,2,3,4,5)
resol : Reolution of Hi-C data
nbprob : Negative binomial probability of chromatin contact
farthest : Searching limit of 1D distance of contact pair
par.noise : Noise removal parameter (default = 1)
network.cluster.method = 'louvain' or 'label.prop'
# --------------------------------------
### 2. Visualization of module
# --------------------------------------
visualizeModule(
<br>
  hicpath = myhic,
  <br>
  HC.object = myhic$ModuleSummary,
  <br>
  moduleNum = 1,
  <br>
  resolution = 5000,
  <br>
  hic.norm = 'VC_SQRT',
  <br>
  heatmap.color.range = c(0,10),
  <br>
  heatmap.color = colorRampPalette(c("white", "red")),
  <br>
  arc.depth = 10,
  <br>
  arc.color = "gray80",
  <br>
  nbnom.param = 0.99,
  <br>
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
  <br>
  highlight.centrality = TRUE,
  <br>
  highlight.cent.col = 'purple',
  <br>
  highlight.node = NULL,
  <br>
  highlight.node.col = NULL,
  <br>
  show.sig.int = TRUE,
  <br>
  netinfo = TRUE
  <br>
)



Author: Sora Yoon (sora.yoon@pennmedicine.upenn.edu)


