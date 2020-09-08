install.packages("BiocManager")
BiocManager::install("KEGGREST")
BiocManager::install("KEGGgraph")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("magick")
library("KEGGREST")
library("KEGGgraph")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("pathfindR")
## demo input file: RA_input
?RA_input
dim(RA_input)
head(RA_input)
## pathway enrichment 
RA_demo <- run_pathfindR(RA_input,
                         iterations = 1,
                         # keeps running time low - default is 10
                         visualize_enriched_terms = FALSE) # needed until bug fixed in next release
head(RA_demo)
## clutser enriched terms
RA_demo_clu <- cluster_enriched_terms(RA_demo2)
## term-gene graph of top 10 terms
term_gene_graph(RA_demo2)
## pathway enrichment can use different networks, different gene sets
RA_demo2 <- run_pathfindR(na.omit(ttest),
                          gene_sets = "Reactome",
                          pin_name_path = "GeneMania",
                          visualize_enriched_terms = FALSE)
head(RA_demo2)




#
# download annotation for a specific chip – in this example the chip is
# Affy Human Genome U133A but there is one for almost any major platform
#
library("hgu133a.db")
# extract map of interest (probeID to GENE SYMBOL)
my.map <- hgu133aSYMBOL
# not all probeID have a mapping (ie an annotation)
mapped.probes <- mappedkeys(my.map)
# get Entrez ID for the Affy ID of interest (ie. first five)
my.first10.symbols <- as.data.frame(my.map[mapped.probes[1:10]])
# inspect result
my.first10.symbols[1:10,]
# check what other maps are available
ls("package:hgu133a.db")

