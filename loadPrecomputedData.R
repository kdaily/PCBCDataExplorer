## Load precomputed data
## Should check in precomputed_data/preCompute.R to make sure there are no
## needed updates

## Load in the hg19 annotations and grouped annotations
hg19_annot_obj <- synGet("syn4943381")
hg19_annot <- readRDS(getFileLocation(hg19_annot_obj))

hg19_grpd_obj <- synGet("syn4943391")
hg19_grpd <- readRDS(getFileLocation(hg19_grpd_obj))

#sample gene list of the user input area
sample_genes_obj <- synGet("syn4943393")
df <- read.delim(getFileLocation(sample_genes_obj), sep="\t")
sample_gene_list <- as.character(unique(df$feature))

# sample_miRNAs_obj <- synGet("syn4943396")
# df <- read.delim(getFileLocation(sample_miRNAs_obj), sep="\t")
# sample_miRNAs <- as.character(unique(df$feature))

sample_miRNAs_obj <- synGet("syn4609631")
df <- read.delim(getFileLocation(sample_miRNAs_obj), sep="\t")
sample_miRNAs <- as.character(unique(df$GeneSymbol))

sample_methyl_obj <- synGet("syn4943397")
df <- read.delim(getFileLocation(sample_methyl_obj), sep="\t")
sample_methyl <- as.character(unique(df$feature))

#get the list siginificant genes from comparative analysis in synapse
flog.info('Reading the precomputed significant gene list')
sigGenesObj <- synGet("syn4484232")
sigGenes <- fread(getFileLocation(sigGenesObj), data.table=FALSE)

# Standardized comparison names and levels
compNames <- synTableQuery("SELECT * FROM syn4483642")@values

sigGenes2 <- sigGenes %>%
  dplyr::rename(feature=GeneSymbol, comparison=Comparison) %>%
  filter(abs(logFC) >= 1.5,
    adj.P.value <= 0.001,
         str_detect(comparison, "^All__")) %>%
  mutate(direction=str_extract(comparison, "__.*"),
         direction=str_replace(direction, ".*__", ""),
         comparison=str_replace(comparison, "__up", ""),
         comparison=str_replace(comparison, "__down", "")) %>%
  left_join(compNames) %>%
  mutate(pretty=str_c(variable1Short, variable2Short, sep=" vs. ")) %>%
  filter(class != "Culture_Conditions") # for now, until the names are cleaned up

sigGenesList <- dlply(sigGenes2, .(pretty),
                      function(x) x$feature)

# #########
# #read the precomputed enriched pathway list
# ########
# flog.info('Reading the precomputed enriched pathway list')
# df_precomputed_enrichedPathways_in_geneLists = readRDS("precomputed_data/precomputed_enrichedPathways_in_geneLists.rds")
# df_precomputed_enrichedPathways_in_geneLists$pathways_with_pvalue =  paste(df_precomputed_enrichedPathways_in_geneLists$pathways,
#                                                                            '#p.adj_',
#                                                                            format.pval(df_precomputed_enrichedPathways_in_geneLists$p.adj,digits=2),
#                                                                            sep='')
# #creating a list of list 
# precomputed_enrichedPathways_in_geneLists = split(df_precomputed_enrichedPathways_in_geneLists$pathways_with_pvalue,
#                                                   df_precomputed_enrichedPathways_in_geneLists$significant_gene_list_name)
# 
# 
# #HACK
# #For each geneList add another PATHWAY TYPE "ALL" which indicates use all the pathways for the shiny SERVER/UI
# # in this case genes in all the enriched pathways will be shown on the heatmap
# precomputed_enrichedPathways_in_geneLists <- lapply(precomputed_enrichedPathways_in_geneLists,
#                                                     function(x) { x[length(x)+1] = 'ALL'; x})
