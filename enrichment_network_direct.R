#Clean workspace and memory ----
rm(list=ls())
gc()

#Set working directory----
gps0=getwd()
gps0=paste(gps0,"/%s",sep="")
rootDir=gps0
setwd(gsub("%s","",rootDir))
set.seed(8)

#Load libraries and data----

library("data.table")
library("matrixStats")
library("plyr")
library("readxl")
library("biomaRt")
library("doParallel")
library("parallel")
source("funs/fun_interactions.R")
library("ReactomePA")
library("org.Mm.eg.db")
library("xlsx")
library("igraph")

## MAIN ## ----

#Set variables ----
input_db_mirna="./data/mmu_MTI.xls"
input_mirna="miR-669c-3p"
input_string_db="./data/10090.protein.links.v11.0.txt"
cut_conf_String=900
ncores=4
output_net1_file="./output/Network_direct/net1.rda"
output_intermediate_file="./output/Network_direct/interm.RData"
output_enr_filex="./output/Network_direct/Unspecific_enrichment.xlsx"
output_enr_file_scan_genes="./output/Network_direct/enr_scan_genes.txt"
output_enr_file_target_genes="./output/Network_direct/enr_target_genes.txt"
MYD88="17874" #MYD88

#Load data ----
load(output_net1_file)

net_targ2scan_e=net_targ2scan_e[!is.na(net_targ2scan_e[,2]),]
gg <- graph.edgelist(as.matrix(net_targ2scan_e), directed=F)
cc = split(V(gg)$name, clusters(gg)$membership)

enr_l=list()
bol=FALSE
mirna_genes=net_targ2scan_e[,1]
scan_genes=net_targ2scan_e[,2]

if(TRUE){
  for(k in 1:length(cc)){
    group=cc[[k]]
    cat("the cluster is big: ", length(group), "\n")
    if((sum(group %in% mirna_genes)>=1 & sum(group %in% scan_genes)>=1)){
      if(sum(MYD88 %in% group)!=0){
        cat("I found MYD in the group: ", k, "\n", "\n")
      }
      genes=unique(group)
      enr=enrichPathway(genes, organism = "mouse", pvalueCutoff = 0.01, pAdjustMethod = "bonferroni", qvalueCutoff = 0.05, minGSSize = 5, maxGSSize = 500, readable = FALSE)
      enr_l[[k]]=as.data.frame(enr)
      
      if(k!=1){bol=TRUE}
      if(nrow(as.data.frame(enr))!=0){
        write.xlsx(as.data.frame(enr), 
                   output_enr_filex, 
                   sheetName=as.character(k),
                   append=bol)
      }
    }else{next;}
  }
}

selected_enr=c(2,6,18,20,22,29)
cc=cc[selected_enr]
res.scan.l=list()
for(k in 1:length(cc)){
  group=cc[[k]]
  scan_candidate=match(group,scan_genes)[!is.na(match(group,scan_genes))]
  scan_candidate=net_targ2scan_name[scan_candidate,2]
  res.scan.l[[k]]=scan_candidate
}

lapply(res.scan.l, write, output_enr_file_scan_genes, append=TRUE, ncolumns=1000)

res.target.l=list()
for(k in 1:length(cc)){
  group=cc[[k]]
  target_candidate=match(group,mirna_genes)[!is.na(match(group,mirna_genes))]
  target_candidate=net_targ2scan_name[target_candidate,1]
  res.target.l[[k]]=target_candidate
}

lapply(res.target.l, write, output_enr_file_target_genes, append=TRUE, ncolumns=1000)
