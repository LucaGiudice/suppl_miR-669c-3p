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
output_net1_file="./output/Network_expanded/net1.rda"
output_intermediate_file="./output/Network_expanded/interm.RData"
output_enr_filex="./output/Network_expanded/Unspecific_enrichment.xlsx"
output_enr_file_scan_genes="./output/Network_expanded/enr_scan_genes.txt"
output_enr_file_target_genes="./output/Network_expanded/enr_target_genes.txt"
MYD88="17874" #MYD88

#Load data ----
load(output_net1_file)

net_targ2scan_e=net_targ2scan_e[!is.na(net_targ2scan_e[,2]),]
gg <- graph.edgelist(as.matrix(net_targ2scan_e), directed=F)
cc = split(V(gg)$name, clusters(gg)$membership)

enr_l=list()
bol=FALSE
load("output/Network_expanded/map.rda")
mirna_genes=map2ens$entrezgene_id
scan_genes=map2ens_scan$entrezgene_id


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


  res.scan.l=list()
  for(k in 1:length(cc)){
    group=cc[[k]]
    scan_candidate=match(group,scan_genes)[!is.na(match(group,scan_genes))]
    scan_candidate=map2ens_scan[scan_candidate,1]
    res.scan.l[[k]]=toupper(scan_candidate)
  }
  
  lapply(res.scan.l, write, output_enr_file_scan_genes, append=TRUE, ncolumns=1000)
  
  res.target.l=list()
  for(k in 1:length(cc)){
    group=cc[[k]]
    target_candidate=match(group,mirna_genes)[!is.na(match(group,mirna_genes))]
    target_candidate=map2ens[target_candidate,1]
    res.target.l[[k]]=toupper(target_candidate)
  }
  
  lapply(res.target.l, write, output_enr_file_target_genes, append=TRUE, ncolumns=1000)
}

selected_enr=c(2,6,18,9,19,25)
cc=cc[selected_enr]
map.colours.l2=list()
map.colours.l=list()
nets.l=list()
for(k in 1:length(cc)){
  group=cc[[k]]
  sub_g=induced_subgraph(gg, group)
  sub_g=as.data.frame(as_edgelist(sub_g, names = TRUE));sub_g$V1=as.character(sub_g$V1);sub_g$V2=as.character(sub_g$V2)
  
  scan_candidate=match(group,scan_genes)
  scan_candidate=toupper(map2ens_scan[scan_candidate,1])
  map=cbind(group,scan_candidate)
  
  map.colours=cbind(scan_candidate[!is.na(scan_candidate)],rep(100,length(scan_candidate[!is.na(scan_candidate)])))
  colnames(map.colours)=c("gene","colour");map.colours=as.data.frame(map.colours)
  map.colours.l[[k]]=map.colours
  
  target_candidate=match(group[is.na(scan_candidate)],mirna_genes)
  target_candidate=toupper(map2ens[target_candidate,1])
  map[is.na(scan_candidate),2]=target_candidate
  
  map.colours2=cbind(target_candidate[!is.na(target_candidate)],rep(50,length(target_candidate[!is.na(target_candidate)])))
  colnames(map.colours2)=c("gene","colour");map.colours2=as.data.frame(map.colours2)
  map.colours.l2[[k]]=map.colours2
  
  sub_g$V1n=mapvalues(sub_g$V1,map[,1],map[,2], warn_missing = FALSE)
  sub_g$V2n=mapvalues(sub_g$V2,map[,1],map[,2], warn_missing = FALSE)
  sub_g=sub_g[,c(3,4)]
  nets.l[[k]]=sub_g
}

cyto.net=rbindlist(nets.l);class(cyto.net)="data.frame"
colour.map=rbindlist(c(map.colours.l,map.colours.l2));class(colour.map)="data.frame"

write.csv(colour.map,file="colour_map.csv",quote = FALSE, row.names = FALSE)
write.csv(cyto.net,file="cyto_map.csv",quote = FALSE, row.names = FALSE)












