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

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("matrixStats"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("readxl"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("doParallel"))
suppressPackageStartupMessages(library("parallel"))
source("funs/fun_interactions.R")

## MAIN ## ----
#http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/mmu_MTI.xls
#https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz
#Set variables ----
input_db_mirna="./data/mmu_MTI.xls"
input_mirna="miR-669c-3p"
input_string_db="./data/10090.protein.links.v11.0.txt"
input_scan_db="./data/TargetScan7.2__miR-669c-3p.predicted_targets.xlsx"
cut_conf_String=900
ncores=4
high_conf=0.85
output_net1_file="./output/Network_expanded/net1.rda"
output_mapping_file="./output/Network_expanded/map.rda"
source("funs/fun_interactions.R")

#Load data ---
cat("Load data \n")
db=read_excel(input_db_mirna);class(db)="data.frame";db=db[,c(2,4)];
scan_db=read_excel(input_scan_db);class(scan_db)="data.frame";scan_db=scan_db[,c(1,11)];
#Find which is the score above there is the 15% of the best scan genes
scan_db[,2]=abs(scan_db[,2]);qs=quantile(scan_db[,2], probs = 0.85);th_score=qs[1];

#Detection of target genes ----
indxs_id=grep(input_mirna,db$miRNA)
mirna_genes=db$`Target Gene`[indxs_id]
cat("number of mirna_genes:", length(mirna_genes), "\n")

#Detection of predict target genes ----
scan_genes=scan_db[scan_db[,2]>th_score,1]
cat("number of scan_genes in the best 15%:", length(scan_genes), "\n")

#Map from genes symbol to ENSMUSP ensemble protein ----
cat("Map from entry name to ENSMUSP ensemble protein \n")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")

#View(listFilters(mouse));entrezgene_id
map2ens <- getBM(filters= "external_gene_name", attributes= c("external_gene_name","ensembl_peptide_id","entrezgene_id"),
                values=mirna_genes,mart= mouse)
map2ens=map2ens[apply(as.matrix(map2ens$ensembl_peptide_id),1,nchar)>0,]
map2ens$ensembl_peptide_id=paste("10090.",map2ens$ensembl_peptide_id,sep="")

map2ens_goal <- getBM(filters= "external_gene_name", attributes= c("external_gene_name","ensembl_peptide_id","entrezgene_id"),
                 values="MyD88",mart= mouse)
map2ens_goal=map2ens_goal[apply(as.matrix(map2ens_goal$ensembl_peptide_id),1,nchar)>0,]
map2ens_goal$ensembl_peptide_id=paste("10090.",map2ens_goal$ensembl_peptide_id,sep="")
map2ens_goal=map2ens_goal[1,]

map2ens_scan <- getBM(filters= "external_gene_name", attributes= c("external_gene_name","ensembl_peptide_id","entrezgene_id"),
                 values=scan_genes,mart= mouse)
map2ens_scan=map2ens_scan[apply(as.matrix(map2ens_scan$ensembl_peptide_id),1,nchar)>0,]
map2ens_scan$ensembl_peptide_id=paste("10090.",map2ens_scan$ensembl_peptide_id,sep="")

save(map2ens,map2ens_goal,map2ens_scan,file=output_mapping_file)

#Creation of the network ----
string_db=fread(input_string_db)
class(string_db)="data.frame";string_db=string_db[string_db$combined_score>=cut_conf_String,];
string_db=string_db[,c(1,2)]

cat("Searching the interaction of interest \n")
genes_interest=map2ens$ensembl_peptide_id
l_Ts_1neigh=find_neigh(genes_interest,string_db,ncores)
net_Ts_1neigh=rbindlist(l_Ts_1neigh);class(net_Ts_1neigh)="data.frame";
Ts_1neigh=net_Ts_1neigh[,2]

genes_interest=map2ens_scan$ensembl_peptide_id
l_Sc_1neigh=find_neigh(genes_interest,string_db,ncores)
net_Sc_1neigh=rbindlist(l_Sc_1neigh);class(l_Sc_1neigh)="data.frame";
Sc_1neigh=net_Sc_1neigh[,2]

#Assemble final network
keep_neigh=c(as.character(intersect(Ts_1neigh,Sc_1neigh)),map2ens$ensembl_peptide_id,map2ens_scan$ensembl_peptide_id)
colnames(net_Sc_1neigh) = colnames(net_Ts_1neigh) = c("main_protein","found_protein")
net_targ2scan=rbind(net_Sc_1neigh,net_Ts_1neigh);class(net_targ2scan)="data.frame"
net_targ2scan=remove_inter(net_targ2scan,ncores)
keep_neigh=match(keep_neigh,net_targ2scan[,2])
cat("number of intersection:",sum(!is.na(keep_neigh)),"\n")
net_targ2scan=net_targ2scan[keep_neigh[!is.na(keep_neigh)],]

#Save point ----
cat("Check and save safe point \n")

#Map fron ENSEMBLE to NAME
cat("Map fron ENSEMBLE to NAME \n")
nodes2map=unique(unlist(net_targ2scan))
nodes2map=gsub("10090.","",nodes2map)
map2name_net <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","entrezgene_id","external_gene_name"),
                      values=nodes2map, mart= mouse)
map2name_net$external_gene_name=toupper(map2name_net$external_gene_name)
map2name_net$ensembl_peptide_id=paste("10090.",map2name_net$ensembl_peptide_id,sep="")

#entrezgene_id
net_targ2scan$Ae <- mapvalues(net_targ2scan$main_protein, from=map2name_net$ensembl_peptide_id, 
                     to=map2name_net$entrezgene_id, warn_missing = FALSE)
net_targ2scan$Be <- mapvalues(net_targ2scan$found_protein, from=map2name_net$ensembl_peptide_id, 
                     to=map2name_net$entrezgene_id, warn_missing = FALSE)
net_targ2scan_e=net_targ2scan[,c(3,4)]

#external_gene_name
net_targ2scan$Aname <- mapvalues(net_targ2scan$main_protein, from=map2name_net$ensembl_peptide_id, 
                     to=map2name_net$external_gene_name, warn_missing = FALSE)
net_targ2scan$Bname <- mapvalues(net_targ2scan$found_protein, from=map2name_net$ensembl_peptide_id, 
                     to=map2name_net$external_gene_name, warn_missing = FALSE)
net_targ2scan_name=net_targ2scan[,c(5,6)]
cat("dimension of the final network:",nrow(net_targ2scan_e), "\n")

cat("Save results \n")
save(net_targ2scan_e,
     net_targ2scan_name,
     net_targ2scan,
     map2name_net,
     file=output_net1_file)

