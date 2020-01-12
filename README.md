# Supplementary repository to replicate the network analysis in the paper: TITLE

The analysis, the plots and the enrichment tables are already in this repository (output directory) to allow immediately the best comprehension of the work made. Enterign more into the analysis, this is used to highlight the fact that: MYD88 is a strong candidate gene to study when miR-669c-3P is under investigation. 

Specifically, MYD88 and the TLR signaling pathway are important factors bound to miR-669c-3P target genes for two big considerations made with the result of this analysis:

- creating a network with only the direct strong interactions between targetscan genes (tsg) and the mirna targets (tg), MYD88 is one of few targetscan genes (MDGA1, FBXW11, IGFBP4, MYD88, FOXO1, CXCR1) showing both a connection to a real target and a statistically significant enrichment of a pathway connected by literature to neuroinflammation. 
- creating a network taking into account the first level of neighbours of each targetscan and real target genes, MYD88 is the only targetscan gene among the best previous ones that enhances its connection to miR-669c-3P. In fact, first: a lot of its neighbours are targetscan genes that create a path to two real targets (TLR6 and ATF3), second: MYD88 and its neighbours increase the statistical significance of TLR signaling pathway associated to miR-669c-3P. P.S. More precisely, the network has direct connections betweenn tsg and tg, undirect connections between a neighbour of a tsg than reaches a tg thanks to its tsg brother. 

********************************
***WORKFLOW OF THE METHOD***
********************************

![Test Image 8](https://raw.githubusercontent.com/LucaGiudice/suppl_miR-669c-3p/master/output/Network_expanded/images/network_analysis_overview.png)

********************************
***INSTALLATION AND USAGE***
********************************
It is possible to download this structered repository, install the dependencies and change the multiple paths in the Rscripts based on your operating system for allowing R to work properly with the and the input files (directories data and funs).


#### By the download of the github repository:
- Our analysis had been developed on R 3.6.1
- Install the following packages in R:
    ```                                 
  install.packages(c("data.table","matrixStats","plyr","readxl","data.table","doParallel","parallel","igraph","xlsx","BiocManager"))
  BiocManager::install(c("biomaRt","org.Mm.eg.db","ReactomePA"))
    ```
- Download and extract this repository
- Replace the dummy StringDB in the data directory with the real file (be carefull, it requires 1GB of disk space)
    ```
    https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz
    ```
- Now you can run the main scripts in the following order:
   - network_direct_analysis.R (it creates the direct network between the real mirna target genes and the targetscan genes)
   - enrichment_network_direct.R (it does the enrichment of the connected components in the "direct network")
   - network_expanded_analysis.R (it creates the expanded network)
   - enrichment_network_expanded.R (it does the enrichment of the connected components in the "expanded network")

********************************
***ADDITIONAL INFORMATION***
********************************

## License
Distributed under the MIT License.

## Contact for issues or support about the usage of the method
Luca Giudice - [github](https://github.com/LucaGiudice/) - luca.giudice@univr.it
