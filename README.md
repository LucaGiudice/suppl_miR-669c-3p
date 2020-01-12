# Supplementary repository to replicate the network analysis in the paper: TITLE

[![Build Status](https://travis-ci.org/networkx/networkx.svg?branch=master)](https://singularity-hub.org/collections/3653)

Propagation of gene expression through HiC-derived networks illuminates active annotated enhancers and novel cell type-specific regulatory regions

********************************
***WORKFLOW OF THE METHOD***
********************************

![Test Image 8](https://raw.githubusercontent.com/LucaGiudice/ncl2/master/github_suppl/workflow.png?token=AHESZ3QXMHRY6R46YY5BABK53ZYXM)

********************************
***INSTALLATION AND USAGE***
********************************
The fastest way to use our method is by Singularity.
Otherwise, it is possible to download this structered repository, install the dependencies and change the multiple paths in the Rscripts based on your operating system for allowing R to work properly with the python softwares (method directory) and the input files (input directories).

#### By Singularity:
- install Singularity following the guide of the official site for your operating system: https://sylabs.io/guides/3.0/user-guide/installation.html#overview
- Open a terminal window
- Paste this command into the terminal:
    ```
    sudo singularity build --sandbox ncl/ shub://LucaGiudice/ncl
    ```
- Now you can access to the sandbox with everything already install to run or modify our method:
    ```
    sudo singularity shell --writable ncl/
    ```
- If you want to replicate our results:
    ```
    cd $HOME/../home/ && bash run_Rs.sh > log_execution.txt
    ```
- If you want to personalize the method, for instance, you want to use your own data:
    ```
    cd $HOME/../home/castle_vG/op1_data_processing
    cd input/ && rm Normalization1.DESeq2 && cp <path_to_your_DESeq2_dataset> Normalization1.DESeq2
    ```
- In case you want to work with our method, the Singularity container provides you Rstudio desktop:
    ```
    rstudio
    ```

![screenshot](https://raw.githubusercontent.com/LucaGiudice/ncl2/master/github_suppl/by_singularity.gif?token=AHESZ3QT725DCCEUEKORJ2C53ZYZM)

#### By the download of the github repository:
- Our method requires bash, anaconda 2 and R 3.4.4
- Install the following packages in R:
    ```                                 
  install.packages(c("pkgload","roxygen2","devtools","matrixStats","data.table","plyr","VertexSort","stringr","foreach","doSNOW","grid","doParallel","gridExtra","scales","svMisc","ggplot2","igraph","bigmemory","biganalytics","caTools"))
  biocLite()
  biocLite(c("Biobase","GenomicRanges","scater,DESeq2"))
    ```
- Then, download this repository
- Install the conda enviroment required to run Random walk and Singular walk:
    ```
    conda create --name rw --file $HOME/<your_path_to_this_repository>/methods/Random_walk/anconda_env/spec-file.txt
    ```
- Last steps:
   - change the python path in p3_2s_propagation/op1_propagation.R based on the position of your python
   - change the python path in p4_sing_propagation/op1_propagation.R based on the position of your python
   - download the following files and put them in their corresponding input directories of the repository:
    ```
    wget --no-check-certificate -O p1_input.zip https://univr-my.sharepoint.com/:u:/g/personal/luca_giudice_univr_it/EdMyCJsgEf5JsKw7GDn9QosB31dUouSV8NzGUrghqGYWtA?download=1
    wget --no-check-certificate -O p2_input.zip https://univr-my.sharepoint.com/:u:/g/personal/luca_giudice_univr_it/EZtv0o4mt0JBrz3iig2Z2D0BCJWuH709Zpj_vep3Vz98cQ?download=1
    wget --no-check-certificate -O p5_input.zip https://univr-my.sharepoint.com/:u:/g/personal/luca_giudice_univr_it/EWDiyK5P7r9BpSHX5_mx9zoBVUIqnzApvymfif-cRB10sg?download=1
    wget --no-check-certificate -O p6_input.zip https://univr-my.sharepoint.com/:u:/g/personal/luca_giudice_univr_it/ERfVEzJTRCJCrtIDHHY-kbgBZ87VdnENZZY0Lbo0M4n29w?download=1
    wget --no-check-certificate -O p8_input.zip https://univr-my.sharepoint.com/:u:/g/personal/luca_giudice_univr_it/Ead5JKvB5Y9CtkCnLNTcmtQBNIHeRsiYUSBWFbgpJn_3SA?download=1
    ```


********************************
***DESCRIPTION OF EACH FILE AND SUBDIRECTORY***
********************************

- ***methods/***
   - Random_walk/: contains the python code to compute a closed form of the random walk proposed by NBS2
   - Singular_walk/: contains the python code to compute the propagation of each element containing information in a cell-type specific expression profile
   - util_dfs_creation.R: takes each propagation done with the application of the pipeline and it puts them in a single compact list saved in Random_walk/output/
   
- ***p1_cell_profiles/***
   - input/: contains the cell-type specific expression dataset 
   - output/: when done, it contains the cell-type specific expression profile of each cell-type in the dataset. 
   - output/shuffled/: when done, it contains the shuffled cell-type specific expression profiles of each cell-type in the dataset (for the validation step)
   - op1_data_processing.R: processes the input dataset and creates a profile (1D vector of expression values) for each cell type. A profile is made computing the average of the expression of every gene in all the cells of the same type.
   - op2_profile_shuffling.R: takes the cell-type specific profiles and shuffles them k times (k defined by the user). 

- ***p2_net_creation/***
   - input/: dataset of chromatin interactions, fragment A | associated gene | fragment B | associated gene | cell-type score
   - output/: when done, it contains the user defined networks (by default: the Monocyte specific network, Tcell specific network, Neutrophil specific network, the Mono U Tcell U Neut network, the overall network of 13 cell types)
   - op1_net_creation.R: processes the dataset and creates the desired networks

- ***p2r_net_rewiring/***
   - output/: when done, it contains the rewired versions of the networks which have been built in the p2 step (for the validation step)
   - op1_net_creation.R: It takes the networks built in the p2 step, it converts them in igraph objects and it applies a rewiring algorithm preserving the starting degree of the nodes (the most restrictive form of rewiring)

- ***p3_2s_propagation/***
   - output/: when done, it contains the propagations of the profiles built in the step p1 with the networks built in the p2 step
   - op1_propagation.R: loads the profiles and the networks. It calls the python Random Walk software to propagate the profiles on the networks. 

- ***p3r_2s_propagation/***
   - output/rewiring/: when done, it contains the propagations of the profiles built in the step p1 with the rewired networks built in the p2r step
   - output/shuffling/: when done, it contains the propagations of the shuffled profiles built in the step p1 with the networks built in the p2 step
   - op1_propagation_rew.R: loads the profiles and the rewired networks. It calls the python Random Walk software to propagate the profiles on the rewired networks.
   - op1_propagation_shuf.R: loads the shuffled profiles and the standard networks. It calls the python Random Walk software to propagate the shuffled profiles on the networks.

- ***p4_sing_propagation/***
   - input/: during the p4 run, temporary files are saved here
   - output/: when done, the resulting propagations from the single element random walk diffusions are save here
   - op1_propagation.R: loads the profiles and the networks. It calls the python Singular_walk software. This decomposes the profiles such that each element is a new distinct profile. Then each element is propagated. In this way, this step gains the propagation result from each element.
   - op2_res_processing.R: loads the propagations made with op1 with respect each specific cell-type, it summarizes the results and it computes the contribution of each propagated element to the information gained to each node in the network.

- ***p5_validation_chr_enh and p6_validation_fantom5/***
   - input/: active cell-type specific enhancers from a database
   - output/: plot showing the comparison of how annotated enhancers have been ranked differently using different cell-type specific expression profiles.
   - op1_chr_enrich.R: load the annotated enhancers, ranks them based on cell-type specific profiles and plots the comparison. 

- ***p7_subg_plot/***
   - input/: contains the result of the p4 step with the contribution of each element in the cell-type specific profile
   - output/: when done, it contains the plot of the subgraph with the fragment of interest and the genes contributing to it
   - op1_plot_subg.R: takes the result of the p4 step, takes a user-defined network, a user-defined fragment of interest and finds the genes which contribute to the information which the fragment has gained with the propagation. Finally, it plots the subgraph.

********************************
***ADDITIONAL INFORMATION***
********************************

## License
Distributed under the MIT License.

## Contact for issues or support about the usage of the method
Luca Giudice - [github](https://github.com/LucaGiudice/) - luca.giudice@univr.it
