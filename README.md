# R-code-S4_Class-protein-clustering-based-on-data-integration-of-corum-and-inparanoid

A machine learning-based approach to identify reliable gold standards for protein complex composition prediction.

![Small world analysis](smallworldgraph (1).png)
---
## Introduction
The project addresses a fundamental issue of using CORUM complexes as reference knowledge in evaluation of protein complex predictions. In plant cell extracts, co-fraction mass spectrometry (CFMS) data and some published literature indicated that fully assembled CORUM complexes rarely exist. Using inaccurate gold standards would result in a wrong validation dataset, misleading prediction models to unreliable predictions. To overcome this issue, we have developed machine learning approaches to identify a refined set of gold standards by integrating the information of CORUM and CFMS. I hope you and your reviewers find the.

 
## Install

1. First recommend to have R >= 4.2.0 


```
$ install.packages("apcluster")
$ install.packages("kohonen")
```
2. once installation of ‘kohonen’ completed, find path in your computer:
     path...\R\win-library\R version...\kohonen\Distances
3. copy and paste the file “wcc3.cpp” into “Distances” folder found in step 2.
4. Set the path in row in the following R code and copy and paste "S4 data integration corum inparanoid clustering.R" and the input data for rice clustering in the folder under the path that just created.


Here is a list of dependent packages:

```
1. library(tempR)
2. library(methods)
3. library(MASS)
4. library(dplyr)
5. library(ggplot2)
6. library(ggrepel)
7. library(openxlsx)
8. library(tibble)
9. library(tidyverse)
10. library(kohonen) # self organizing map
11. library(apcluster) # affinity propagation
```

---

## Run
Here the R commands that you need to run:
```
rm(list = ls()) 
library(kohonen)  
library(apcluster) 
library(tempR)
library(methods)
library(MASS)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(tibble)
library(tidyverse)
require(Rcpp) # help to integrate R and C++ via R functions
sourceCpp(system.file("Distances", "wcc3.cpp", package = "kohonen",  mustWork = TRUE)) # expose C++ function wcc3 to R level

getwd()
workingDir = "create your path"
setwd(workingDir)
source("S4 data integration corum inparanoid clustering.R")
options(stringsAsFactors = FALSE)

Protein.Complexes.Gold.Standard.Prediction.Rice <- new("Protein.Complexes.Gold.Standard.Prediction",
                                                       mammalPlantSpecies=c("human", "rice"),
                                                       mammalCORUM.Data=get(load(file="Data_CORUM_human.RData")),
                                                       mammalPlantOrtholog.Data=read.table("SQLtable.Human.fasta-Rice.FA", header = T),
                                                       mammalCORUM.StoichiometryData=get(load(file="CORUM_human_stoichiometry.RData")),
                                                       plantProtnames.Data=get(load(file="rice_protnames_data.RData")),
                                                       plantMmono.Data=get(load(file="rice_mmono_data.RData")),
                                                       peakProfileMatrixB1=get(load("Matrix_SEC_rice_1.RData")),
                                                       peakProfileMatrixB2=get(load("Matrix_SEC_rice_2.RData")),
                                                       reproduciblePeakFeatures.Data=get(load(file="Rice_Reproducible_SEC_data_plant.RData")),
                                                       singlePeakProts=get(load(file="Rice_prots_reproducible_single_peak.RData")),
                                                       multiplePeakProts=get(load(file="Rice_protsreproducible_multiple_peak.RData")),
                                                       mostCommonProtIDExample="LOC_Os03g51600.1",
                                                       genomeCoverageCutoff=0.66666,
                                                       RappCutoffNonmono = 1.6,
                                                       RappCutoffChoosingThreshold=1,
                                                       parametersSOM=list(num_of_multiple_of_sample_size=500, topo_shape="rectangular", toroidal_choice=F,  
                                                                          user_weight=c(0.35, 0.35, 0.15, 0.15), 
                                                                          dist_measure=c("WCCd3","WCCd3", "euclidean", "euclidean")),
                                                       parametersAPC=list(q_1=0.95,q_2=0.85,q_3=0.8,maxits=50000, convits=10000, lam=0.9,detail=F),
                                                       crossCorrelationRadius=3,
                                                       exogenousClusteringResult=get(load(file="Clustering_resultRice.RData")),
                                                       targetCluster="X1000")

clustering_rice <- ProteinComplexGoldStandardPredictionWithinCORUMOrthologousComplexes(Protein.Complexes.Gold.Standard.Prediction.Rice)
```



---

## Parameter Definition


    * If you want to run EPPC with SPIFFED scores, then you can set this parameter to "<b>`-s`  `11101001`</b>". (* note that there are 8 characters in the string). In this example, it will use Mutual Information, Bayes Correlation, Euclidean Distance, Jaccard Score and Apex Score. To specify the correlation scores to use:

2. <b>`input_directory`</b>: This parameter stores the input directory where you store your elution profile file. It is recommended to use the abosulte path instead of relative path.


3.  (`-c` <b>`gold_standard_file_path`</b>) or (`--cluster` <b>`gold_standard_file_path`</b>): This parameter stores the path to the gold standard file that you curated.

4. <b>`output_directory`</b>: This parameter stores the path to the ouput directory. Make sure that you've already created the directory before running the command. It is recommended to use the abosulte path instead of relative path.

5. (`-o` <b>`output_filename_prefix`</b>) or (`--output_prefix` <b>`output_filename_prefix`</b>): You can specify a prefix name for all the output files. The default is "Out"

6. (`-M` <b>`training_method`</b>) or (`--classifier` <b>`training_method`</b>): This parameter specifies what kind of classifier that you use. Possible options include <b>`RF`</b>, <b>`CNN`</b>, <b>`LS`</b>. Note that <b>`RF`</b> must comes with selected SPIFFED scores like "<b>`-s`  `11101001`</b>" instead of raw elution profile ("<b>`-s` `000000001`</b>"). <b>`CNN`</b> and <b>`LS`</b> must come with raw elution profile ("<b>`-s` `000000001`</b>").

7. (`-n` <b>`number_of_cores`</b>) or (`--num_cores` <b>`number_of_cores`</b>): You need to specify the number of cores used to run EPPC, the default number is 1. Assume you want to use six cores to run SPIFFED, you can set "<b>`-n` `6`</b>"

8. `--LEARNING_SELECTION` <b>`learning method selection`</b>: This parameter specifies whether you want to use <b>supervised learning</b> or <b>semi-supervised learning</b>. If you want to run with <b>supervised learning</b>, then set "<b>`--LEARNING_SELECTION` `sl`</b>" (Your <b>`training_method`</b> can be <b>`RF`</b> or <b>`CNN`</b>); if you want to run with <b>semi-supervised learning</b>, then set "<b>`--LEARNING_SELECTION` `ssl`</b>" (Your <b>`training_method`</b> can be <b>`CNN`</b> or <b>`LS`</b>).

9. `--K_D_TRAIN` <b>`fold_or_direct_training`</b>: Set <b>`d`</b> to directly train the model; set <b>`k`</b> to run with k-fold training. (options: <b>`d`</b> and <b>`k`</b>; default: <b>`d`</b>)

10. `--FOLD_NUM` <b>`number_of_folds`</b>: If you set `--K_D_TRAIN` <b>`k`</b>, then this parameter stores how many folds you are going to evaluate your mode. Note that this parameter must be bigger than `2`. (default: <b>`5`</b>)

11. `--TRAIN_TEST_RATIO` <b>`testing_data_ratio`</b>: This parameter stores the ratio of testing data to all data. (default: <b>`0.3`</b>)

12. `--POS_NEG_RATIO` <b>`negative_PPIs_ratio`</b>: This parameter stores the ratio of negative PPIs to positive PPIs. (default: <b>`1`</b>)

13. `--NUM_EP` <b>`number_of_elution_profiles`</b>: This parameter stores the number of elution profiles inside each PPI. (default: <b>`2`</b>)

14. `--NUM_FRC` <b>`number_of_fractions`</b>: This parameter stores the number of fractions in the elution profile file. (default: <b>`27`</b>)

15. `--CNN_ENSEMBLE` <b>`number_of_fractions`</b>: This parameter is a boolean value. If it's `0`, users need to provide one elution profile; if it's `1`, users need to provide multiple elution profiles.


---

## SPIFFED Command Examples

To run SPIFFED:
    
> `python ./main.py -s 000000001 /ccb/salz3/kh.chao/SPIFFED/input/EPIC_DATA/beadsALF -c /ccb/salz3/kh.chao/SPIFFED/input/EPIC_DATA/Worm_reference_complexes.txt /ccb/salz3/kh.chao/SPIFFED/output/EPIC_DATA/beadsALF/TEST/CNN_SL/FOLDS/beadsALF__K_D__k__CNN_SL__fold_number_5__negative_ratio_5__test_ratio_30 -o TEST -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN k --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.3 --POS_NEG_RATIO 5 --CNN_ENSEMBLE 0
`

---

## Multiple

To run SPIFFED with ensemble model:
> `python ./main.py -s 000000001 /home/kuan-hao/SPIFFED/input/OUR_DATA/intensity_HML_ensemble/ -c /home/kuan-hao/SPIFFED/input/OUR_DATA/gold_standard.tsv /home/kuan-hao/SPIFFED/output/SELF_DATA/intensity_HML_ensemble__negative_ratio_5/ -o out -M CNN -n 10 -m EXP -f STRING --LEARNING_SELECTION sl --K_D_TRAIN d --FOLD_NUM 5 --TRAIN_TEST_RATIO 0.7 --POS_NEG_RATIO 5 --NUM_EP 2 --NUM_FRC 27 --CNN_ENSEMBLE 1`
