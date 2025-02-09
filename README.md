# sRNARFTarget: A machine learning-based approach for fast sRNA Target Prediction #
  
  ## Introduction

This repository contains the nextflow pipeline sRNARFTarget to obtain transcriptome-wide bacterial sRNA target predictions, and all the code, data, results and supplementary files related to the sRNARFTarget's manuscript. In the text below, we provide instructions to run sRNARFTarget.
    
  ## Requirements
  
**Requirements for running sRNARFTarget with our Python 3.8 docker container:
1. [Nextflow v21.04.1](https://www.nextflow.io/)
2. [Docker](https://docs.docker.com/)

**Requirements for running sRNARFTarget with Python installed locally:
  1. [Nextflow v21.04.1](https://www.nextflow.io/)
  2. Python3 (tested on 3.8.10)
  3. Python modules: [pickle](https://docs.python.org/3/library/pickle.html), [biopython v1.79](https://biopython.org/), [pprint](https://docs.python.org/3/library/pprint.html), [scikit-bio v0.5.6](http://scikit-bio.org/), [itertools](https://docs.python.org/3/library/itertools.html), [scikit-learn v0.24.1](https://scikit-learn.org/stable/), [pandas v1.2.1](https://pandas.pydata.org/), [numpy v1.19.5](https://numpy.org/).
        
The following modules are additionally required for running sRNARFTarget_SHAP and sRNARFTarget_CP
  1. [shap v0.39](https://pypi.org/project/shap/), [pyCeterisParibus v0.5.2](https://github.com/ModelOriented/pyCeterisParibus) and matplotlib v3.3.4
  
These modules can be installed using pip.

  ## Instructions to run sRNARFTarget
  
  1. Clone the repository. For instructions about how to clone GitHub repositories see [this](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository-from-github/cloning-a-repository).
  2. Create the fasta files with the sRNA and mRNA nucleotide (i.e, A,C,G,T) sequences in the folder containing the sRNARFTarget.nf script (referred from now on as sRNARFTarget folder/directory).
  3. Go to the sRNARFTarget folder so that it is the current working directory.
  4. OPTION A (with our [Python 3.8 docker container](https://hub.docker.com/r/penacastillolab/python38env), recommended). After pulling the docker container, type the command below to run sRNARFTarget replacing sRNA.fasta (--s parameter) and mRNA.fasta (--m parameter) with the corresponding filenames of the fasta files containing the sRNAs and mRNAs sequences, respectively. Both files should be located in the sRNARFTarget directory.
   ```
    nextflow run sRNARFTarget.nf --s sRNA.fasta --m mRNA.fasta -with-docker penacastillolab/python38env
   ```
   4. OPTION B (with Python installed locally). Type the command below to run sRNARFTarget replacing sRNA.fasta (--s parameter) and mRNA.fasta (--m parameter) with the corresponding filenames of the fasta files containing the sRNAs and mRNAs sequences, respectively. Both files should be located in the sRNARFTarget directory.
   ```
    nextflow run sRNARFTarget.nf --s sRNA.fasta --m mRNA.fasta
   ```
   
  ## Creation of all possible sRNA-mRNA pairs
  
sRNARFTarget creates all possible pairs from the input sRNA and mRNA sequences. Each sRNA is paired with all mRNAs. For example, if the input sRNA file has 5 sRNA sequences and mRNA file has 9 mRNA sequences, then it will create 45 sRNA-mRNA pairs, 9 pairs per sRNA.
 
  ## sRNARFTarget Results
    
  1. On your terminal, you should see something like this after sRNARFTarget's execution:
  ```
  N E X T F L O W  ~  version 21.04.1
Launching `sRNARFTarget.nf` [gloomy_easley] - revision: 273666007b
executor >  local (5)
[d0/aab5ef] process > createAllPossiblePairs              [100%] 1 of 1 ✔
[08/c1a4ba] process > getsRNATrinucleotidesFrequncies (1) [100%] 1 of 1 ✔
[c8/e2134c] process > getmRNATrinucleotidesFrequncies (1) [100%] 1 of 1 ✔
[eb/0b3d90] process > runRandomForestModel (1)            [100%] 1 of 1 ✔
[c8/b7f154] process > generateSortedResultFile (1)        [100%] 1 of 1 ✔

Pipeline execution summary
---------------------------
Run as : nextflow run sRNARFTarget.nf --s Multocida_sRNA_gcvb.fasta --m Multocida_mRNA.fasta
Completed at: 2021-06-07T14:28:19.852-02:30
Duration : 40.7s
Success : true
workDir : Afolder/sRNARFTarget/work
exit status : 0
```
  2. sRNARFTarget's output files are saved in the folder 'sRNARFTargetResult' which is created in the working directory. This folder will contain two files: Prediction\_probabilities.csv and FeatureFile.csv.
  * Prediction_probabilities.csv: this is the main result file and contains results sorted by predicted interaction probability from high to low, rounded to five decimals. It contains three columns, sRNA_ID, mRNA_ID and Prediction\_Probability. Here are some lines of a Prediction\_probabilities.csv file generated:
  ```
sRNA_ID mRNA_ID Prediction_Probability
gcvb    PM0494(+)       0.57444
gcvb    PM_RS03970(-)   0.55257
gcvb    PM_RS00560(-)   0.55193
gcvb    PM_RS06810(-)   0.54968
gcvb    PM_RS02870(+)   0.54926
gcvb    PM_RS00565(-)   0.54756
  ```
  * FeatureFile.csv: this file contains features for all the sRNA-mRNA pairs. This file consists of 66 columns. The first two columns are sRNA_ID and mRNA_ID. The remaining 64 columns are corresponding trinucleotide frequency difference of sRNA-mRNA pairs. This file is later used by sRNARFTarget interpretability scripts.

  ## Interpretation of sRNARFTarget Predictions 
  
  1. We created two python scripts for understanding the predictions generated by sRNARFTarget: sRNARFTarget_SHAP.py and sRNARFTarget_CP.py
  2. You need to run sRNARFTarget first so that the Prediction\_probabilities.csv and FeatureFile.csv files are generated.
     
   #### Instructions to run sRNARFTarget_SHAP
   
   1. Choose an sRNA-mRNA pair of interest from Prediction\_probabilities.csv file.
   2. OPTION A (with our [Python 3.8 docker container](https://hub.docker.com/r/penacastillolab/python38env)). After pulling the docker container, run the docker container to execute sRNARFTarget_SHAP.py as shown below.     The -v command indicates to docker that the folder "/ABSOLUTE_PATH_TO/sRNARFTarget" will be referred to as "/data" in the command.

   ```
   docker run -i -v /ABSOLUTE_PATH_TO/sRNARFTarget:/data --rm penacastillolab/python38env python /data/sRNARFTarget_SHAP.py '/data/sRNARFTargetResult/' 'gcvb' 'PM0494(+)'
   ```
   
    
   2. OPTION B (with Python installed locally). Run sRNARFTarget_SHAP using the below command.
   ```
    python sRNARFTarget_SHAP.py 'PATH_TO_FeatureTable.csv' 'sRNA_ID' 'mRNA_ID'
   ``` 
      
   Example usage: python sRNARFTarget_SHAP.py 'sRNARFTargetResult/' 'omrA' 'ompT'
   
   Make sure to use single quotations around the IDs and write the IDs exactly as they appear in the Prediction probabilities.csv.
   
   3. sRNARFTarget_SHAP will create a decisionPlot.pdf file and ForcePlot.html file
     
   #### Instructions to run sRNARFTarget_CP
   
   1. For the same sRNA-mRNA pair that was chosen to run sRNARFTarget_SHAP, choose a feature/variable by looking at the plots generated by sRNARFTarget_SHAP or any one variable of interest.
  2. OPTION A (with our [Python 3.8 docker container](https://hub.docker.com/r/penacastillolab/python38env)). After pulling the docker container, run the docker container to execute sRNARFTarget_CP as shown below.   The -v command indicates to docker that the folder "/ABSOLUTE_PATH_TO/sRNARFTarget" will be referred to as "/data" in the command.

  ```
  docker run -i -v /ABSOLUTE_PATH_TO/sRNARFTarget:/data --rm penacastillolab/python38env python /data/sRNARFTarget_CP.py '/data/sRNARFTargetResult/' 'gcvb' 'PM0494(+)' 'TTA'
  ```
    
  2. OPTION B (with Python installed locally). Run sRNARFTarget_CP using the below command.
  ```
  python sRNARFTarget_CP.py 'PATH_TO_FeatureTable.csv' 'sRNA_ID' 'mRNA_ID' 'feature_name'
  ```
         
  Example usage: python sRNARFTarget_CP.py 'sRNARFTargetResult/' 'omrA' 'ompT' 'GCG'
  
  Make sure to use single quotations around each parameter, and write the sRNA and mRNA ID exactly as they appear in the Prediction probabilities.csv file.
  
  3. sRNARFTarget_CP.py will create a directory called \_plots_files and will open the file plots0.html automatically in the default web browser. This file contains an interactive plot showing the predicted interaction probability as a function of the value of the feature provided.
      
   ## Notes
1. sRNARFTarget can be run for any number of sRNAs and mRNAs at a time.
2. sRNARFTarget_SHAP program can only be run for one sRNA-mRNA pair at a time.
3. sRNARFTarget_CP can only be executed for single sRNA-mRNA pair and a single feature/variable at a time.
4. Make sure the folder provided to  sRNARFTarget_CP.py and sRNARFTarget_SHAP.py as their first argument contains the files Prediction_probabilities.csv and FeatureFile.csv generated by sRNARFTarget.nf.

   ## Citation
   If you use this software please cite:
   
Kratika Naskulwar & Lourdes Peña-Castillo (2022) sRNARFTarget: a fast machine-learning-based approach for transcriptome-wide sRNA target prediction, RNA Biology, 19:1, 44-54, [DOI: 10.1080/15476286.2021.2012058](https://doi.org/10.1080/15476286.2021.2012058)

