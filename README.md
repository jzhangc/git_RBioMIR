# RBioMIR
miRNA-seq analysis package

To cite in publication
  
    Zhang J, Hadj-Moussa H, Storey KB. 2016. Current progress of high-throughput microRNA differential expression analysis and random forest gene selection for model and non-model systems: an R implementation. J Integr Bioinform. 13: 306.
    

NOTE:

  - Currently the package is for Unix-like system only (e.g., Linux or macOS).
  - (for testing only) For novel miRNA discovery, be sure to copy training.smir and test.smir files to your working directory, and install pcregrep and parallel commands.


Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        source("https://bioconductor.org/biocLite.R")
      
        biocLite()
    
  - Install stable release
  
        devtools::install_github("jzhangc/git_RBioMIR/RBioMIR", repos = BiocInstaller::biocinstallRepos())    

  - Install development build
  
        devtools::install_github("jzhangc/git_RBioArray/RBioArray", repos = BiocInstaller::biocinstallRepos(), ref = "beta")  


Update log

    0.1.8 (Feature preview)
    (ICEBOX)
      - mirProcess() re-written
      - A preview of UNIX Shell program mirna_processing added
      - Other updates
        - Information page updated
