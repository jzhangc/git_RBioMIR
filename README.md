# RBioMIR
miRNA-seq analysis package

To cite in publication
  
    Zhang J, Hadj-Moussa H, Storey KB. 2016. Current progress of high-throughput microRNA differential expression analysis and random forest gene selection for model and non-model systems: an R implementation. J Integr Bioinform. 13: 306.
    

NOTE:
  - (For testing only) For novel miRNA discovery, be sure to copy training.smir and test.smir files to your working directory, and install pcregrep and parallel commands.
  - Small RNA RNAseq file processing program mirna_processing.sh only supports UNIX/UNIX-like systems. 

Installation:

  - Install devtools
  
        install.packages("devtools")
    
  - Install bioconductor
  
        if (!requireNamespace("BiocManager"))
            install.packages("BiocManager")
            
        BiocManager::install()
    
  - Install stable release
  
        devtools::install_github("jzhangc/git_RBioMIR/RBioMIR", repos = BiocManager::repositories())    

  - Install development build
  
        devtools::install_github("jzhangc/git_RBioMIR/RBioMIR", repos = BiocManager::repositories(), ref = "beta")  
        

Update log

    0.2.7 (Mar.6.2022)
      - mircount(s4) and mir_count now have a sample_groups_var_name item

    0.2.6 (Jan.17.2022)
      - Small typo fixes

    0.2.5 (Dec.23.2021)
      - S4 class "mircount" added as an improved miRNA read class over the S3 class "mir_count"

    0.2.4 (Dec.7.2021)
      - New items added to the "mir_count" object
        - working_gene_annot_var_name: for the object derived from mirProcess(), the value for this is always "mirna"
        - target_annotation_file_processed

      - mirDeepProcess() function added to process conserved miRNA count resutls from miRDeep2
      - The "mir_count" object updated with "genes_complete_annotation" data frame
        - The print function for "mir_count" updated accordingly
      - The "mir_count" object's "raw_read_count" now uses gene names as row names

    0.2.3 (Jan.17.2019)
      - New bioconductor installation instructions added
      
      
    0.2.2
      - "files_processed" added to "mir_count" object
      - A bug fixed for mirProcess function where the function will improperly check target.annot.file even if it's set at NULL
    
    
    0.2.1
      - Adjustments made to "mir_count" objects
        - The object now outputs genes and counts separately
        - The object now outputs library size for the samples
        
      - Bug fixes


    0.2.0
      - New functions
        - S3 print function for "mir_count"
          
      - Updates to exisiting functions:
        - mirProcess() re-written, with:
          - simplied code base
          - foreach implementation of parallel computing
          - compatibility for UNIX/UNIX-like/Windows systems
          - non-mandatory naming convention for input .txt read count files
          - ability to take annotation file
          - output in S3 object "mir_count"
              
      - A preview of UNIX Shell program mirna_processing added
      - Other updates
          - Information page updated
