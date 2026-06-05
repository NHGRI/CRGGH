<img align="right" width="200" height="200" src="https://github.com/user-attachments/assets/5361ef5a-1540-4e73-b555-75a4441ebb3c"> <br>
# Skeletal Muscle Single Nucleus RNA Sequencing in Type 2 Diabetes Shows Transcriptional Heterogeneity and Differential Trait Enrichment 

---------
## Authors:
Elisabeth F. Heuston1, Ayo P. Doumatey1, Lin Lei1, Jie Zhou1, Faiza Naz2, Shamima Islam2, Clement A. Adebamowo3, Stefania Dell’Orso2, Adebowale A. Adeyemo1, Charles N. Rotimi1

1. Center for Research on Genomics & Global Health, National Human Genome Research Institute, National Institutes of Health Bethesda, MD, USA; 
2. Genomic Technology Section, National Institute of Arthritis and Musculoskeletal and Skin Diseases, National Institutes of Health, Bethesda, MD, USA; 
3. University of Maryland School of Medicine, Baltimore, MD, USA

*Corresponding author:*  
>Charles N. Rotimi, PhD  
Center for Research on Genomics & Global Health,   
National Human Genome Research Institute,   
National Institutes of Health,   
Bethesda, MD 20892  
rotimic@mail.nih.gov    

## Abstract
Single cell studies have yielded new insights into the molecular pathophysiology of type 2 diabetes (T2D) with most focusing on pancreatic islet cell populations. To investigate how other tissues that are important in glucose metabolism may contribute to this growing understanding, we undertook single nucleus RNA sequencing of human skeletal muscle obtained from West Africans using the 10X Genomics platform. By analyzing 91,454 high quality nuclei we identified 21 transcriptional clusters, of which ten clusters (including myocyte progenitor, endothelial, smooth muscle, and neuronal cell nuclei) showed differential nuclei abundance between individuals with T2D (T2D+) compared to those without diabetes (T2D-). Four myocyte clusters displayed an upregulated RNA profile associated with Mitochondrial Dysfunction, Sirtuin Signaling Pathway, and Granzyme A Signaling as well as a downregulated RNA profile associated with Oxidative Phosphorylation, Respiratory Electron Transport and Neutrophil Extracellular Trap Signaling. Pseudotime analysis revealed that T2D+ nuclei showed different trajectory patterns in comparison with T2D- nuclei, implying cellular developmental differences associated with T2D and consistent with observed differential abundance in myocyte progenitor sub-clusters between the two groups. Transcriptional profiles of differentially abundant clusters revealed enrichment for genome-wide significant loci from GWAS of T2D and related traits, with varying patterns of cluster-GWAS enrichment profiles. Our findings of differential abundance, differential gene expression profiles and differences in cell trajectories highlight the T2D-associated molecular and pathway changes observed in skeletal muscle at single cell resolution. 

--------

## Repo structure

```
├── RFunctions/                        # helper functions
│   ├── ColorPallete.R                 # standard color palette
│   ├── AssignMetadata.R               # assign metadata during seurat object creation
│   ├── PercentVariance.R              # calculate percent variance to determine dimensionality
│   ├── runDoubletFinder.R             # calls DoubletFinder
│   ├── xlsx.tablefy.R                 # generates sorted and colored excel tables from `FindAllMarkers()`
│   ├── enrichR_fucntions.R            # call and process enrichR analysis data
│   ├── ipa.networks.R                 # subset list of IPA networks
│   └── dge_listToTable.R              # convert list to table
│   └── custom_radarchart-fmsb.R       # based on fmsb code to modify radarcharts
├── PrimaryAnalysis/                   # processing data and generating objects
│   ├── cellranger_swarm_generator.sh  # demultiplexing and alignment
│   ├── cellbender_swarm_generator.sh  # ambient RNA removal
│   ├── setup.R                        # called at beginning of scripts
│   ├── preprocess_seurat.R            # original processing pipeline (includes all data)
│   ├── skmus_seurat.R                 # secondary processing pipeline to remove debris
│   └── monocle_trajectory_analysis.R  # generate CellDataSet object
├── SecondaryAnalysis/                 # work from objects generated in PrimaryAnalysis
│   └── DifferentialAbundanceTesting.R # perform differential abundance testing
│   └── enrichR.R                      # run GWAS_Catalog_2025 enrichment on `FindAllMarkers()`
│   └── pseudobulk.R                   # run `AggregateExpression()` to compare clusters and T2D vs no T2D per cluster
│   └── garnett.R                      # commands to perform cell type annotation with Garnett
│   └── IPA.R                          # command to generate gene lists for IPA input and to reformat IPA output
│   └── Figures.R                      # generate figures in the manuscript
└── README.md                          # this file
```
## Notes
* All file paths have been removed or replaced with relative paths
* No identifying information is contained within these files
--------
## Contact
For questions regarding this repository, please contact:  
>Elisabeth F. Heuston, PhD  
*Center for Research on Genomics & Global Health, NHGRI, NIH  
elisabeth.heuston@nih.gov*  
