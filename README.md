# Psilocybin.Female.Mouse.mPFC
This repo contains the code and data to reproduce analysis and figures for the manuscript titled "Single-cell sequencing reveals sustained cell-type specific plasticity in mouse medial prefrontal cortex".

Single cell RNA-Seq (scRNA-Seq) related analysis are contained in Code/Analysis folder.
Visualization of processed scRNA-Seq data and other results presented in the figures are contained in Code/Figures folder using preprocessed data in the Data folder.

## Platforms and software versions
R scripts were performed on Compute Canada Béluga cluster. The following line was used to load required modules:
```
module load  StdEnv/2020 geos/3.9.1 gcc r/4.2.1 hdf5
```
The R package versions can be found in the R_sessioninfo.txt
Analysis with CellBender was performed on Compute Canada Béluga cluster with Python virtual enrivonment.
Other Python scripts were executed using WinPython with package versions in WinPython.txt

Packages and dependencies required to extract the gene modules using supervised latent factor modeling (and all downstream analyses shown in Fig 3 and Extended Figures 6–8) and listed in ```python.package.version.txt``` . The conda environment can be recreated with ```psilocybin_mouse_pls.yml```.


## Analysis Pipeline
Single cell RNA-Seq analysis and R related visualization were performed using Compute Canada Béluga cluster.

Supervised latent factor modeling analyses, and mapping to the [Imputed MERFISH spatial transcriptomics of a single adult mouse brain](https://alleninstitute.github.io/abc_atlas_access/descriptions/MERFISH-C57BL6J-638850-imputed.html) were performed in Python 3.13.2 on a private server hosted by the McConnell Brain Imaging Centre (Montreal Neurological Institute).



### Preprocess
CellRanger 7.0.0 was used to generate count matrix and cell calling with pre-compiled mouse reference.

CellBender 0.3.0 was used to remove background RNA. Related package versions can be found in cellbender.package.version.txt.

Cells were filtered on the number of genes (greater than 800 and within two standard deviations from average), the number of total reads (within two standard deviations from average) and mitochondrial content (lower than 10%). In each sample, only genes detected in more than 3 cells were retained. After quality control, two saline-injected samples were pooled due to low cell number.

Cell type labeling was performed using Scibet, SingleR and SingleCellNet, with previously published data sets (cells collected from saline treated P60 animals published by Bhattacherjee 2019 & Allen Mouse Brain Atlas with region label ILA;PL;ORB) as references. A consensus annotation was assigned when at least two methods agreed. Cells were only assigned as neuron if labeled as neuron by both references. Neuron subtypes were labeled from the Allen Mouse Brain Atlas. The expression of cell type and neuron subtype markers was examined to confirm the cell type assignment.

The dataset was scaled to 10,000 reads per cell and log transformed. Samples were integrated by merging samples without batch correction. Variable features were identified with the variance-stabilizing transform method. Dimension reduction was performed with the 2000 most variable features and the top 30 dimensions were used for clustering and projection to two dimensions with Uniform Manifold Approximation and Projection.

### Cell Prioritization with Augur
Augur was used to measure the response level of each cell type / subtype to psilocybin injection.
A parameter sweep with vquant={0.5,0.7,0.9}, ntree={50,100,150}, mtry={10,15,20} was performed to determine the best combination.
Lower vquant in general produced lowest scores and takes more computational time. Higher ntree and mtry produced higher scores.
The combination of vquant=0.9, ntree=150, mtry=20 produced best score and the improvement was plateaued comparing to the other combinations. Further analysis were performed with vquant=0.9, ntree=150, mtry=20.
The entire dataset was used to establish the overall response. The cell condition labels were then shuffled (with Augur "permute" option) to determine whether the Augur score was significantly different above chance.
Neurons expressing serotonin receptors were then examined by subsetting the neurons according to their serotonin receptor expression.


### Cell-Cell Communication with CellChat and NeuronChat
CellChat was used to estimate cell-cell communication in either control group or psilocybin group. The numbers of interactions identified were then compared. Most difference in cell-cell communication were found in neurons. Thus, we used NeuronChat to further identify cell-cell communications in the neurons in a similar manner, but with a model and database optimized for neurons.

### Differential Gene Expression with Libra
Libra is run with the following setting:
```
de_family = 'mixedmodel', de_method = 'negbinom_offset', de_type = 'Wald'
```
The mixed model with Wald test provides a good compromise between the stat power of single cell methods and the sample variation of pseudobulk methods. However many this method produced genes with close to zero Log2 Fold Changes have a stong p value, therefore a filtering of L2FC>0.1 is performed using python.

### Gene Set Enrichment with clusterProfiler
Gene symbols were first converted to Entrez ID using DAVID conversion tool at 
https://david.ncifcrf.gov/conversion.jsp
Genes without converted Entrez ID were discarded.
GSE analysis was then performed using clusterProfiler.
The ancestors of the enriched Gene Ontology terms were obtained using EMBL-EBI's API at 
https://www.ebi.ac.uk/QuickGO/api/index.html
and aggregated for the top 50 of either positively or negatively enriched terms in each cell type or subtype. Representative ancestor terms were manually selected for the heatmap.
 

## Visualization
Seurat related figures such as UMAP and violin plots are performed on Compute Canada Béluga.
Python scripts were executed using WinPython with package versions in WinPython.txt

## Biorender
Portions of Figures 1, 3 and Supp. Figures S4 were created with Biorender. These graphics were not included in this repo.
