#Cellranger-7.0.0
#Using locally installed cellranger as the server version is outdated
#Using CellRanger pre-compiled Mouse reference refdata-gex-mm10-2020-A


cellranger count --id=$sample \
                 --transcriptome=<path\to\CellRanger\mouse\reference> \
                 --fastqs=<\path\to\fastq> \
                 --sample=$sample \
                 --localcores=40 \
                 --localmem=185
