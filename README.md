# individual-project
This is github page includes the codes I used to implement my individual project. In general, the codes implements the following functions:
fastaqc(https://github.com/s-andrews/FastQC) is ued for quality control

Getorganelle(https://github.com/Kinggerm/GetOrganelle): The chloroplast genome of Cochlearia excelsa was assembled for reference

quast（https://github.com/ablab/quast）is a genome assembly evaluation tool

Minimap2(https://github.com/lh3/minimap2) is used for alignment betweeen reference genome and query_genome

samtools(https://github.com/samtools) is used to convert gene sequence formats

bcftools（https://github.com/samtools/bcftools) is ued for SNPs calling

vcftools（https://github.com/vcftools/vcftools）extracts SNPs and generate FASTA file

snp-site (https://github.com/sanger-pathogens/snp-sites) is used for finding SNP sites from a multi-FASTA alignment file

RAxML (https://github.com/amkozlov/raxml-ng) is the tools used for phylogenetic tree inference



all these bioinformatics tools mentioned above could be installed by conda use following codes:
eg. conda install -c bioconda bcftools
conda install -c bioconda samtools
conda install -c bioconda getorganelle
etc


