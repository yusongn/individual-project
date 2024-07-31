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
all these bioinformatics tools mentioned above could be installed by conda use following codes
eg. conda install -c bioconda bcftools
conda install -c bioconda samtools
conda install -c bioconda getorganelle
etc


##quality_control
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=100g
#SBATCH --time=6-08:08:08
#SBATCH --output=fastqc_out_%j.out
#SBATCH `--error=fastqc_error_%j.err```

fastqc /gpfs01/home/mbzra2/yusong/ref_excelsa/* -o output_folder

##getorganelle
#!/bin/bash
#SBATCG --job-name=ad_getorg
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=300g
#SBATCH --time=6-08:08:08
#SBATCH --output=ad_getorg_pt_output_%j.out
#SBATCH --error=ad_getorg_pt_error_%j.err

get_organelle_from_reads.py -1 /gpfs01/home/mbzra2/yusong/ref_excelsa/Cochlearia_excelsa_linked_no_idx_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq -2 /gpfs0$

##quast
#!/bin/bash
#SBATCG --job-name=quast
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=300g
#SBATCH --time=6-08:08:08
#SBATCH --output=quast_pt_output_%j.out
#SBATCH --error=quast_pt_error_%j.err

#run quast for assembled chloroplast genome
quast.py -o quast_pt_output /gpfs01/home/mbxyn1/ad_Org_pt_Assembly/embplant_pt.K105.complete.graph1.1.path_sequence.fasta

script used for genome alighment and SNPs calling
#This script run the 107 query_folder under "trim" file on HPC

#!/bin/bash
#SBATCH --job-name=pt_alignment
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=1400g
#SBATCH --time=6-08:08:08
#SBATCH --output=pt_align_out_%j.out
#SBATCH --error=pt_align_error_%j.err

set -e

Set the reference sequence file path
reference="/gpfs01/home/mbxyn1/ad_Org_pt_Assembly/embplant_pt.K105.complete.graph1.1.path_sequence.fasta"

Create final and intermediate output directories
mkdir -p "final_output"
mkdir -p "output"

Function to shorten contig names in VCF files
shorten_contig_names() {
input_vcf=$1
output_vcf=$2

bcftools annotate --rename-chrs <(
    grep "^>" "$reference" | awk -F'[>| ]' '{print $2"\t"$2}' | nl -v 1 | awk '{print $2"\t"NR}' 
) "$input_vcf" -Oz -o "$output_vcf"

tabix -p vcf "$output_vcf"
}

Process all query folders
for query_folder in /gpfs01/home/mbzra2/yusong/trim/*; do
if [ -d "$query_folder" ]; then
output_folder=$(basename "$query_folder")
mkdir -p "output/${output_folder}"

    echo "Processing $query_folder..."

    # Process each FASTQ file
    for fastq_file in ${query_folder}/*.fastq.gz; do
        file_name=$(basename "$fastq_file" .fastq.gz)
        
        # Extracting chloroplast genome, Aligning with reference genome using Minimap2 and converting SAM to BAM, sorting and indexing BAM files.
        minimap2 -ax map-ont "$reference" "$fastq_file" | \
        samtools view -bS -F 4 - | \
        samtools sort -o "output/${output_folder}/${file_name}.sorted.bam"
        
        samtools index "output/${output_folder}/${file_name}.sorted.bam"
        
        # SNP（Single  Nucleotide  Polymorphisms） calling using bcftools
        bcftools mpileup -Ou -f "$reference" "output/${output_folder}/${file_name}.sorted.bam" | \
        bcftools call -mv -Oz -o "output/${output_folder}/${file_name}.vcf.gz"
        
        # Index the VCF file using tabix
        tabix -p vcf "output/${output_folder}/${file_name}.vcf.gz"

        # Filter the VCF file for high quality SNPs
        bcftools filter -i 'QUAL>30 && DP>10' "output/${output_folder}/${file_name}.vcf.gz" -Oz -o "output/${output_folder}/${file_name}.filtered.vcf.gz"

        # Index the filtered VCF file using tabix
        tabix -p vcf "output/${output_folder}/${file_name}.filtered.vcf.gz"

        # Shorten contig names in the filtered VCF file
        shorten_contig_names "output/${output_folder}/${file_name}.filtered.vcf.gz" "output/${output_folder}/${file_name}.filtered.short.vcf.gz"

        # Delete the BAM file and its index file to save disk space
        rm "output/${output_folder}/${file_name}.sorted.bam"
        rm "output/${output_folder}/${file_name}.sorted.bam.bai"
        rm "output/${output_folder}/${file_name}.vcf.gz"
        rm "output/${output_folder}/${file_name}.filtered.vcf.gz"
    done

    # Combine all filtered VCF files within the same folder
    bcftools merge --force-samples "output/${output_folder}"/*.filtered.short.vcf.gz -Oz -o "output/${output_folder}/combined.filtered.short.vcf.gz"

    # Index the combined VCF file using tabix
    tabix -p vcf "output/${output_folder}/combined.filtered.short.vcf.gz"


    echo "$query_folder processing complete."
fi
done

Combine all combined filtered VCF files from all folders
bcftools merge --force-samples "output"/*/combined.filtered.short.vcf.gz -Oz -o "final_output/all_combined.filtered.short.vcf.gz"

Index the final combined VCF file using tabix
tabix -p vcf "final_output/all_combined.filtered.short.vcf.gz"

Extract SNPs and generate FASTA file using vcftools
vcftools --gzvcf "final_output/all_combined.filtered.short.vcf.gz" --out "final_output/snps" --maf 0.05 --recode --recode-INFO-all
vcftools --vcf "final_output/snps.recode.vcf" --out "final_output/snps" --extract-FORMAT-info GT

Generate MSA FASTA file using snp-sites
snp-sites -v "final_output/snps.recode.vcf" -o "final_output/alignment.fasta"

#In addition, use similar method on outgroup ionopsiduiom data, still using chloroplast genome assembled before as reference,obtain the SNPs（Single Nucleotide Polymorphisms） multi-sequence alignment data
#!/bin/bash
#SBATCH --job-name=outgroup_alignment
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=1000g
#SBATCH --time=6-08:08:08
#SBATCH --output=pt_outgroup_align_out_%j.out
#SBATCH --error=pt_outgroup_align_error_%j.err

set -e # Exit on error

Set the reference sequence file path
reference="/gpfs01/home/mbxyn1/ad_Org_pt_Assembly/embplant_pt.K105.complete.graph1.1.path_sequence.fasta"

Create final and intermediate output directories
mkdir -p "final_output"
mkdir -p "output"

Function to shorten contig names in VCF files
shorten_contig_names() {
input_vcf=$1
output_vcf=$2

echo "Shortening contig names for $input_vcf..."
bcftools annotate --rename-chrs <(
    grep "^>" "$reference" | awk -F'[>| ]' '{print $2"\t"$2}' | nl -v 1 | awk '{print $2"\t"NR}' 
) "$input_vcf" -Oz -o "$output_vcf"

tabix -p vcf "$output_vcf"
echo "Finished shortening contig names."
}

Process all query folders
for query_folder in /gpfs01/home/mbzra2/yusong/ionopsiduiom/*; do
if [ -d "$query_folder" ]; then
output_folder=$(basename "$query_folder")
mkdir -p "output/${output_folder}"

    echo "Processing $query_folder..."

    # Process each FASTQ file
    for fastq_file in ${query_folder}/*.fq.gz; do
        file_name=$(basename "$fastq_file" .fq.gz)
        
        echo "Aligning and processing $fastq_file..."
        # Aligning using Minimap2 and converting SAM to BAM, sorting and indexing BAM files.
        bwa mem -t 10 "$reference" "$fastq_file" | \
        samtools view -bS -F 4 - | \
        samtools sort -o "output/${output_folder}/${file_name}.sorted.bam"
        
        samtools index "output/${output_folder}/${file_name}.sorted.bam"
        
        echo "Calling SNPs for $fastq_file..."
        # SNP calling using bcftools with annotation flags
        bcftools mpileup -Ou -f "$reference" "output/${output_folder}/${file_name}.sorted.bam" | \
        bcftools call -vm -Oz -o "output/${output_folder}/${file_name}.vcf.gz" 
        
        # Index the VCF file using tabix
        tabix -p vcf "output/${output_folder}/${file_name}.vcf.gz"

        echo "Filtering SNPs for quality..."

          # Filter the VCF file for high quality SNPs
        bcftools filter -i 'QUAL>20 && DP>5' "output/${output_folder}/${file_name}.vcf.gz" -Oz -o "output/${output_folder}/${file_name}.filtered.vcf.gz"

        # Index the filtered VCF file using tabix
        tabix -p vcf "output/${output_folder}/${file_name}.filtered.vcf.gz"

        # Shorten contig names in the filtered VCF file
        shorten_contig_names "output/${output_folder}/${file_name}.filtered.vcf.gz" "output/${output_folder}/${file_name}.filtered.short.vcf.gz"

        # Delete the BAM file and its index file to save disk space
        rm "output/${output_folder}/${file_name}.sorted.bam"
        rm "output/${output_folder}/${file_name}.sorted.bam.bai"
        rm "output/${output_folder}/${file_name}.vcf.gz"
        rm "output/${output_folder}/${file_name}.filtered.vcf.gz"
    done

    # Combine all filtered VCF files within the same folder
    bcftools merge --force-samples "output/${output_folder}"/*.filtered.short.vcf.gz -Oz -o "output/${output_folder}/combined.filtered.short.vcf.gz"

    
    # Index the combined VCF file using tabix
    tabix -p vcf "output/${output_folder}/combined.filtered.short.vcf.gz"

    echo "$query_folder processing complete."
fi
done

Combine all combined filtered VCF files from all folders
bcftools merge --force-samples "output"/*/combined.filtered.short.vcf.gz -Oz -o "final_output/all_combined.filtered.short.vcf.gz"

Index the final combined VCF file using tabix
tabix -p vcf "final_output/all_combined.filtered.short.vcf.gz"

Extract SNPs and generate FASTA file using vcftools
vcftools --gzvcf "final_output/all_combined.filtered.short.vcf.gz" --out "final_output/snps" --maf 0.01 --recode --recode-INFO-all

Validate if the VCF file has any sites before extracting FORMAT info
if [ -s "final_output/snps.recode.vcf" ]; then
vcftools --vcf "final_output/snps.recode.vcf" --out "final_output/snps" --extract-FORMAT-info GT

# Generate MSA FASTA file using snp-sites
snp-sites -v "final_output/snps.recode.vcf" -o "final_output/alignment.fasta"
else
echo "No SNPs found after filtering. Check the VCF files and filtering criteria."
fi

##Then, we have two multi-sequence alignemnt fasta file. Merging them and we can run RAxML

#!/bin/bash
#SBATCG --job-name=pt_phyloTree
#SBATCH --partition=hmemq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=1000g
#SBATCH --time=6-08:08:08
#SBATCH --output=pt_phylotree_out_%j.out
#SBATCH --error=pt_phylotree_error_%j.err

Phylogentic tree construction using RAxML
using GTR model in this case for high performance
cat /gpfs01/home/mbxyn1/Minimap2/final_output/alignment.fasta /gpfs01/home/mbxyn1/Minimap2/pt_outgroup/final_output/alignment.fasta > merged.fasta

Build phylogenetic tree using RAxML
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s "/gpfs01/home/mbxyn1/phylogenetic_Analysis/merged.fasta" -n tree
