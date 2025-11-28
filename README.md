# bottleneck-passeriformes
Repository containing scripts for quality control analysis reads, mapping, SNPs identification, and population analysis.
-----------------------------------------------------------------------------------------------
# QUALITY CONTROL OF READS
fastqc -t 30 seq.fastq.gz

# Comparative control between sequences of the same sample
multiqc .  # inside the folder with all sequences

# Filtering the reads
for file in *.fastq.gz do
	output_file = "filtrado_q20_${file}"
	log_file = "log_{$file%fastq.gz}.txt
	cutadapt -q 20 -o "output_file" "${file}" 2>&1 -t 30
	echo "Process finished for ${file}. Log salvo em $log_file"
	done
	echo "All the files were processed."
  
------------------------------------------------------------------------------------------------
# MAPPING AND ALIGNEMENT

# Creating index
bwa index /adress/of/genome/reference.fasta

# Combining the R1 sequences of the same individual, with different runs
cat seq_L003_R1.fastq.gz seq_L004_r1.fastq.gz > seq_R1.fastq.gz # Do the same with R2

# Mapping the reference reads and creating sam file
bwa mem -t 64 /adress/of/ref/genome/with/index.fasta
seq_R1.fastq.gz seq_R2.fastq.gz | gzip -3 > seq.sam.gz

# Converting sam file in bam file
samtools view -h -b -@ 64 sample.sam > sample.bam

# Removing reads not paired
samtools view -b -F 4 -@ 64 sample.bam > sample_unflags.bam

# Classifying the file
samtools sort -@ 64 sample_unflags.bam -o sample_sorted.bam 

# Indexing  for individuals that have the same position in reads
samtools index sample_sorted.bam 

# Quality control of Alignment
samtools flagstat sample_sorted.bam

# Removing duplicates
  # Step 1
  java -jar /home/cmarinho/picard.jar MarkDuplicates \
    I=/adress/sample_sorted.bam \
    O=/adress/sample_withoutdupl.bam \
    M=/adress/sample_marked_dup_metrics.txt \
    REMOVE_DUPLICATES=true
    
  # if there is any error
  samtools view -H /endereco/de/sample_sorted.bam | grep '@RG'
  # if it does not return nothing, realize the following step:
  java -jar /home/cmarinho/picard.jar AddOrReplaceReadGroups \
    I=/adress/sample_sorted.bam \
    O=/adress/sample_sorted_RG.bam \
    RGID=sample_name \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=sample_name
      # After, return to step 1, with sample_sorted_RG.bam
      # About file.metrics.txt --> if number of duplicates <= 20%, it is acceptable; between 30-40% is moderate; >50% is high.
      # don't forget to index the final file
      samtools index sample_sorted_RG.bam

-------------------------------------------------------------------------------------------------
# Identifying SNPS

# Indexing the reference genome # it is done only one time
samtools faidx /adress/refseq.fasta

# Creating a reference dictionary #also, only one time
gatk CreatingSequenceDictionary -R refseq.fasta

# Happlotype calling
gatk --java-options "-Xmx100g" HaplotypeCaller \
    -R reference.fasta \
    -I sample.bam \
    -O sample.vcf \
    -ERC GVCF
      # This step spends a long time, around 5-9 days with my data

# Combining all gvcf file in only one file
gatk --java-options "-Xmx100g" CombineGVCFs -R fererence.fasta -V individuo1.gvcf -V individuo2.gvcf -V individuo3.gvcf -O output_combined.gvcf
	# através do gatk, comando Combine GVCFs, todos os arquivos dos indivíduos se unirão em um único. ps.: "-Xmx100g" --> não ultrapassar 100Gb de RAM

  # To validate the GVCFs:
  gatk ValidateVariants -V file.gvcf
  # If there is any error:
  gatk --java-options "-Xmx100g" LeftAlignAndTrimVariants \
  -R /adress/refseq.fna \
  -V /adress/sample.gvcf \
  -O /adress/sample_fixed.gvcf \
  > /adress/sample_log.txt \
  2> /adress/sample_errors.txt

# Genotyping
gatk --java-options "-Xmx100g" GenotypeGVCFs -R reference.fasta -V input_combined.gvcf -O output.vcf

# Compaction and indexing vcf
bgzip output.vcf

# Converting vcf in bcf file
bcftools view input.vcf.gz -Ob -o output.bcf.gz

# Statistics of bcf file
bcftools stats -s - input.bcf.gz > stats.txt

# Hardfilter
bcftools filter -e "QD<2 || MQRankSum<-12.5 || FS>60 || SOR>3 || ReadPosRankSum<-8 || QUAL<20 || MQ<20 || MAF<0.05 || MEAN(FORMAT/DP)<0.8 || MEAN(FORMAT/DP)>50 || F_MISSING>0.5" --SnpGap 10 input.bcf.gz -Ob -o parcial.bcf.gz

# Extracting bialelic SNPS
bcftools view -c1 -v snps -m2 -M2 PARCIAL.bcf.gz -Ov -o PARCIAL.vcf

# Excluding SNPS which desviates significativally of HWE (p<0.001)
vcftools --gzvcf parcial.vcf.gz \
  --hwe 0.001 \
  --recode --recode-INFO-all \
  --out hwe_parcial

  # Counting the quantity of SNPS before and after HWE filter:
  bcftools view -v snps -H parcial.vcf.gz | wc -l # before
  grep -v "^#" hwe_parcial.recode.vcf | wc -l # after

# Removing bariants with very low frequency
vcfutils.pl varFilter -1 0.000001 hwe_parcial.recode.vcf > species_final.vcf

# Creating a file with statistics about SNPS, Indels, allelic frequencies, genotypes and coverage
bcftools stats -s - species_final.vcf > stats_hardfilters.txt

# Average mapping coverage per sample
samtools coverage -q 20 - Q 20 file.bam

-----------------------------------------------------------------------------------------------
# Performing Principal Component Analysis (PCA)

# Installing PLINK
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip -O plink.zip
unzip plink.zip -d .

# Adding to Path
echo 'export PATH=$HOME/plink:$PATH' >> ~/.bashrc
source ~/.bashrc

# Creating files for PCA analysis using vcf file with no pattern chromossomes
plink --vcf species_final.vcf --pca --out chiroxiphia --mind 1.0 --allow-extra-chr

# <ON ENVIRONMENT>
# Loading libraries
library(ggplot2)
library(stringr)
library(viridis)

# 1. Load and process engenvalues
Val <- read.table("/home/cmarinho/Chiroxiphia_caudata_genome/hwe/chiroxiphia.eigenval", header = FALSE)
colnames(Val) <- c("eigenvalue")  
Val$PC <- 1:nrow(Val)  
Val$percent <- (Val$eigenvalue / sum(Val$eigenvalue)) * 100  

# 2. Plotting eigenvalues (Scree Plot)
ggplot(Val[1:10,], aes(x = PC, y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
  labs(y= "% Variance", x = "Principal Component", title = "Scree Plot") +
  scale_y_continuous(limits = c(0, max(Val$percent))) + 
  theme_minimal()

# 3. Load and process eigenvectors
Vec <- read.table("/home/cmarinho/Chiroxiphia_caudata_genome/hwe/chiroxiphia.eigenvec", header = FALSE)  
colnames(Vec) <- c("POP", "Sample", paste0("PCA", 1:(ncol(Vec)-2)))  
Vec$POP <- as.factor(Vec$POP)

# Choose the dynamic color pallet
num_pop <- length(unique(Vec$POP))
color_palette <- viridis(num_pop)

# 4. Plot PCA1 vs. PCA2
p1 <- ggplot(Vec, aes(x= PCA1, y= PCA2, color = POP)) +  
  geom_point(size=4, alpha = 0.8) +  
  theme_classic() +  
  scale_color_manual(values = color_palette) +
  labs(x="PC1", y="PC2", title="PCA - PC1 vs PC2") +
  theme(legend.position="right")
print(p1)

# 5. Plot PCA2 vs. PCA3
p2 <- ggplot(Vec, aes(x= PCA2, y= PCA3, color = POP)) +
  geom_point(size=4, alpha = 0.8) +
  theme_classic() +
  scale_color_manual(values = color_palette) +
  labs(x="PC2", y="PC3", title="PCA - PC2 vs PC3") +
  theme(legend.position="right")
print(p2)

# 6. Calculate % explained variation
pve <- data.frame(PC = 1:nrow(Val), pve = Val$percent)

------------------------------------------------------------------------------------------------
# CALCULATING HETEROZIGOSITY PER INVIDIVUAL AFTER HARDFILTER FOR SNPS VIA PLINK

# 1) Convert the VCF to binary:
plink --vcf chiroxiphia_final.vcf.gz --allow-extra-chr --make-bed --out chiroxiphia

# 2) Calculate heterozygosity per individual from the newly generated file:
plink --bfile chiroxiphia --double-id --het --out het_stats

  # het_stats.log, .het and .nosex were generated

# Evaluating the .het: heterozygosity of animals post-accident greater than pre-accident (excess of heterozygotes) (F<0). BUT N(NM) (very variable coverage difference - may bias results).
#Redoing the --het in sites common to all and with minimum QC:

# 1) VCF with biallelic SNPs A/C/G/T and useful tags (AF):
bcftools view -m2 -M2 -v snps chiroxiphia_final.vcf.gz \
| bcftools +fill-tags -Ou -- -t AF \
| bcftools view -Oz -o chiroxiphia.snps.vcf.gz
tabix -p vcf chiroxiphia.snps.vcf.gz

# 2) Run PLINK requiring everyone genotyped at the site (--geno 0) and accepting scaffolds:
plink --vcf chiroxiphia.snps.vcf.gz --allow-extra-chr \ 
--snps-only just-acgt --biallelic-only strict\ 
--geno 0 --double-id --het --out het_allpresent 

# 3) Then separate by group (same set of sites) and recompute
# ps.: pre.txt and pos.txt are the FID and FII columns respectively with the data present in # het_allpresent.het for each pre- and post-accident sample.

# pre.txt and pos.txt with "FID IID" per line
plink --vcf chiroxiphia.snps.vcf.gz --allow-extra-chr \ 
--snps-only just-acgt --biallelic-only strict --geno 0 \ 
--keep pre.txt --double-id --het --out het_pre_comm

plink --vcf chiroxiphia.snps.vcf.gz --allow-extra-chr \ 
--snps-only just-acgt --biallelic-only strict --geno 0 \ 
--keep pos.txt --double-id --het --out het_pos_comm
Checking averages:
awk 'NR>1{n++;s+=$6}END{print "mean_F_pre_common=",s/n}' het_pre_comm.het
awk 'NR>1{n++;s+=$6}END{print "mean_F_pos_common=",s/n}' het_pos_comm.het

## <ON ENBIRONMENT>
## Mann-Whitney test to compare the heterozigosity difference found in the groups

# 1. Load the data

# Load the "Pre" file (ignoring the first header line)
data_pre <- read.table("het_pre_comm.het", skip = 1, header = FALSE)

# Load the "Post" file (ignoring the first header line)
data_pos <- read.table("het_pos_comm.het", skip = 1, header = FALSE)


# 2. Extract the F-statistics column (column 6)

# Column 6 (FIS) is extracted and stored in separate vectors
fis_pre <- data_pre[, 6]
fis_pos <- data_pos[, 6]

# Optional: Remove any 'Not a Number' (NaN) or 'Inf' values

# that may have been generated in the F_IS calculation for some loci
fis_pre <- fis_pre[!is.nan(fis_pre) & is.finite(fis_pre)]
fis_pos <- fis_pos[!is.nan(fis_pos) & is.finite(fis_pos)]

# 3. Running the Mann-Whitney Test

result_mann_whitney <- wilcox.test(fis_pre, fis_pos,
alternative = "two.sided")
print(result_mann_whitney)

------------------------------------------------------------------------------------------------
# Calculating pi (nucleotide diversity) by windows (100kb):

# pre (mfv)
vcftools --gzvcf chiroxiphia.snps.vcf.gz \
--keep pre.txt \
--window-pi 100000 \
--out pi_pre

# post (mar)
vcftools --gzvcf chiroxiphia.snps.vcf.gz \
--keep pos.txt \
--window-pi 100000 \
--out pi_pos

# Main output: *.windowned.pi (columns: CHR START END N_VARIANTS PI)
# <ON ENVIRONMENT>
# Summarizing (mean/median) and number of windows:

# pre
awk 'NR>1{a[++n]=$5; s+=$5} END{
asort(a); med=(n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2); 
printf "pre nwin=%d mean_pi=%.6g median_pi=%.6g\n", n, s/n, med
}' pi_pre.windowed.pi

# post
awk 'NR>1{a[++n]=$5; s+=$5} END{ 
assort(a); med=(n%2?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2); 
printf "pos nwin=%d mean_pi=%.6g median_pi=%.6g\n", n, s/n, med
}' pi_pos.windowed.pi

## Recommended: compare only the windows in common (same CHR: START:END)
# 1) Generate key CHR:START:END and PI
awk 'NR>1{printf "%s:%s:%s\t%f\n",$1,$2,$3,$5}' pi_pre.windowed.pi | sort -k1,1 > pre.pi.tsv
awk 'NR>1{printf "%s:%s:%s\t%f\n",$1,$2,$3,$5}' pi_pos.windowed.pi | sort -k1,1 > pos.pi.tsv

#2) Intersection of windows and calculation of paired means/medians
join -t $'\t' -j 1 pre.pi.tsv pos.pi.tsv > pi_common.tsv
# pi_common.tsv: KEY \t PI_pre \t PI_pos

#3) Quick Summary in Common Windows
awk '{ 
n++; pre+= $2; pos+= $3; diff+= ($3-$2); 
a[n]=$2; b[n]=$3; d[n]=($3-$2)
}
END{ 
# median helper 
function median(arr, n, i){asort(arr); return (n%2?arr[(n+1)/2]:(arr[n/2]+arr[n/2+1])/2)} 
printf "COMMON nwin=%d mean_pre=%.6g mean_pos=%.6g mean_delta(pos-pre)=%.6g\n", n, pre/n, pos/n, diff/n; 
printf "median_pre=%.6g median_pos=%.6g median_delta=%.6g\n", median(a,n), median(b,n), median(d,n);
}' pi_common.tsv

# Mean and mean delta (without median):
awk '{n++; pre+=$2; pos+=$3; diff+=($3-$2)} 
END{printf "COMMON nwin=%d mean_pre=%.6g mean_pos=%.6g mean_delta=%.6g\n", n, pre/n, pos/n, diff/n}' \ 
pi_common.tsv

# files with values ​​(pre, post, delta):
awk '{print $2}' pi_common.tsv | sort -n > pre_vals.txt
awk '{print $3}' pi_common.tsv | sort -n > pos_vals.txt
awk '{print $3-$2}' pi_common.tsv | sort -n > delta_vals.txt

# median function:
median() { awk '{a[NR]=$1} END{ if (NR%2) print a[(NR+1)/2]; else print (a[int(NR/2)] + a[int(NR/2)+1]) / 2 }' "$1"; }
printf "median_pre="; median pre_vals.txt
printf "median_pos="; median pos_vals.txt
printf "median_delta="; median delta_vals.txt

# fraction of windows with smaller pi in the post:
awk '{if($3<$2) c++} END{printf "frac(pos<pre)=%.3f (%d/%d)\n", c/NR, c, NR}' pi_common.tsv
------------------------------------------------------------------------------------------------
# Calculating Runs of Homozigosity (ROH) via PLINK 

# All individuals together
plink --vcf chiroxiphia.snps.vcf.gz \
--allow-extra-chr --double-id \
--snps-only just-acgt --biallelic-only strict \
--homozyg --out roh_all

# Generates roh_all.hom (list of each ROH segment (row = one ROH) - use to see long ROHs, positions, etc.)

roh_all.hom.indiv (summary per individual. Kb = sum of the length of all ROHs for that individual) - main for comparison
roh_all.hom.summary (global statistics)
roh_all.hom.log (execution log)
roh_all.nosex (warning of undefined sex)

# To view the data:

column -t roh_all.hom.indiv | less -S

#Result: Everything zeroed out. No segment long enough with standard parameters was found - common in scaffold genomes (short contigs break long ROHs). It's not an error, it's that the thresholds are too high.

# moderate set, no additional filters

plink --vcf chiroxiphia.snps.vcf.gz \
--allow-extra-chr --double-id \
--snps-only just-acgt --biallelic-only strict \
--homozyg \
--homozyg-snp 50 \
--homozyg-kb 500 \
--homozyg-density 50 \
--homozyg-gap 1000 \
--homozyg-window-snp 50 \
--homozyg-window-het 1 \
--homozyg-window-missing 1 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_mod 

#Result: everything reset- 

# Permissive set (for short scaff): 
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--homozyg\ 
--homozyg-snp 30\ 
--homozyg-kb 300 \ 
--homozyg-density 30 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 50 \ 
--homozyg-window-het 1 \ 
--homozyg-window-missing 1 \ 
--homozyg-window-threshold 0.05\ 
--out roh_relaxed 

#Result: 
#summary by group: 
# PRE (mfv*) 
awk '$2 ~ /^mfv/ && NR>1 {n++; nroh+=$4; sroh+=$5} END{printf "PRE n=%d mean_nROH=%.2f mean_SROH(Mb)=%.3f\n", n, nroh/n, (sroh/1000)/n}' roh_relaxed.hom.indiv 

# POS (mar*) 
awk '$2 ~ /^mar/ && NR>1 {n++; nroh+=$4; sroh+=$5} END{printf "POS n=%d mean_nROH=%.2f mean_SROH(Mb)=%.3f\n", n, nroh/n, (sroh/1000)/n}' roh_relaxed.hom.indiv

# Separating ROH into pre and post
cd /home/cmarinho/Chiroxiphia_caudata_genome/hwe

# PRE
awk 'NR==FNR{pre[$1]=1; next} 
FNR==1{hdr=$0; print hdr > "roh_pre.hom"; next} 
($2 in pre){print >> "roh_pre.hom"}' pre.txt roh_relaxed.hom

# POS
awk 'NR==FNR{pos[$1]=1; next} 
FNR==1{hdr=$0; print hdr > "roh_pos.hom"; next}
($2 in pos){print >> "roh_pos.hom"}' pos.txt roh_relaxed.hom

# conference
echo "PRE rows: $(($(wc -l < ​​roh_pre.hom)-1)) | POS rows: $(($(wc -l < ​​roh_pos.hom)-1))"
echo "PRE unique IIDs:"; awk 'NR>1{print $2}' roh_pre.hom | sort -u | wc -l
echo "POS unique IIDs:"; awk 'NR>1{print $2}' roh_pos.hom | sort -u | wc -l

# Normalizing pre and pos lists
# Creates lists with only the LAST column (IID), without header and without empty rows
awk 'NR>1 && $NF!="IID"{print $NF}' pre.txt | sed '/^[[:space:]]*$/d' | sort -u > pre_iid.txt
awk 'NR>1 && $NF!="IID"{print $NF}' pos.txt | sed '/^[[:space:]]*$/d' | sort -u > pos_iid.txt

###Calculating FROH via console
setwd("/home/cmarinho/Chiroxiphia_caudata_genome/hwe")
library(dplyr); library(readr); library(stringr); library(tidyr)

genome_len_bp <- 1089631598
thresholds <- c(5e5, 1e6, 2e6, 5e6)

read_plink_hom <- function(path){ 
df <- read.table(path, header=TRUE, sep="", stringsAsFactors=FALSE, check.names=FALSE) 
nms <- names(df) 
col_iid <- if ("IID" %in% nms) "IID" else if ("Id" %in% nms) "Id" else nms[2] 
col_kb <- nms[tolower(nms)=="kb"][1] 
if (is.na(col_kb)) stop("KB column not found in: ", path) 
tibble(sample = as.character(df[[col_iid]]), 
length_bp = as.numeric(df[[col_kb]])*1000)
}

pre_ids <- read_lines("pre_iid.txt") |> trimws()
pos_ids <- read_lines("pos_iid.txt") |> trimws()

seg_pre <- if (file.exists("roh_pre.hom")) read_plink_hom("roh_pre.hom") else tibble(sample=character(), length_bp=numeric())
seg_pos <- if (file.exists("roh_pos.hom")) read_plink_hom("roh_pos.hom") else tibble(sample=character(), length_bp=numeric())

seg <- bind_rows(seg_pre %>% mutate(group="PRE"), 
seg_pos %>% mutate(group="POS"))

per_ind <- sec %>% 
group_by(sample, group) %>% 
summarise(roh_total_bp = sum(length_bp, na.rm=TRUE), .groups="drop")

per_ind <- bind_rows( 
per_ind, 
tibble(sample = setdiff(pre_ids, per_ind$sample), group="PRE", roh_total_bp=0), 
tibble(sample = setdiff(pos_ids, per_ind$sample), group="POS", roh_total_bp=0)
)

per_ind <- per_ind %>% 
mutate(genome_len_bp = genome_len_bp, 
FROH_total = roh_total_bp / genome_len_bp)

for (thr in thresholds) { 
tmp <- sec %>% filter(length_bp >= thr) %>% 
group_by(sample) %>% summarise(val = sum(length_bp, na.rm=TRUE), .groups="drop") 
nm_bp <- paste0("roh_ge_", format(thr, scientific=FALSE), "_bp") 
nm_fr <- paste0("FROH_ge_", format(thr, scientific=FALSE)) 
per_ind <- per_ind %>% 
left_join(tmp, by="sample") %>% 
mutate(val = replace_na(val, 0), 
"{nm_bp}" := val, 
"{nm_fr}" := val/genome_len_bp) %>% 
select(-val)
}

per_ind <- per_ind %>% arrange(group, sample)
write_csv(per_ind, "froh_per_individual.csv")

num_cols <- names(per_ind)[str_detect(names(per_ind), "^FROH|_bp$")]
group_summary <- per_ind %>% 
group_by(group) %>% 
summarize( 
N = dplyr::n(), 
across(all_of(num_cols), 
list(mean=~mean(.x, na.rm=TRUE), 
sd=~ifelse(N>1, sd(.x, na.rm=TRUE), NA_real_), 
median=~median(.x, na.rm=TRUE)), 
.names="{.col}_{.fn}"),
.groups="drop"

) %>% arrange(group)

write_csv(group_summary, "froh_group_summary.csv")
cat("OK: froh_per_individual.csv and froh_group_summary.csv\n")

# Comparison between pre and post + graph
library(ggplot2)

# Uses the file you already generated (or the new one)
df <- read_csv("froh_per_individual.csv", show_col_types = FALSE)

# Non-parametric test (attention: very small n → low power)
w <- wilcox.test(FROH_total ~ group, data=df, exact=FALSE)

w

# Effect (Cliff's delta) — simple estimate
cliffs_delta <- function(x, y){

# Returns delta in [-1,1] 
sum(outer(x, y, FUN = function(a,b) sign(a-b))) / (length(x)*length(y))
}
delta <- with(df, cliffs_delta(FROH_total[group=="POS"], FROH_total[group=="PRE"]))
delta

# boxplot + points
ggplot(df, aes(group, FROH_total)) + 
geom_boxplot(outlier.shape = NA) + 
geom_jitter(width = 0.1, height = 0, size=2) + 
labs(title="FROH_total per group", x="", y="FROH_total") + 
theme_minimal()

# Running for larger sizes (between 5 and 8 MB - to catch recent events)

plink --vcf chiroxiphia.snps.vcf.gz \
--allow-extra-chr --double-id \
--snps-only just-acgt --biallelic-only strict \
--homozyg \
--homozyg-snp 100 \
--homozyg-kb 8000 \
--homozyg-density 50 \
--homozyg-gap 1000 \
--homozyg-window-snp 100 \
--homozyg-window-het 1 \
--homozyg-window-missing 5 \
--homozyg-window-threshold 0.05 \
--out roh_pos_8mb

#Result: zero for all

# Running for size from 4mb
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pos.txt \ 
--homozyg\ 
--homozyg-kb 4000\ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 1 \ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pos_4mb 
#Result: zero for everyone

# Rotating to size a from 2mb
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pos.txt \ 
--homozyg\ 
--homozyg-kb 2000\ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pos_2mbv2 
#Zero

# Running with relaxadp filter, 2mb and het=2
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pos.txt \ 
--homozyg\ 
--homozyg-kb 2000\ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pos_2mb_het2_0810 
#Result: zero

# Running with Het=2
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pos.txt \ 
--homozyg\ 
--homozyg-kb 2000\ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pos_2mb_het2_08102 
#Result: zero

#The strategy now Use the Het=2 filter on all tracks and compare the FROH at 1 Mb and 500 kb between the "pre" and "post" groups.

plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pos.txt \ 
--homozyg\ 
--homozyg-kb 1000\ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pos_1mb_het2 
#Result: 1.56mb for Mar051

#Running for Pre
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pre.txt\ 
--homozyg\ 
--homozyg-kb 1000\ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pre_1mb_het2 
#Result: zero

#Running from 500kb for pre
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pre.txt\ 
--homozyg\ 
--homozyg-kb 500\ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pre_500kb_het2 
#Result: 2.9mb for one of the samples- ancestral inbreeding

#running from 1mb to pre
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pre.txt\ 
--homozyg\ 
--homozyg-kb 1000\ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pre_1mb_het2 
#Result: zero

# Running the total FROH for the POST group using the filter Het=2
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pos.txt \ 
--homozyg\ 
--homozyg-kb 300 \ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \ 
--homozyg-window-threshold 0.05 \ 
--out roh_pos_total_het2 
#Result: 1 ROH for Mar 051 (1,559mb) 
# Running the total FROH for the pre group using the Het=2 filter 
plink --vcf chiroxiphia.snps.vcf.gz \ 
--allow-extra-chr --double-id \ 
--snps-only just-acgt --biallelic-only strict\ 
--keep pre.txt\ 
--homozyg\ 
--homozyg-kb 300 \ 
--homozyg-snp 100 \ 
--homozyg-density 50 \ 
--homozyg-gap 1000\ 
--homozyg-window-snp 100 \ 
--homozyg-window-het 2\ 
--homozyg-window-missing 5 \
--homozyg-window-threshold 0.05 \

--out roh_pre_total_het2

#Result: ALWAYS IN HOM.INDIV (2 ROH of 851.3kb for mfv4126, 1 of 376kb for mfv4127 and 11 of 5206kb for mfv4129.
------------------------------------------------------------------------------------------------
# Generating SFS with Hardfilter via Angsd

# Note: angsd is located in /home/cmarinho/angsd/angsd

# Note: realSFS is located in angsd/misc/

# Input files: bamlist_pre.txt and bamlist_pos.txt (located in /home/cmarinho/Chiroxiphia_caudata_genome/hwe/angsd_results)

# 1) Generating SAF by group

# 1.1. from the pre-accident group

/home/cmarinho/angsd/angsd/angsd -b /home/cmarinho/Chiroxiphia_caudata_genome/hwe/angsd_results/bamlist_pre.txt \ 
-ref "/home/cmarinho/chiroxiphia_refseq/GCF_009829145.1_bChiLan1.pri_genomic.fna" \ 
- anc "/home/cmarinho/chiroxiphia_refseq/GCF_009829145.1_bChiLan1.pri_genomic.fna" 
-out pre\ 
-uniqueOnly 1 \ 
-remove_bads 1 \ 
-only_proper_pairs 1 \ 
-baq 1\ 
-C 50 \ 
-minMapQ 20 \ 
-minQ 20 \ 
-GL 1 \ 
-doCounts 1\ 
-doMajorMinor 1 \ 
-doMaf 1 \ 
-minMaf 0.05 \ 
-SNP_pval 1e-6 \ 
-skipTriallelic 1\ 
-minInd 1 \ 
-doSaf 1 

# NOTE: 
# -minMapQ 20 \ # ~ MQ<20 
# -minQ 20\# base quality 
# -GL 1 \ # likelihoods SAMtools; You can use -GL 2 (GATK) if you prefer
# -minMaf 0.05 \ # MAF<0.05
# -SNP_pval 1e-6 \ # restricts to SNPs (approx. of your -v snps)
# -skipTriallelic 1 \ # biallelics (approx. -m2 -M2)
# -minInd 1 # F_MISSING>0.5 ⇒ requires ≥50% present

# 1.2. from the post-accident group
/home/cmarinho/angsd/angsd/angsd -b /home/cmarinho/Chiroxiphia_caudata_genome/hwe/angsd_results/bamlist_pos.txt \
-ref "/home/cmarinho/chiroxiphia_refseq/GCF_009829145.1_bChiLan1.pri_genomic.fna" \
-anc "/home/cmarinho/chiroxiphia_refseq/GCF_009829145.1_bChiLan1.pri_genomic.fna" \
-out pos \
-uniqueOnly 1 \
-remove_bads 1 \
-only_proper_pairs 1 \
-baq 1 \
-C 50 \
-minMapQ 20 \
-minQ 20 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-minMaf 0.05 \
-SNP_pval 1e-6 \
-skipTriallelic 1 \
-minInd 1 \
-doSaf 1

# ps.: It took around 5-6 hours per group

# 2) Generating SFS per group
/home/cmarinho/angsd/angsd/misc/realSFS pre.saf.idx -P 8 > pre.sfs
/home/cmarinho/angsd/angsd/misc/realSFS pos.saf.idx -P 8 > pos.sfs

# 3) Generating graph and statistical test

# realSFS files (1D)
# reads the SFS distributions
pre <- scan("/home/cmarinho/Chiroxiphia_caudata_genome/hwe/angsd_results/pre.sfs")
pos <- scan("/home/cmarinho/Chiroxiphia_caudata_genome/hwe/angsd_results/pos.sfs")

# normalizes to proportion
pre <- pre / sum(pre)
pos <- pos / sum(pos)

# defines the k-axis for each one
k_pre <- 0:(length(pre) - 1)
k_pos <- 0:(length(pos) - 1)

# finds the largest k
k_max <- max(k_pre, k_pos)

# creates a vector from 0 to the maximum
plot(0:k_max, c(pre, rep(NA, k_max + 1 - length(pre))),

type = "l", lwd = 2, col = "blue",
xlab = "Derived allele count (k)", ylab = "Proportion (SFS)",
main = "SFS per population", ylim = c(0, max(pre, pos) * 1.1))

# Adds the post-accident line
lines(k_pos, pos, lwd = 2, col = "red")
legend("topright", c("Pre-accident", "Post-accident"),
lwd = 2, col = c("blue", "red"))

# Difference test (K-S) on normalized vectors
ks <- ks.test(pre, pos)
print(ks)

-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Estimating demographic history modeling with Dadi
<IN PREP>
