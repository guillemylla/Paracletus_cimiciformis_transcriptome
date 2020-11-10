
# *Paracletus cimiciformis* transcriptome

- Authors: Miquel Barberà & Guillem Ylla, 2020
- Citation ""


## Read filtering

- Before assembling the reads, we filter out those reads that belong to other organisms such as bacteria, human and virus. The rest is supposed be mainly aphid reads.
 To do so, we will remove reads that map against reference genomes from organisms we are not interested.


### Install software

- Trinity
- Salmon
- Samtools
- BUSCO
- SortMeRNA

```

# Install latest(ish) versions of Trinity
# it will require bowtie2, jellyfish and salmon
cd ~/projects/paracletus/user_installed_software/
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.1/bowtie2-2.4.1-source.zip -O bowtie -O bowtie2-2.4.1-source.zip
unzip bowtie2-2.4.1-source.zip
cd bowtie2-2.4.1/
make
echo export PATH=\$PATH:`pwd`\ >> ~/.bashrc && source ~/.bashrc

# salmon
cd ~/projects/paracletus/user_installed_software/
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz
tar zxvf salmon-1.1.0_linux_x86_64.tar.gz
cd salmon-latest_linux_x86_64_1.1.0
echo export PATH=\$PATH:`pwd`\ >> ~/.bashrc && source ~/.bashrc

# install Samtools
cd ~/projects/paracletus/user_installed_software/
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar -vxjf samtools-1.10.tar.bz2
cd samtools-1.10/
make
echo export PATH=\$PATH:`pwd`\ >> ~/.bashrc && source ~/.bashrc
PATH=$PATH:~/augustus/bin:~/augustus/scripts

# Trinity
cd ~/projects/paracletus/user_installed_software/
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.9.1/trinityrnaseq-v2.9.1.FULL_with_extendedTestData.tar.gz -O Trinity-v2.9.1.tar.gz
tar zxvf Trinity-v2.9.1.tar.gz; rm Trinity-v2.9.1.tar.gz
cd trinityrnaseq-v2.9.1
# if required cmake → sudo apt install cmake
make
make plugins
export TRINITY_HOME=/home/quelo/projects/paracletus/user_installed_software/trinityrnaseq-v2.9.1/
echo export PATH=\$PATH:`pwd`\ >> ~/.bashrc && source ~/.bashrc

# install BUSCO
cd ~/projects/paracletus/user_installed_software/
git clone https://gitlab.com/ezlab/busco.git
cd busco/
# busco code and config file
python3 setup.py install --user
# configure config.ini file for BUSCO
python3 /home/quelo/projects/paracletus/user_installed_software/busco/scripts/busco_configurator.py config/config.ini config/myconfig.ini
export BUSCO_CONFIG_FILE="/home/quelo/projects/paracletus/user_installed_software/busco/config/myconfig.ini"


# SortMeRNA
cd ~/projects/paracletus/user_installed_software/
wget https://github.com/biocore/sortmerna/releases/download/v4.2.0/sortmerna-4.2.0-Linux.sh
mkdir sortmerna
bash sortmerna-4.2.0-Linux.sh --skip-license --prefix=sortmerna
export PATH=$PATH:/home/quelo/projects/paracletus/user_installed_software/sortmerna/bin:$PATH


wget http://bioinfo.lifl.fr/RNA/sortmerna/code/sortmerna-2.1-linux-64-multithread.tar.gz
tar -xvf sortmerna-2.1-linux-64-multithread.tar.gz sortmerna-2.1b/
cd sortmerna-2.1b/
echo export PATH=\$PATH:`pwd`\ >> ~/.bashrc && source ~/.bashrc
cd scripts
echo export PATH=\$PATH:`pwd`\ >> ~/.bashrc && source ~/.bashrc
# index all included databases
indexdb_rna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:\./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:\./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:\./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:\./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:\./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:\./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db


```

### Download genomes of putative contaminants

```

# download tables of refseq databases from "https://www.ncbi.nlm.nih.gov/genome" and select/filter/search groups of genomes of interest
cd databases
# most common endosymbionts in aphids are the obligate Buchnera aphidicola, and the facultative Serratia symbiotica, Regiella insecticola, Hamiltonella defensa, Rickettsiella viridis, Candidatus Sodalis sp. SoCistrobi and Candidatus Fukatsuia symbiotica. Select those that are hosted by aphids.

# Buchnera
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/237/295/GCF_005237295.1_ASM523729v1/GCF_005237295.1_ASM523729v1_genomic.fna.gz -P endosymbionts/buchnera
# Buchnera baizongia (NC_004545.1, NC_004555.1)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/725/GCF_000007725.1_ASM772v1/GCF_000007725.1_ASM772v1_genomic.fna.gz -P endosymbionts/buchnera
# Hamiltonella defensa (from A. pisum, it is also hosted by other insects)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/021/705/GCF_000021705.1_ASM2170v1/GCF_000021705.1_ASM2170v1_genomic.fna.gz -P endosymbionts/hamiltonella
# Serratia
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/370/165/GCF_008370165.1_ASM837016v1/GCF_008370165.1_ASM837016v1_genomic.fna.gz -P endosymbionts/serratia
# Regiella
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/143/625/GCF_000143625.1_ASM14362v1/GCF_000143625.1_ASM14362v1_genomic.fna.gz -P endosymbionts/regiella
# Rickettsia


wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/966/755/GCF_003966755.1_ASM396675v1/GCF_003966755.1_ASM396675v1_genomic.fna.gz -P endosymbionts/rickettsia
# Sodalis
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/143/145/GCF_900143145.1_SoCistrobi_v1.0/GCF_900143145.1_SoCistrobi_v1.0_genomic.fna.gz -P endosymbionts/sodalis
# Fukatsuia
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/122/425/GCF_003122425.1_ASM312242v1/GCF_003122425.1_ASM312242v1_genomic.fna.gz -P endosymbionts/fukatsuia

# extract all databases and combine into one endosymbiont database
gunzip -r endosymbionts/
# get number of contigs in each file recursively → 582 in this case
#find endosymbionts/ -name '*.fna' -exec grep -c '>' {} \;
# merge all endosymbionts contigs in a single database
find endosymbionts/ -name '*.fna' -exec cat {} \; > endosymbionts.fna


# mitochondrial aphid genomes
mkdir mitochondrion
# Acyrthosiphon pisum
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/508/785/GCF_005508785.1_pea_aphid_22Mar2018_4r6ur/GCF_005508785.1_pea_aphid_22Mar2018_4r6ur_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz -P mitochondrion; mv mitochondrion/chrMT.fna.gz mitochondrion/ApichrMT.fna.gz
# Aphis gossypii
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/010/815/GCF_004010815.1_ASM401081v1/GCF_004010815.1_ASM401081v1_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz -P mitochondrion; mv mitochondrion/chrMT.fna.gz mitochondrion/AgochrMT.fna.gz
# Diuraphis noxia
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/186/385/GCF_001186385.1_Dnoxia_1.0/GCF_001186385.1_Dnoxia_1.0_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz -P mitochondrion; mv mitochondrion/chrMT.fna.gz mitochondrion/DnochrMT.fna.gz
# Myzus persicae
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/856/785/GCF_001856785.1_MPER_G0061.0/GCF_001856785.1_MPER_G0061.0_assembly_structure/non-nuclear/assembled_chromosomes/FASTA/chrMT.fna.gz -P mitochondrion; mv mitochondrion/chrMT.fna.gz mitochondrion/MpechrMT.fna.gz
# extract all databases and combine into one endosymbiont database
gunzip -r mitochondrion/
find mitochondrion/ -name '*.fna' -exec cat {} \; > aphidsMT.fna


# download human database
wget https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/reference/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -P human/

# download C. elegans database
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz -P celegans/

# download A. thaliana database (best annotated plant) and Triticum aestivum (closest to Paracletus poacea feeding)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz -P plant/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/220/415/GCA_002220415.3_Triticum_4.0/GCA_002220415.3_Triticum_4.0_genomic.fna.gz -P plant/

# download S. cerevisiae database
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz -P yeast/

# download E. coli databases (2 reference genomes)
wget https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz -P ecoli/
wget https://ftp.ncbi.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz -P ecoli

# download virus (refseq AND invertebrate host) databases. 1905 genomes in total.
inv_virus_list=`cat refseq_reference_tables/inv_virus.list`
for f in $inv_virus_list;
do wget $f/*_genomic.fna.gz -P virus/;
do rm !(*_DNA.fasta.gz)
done

# extract and combine all virus databases in one
gunzip -r virus/
cd virus
cat *.fna > ../virus.fna
cd ..

# combine all "contamination" databases in one, including human, C. elegans, plants, S. cerevisiae, E.coli and viruses.
# moved all databases to /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/databases/contamination and concatenate them into contamination_all.fna


#---------------------------------------------------------------------------------------------------------------------------#
# MAPPING CONTAMINATION
#---------------------------------------------------------------------------------------------------------------------------#
# map reads to each of the databases with bowtie2
# combined databases "contamination.fna", "aphidsMT.fna" and "endosymbionts.fna" are all in databases folder
# for now, only use contamination.fna
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads
mkdir mapping
mkdir mapping/contamination
# move databases to correspondig mapping folder
#mv databases/contamination.fna mapping/contamination/
# index each database
#cd mapping/
#bowtie2-build contamination/contamination.fna contamination/contamination
# indexing of contamination.fna took >36 h. Index each contamination separated.
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/databases/contamination
cat GCF_000005845.2_ASM584v2_genomic.fna GCF_000008865.2_ASM886v2_genomic.fna  > E.coli_genomic.fna
# index with bowtie2-build
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/contamination
mkdir triticum arabidopsis ecoli human worm virus yeast
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/
bowtie2-build databases/contamination/GCA_002220415.3_Triticum_4.0_genomic.fna mapping/contamination/triticum/Triticum
bowtie2-build databases/contamination/GCF_000001735.4_TAIR10.1_genomic.fna mapping/contamination/arabidopsis/arabidopsis
bowtie2-build databases/contamination/E.coli_genomic.fna mapping/contamination/ecoli/ecoli
bowtie2-build databases/contamination/GCF_000001405.39_GRCh38.p13_genomic.fna mapping/contamination/human/human
bowtie2-build databases/contamination/GCF_000002985.6_WBcel235_genomic.fna mapping/contamination/worm/worm
bowtie2-build databases/contamination/virus.fna mapping/contamination/virus/virus
bowtie2-build databases/contamination/GCF_000146045.2_R64_genomic.fna mapping/contamination/yeast/yeast

# map with bowtie2 (adapted from http://www.metagenomics.wiki/tools/short-read/remove-host-sequences)

# bowtie2 mapping against HUMAN sequences database, keep both mapped and unmapped reads (paired-end reads)
cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/human/human -p 14 -q --rf -1 ../../R1_pairedout.fastq -2 ../../R2_pairedout.fastq 2>contamination/human/align_stats.txt | samtools view -@14 -bS -> contamination/human/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 14 -b -f 12 -F 256 contamination/human/contamination_mapped_and_unmapped.bam > contamination/human/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 14 -n contamination/human/SAMPLE_bothEndsUnmapped.bam -o contamination/human/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i contamination/human/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont1_R1.fastq -fq2 ../reads/decont1_R2.fastq

# bowtie2 mapping against E. coli sequences databases, keep both mapped and unmapped reads (paired-end reads)
#cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/ecoli/ecoli -p 14 -q --rf -1 ../reads/decont1_R1.fastq -2 ../reads/decont1_R2.fastq 2>contamination/ecoli/align_stats.txt | samtools view -@14 -bS -> contamination/ecoli/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 14 -b -f 12 -F 256 contamination/ecoli/contamination_mapped_and_unmapped.bam > contamination/ecoli/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 14 -n contamination/ecoli/SAMPLE_bothEndsUnmapped.bam -o contamination/ecoli/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i contamination/ecoli/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont2_R1.fastq -fq2 ../reads/decont2_R2.fastq

# bowtie2 mapping against YEAST sequences databases, keep both mapped and unmapped reads (paired-end reads)
#cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/yeast/yeast -p 14 -q --rf -1 ../reads/decont2_R1.fastq -2 ../reads/decont2_R2.fastq 2>contamination/yeast/align_stats.txt | samtools view -@14 -bS -> contamination/yeast/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 14 -b -f 12 -F 256 contamination/yeast/contamination_mapped_and_unmapped.bam > contamination/yeast/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 14 -n contamination/yeast/SAMPLE_bothEndsUnmapped.bam -o contamination/yeast/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i contamination/yeast/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont3_R1.fastq -fq2 ../reads/decont3_R2.fastq

# bowtie2 mapping against WORM sequences databases, keep both mapped and unmapped reads (paired-end reads)
#cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/worm/worm -p 14 -q --rf -1 ../reads/decont3_R1.fastq -2 ../reads/decont3_R2.fastq 2>contamination/worm/align_stats.txt | samtools view -@14 -bS -> contamination/worm/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 14 -b -f 12 -F 256 contamination/worm/contamination_mapped_and_unmapped.bam > contamination/worm/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 14 -n contamination/worm/SAMPLE_bothEndsUnmapped.bam -o contamination/worm/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i contamination/worm/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont4_R1.fastq -fq2 ../reads/decont4_R2.fastq

# bowtie2 mapping against ARABIDOPSIS sequences databases, keep both mapped and unmapped reads (paired-end reads)
#cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/arabidopsis/arabidopsis -p 14 -q --rf -1 ../reads/decont4_R1.fastq -2 ../reads/decont4_R2.fastq 2>contamination/arabidopsis/align_stats.txt | samtools view -@14 -bS -> contamination/arabidopsis/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 14 -b -f 12 -F 256 contamination/arabidopsis/contamination_mapped_and_unmapped.bam > contamination/arabidopsis/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 14 -n contamination/arabidopsis/SAMPLE_bothEndsUnmapped.bam -o contamination/arabidopsis/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i contamination/arabidopsis/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont5_R1.fastq -fq2 ../reads/decont5_R2.fastq

# bowtie2 mapping against TRITICUM sequences databases, keep both mapped and unmapped reads (paired-end reads)
#cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/triticum/triticum -p 14 -q --rf -1 ../reads/decont5_R1.fastq -2 ../reads/decont5_R2.fastq 2>contamination/triticum/align_stats.txt | samtools view -@14 -bS -> contamination/triticum/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 14 -b -f 12 -F 256 contamination/triticum/contamination_mapped_and_unmapped.bam > contamination/triticum/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 14 -n contamination/triticum/SAMPLE_bothEndsUnmapped.bam -o contamination/triticum/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i contamination/triticum/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont6_R1.fastq -fq2 ../reads/decont6_R2.fastq

# bowtie2 mapping against VIRUS sequences databases, keep both mapped and unmapped reads (paired-end reads)
cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/virus/virus -p 14 -q --rf -1 ../reads/decont6_R1.fastq -2 ../reads/decont6_R2.fastq 2>contamination/virus/align_stats.txt | samtools view -@14 -bS -> contamination/virus/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 14 -b -f 12 -F 256 contamination/virus/contamination_mapped_and_unmapped.bam > contamination/virus/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 14 -n contamination/virus/SAMPLE_bothEndsUnmapped.bam -o contamination/virus/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i contamination/virus/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont7_R1.fastq -fq2 ../reads/decont7_R2.fastq

# check number of reads removed on each step with echo $(cat decont1_R1.fastq|wc -l)/4|bc


### FURTHER REMOVE CONTAMINATION WITH RELAXED BOWTIE2 PARAMETERS ###

# concatenate all "contamination" databases into a single file (74040 fasta entries)
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/contamination

find -name '*.fna' -exec cat {} \; > contamination.fna
# index database
bowtie2-build --large-index --threads 15 contamination.fna contamination

# map with bowtie2 relaxed parameters
cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/contamination -p 14 --mp 4 --rdg 4,2 --rfg 4,2 -q --rf -1 ../reads/decont7_R1.fastq -2 ../reads/decont7_R2.fastq 2>contamination/align_stats.txt | samtools view -@15 -bS -> contamination/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 15 -b -f 12 -F 256 contamination/contamination_mapped_and_unmapped.bam > contamination/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 15 -n contamination/SAMPLE_bothEndsUnmapped.bam -o contamination/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i contamination/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont8_R1.fastq -fq2 ../reads/decont8_R2.fastq

# ENDOSYMBIONTS
# bowtie2 mapping against ENDOSYMBIONTS sequences databases, keep both mapped and unmapped reads (paired-end reads)
cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x endosymbionts/endosymbionts -p 15 -q --rf -1 ../reads/decont8_R1.fastq -2 ../reads/decont8_R2.fastq 2>endosymbionts/align_stats.txt | samtools view -@15 -bS -> endosymbionts/endosymbionts_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 15 -b -f 12 -F 256 endosymbionts/endosymbionts_mapped_and_unmapped.bam > endosymbionts/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -@ 15 -n endosymbionts/SAMPLE_bothEndsUnmapped.bam -o endosymbionts/SAMPLE_bothEndsUnmapped_sorted.bam
bedtools bamtofastq -i endosymbionts/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../reads/decont9_R1.fastq -fq2 ../reads/decont9_R2.fastq


# MITOCHONDRIAL reads were not removed





#---------------------------#
# SortMeRNA to remove rRNA
#---------------------------#
# remove rRNA from filtered contaminations reads. Use all databases
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/reads/sortmerna
# use script to merge R1 and R2 fastq files
merge-paired-reads.sh ../decont9_R1.fastq ../decont9_R2.fastq decont.merged.fastq
# add environmental variables as paths to database and indexes
export SORTMERNA_DB=$HOME/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases
export SORTMERNA_INDEX=$HOME/projects/paracletus/user_installed_software/sortmerna-2.1b/index
/home/quelo/projects/paracletus/cleaned_trimmed_reads/reads
# filter rRNA reads (with 20GB RAM). Only removes paired reads that have both aligned to rRNA
sortmerna --ref $SORTMERNA_DB/rfam-5s-database-id98.fasta,$SORTMERNA_INDEX/rfam-5s-db:\
$SORTMERNA_DB/rfam-5.8s-database-id98.fasta,$SORTMERNA_INDEX/rfam-5.8s-db:\
$SORTMERNA_DB/silva-arc-16s-id95.fasta,$SORTMERNA_INDEX/silva-arc-16s-db:\
$SORTMERNA_DB/silva-arc-23s-id98.fasta,$SORTMERNA_INDEX/silva-arc-23s-db:\
$SORTMERNA_DB/silva-bac-16s-id90.fasta,$SORTMERNA_INDEX/silva-bac-16s-db:\
$SORTMERNA_DB/silva-bac-23s-id98.fasta,$SORTMERNA_INDEX/silva-bac-23s-db:\
$SORTMERNA_DB/silva-euk-18s-id95.fasta,$SORTMERNA_INDEX/silva-euk-18s-db:\
$SORTMERNA_DB/silva-euk-28s-id98.fasta,$SORTMERNA_INDEX/silva-euk-28s \
--reads decont.merged.fastq -a 15 -m 20480 --paired_in --sam --num_alignments 1 --fastx --aligned decont.rRNA --other decont.norRNA --log -v

# split output fastq into R1 and R2
unmerge-paired-reads.sh decont.norRNA.fastq decont.norRNA.R1.fastq decont.norRNA.R2.fastq

echo $(cat decont.norRNA.R1.fastq|wc -l)/4|bc
echo $(cat decont.norRNA.R2.fastq|wc -l)/4|bc



#---------------------------------------------------------------------------------------------------------------------------#
# ASSEMBLY 3: removed contamination + relaxed decontamination + remove endosymbionts + removed rRNA (kept mirochondrial)
#---------------------------------------------------------------------------------------------------------------------------#

# running assembly
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads
mkdir assemblies/trinity3
Trinity --seqType fq --max_memory 60G --SS_lib_type RF --min_kmer_cov 2 --left reads/sortmerna/decont.norRNA.R1.fastq --right reads/sortmerna/decont.norRNA.R2.fastq --CPU 15 --output assemblies/trinity3/ 1>assemblies/trinity3/Trinity.log 2>assemblies/trinity3/Trinity.err

# get Trinity basic statistics
TrinityStats.pl assemblies/trinity3/Trinity.fasta > assemblies/trinity3/Stats.trinity3.txt

# use BUSCO to assess completeness of the assembly against arthropoda
#export BUSCO_CONFIG_FILE="/home/quelo/projects/paracletus/user_installed_software/busco/config/myconfig.ini"
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assemblies/trinity3/
mkdir BUSCO
cd BUSCO
busco -c 15 -m transcriptome -i ../Trinity.fasta -o arthr -l arthropoda_odb10 1>busco_arth.log 2>busco_arth.err
busco -c 15 -m transcriptome -i ../Trinity.fasta -o insect -l insecta_odb10 1>busco_insect.log 2>busco_insect.err
busco -c 15 -m transcriptome -i ../Trinity.fasta -o bact -l bacteria_odb10 1>busco_bact.log 2>busco_bact.err


# backmapping
mkdir backmap
cd backmap
bowtie2-build --threads 15 ../Trinity.fasta Trinity.fasta
bowtie2 --local --no-unal -x Trinity.fasta -p 15 \
      -q --rf -1 ../../../reads/sortmerna/decont.norRNA.R1.fastq -2 ../../../reads/sortmerna/decont.norRNA.R2.fastq 2>align_stats.txt\
      | samtools view -b -@ 15 | samtools sort -@ 15 -o bowtie2.bam







#######################################################
# Quantification of transcripts with bowtie2 and RSEM #
#######################################################

cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assemblies/trinity3/
mkdir quantification
cd quantification

# create a sample read file. NOTE! Directories must be in reference to actual working directory.
cat > reads.txt << "EOF"
R_morph	R_morph_rep1	../../../../reads/1V_S6_.R1_pairedout.fq	../../../../reads/1V_S6_.R2_pairedout.fq
R_morph	R_morph_rep2	../../../../reads/2V_S3_.R1_pairedout.fq	../../../../reads/2V_S3_.R2_pairedout.fq
R_morph	R_morph_rep3	../../../../reads/3V_S2_.R1_pairedout.fq	../../../../reads/3V_S2_.R2_pairedout.fq
F_morph	F_morph_rep1	../../../../reads/7B_S5_.R1_pairedout.fq	../../../../reads/7B_S5_.R2_pairedout.fq
F_morph	F_morph_rep2	../../../../reads/8B_S1_.R1_pairedout.fq	../../../../reads/8B_S1_.R2_pairedout.fq
F_morph	F_morph_rep3	../../../../reads/10B_S4_.R1_pairedout.fq	../../../../reads/10B_S4_.R2_pairedout.fq
EOF

#export READS=/home/quelo/projects/paracletus/cleaned_trimmed_reads/reads


# generate count matrices with RSEM. First align the original (the cleaned and filtered by quality) rna-seq reads back against the Trinity transcripts, then run RSEM to estimate the number of rna-seq fragments that map to each contig
# The Trinity script will run the Bowtie2 aligner to align reads to the Trinity transcripts, and RSEM will then evaluate those alignments to estimate expression values separately for each sample.
align_and_estimate_abundance.pl --seqType fq --samples_file reads.txt --transcripts ../paracletus_contigs.fa --est_method RSEM --aln_method bowtie2 --output_dir rsem_outdir --SS_lib_type RF --thread_count $NPROC --trinity_mode --prep_reference

# we have a table with counts and normalized expression for each sample
# now create two matrices: one with the estimated counts and another containing TMP normalized expression values
# create a list of the quant.sf files for isoforms
find R_morph_* F_morph_* -name "RSEM.isoforms.results" | tee quant_files.isoforms.list

# Generate Gene/Transcript relationships. For Trinity assemblies is done with the following transcript (newer versions produce this file automatically by --trinity_mode)
#get_Trinity_gene_to_trans_map.pl ../paracletus_contigs.fa > Trinity.fasta.gene_trans_map

# generate the count and expression matrices for both the transcript isoforms and sepearate files for 'gene's
abundance_estimates_to_matrix.pl --est_method RSEM --quant_files quant_files.isoforms.list --gene_trans_map ../paracletus_contigs.fa.gene_trans_map --out_prefix paracletus_contigs.rsem --name_sample_by_basedir

#among the results there is a counts matrix and a normalized expression matrix (TPM) for both isoforms and genes. The results have been also cross-sample normalized (TMM).


### Contig Ex90N50 Statistic and Ex90 Gene Count ###
#check quality of assembly (https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats#contig-ex90n50-statistic-and-ex90-gene-count)
export TRINITY_HOME=$HOME/projects/paracletus/user_installed_software/trinityrnaseq-v2.9.1
mkdir ExN50.stats
#$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl assembly3.rsem.isoform.TPM.not_cross_norm ../Trinity.fasta | tee ExN50.stats/ExN50.isoforms.stats
#${TRINITY_HOME}/util/misc/plot_ExN50_statistic.Rscript ExN50.stats/ExN50.isoforms.stats
#xpdf ExN50.stats/ExN50.stats.plot.pdf
#ExN50 referes to the N50 but for the top 100 more expressed genes. This is more menaningful in transcriptomes. Our assembly peaks around 85% (optimum is 90%).



### Counting Numbers of Expressed Transcripts or Genes ###
mkdir counting
$TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
	assembly3.rsem.isoform.TPM.not_cross_norm | tee counting/trans_matrix.TPM.not_cross_norm.counts_by_min_TPM

$TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
	assembly3.rsem.gene.TPM.not_cross_norm | tee counting/gene_matrix.TPM.not_cross_norm.counts_by_min_TPM

# plot the numer of genes with R
R
data = read.table("counting/gene_matrix.TPM.not_cross_norm.counts_by_min_TPM", header=T)
plot(data, xlim=c(-100,0), ylim=c(0,100000), t='b')
data2 = read.table("counting/trans_matrix.TPM.not_cross_norm.counts_by_min_TPM", header=T)
plot(data2, xlim=c(-100,0), ylim=c(0,100000), t='b')
q()


########################################################
### TRANSCRIPTS FILTERING BASED ON EXPRESSION VALUES ###
########################################################

/home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assemblies/trinity3/quantification
# attempt to reduce redundancy in the transcriptome by retaining the most expressed isoform only
mkdir ../reduced_asssemblies/gene_assembly
filter_low_expr_transcripts.pl \
	--matrix assembly3.rsem.isoform.TPM.not_cross_norm \
	--transcripts ../Trinity.fasta \
	--highest_iso_only \
	--trinity_mode \
	> ../reduced_asssemblies/gene_assembly/a3_gene.fa

# Retained 319789 / 472478 = 67.68% of total transcripts. NOTE! BUSCO_insecta down to 69.1%
#---------------------------------------------------------------------------------------------#

# attempt to reduce redundancy in the transcriptome by retaining transcripts with TPM >0.5
mkdir ../reduced_asssemblies ../reduced_asssemblies/exp0.5
filter_low_expr_transcripts.pl \
	--matrix assembly3.rsem.isoform.TPM.not_cross_norm \
	--transcripts ../Trinity.fasta \
	--min_expr_any 0.5 \
	--trinity_mode \
	> ../reduced_asssemblies/exp0.5/a3_exp0.5.fa

# Retained 361021 / 472478 = 76.41% of total transcripts. NOTE! BUSCOinsecta down to XXX%
#---------------------------------------------------------------------------------------------#

# attempt to reduce redundancy in the transcriptome by retaining transcripts with TPM >1
mkdir ../reduced_asssemblies/exp1
filter_low_expr_transcripts.pl \
	--matrix assembly3.rsem.isoform.TPM.not_cross_norm \
	--transcripts ../Trinity.fasta \
	--min_expr_any 1 \
	--trinity_mode \
	> ../reduced_asssemblies/exp1/a3_exp1.fa

# Retained 201333 / 472478 = 42.61% of total transcripts. NOTE! BUSCO_insecta down to 93%
#---------------------------------------------------------------------------------------------#

# attempt to reduce redundancy in the transcriptome by retaining transcripts with TPM >1.5
mkdir ../reduced_asssemblies/exp1.5
filter_low_expr_transcripts.pl \
	--matrix assembly3.rsem.isoform.TPM.not_cross_norm \
	--transcripts ../Trinity.fasta \
	--min_expr_any 1.5 \
	--trinity_mode \
	> ../reduced_asssemblies/exp1.5/a3_exp1.5.fa

# Retained 119723 / 472478 = 25.34% of total transcripts. NOTE! BUSCO_insecta down to 92.1%
#---------------------------------------------------------------------------------------------#

# attempt to reduce redundancy in the transcriptome by retaining transcripts with TPM >1.7
mkdir ../reduced_asssemblies/exp1.7
filter_low_expr_transcripts.pl \
	--matrix assembly3.rsem.isoform.TPM.not_cross_norm \
	--transcripts ../Trinity.fasta \
	--min_expr_any 1.7 \
	--trinity_mode \
	> ../reduced_asssemblies/exp1.7/a3_exp1.7.fa

# Retained 100685 / 472478 = 21.31% of total transcripts. NOTE! BUSCOinsecta down to 91.6%
#---------------------------------------------------------------------------------------------#

# attempt to reduce redundancy in the transcriptome by retaining transcripts with TPM >2
mkdir ../reduced_asssemblies/exp2
filter_low_expr_transcripts.pl \
	--matrix assembly3.rsem.isoform.TPM.not_cross_norm \
	--transcripts ../Trinity.fasta \
	--min_expr_any 2 \
	--trinity_mode \
	> ../reduced_asssemblies/exp2/a3_exp2.fa

# Retained 79025 / 472478 = 16.73% of total transcripts. NOTE! BUSCOinsecta down to 90.4%
#---------------------------------------------------------------------------------------------#












































###############################################
# DE with Trinity scripts using DESeq2 method #
###############################################


#DEseq2 is used to identify DE transcripts
run_DE_analysis.pl \
      --matrix Trinity5.rsem.isoform.counts.matrix \
      --samples_file ../reads.txt \
      --method DESeq2 \
      --output DESeq2_trans


#now the analysis at the gene level
run_DE_analysis.pl \
      --matrix Trinity5.rsem.gene.counts.matrix \
      --samples_file ../reads.txt \
      --method DESeq2 \
      --output DESeq2_gene

--------------------------------------------------------------------------
#Extracting differentially expressed transcripts and generating heatmaps #
--------------------------------------------------------------------------

cd DESeq2_trans
#Extract those differentially expressed (DE) transcripts that are at least 4-fold differentially expressed at a significance of <= 0.001 in any of the pairwise sample comparisons
analyze_diff_expr.pl \
      --matrix ../Trinity5.rsem.isoform.TMM.EXPR.matrix \
      --samples ../../reads.txt \
      -P 1e-3 -C 2 

#repeat at the gene level

cd DESeq2_gene
#Extract those differentially expressed (DE) transcripts that are at least 4-fold differentially expressed at a significance of <= 0.001 in any of the pairwise sample comparisons
analyze_diff_expr.pl \
      --matrix ../Trinity5.rsem.gene.TMM.EXPR.matrix \
      --samples ../../reads.txt \
      -P 1e-3 -C 2 

------------------------------------------------------------------------------
#Extract transcript clusters by expression profile by cutting the dendrogram #
------------------------------------------------------------------------------

#Extract clusters of transcripts with similar expression profiles by cutting the transcript cluster dendrogram at a given percent of its height (ex. 60%)
define_clusters_by_cutting_tree.pl \
       --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData

#for transcripts
cd ../DEseq2_trans

define_clusters_by_cutting_tree.pl \
       --Ptree 60 -R diffExpr.P1e-3_C2.matrix.RData

































































































































#---------------------------------------------------------------------------------------------------------------------------#
# ASSEMBLY 2: removed contamination and rRNA
#---------------------------------------------------------------------------------------------------------------------------#

#---------------------------#
# SortMeRNA to remove rRNA
#---------------------------#
# remove rRNA from filtered contaminations reads. Use all databases (?)
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/reads
mkdir sortmerna
cd sortmerna
# use script to merge R1 and R2 fastq files
merge-paired-reads.sh ../decont7_R1.fastq ../decont7_R2.fastq decont.merged.fastq
# add environmental variables as paths to database and indexes
export SORTMERNA_DB=$HOME/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases
export SORTMERNA_INDEX=$HOME/projects/paracletus/user_installed_software/sortmerna-2.1b/index

# filter rRNA reads (with 20GB RAM). Only removes paired reads that have both aligned to rRNA
sortmerna --ref $SORTMERNA_DB/rfam-5s-database-id98.fasta,$SORTMERNA_INDEX/rfam-5s-db:\
$SORTMERNA_DB/rfam-5.8s-database-id98.fasta,$SORTMERNA_INDEX/rfam-5.8s-db:\
$SORTMERNA_DB/silva-arc-16s-id95.fasta,$SORTMERNA_INDEX/silva-arc-16s-db:\
$SORTMERNA_DB/silva-arc-23s-id98.fasta,$SORTMERNA_INDEX/silva-arc-23s-db:\
$SORTMERNA_DB/silva-bac-16s-id90.fasta,$SORTMERNA_INDEX/silva-bac-16s-db:\
$SORTMERNA_DB/silva-bac-23s-id98.fasta,$SORTMERNA_INDEX/silva-bac-23s-db:\
$SORTMERNA_DB/silva-euk-18s-id95.fasta,$SORTMERNA_INDEX/silva-euk-18s-db:\
$SORTMERNA_DB/silva-euk-28s-id98.fasta,$SORTMERNA_INDEX/silva-euk-28s \
--reads decont.merged.fastq -a 14 -m 20480 --paired_out --sam --num_alignments 1 --fastx --aligned decont.rRNA --other decont.norRNA --log -v

# split output fastq into R1 and R2
unmerge-paired-reads.sh decont.norRNA.fastq decont.norRNA.R1.fastq decont.norRNA.R2.fastq

echo $(cat decont.norRNA.R1.fastq|wc -l)/4|bc
echo $(cat decont.norRNA.R2.fastq|wc -l)/4|bc


# running assembly with decontaminated and rRNA removed
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads
mkdir assemblies/trinity_decont_norRNA
Trinity --seqType fq --max_memory 60G --SS_lib_type RF --min_kmer_cov 2 --left reads/sortmerna/decont.norRNA.R1.fastq --right reads/sortmerna/decont.norRNA.R2.fastq --CPU 14 --output assemblies/trinity_decont_norRNA/ 1>assemblies/trinity_decont_norRNA/Trinity.log 2>assemblies/trinity_decont_norRNA/Trinity.err

# get Trinity basic statistics
TrinityStats.pl Trinity.fasta > Stats.decont.norRNA.txt

# use BUSCO to assess completeness of the assembly against arthropoda
#export BUSCO_CONFIG_FILE="/home/quelo/projects/paracletus/user_installed_software/busco/config/myconfig.ini"
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assemblies/trinity_decont_norRNA
mkdir BUSCO
cd BUSCO
busco -c 14 -m transcriptome -i ../Trinity.fasta -o arthr -l arthropoda_odb10 1>busco_arth.log 2>busco_arth.err
busco -c 14 -m transcriptome -i ../Trinity.fasta -o insect -l insecta_odb10 1>busco_insect.log 2>busco_insect.err
busco -c 14 -m transcriptome -i ../Trinity.fasta -o bact -l bacteria_odb10 1>busco_bact.log 2>busco_bact.err

# cd-hit remove redundancy
#cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assembly/trinity_decont_norRNA
#mkdir cdhit cdhit/0.95
#cd cdhit/0.95
#cd-hit-est -i ../../Trinity.fasta -c 0.95 -o cdhit0.95.fasta -c 0.95 -n 8 -p 1 -g 1 -M 81920 -T 14 -d 40 1>cd-hit-est_0.95.log 2>cd-hit-est_0.95.err
# reduced the assembly to 411k contigs. Still too much!


















































#find . -maxdepth 1 -type d \( -path ./endosymbionts -o -path ./virus -o -path ./refseq_reference_tables -o -path ./mitochondrion \) -prune -o -print | find . -name '*.fna' -exec cat {} \; > contamination2.fna
#cat contamination2.fna virus.fna > contamination.fna



#---------------------------------------------------------------------------------------------------------------------------#
# MAPPING
# map reads to each of the databases with bowtie2
# combined databases "contamination.fna", "aphidsMT.fna" and "endosymbionts.fna" are all in databases folder
# first, the mitochondrial reads will be set apart. Then, the endosymbionts will also be set apart. Finally, the reads mapping to contamination will be discarded.
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads
mkdir mapping
mkdir mapping/contamination mapping/endosymbionts mapping/mitochondrion
# move databases to correspondig mapping folder
mv databases/endosymbionts.fna mapping/endosymbionts/
mv databases/contamination.fna mapping/contamination/
mv databases/aphidsMT.fna mapping/mitochondrion/
#index each database
cd mapping/
bowtie2-build contamination/contamination.fna contamination/contamination
bowtie2-build endosymbionts/endosymbionts.fna endosymbionts/endosymbionts
bowtie2-build mitochondrion/aphidsMT.fna mitochondrion/aphidsMT

# map with bowtie2 (adapted from http://www.metagenomics.wiki/tools/short-read/remove-host-sequences)

# bowtie2 mapping against mitochondrial sequences database, keep both mapped and unmapped reads (paired-end reads)
bowtie2 -x mitochondrion/aphidsMT -p 14 -q --rf -1 ../../R1_pairedout.fastq -2 ../../R2_pairedout.fastq 2>mitochondrion/align_stats.txt | samtools view -@14 -bS -> mitochondrion/MT_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -b -f 12 -F 256 MT_mapped_and_unmapped.bam > removed_MT_bothEndsUnmapped.bam
# SAMtools SAM-flag filter: get mapped pairs (both ends mapped)
samtools view -b -f 2 -F 256 MT_mapped_and_unmapped.bam > MT_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -n removed_MT_bothEndsUnmapped.bam removed_MT_bothEndsUnmapped_sorted
bedtools bamtofastq -i removed_MT_bothEndsUnmapped_sorted.bam -fq removed_MT_R1.fastq -fq2 removed_MT_R2.fastq


# bowtie2 mapping against endosymbionts sequences database, keep both mapped and unmapped reads (paired-end reads)
bowtie2 -x endosymbionts/endosymbionts -p 14 -q --rf -1 ../../R1_pairedout.fastq -2 ../../R2_pairedout.fastq 2>endosymbionts/align_stats.txt | samtools view -@14 -bS -> endosymbionts/endosym_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -b -f 12 -F 256 virus_mapped_and_unmapped.bam > notvirus_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -n notvirus_bothEndsUnmapped.bam notvirus_bothEndsUnmapped_sorted
bedtools bamtofastq -i notvirus_bothEndsUnmapped_sorted.bam -fq virus_removed_r1.fastq -fq2 virus_removed_r2.fastq

# map with bowtie2 (adapted from http://www.metagenomics.wiki/tools/short-read/remove-host-sequences)
# bowtie2 mapping against contamination sequences database, keep both mapped and unmapped reads (paired-end reads)
cd /projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping/
bowtie2 -x contamination/contamination -p 14 -q --rf -1 ../../R1_pairedout.fastq -2 ../../R2_pairedout.fastq 2>contamination/align_stats.txt | samtools view -@14 -bS -> contamination/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -b -f 12 -F 256 contamination/contamination_mapped_and_unmapped.bam > contamination/SAMPLE_bothEndsUnmapped.bam
# split paired-end reads into separated fastq files .._r1 .._r2
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -n contamination/SAMPLE_bothEndsUnmapped.bam contamination/SAMPLE_bothEndsUnmapped_sorted
bedtools bamtofastq -i contamination/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../filtered_reads/cont_removed_R1.fastq -fq2 ../filtered_reads/cont_removed_R2.fastq

# running assembly on reads with filtered out contamination (not endosymbionts, nor mitochondrial)
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads
mkdir assemblies assemblies/Trinity_removed_cont/
Trinity --seqType fq --max_memory 60G --SS_lib_type RF --min_kmer_cov 2 --left reads/cont_removed_R1.fastq --right reads/cont_removed_R2.fastq --CPU 14 --output assemblies/Trinity_removed_cont/ 1>assemblies/Trinity_removed_cont/Trinity.log 2>assemblies/Trinity_removed_cont/Trinity.err

# get Trinity basic statistics
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assembly/Trinity_removed_cont
mkdir TrinityStats
cd TrinityStats
TrinityStats.pl ../Trinity.fasta > Stats.txt


# use BUSCO to assess completeness of the assembly against arthropoda
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assemblies/Trinity_removed_cont
mkdir BUSCO
cd BUSCO
busco -c 14 -m transcriptome -i ../Trinity.fasta -o arthr -l arthropoda_odb10 1>busco_arth.log 2>busco_arth.err
busco -c 14 -m transcriptome -i ../Trinity.fasta -o bact -l bacteria_odb10 1>busco_bact.log 2>busco_bact.err
busco -c 14 -m transcriptome -i ../Trinity.fasta -o insect -l insecta_odb10 1>busco_insect.log 2>busco_insect.err





# remove rRNA from filtered contaminations reads. Use all databases (?)
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/reads
mkdir sortmerna
cd sortmerna
# use script to merge R1 and R2 fastq files
merge-paired-reads.sh ../cont_removed_R1.fastq ../cont_removed_R2.fastq cont_removed_merged_reads.fastq
# add environmental variables as paths to database and indexes
export SORTMERNA_DB=$HOME/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases
export SORTMERNA_INDEX=$HOME/projects/paracletus/user_installed_software/sortmerna-2.1b/index

# filter rRNA reads (with 20GB RAM)
#sortmerna --ref /home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-bac-16s-db:\/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-bac-23s-db:\/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-arc-16s-db:\/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-arc-23s-db:\/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-euk-18s-db:\/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-euk-28s:\/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/rfam-5s-db:\/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/rfam-5.8s-db\ --reads cont_removed_merged_reads.fastq -m 20480 --paired_out --sam --num_alignments 1 --fastx --aligned cont_removed_rRNA --other cont_removed_non_rRNA --log -v







sortmerna --ref /home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-euk-18s-db --reads cont_removed_merged_reads.fastq -a 14 -m 20480 --paired_out --sam --num_alignments 1 --fastx --aligned cont_removed_rRNA --other cont_removed_non_rRNA --log -v

sortmerna --ref /home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-euk-28s --reads cont_removed_non_rRNA.fastq -a 14 -m 20480 --paired_out --sam --num_alignments 1 --fastx --aligned cont_removed_rRNA2 --other cont_removed_non_rRNA2 --log -v

sortmerna --ref /home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-bac-16s-db:/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-bac-23s-db:/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-arc-23s-db:/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/silva-arc-16s-db:/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/rfam-5s-db:/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,/home/quelo/projects/paracletus/user_installed_software/sortmerna-2.1b/index/rfam-5.8s-db --reads cont_removed_non_rRNA2.fastq -a 14 -m 20480 --paired_out --sam --num_alignments 1 --fastx --aligned cont_removed_rRNA3 --other cont_removed_non_rRNA3 --log -v

# count number of fastq reads non rRNA → 131262701, i.e. ~ 3.75M reads removed.
echo $(cat cont_removed_rRNA3.fastq|wc -l)/8|bc

# split reads into R1 and R2 fastq files
unmerge-paired-reads.sh cont_removed_non_rRNA3.fastq nonrRNA_decont_R1.fastq nonrRNA_decont_R2.fastq

echo $(cat nonrRNA_decont_R1.fastq|wc -l)/4|bc
echo $(cat nonrRNA_decont_R2.fastq|wc -l)/4|bc

# running assembly with decontaminated and rRNA removed

cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads
mkdir assemblies/Trinity_removed_cont_nonrRNA/
Trinity --seqType fq --max_memory 60G --SS_lib_type RF --min_kmer_cov 2 --left reads/cont_removed_R1.fastq --right reads/cont_removed_R2.fastq --CPU 14 --output assemblies/Trinity_removed_cont_nonrRNA/ 1>assemblies/Trinity_removed_cont_nonrRNA/Trinity.log 2>assemblies/Trinity_removed_cont_nonrRNA/Trinity.err

# use BUSCO to assess completeness of the assembly against arthropoda
#export BUSCO_CONFIG_FILE="/home/quelo/projects/paracletus/user_installed_software/busco/config/myconfig.ini"
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assemblies/Trinity_removed_cont_nonRNA
mkdir BUSCO
cd BUSCO
busco -c 14 -m transcriptome -i ../Trinity.fasta -o arthr -l arthropoda_odb10 1>busco_arth.log 2>busco_arth.err

busco -c 14 -m transcriptome -i ../Trinity.fasta -o bact -l bacteria_odb10 1>busco_bact.log 2>busco_bact.err

# get Trinity basic statistics
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assembly/Trinity_removed_cont_nonrRNA
mkdir TrinityStats
cd TrinityStats
TrinityStats.pl ../Trinity.fasta


# cd-hit remove redundancy
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assembly
mkdir cdhit cdhit/0.95
cd cdhit/0.95
cd-hit-est -i ../../Trinity_removed_cont_nonrRNA/Trinity.fasta -c 0.95 -o decont_nonrRNA_cdhit0.95.fasta -c 0.95 -n 8 -p 1 -g 1 -M 81920 -T 14 -d 40 1>cd-hit-est_0.95.log 2>cd-hit-est_0.95.err
# reduced the assembly to 411k contigs. Still too much!


#-----------------------------------------------------------------------------------------------------------------------#

# apply decontamination with less astringent bowtie2
 # bowtie2 mapping against contamination sequences database, keep both mapped and unmapped reads (paired-end reads)
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/mapping
mkdir contamination2
# copied again contamination.fna to contamination2 folder because I thought there was some problem with the old index
bowtie2-build contamination2/contamination.fna contamination2/contamination

# bowtie2 mapping against contamination sequences database, keep both mapped and unmapped reads (paired-end reads)
bowtie2 -x contamination2/contamination -p 14 -q --rf --mp 4 --rdg 4,2 --rfg 4,2 -1 ../../R1_pairedout.fastq -2 ../../R2_pairedout.fastq 2>contamination2/align_stats.txt | samtools view -@14 -bS -> contamination2/contamination_mapped_and_unmapped.bam
bowtie2 -x contamination/contamination -p 14 -q --rf -1 ../../R1_pairedout.fastq -2 ../../R2_pairedout.fastq 2>contamination/align_stats.txt | samtools view -@14 -bS -> contamination/contamination_mapped_and_unmapped.bam 
# filter required unmapped reads
# SAMtools SAM-flag filter: get unmapped pairs (both ends unmapped)
samtools view -@ 14 -b -f 12 -F 256 contamination2/contamination_mapped_and_unmapped.bam > contamination2/SAMPLE_bothEndsUnmapped.bam
# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort -n contamination/SAMPLE_bothEndsUnmapped.bam contamination/SAMPLE_bothEndsUnmapped_sorted
# split paired-end reads into separated fastq files .._r1 .._r2
bedtools bamtofastq -i contamination2/SAMPLE_bothEndsUnmapped_sorted.bam -fq ../filtered_reads/cont2_removed_R1.fastq -fq2 ../filtered_reads/cont2_removed_R2.fastq

# running assembly on reads with filtered out contamination2 (not endosymbionts, nor mitochondrial)
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads
mkdir assemblies assemblies/Trinity_removed_cont2/
Trinity --seqType fq --max_memory 60G --SS_lib_type RF --min_kmer_cov 2 --left reads/cont2_removed_R1.fastq --right reads/cont2_removed_R2.fastq --CPU 14 --output assemblies/Trinity_removed_cont2/ 1>assemblies/Trinity_removed_cont2/Trinity.log 2>assemblies/Trinity_removed_cont2/Trinity.err

# get Trinity basic statistics
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assembly/Trinity_removed_cont2
mkdir TrinityStats
cd TrinityStats
TrinityStats.pl ../Trinity.fasta > Stats.txt


# use BUSCO to assess completeness of the assembly against arthropoda
cd /home/quelo/projects/paracletus/cleaned_trimmed_reads/filtered_reads/assemblies/Trinity_removed_cont2
mkdir BUSCO
cd BUSCO
busco -c 14 -m transcriptome -i ../Trinity.fasta -o arthr -l arthropoda_odb10 1>busco_arth.log 2>busco_arth.err

busco -c 14 -m transcriptome -i ../Trinity.fasta -o bact -l bacteria_odb10 1>busco_bact.log 2>busco_bact.err

busco -c 14 -m transcriptome -i ../Trinity.fasta -o insect -l insecta_odb10 1>busco_insect.log 2>busco_insect.err


### DUBTES ###

grep "@" R1_pairedout.fastq | sed 's/\ 1//g'  > h1.txt
grep "@" R2_pairedout.fastq | sed 's/\ 2//g'  > h2.txt
diff h1.txt h2.txt > only_in_1.txt



grep "@" decont7_R1.fastq | sed 's/\/1//g' > h1.txt; grep "@" decont7_R2.fastq | sed 's/\/2//g' > h2.txt
diff h1.txt h2.txt > only_in_1.txt
diff h2.txt h1.txt > only_in_2.txt



#error relaxed decont → apparently, in contamination.fna there are duplicated fasta entries for these five sequences  of ecoli.
[W::sam_hdr_create] Duplicated sequence 'NC_000913.3'
[W::sam_hdr_create] Duplicated sequence 'NC_002695.2'
[W::sam_hdr_create] Duplicated sequence 'NC_002128.1'
[W::sam_hdr_create] Duplicated sequence 'NC_002127.1'
NC_000913.3
```
