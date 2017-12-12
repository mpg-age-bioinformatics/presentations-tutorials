#!/bin/bash

mkdir -p ~/tuxedo-update/slurm_logs
LOGS=$(readlink -f ~/tuxedo-update/slurm_logs)/

module load shifter
SHIFTER="shifter --image=mpgagebioinformatics/bioinformatics_software:v1.1.1"

# ************************************************************************
# download raw data
mkdir -p ~/tuxedo-update/raw_data
cd ~/tuxedo-update/raw_data

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32038/suppl/GSE32038%5Fsimulated%5Ffastq%5Ffiles%2Etar%2Egz
tar -zxvf GSE32038_simulated_fastq_files.tar.gz

# ************************************************************************

# download reference genome and gene models

mkdir -p ~/tuxedo-update/references
cd ~/tuxedo-update/references

wget ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/dna//Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-90/gtf/drosophila_melanogaster//Drosophila_melanogaster.BDGP6.90.gtf.gz

unpigz Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
unpigz Drosophila_melanogaster.BDGP6.90.gtf.gz

# ************************************************************************

# index reference genome

cd ~/tuxedo-update/references

${SHIFTER} << SHI
#!/bin/bash
source ~/.bashrc
module load hisat
hisat2-build Drosophila_melanogaster.BDGP6.dna.toplevel.fa Drosophila_melanogaster.BDGP6.dna.toplevel.fa
SHI

# ************************************************************************

# relabel samples with added metadata

cd ~/tuxedo-update/raw_data

ln -s GSM794483_C1_R1_1.fq.gz S_001-F_TuUp-L____C1-___-____-REP_1-READ_1.fastq.gz
ln -s GSM794483_C1_R1_2.fq.gz S_001-F_TuUp-L____C1-___-____-REP_1-READ_2.fastq.gz
ln -s GSM794484_C1_R2_1.fq.gz S_002-F_TuUp-L____C1-___-____-REP_2-READ_1.fastq.gz
ln -s GSM794484_C1_R2_2.fq.gz S_002-F_TuUp-L____C1-___-____-REP_2-READ_2.fastq.gz
ln -s GSM794485_C1_R3_1.fq.gz S_003-F_TuUp-L____C1-___-____-REP_3-READ_1.fastq.gz
ln -s GSM794485_C1_R3_2.fq.gz S_003-F_TuUp-L____C1-___-____-REP_3-READ_2.fastq.gz

ln -s GSM794486_C2_R1_1.fq.gz S_004-F_TuUp-L____C2-___-____-REP_1-READ_1.fastq.gz
ln -s GSM794486_C2_R1_2.fq.gz S_004-F_TuUp-L____C2-___-____-REP_1-READ_2.fastq.gz
ln -s GSM794487_C2_R2_1.fq.gz S_005-F_TuUp-L____C2-___-____-REP_2-READ_1.fastq.gz
ln -s GSM794487_C2_R2_2.fq.gz S_005-F_TuUp-L____C2-___-____-REP_2-READ_2.fastq.gz
ln -s GSM794488_C2_R3_1.fq.gz S_006-F_TuUp-L____C2-___-____-REP_3-READ_1.fastq.gz
ln -s GSM794488_C2_R3_2.fq.gz S_006-F_TuUp-L____C2-___-____-REP_3-READ_2.fastq.gz

# ************************************************************************

# Quality control raw data

mkdir -p ~/tuxedo-update/fastqc_output
cd ~/tuxedo-update/raw_data

for file in $(ls S_*.fastq.gz);
  do sbatch << EOF
#!/bin/bash
#SBATCH -p himem,hugemem,blade,dontuseme
#SBATCH -o ${LOGS}${file%.fastq.gz}.fastqc.out
#SBATCH -c 2

${SHIFTER} << SHI
#!/bin/bash
source ~/.bashrc
module load fastqc

fastqc -o ../fastqc_output ${file}

SHI
EOF
done


# ************************************************************************

# Mapping & Assembly & Quantification

mkdir -p ~/tuxedo-update/hisat2_output
mkdir -p ~/tuxedo-update/stringtie_output

cd ~/tuxedo-update/raw_data

R=$(readlink -f ../references/Drosophila_melanogaster.BDGP6.dna.toplevel.fa)

IDS=""
for file in $(ls S_*READ_1.fastq.gz);
  do ID=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH -p himem,hugemem,blade,dontuseme
#SBATCH -o ${LOGS}${file%-READ_1.fastq.gz}.hisat.stringtie.out
#SBATCH -c 8

${SHIFTER} << SHI
#!/bin/bash
source ~/.bashrc

# mapping with hisat2
module load hisat

hisat2 -p 8 --dta-cufflinks \
--met-file ../hisat2_output/${file%-READ_1.fastq.gz}.stats \
-x ${R} -1 ${file} -2 ${file%1.fastq.gz}2.fastq.gz \
-S ../hisat2_output/${file%-READ_1.fastq.gz}.sam

# compress and sort alignments
module load samtools

cd ~/tuxedo-update/hisat2_output
samtools view -@ 8 -bhS -F 4 ${file%-READ_1.fastq.gz}.sam | samtools sort -@ 8 -o ${file%-READ_1.fastq.gz}.bam -

# assemly and quantification with stringtie
module load stringtie

stringtie ${file%-READ_1.fastq.gz}.bam \
-o ../stringtie_output/${file%-READ_1.fastq.gz}.gtf \
-p 18 -f 0.99 \
-G ../references/Drosophila_melanogaster.BDGP6.90.gtf \
-C ../stringtie_output/${file%-READ_1.fastq.gz}_full_cov.gtf \
-b ../stringtie_output/${file%-READ_1.fastq.gz}

SHI
EOF
)

IDS=${IDS}:${ID}
done

# ************************************************************************

# wait for jobs to finish

echo "Waiting for jobs to finish" ${IDS}
srun -p himem,hugemem,blade,dontuseme -d afterok${IDS} echo "jobs completed"

# ************************************************************************

# Merging assemblies

mkdir -p ~/tuxedo-update/cuffmerge_output
cd ~/tuxedo-update/stringtie_output

rm -rf assemblies.txt
for file in $(ls *.gtf | grep -v _full_cov);
  do readlink -f ${file} >> assemblies.txt
done

# bug: S_002-F_TuUp-L____C1-___-____-REP_2.gtf shows 2x the transcript FBtr0346424

cat S_002-F_TuUp-L____C1-___-____-REP_2.gtf | grep -v FBtr0346424 > S_002-F_TuUp-L____C1-___-____-REP_2.gtf_
mv S_002-F_TuUp-L____C1-___-____-REP_2.gtf_ S_002-F_TuUp-L____C1-___-____-REP_2.gtf

${SHIFTER} << SHI
#!/bin/bash
source ~/.bashrc
module load cufflinks
cuffmerge -o ../cuffmerge_output --min-isoform-fraction 1.0 \
-g ../references/Drosophila_melanogaster.BDGP6.90.gtf \
-s ../references/Drosophila_melanogaster.BDGP6.dna.toplevel.fa \
assemblies.txt
SHI

# ************************************************************************

# Post assembly merge and quantification

mkdir -p ~/tuxedo-update/cuffquant_output
cd ~/tuxedo-update/hisat2_output

IDS=""
for file in $(ls *.bam);
  do ID=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH -p himem,hugemem,blade,dontuseme
#SBATCH -o ${LOGS}${file%.bam}.cuffquant.out
#SBATCH -c 8

${SHIFTER} << SHI
#!/bin/bash
source ~/.bashrc

module load cufflinks

cuffquant -p 8 --library-type fr-unstranded -o ../cuffquant_output/${file%.bam} \
../cuffmerge_output/merged.gtf ${file}

SHI
EOF
)

IDS=${IDS}:${ID}

done

# ************************************************************************

# wait for jobs to finish

echo "Waiting for jobs to finish"
srun -p himem,hugemem,blade,dontuseme -d afterok${IDS} echo "jobs completed"

# ************************************************************************

# Differential analysis

mkdir -p ~/tuxedo-update/cuffdiff_output
cd ~/tuxedo-update/cuffquant_output

IDS=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH -p himem,hugemem,blade,dontuseme
#SBATCH -o ${LOGS}cuffdiff.out
#SBATCH -c 8

${SHIFTER} << SHI
#!/bin/bash
source ~/.bashrc

module load cufflinks

cuffdiff -p8 --library-type fr-unstranded \
-L C1,C2 \
-o ../cuffdiff_output \
--dispersion-method per-condition \
../cuffmerge_output/merged.gtf \
S_001-F_TuUp-L____C1-___-____-REP_1/abundances.cxb,S_002-F_TuUp-L____C1-___-____-REP_2/abundances.cxb,S_003-F_TuUp-L____C1-___-____-REP_3/abundances.cxb \
S_004-F_TuUp-L____C2-___-____-REP_1/abundances.cxb,S_005-F_TuUp-L____C2-___-____-REP_2/abundances.cxb,S_006-F_TuUp-L____C2-___-____-REP_3/abundances.cxb
SHI
EOF
)

# ************************************************************************

# wait for jobs to finish

echo "Waiting for jobs to finish:"${IDS}
srun -p himem,hugemem,blade,dontuseme -d afterok:${IDS} echo "jobs completed"

# ************************************************************************

# QC and transcriptome wide plots

mkdir -p ~/tuxedo-update/cummeRbund_output
cd ~/tuxedo-update/

INPUT=$(readlink -f cuffdiff_output)
OUTPUT=$(readlink -f cummeRbund_output)

wget https://raw.githubusercontent.com/mpg-age-bioinformatics/htseq-tools/master/QC.R
chmod +x QC.R

sbatch << EOF
#!/bin/bash
#SBATCH -p himem,hugemem,blade,dontuseme
#SBATCH -o ${LOGS}cummeRbund.out

${SHIFTER} << SHI
#!/bin/bash
source ~/.bashrc

module load rlang
rm -rf ${INPUT}cuffData.db
./QC.R ${INPUT} ${OUTPUT}

SHI
EOF

# ************************************************************************

# Annotation, enrichment and network building

mkdir -p ~/tuxedo-update/adiff_output
cd ~/tuxedo-update

sbatch << EOF
#!/bin/bash
#SBATCH -p himem,hugemem,blade,dontuseme
#SBATCH -o ${LOGS}adiff.out

${SHIFTER} << SHI
#!/bin/bash
source ~/.bashrc

module load python

export PATH=~/.software_container/.python/2.7.13/bin:$PATH

aDiff -D -i cuffdiff_output -o adiff_output \
-G references/Drosophila_melanogaster.BDGP6.90.gtf \
-C cuffmerge_output/merged.gtf \
--dataset dmelanogaster_gene_ensembl \
--filter flybase_gene_id \
--outputBiotypes 'flybase_gene_id gene_biotype' \
--outputGoterms 'flybase_gene_id go_id name_1006' \
--DAVIDid FLYBASE_GENE_ID \
--DAVIDuser jorge.boucas@age.mpg.de \
--organismtag DMEL \
--species='drosophila melanogaster' \
--cytoscape_host=localhost \
--cytoscape_port=1234

SHI
EOF
