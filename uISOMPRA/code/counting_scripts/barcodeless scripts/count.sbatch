#!/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-144%20
#SBATCH --output=./Slurm_Outputs/%A_%a.out
#SBATCH --error=./Slurm_Outputs/%A_%a.err

#We will make the folders to drop the files into
mkdir -p "./Counts"
mkdir -p "./Alignment_Summaries"
mkdir -p "./Slurm_Outputs"

#Load the tools
eval $( spack load --sh bowtie2@2.4.2 )
eval $( spack load --sh samtools@1.13 )

#These eval lines load the spack packages for bowtie2 and samtools.  The SBATCH comments above set the number of CPUs and memory per cpu, and the array sets how many items there (1-144),
#and how many to process at one time to not be a dick (20).

#We will make the folders to drop the files into
mkdir -p "./Counts"
mkdir -p "./Alignment_Summaries"

#This line takes the lookup file and designates that the first entry is the read1 file, the second is the read2 file, and the third is what we're calling the output (matched to the actual sample).
read read1	read2	sample < <( sed -n ${SLURM_ARRAY_TASK_ID}p lookup.txt )

#This specifies run bowtie with 8 processors, using the reference genome we created earlier with bowtie2-build on the reference_genome.fasta file, with -1 and -2 specifying read 1 and read 2 from the
#lookup table, we want very sensitive matching, with local alignments (allows end trimming of read), and no mixed or discordant reads (both reads in a position wise manner on the chip
#have to agree and align).  We then save the output as a SAM file and the alignment summary separately.
bowtie2 -p 12 -x reference -1 ./Reads/${read1} -2 ./Reads/${read2} --very-sensitive-local --no-mixed --no-discordant 1>${sample}.SAM 2>./Alignment_Summaries/${sample}.alignment_summary.txt &&

# This next line of samtools arguments is important.  It is designed, using awk, to extract the SAM file header, and any reads that aligned perfectly to the start of the reference gene.
# It does so by searching for a 1 in the fourth position, which is a base 1 position value of where the read started aligning to the reference sequence after end clipping.
# Reads that did not align to the start of the reference gene will be thrown out.
# It also searches for a NM:i:0 in the eighteenth position.  This indicates no mismatches between the reference gene and the read.
# Filtering on these two features ensures ONLY exact matches come back, and gives the same result as grepping for the sequences.

samtools view -h ${sample}.SAM | awk '($1 ~ /^@/ || ($4==1 && $18=="NM:i:0"))' > ${sample}.exact.SAM &&

#Now we convert to BAM, sort, index, and get counts

samtools view -bS -o ${sample}.BAM ${sample}.exact.SAM &&
samtools sort ${sample}.BAM -@ 12 -o ${sample}.sorted.BAM &&
samtools index -@ 12 ${sample}.sorted.BAM &&
samtools idxstats -@ 12 ${sample}.sorted.BAM > ./Counts/${sample}.counts.txt &&
rm ${sample}.SAM ${sample}.exact.SAM ${sample}.BAM ${sample}.sorted.BAM ${sample}.sorted.BAM.bai
