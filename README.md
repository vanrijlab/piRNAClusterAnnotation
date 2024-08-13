# Annotation of piRNA clusters in non-model insect species
## Purpose of the pipeline
The piRNA pathway is crucial for protecting metazoan genomes from transposable elements, with piRNAs in *Drosophila melanogaster* primarily originating from defective transposons in piRNA clusters. Existing tools are optimized for annotating piRNA clusters in model organisms, and they often rely on assumptions that may not apply to non-model insects, where the piRNA pathway is less understood. We therefore implemented a simple annotation approach to determine piRNA clusters in non-model species that utilizes very little assumptions on the biology of the piRNA pathway. We validated it for *Aedes* mosquitoes, but it should be versatile enough to work with a variety of insect species.
We hope that it will be useful in exploring the piRNA pathway in different non-model insect species but we do not provide any guarantee that the pipeline provides provide any guarantee that the pipeline provides meaningful results for your application.  

<img src=https://github.com/user-attachments/assets/10bedb2f-72de-40e7-b204-96cfa4400f52 width=800>
  
  
## Installation
1.	Clone the repository
```bash
git clone https://github.com/vanrijlab/piRNAClusterAnnotation.git
cd ./piRNAClusterAnnotation
```
2.	Make the scripts executable with 
```bash
chmod +x ExtractpiRNAs.sh
chmod +x AnnotatepiRNAclusters.sh
```

## Requirements
Install [Bedtools software](https://bedtools.readthedocs.io/en/latest/) from the Quinlan lab before using the piRNA cluster annotation pipeline. Additionally, you will need the following files:
- Two BAM files containing small RNA reads mapped to the genome of interest. One file needs to contain all small RNA reads that map at least one time to the genome, while the other should contain the same data, however, reads that map ambigously (more than once) to the genome removed (For an example of how to map small RNA reads to a genome see the section [Example](https://github.com/vanrijlab/piRNAClusterAnnotation/edit/main/README.md#example) ) below.
- A tab-delimited chromosome file defining the chromosome lengths. This file requires two columns indicating the chromosome/ contig name (column 1) and the total length of the chromosome/ contig (column 2). Make sure that these are the same as in your BAM file. (See section [Example](https://github.com/vanrijlab/piRNAClusterAnnotation/edit/main/README.md#example) )

## Usage

### Extracting piRNA 5’ positions
The cluster annotation requires BED files that only contain the 5’ positions of piRNA-sized reads. To generate those, run:
```bash
bash ExtractpiRNAs.sh -i mappedsequences.bam [options]
```
Note that this step has to be performed for both BAM files (multimappers and uniquely mapping reads). The output is a BED file containing the 5’ positions of all reads within the specified size range.  
The command line options are:
- -i: BAM file containing small RNA reads
- -m: Minimal size of piRNA reads (default: 25)
- -M: Maximum size of piRNA reads (default: 30)
- -b: This option will keep an unfiltered BED file containing all small RNA reads
- -o: Folder in which the BED file will be stored. If not specified, the current working directory will be used.
- -x: Run script in debugging-mode. In that case the script will print the code that is executed.  

### Annotation of clusters
The script *AnnotatepiRNAclusters.sh* will perform the actual cluster annotation. Run this with:
```bash
bash AnnotatepiRNAclusters.sh \
	-i multimappers.piRNAs.bed \
	-u uniquemapper.piRNAs.bed \
	-g chromfile.tab \
	[options]
```
The command line options are:
- -i multimappers: a BED file containing all piRNA reads including single and multimapping reads (required)
- -u unique mappers: a BED file containing only piRNA reads that map unambiguously (required)
- -g chromosome file: a tab-delimited file containing the size of all chromosomes or contigs (required)
- -a: Minimal number of piRNAs that need to map to a window (default: 5)
- -d: Maximum distance between two windows to be merged (default: 5000)
- -c: Minimal number of uniquely mapping piRNAs that need to map to a cluster (default: 5)
- -p: Minimal number of unique piRNA positions per cluster (default: 5)
- -s: Minimal size of a cluster (default: 1000)
- -m: Minimal piRNA density per kb of cluster (default: 10)
- -o: Folder in which results will be stored. If not specified, the current working directory will be used. Note: if the specified folder is not empty, files might be overwritten!
- -l: Run script in loop-mode. The script will then add an additional output containing summary statistics of the clusters. This is useful for testing a series of cut-off values for different parameters and comparing the performance.
- -x: Run script in debugging-mode. In that case the script will print the code that is executed and will not delete temporary files.  

## Output
The pipeline will output two files, or three files if run in loop-mode:
- *log file*: short summary of the run
- *cluster file*: all annotated piRNA clusters:   
	   
	|  Column	 | Column name					 | Description	 |  
	|  ------	 | ------						 | ------			 |           
	|  1		 | Chrom						 | Chromosome or scaffold name	 |  
	|  2		 | Start						 | Start coordinate of cluster (zero-based coordinate system)	 |  
	|  3 		 | End 							 | End coordinate of cluster	 |  
	|  4		 | Cluster_length				 | Length of the cluster [bp]	 |  
	|  5		 | Total_piRNAs					 | Number of piRNA reads (not normalized) mapping to a cluster	 |  
	|  6		 | piRNA_density				 | PiRNAs per million piRNA reads per kb mapping to a cluster	 |  
	|  7		 | piRNAs_per_million			 | Number of piRNAs per million piRNAs reads mapping to a cluster	 |  
	|  8		 | Unique_piRNAs_per_million	 | Number of unambiguously mapping piRNAs per million piRNA reads mapping to a cluster	 |  
	|  9		 | Unique_piRNA_positions		 | number of positions at which  piRNAs map unambiguously	 |  
	|  10		 | Perc_plus					 | Fraction of piRNAs in a cluster mapping to the positive strand [%]	 |  
	|  11		 | Perc_min						 | Fraction of piRNAs in a cluster mapping to negative strand [%]	 |  
- *stat file*: summary statistics of the annotation (in loop-mode only):    
  
	| Column	| Column name	| Description	|  
	| ------		| ------			| ------		|  
	| 1		| Min_piRNAs	| Threshold value used for the annotation	|   
	| 2	| Max_dist	|	Threshold value used for the annotation	|  	
	| 3	| Min_unique_piRNAs	| Threshold value used for the annotation	|  	
	| 4	| Min_unique_pos		| Threshold value used for the annotation	|  
	| 5	| Min_size	| Threshold value used for the annotation	|  	
	| 6	| Min_density 	| Threshold value used for the annotation	|  	
	| 7	| No_clusters	|	Number of annotated clusters	|  
	| 8	| Av_length	| Average length of all clusters	|  
	| 9	| Frac_genome	| Fraction of the genome that is covered by piRNA clusters [%]	|  
	| 10	| Frac_piRNAs_clusters	| Fraction of piRNAs that map within annotated clusters [%]	|  
  
    
	
After annotation of the piRNA clusters, the clusters can be visually inspected using IGV or the UCSC Genome browser. Here is an example of a piRNA cluster in *Ae. albopictus* (SWKY01000101:21730-118879):  
  
  
![piRNACluster_example](https://github.com/vanrijlab/General/assets/29331754/df1d78b8-446b-466e-a2ca-98fc913adc4b)

## Example 

#### Mapping
This example is based on small RNA sequencing data from female ovaries of the Asian tiger mosquito *Aedes albopictus*. The raw data is available at the NCBI Sequence Read Archive under the accession number SRR12512270. For this dataset, the adapters are already clipped.  We will map the clipped reads with bowtie to the *Aedes albopictus* [reference sequence]( https://vectorbase.org/common/downloads/release-68/AalbopictusFoshanFPA/fasta/data/VectorBase-68_AalbopictusFoshanFPA_Genome.fasta) available at [Vectorbase]( https://vectorbase.org/vectorbase/app/downloads) (It is possible to use other software as well). We do not allow mismatches to restrict the analysis to high-confidence reads only. However, if your small RNA dataset is not derived from the reference strain the reference annotation is based on, or if a reference genome for your species of interest is not available and you need to map to a closely related species, you might need to increase the number of mismatches that are allowed.  
The mapping step needs to be performed twice:
```bash
# First mapping: considering all mapped reads that map without mismatches
# It is important to randomly distribute multimapping reads with the option --best --strata -M 1 
# Mismatches to the reference genome are not allowed (-v 0)

bowtie genome_index clipped_reads.fastq.gz \
	--best --strata \
	-v 0 -M 1 -S | \
	# convert to BAM file
	samtools view -Sb -F 4 - | \
	# sort BAM file
	samtools sort - -o multimappers.bam
	
# Second mapping: discard reads that map more than one time with option -m 1

bowtie genome_index clipped_reads.fastq.gz \
	-v 0 -m 1 -S | \
	samtools view -Sb -F 4 - | \
	samtools sort - -o uniquemappers.bam 
```
Since these files are large, we will continue only with reads that map to the scaffold SWKY01000014 which corresponds to a region on chromosome 3. These two BAM files are provided as example data and can be downloaded from the ExampleData folder on this repository. 

#### Extracting piRNA reads
The next step is to isolate piRNA-sized reads, and to extract their 5’ end position. However, it can be useful to first remove all reads that potentially contaminate the fraction of genuine piRNAs. For that, first get the genomic positions of all tRNA, rRNA and miRNA genes from a GFF3 gene annotation file with: 

```bash
awk '$3 ~ /rRNA|tRNA|miRNA/' annotationfile.gff > unwanted.gff
```

and then remove all reads with:

```bash
bedtools intersect -a multimappers.bam -b unwanted.gff -v > multimappers_filtered.bam 
```

For the sake of this example, however, we just skip this step and continue with the extraction of piRNA-sized reads from both BAM files. For *Ae. albopictus* we will use 25-30nt.

```bash
# First for all multimapping reads

bash ExtractpiRNAs.sh -i multimappers.bam -m 25 -M 30

# Then run again for all uniquely mapping reads

bash ExtractpiRNAs.sh -i uniquemappers.bam -m 25 -M 30 
```

The resulting files will be called multimappers.piRNAs.bed and uniquemappers.piRNAs.bed and serve as input for the cluster annotation pipeline.  

Note that the size of piRNAs can vary between species. For example, *D. melanogaster* produces shorter piRNAs than *Ae. albopictus* (average of 25nt vs 28nt). If you are not sure what size range to select for your species, it is a good idea to create a simple size profile of all your reads and use this to guide your decision. For that simply run:
```bash
samtools view multimappers.bam | cut -f10 | awk '{ print length }' | sort | uniq -c > sizeprofile.tab
```
and plot the data as bar chart: 
  
<img src="https://github.com/vanrijlab/General/assets/29331754/18da936e-719f-4306-8fc1-9c2709170457" width =500>

As you can see in this example, especially in the ovary of *Ae. albopictus* a peak of longer small RNAs corresponding to piRNAs can be distinguished from miRNAs and siRNAs. Most of those reads fall in the range of ~25-30nt. However, keep in mind that piRNA reads might not be as abundant or clearly distinguishable in every tissue (as this is the case for small RNAs from carcasses of *Ae. albopictus*).  
 
   
  
#### Generating a chromosome file  
In addition to the piRNA BED files, we require a chromosome file that indicates the lengths of all chromosomes (or in this case scaffolds) for the cluster annotation. You can easily create it from one of the BAM files with:

```bash
samtools view -H multimappers.bam | \
	grep @SQ | \
	sed 's/@SQ\tSN:\|LN://g' > chromfile.tab
```
#### Annotating piRNA clusters
Let’s say we want to annotate piRNA clusters using the default values. Then simply use:
```bash
bash AnnotatepiRNAclusters.sh \
	-i multimappers.piRNAs.bed \
	-u uniquemappers.piRNAs.bed \
	-g chromfile.tab
```
With the provided example datasets, around 90% of all piRNAs should be located in 128 annotated piRNA clusters. If we had used the entire dataset, that fraction would be lower though (~43% of piRNAs mapping to clusters).
#### Determining optimal threshold values
However, what if the optimal threshold values are not known *a priori*? Then run the annotation pipeline within a ```for``` loop for testing multiple different values (or even all combinations of values) and extract summary statistics that can guide the decision. Keep in mind that this can take a long time to run, depending on how large your files are and how many different combinations of threshold values you want to test.  
Run the pipeline using the loop-mode:

```bash
# Testing multiple values for the minimal amount of piRNAs that need to map to a window (1, 5, 10, 50 and 100)
# and maximal distance between windows (5000 and 10000 bp)
# while all other thresholds values will remain default

for min_piRNAs in 1 5 10 50 100; do
	for dist in 1000 5000 10000; do
		bash AnnotatepiRNAclusters.sh \
			-i multimappers.piRNAs.bed \
			-u uniquemappers.piRNAs.bed \
			-g chromfile.tab \
			-a "${min_piRNAs}" \
			-d "${dist}" 
	done
done
```
The resulting stat file can be used to select for the best combination of thresholds. A general assumption based on studies in fruit flies, silkworm, zebrafish, and mice is that piRNA clusters are, in general, large regions that produce a considerable part of the total piRNA pool. Therefore, the threshold values should be optimized so that a substantial fraction of piRNAs map to large clusters (*Fract_piRNAs_clusters* and *Av_length* should be maximized). At the same time, the fraction of the genome that is covered by clusters should be relatively low in order to keep the clusters compact (*Frac_genome* should be minimized). The decisions on threshold values will have to strike a fine balance between completeness and accuracy of the cluster annotation and the optimal conditions may additionally differ depending on the aim of your study. A visual example can be found in the supplementary materials of [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02141-w).   
## Publications

### Description of the piRNA cluster annotation

- Halbach, R & van Rij, RP. Annotation of piRNA source loci in the genome of non-model insects. Methods in Molecular Biology (in print)

### Related publications

- Palatini, U., et al. Improved reference genome of the arboviral vector *Aedes albopictus*. Genome Biol 21, 215 (2020). https://doi.org/10.1186/s13059-020-02141-w
- Crava, C.M., et al. Population genomics in the arboviral vector *Aedes aegypti* reveals the genomic architecture and evolution of endogenous viral elements. Mol Ecol, 30: 1594-1611 (2021). https://doi.org/10.1111/mec.15798
  
## Contact

If you encounter any issues or have questions, please open an issue on GitHub or contact us at rebecca [dot] halbach [at] radboudumc [dot] nl.


