#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V

set -e


# Set working directory
date=$(date)
workdir=$(pwd)


# Function to display usage
usage() 
{
cat << EOF

Summary: This script annotates piRNA clusters based on read density on a given genome. 

Usage: 
	bash $(basename "$0") [options] -i multimappers -u unique mappers -g chromosome file
	
	-i multimappers: a bed file containing all piRNA reads including single and multimapping reads
	-u unique mappers: a bed file containing only piRNA reads that map unambigously
	-g chromosome file: a tab-delimited file containing the size of all chromosomes (or contigs)
	
Options:
	-a: Minimal number of piRNAs that need to map to a window (default: 10)
	-d: Maximum distance between two windows to be merged (default: 5000)
	-c: Minimal number of uniquely mapping piRNAs that need to map to a cluster (default: 5)
	-p: Minimal number of unique piRNA positions per cluster (default: 5)
	-s: Minimal size of a cluster (default: 1000)
	-m: Minimal piRNA density per kb of cluster (default: 10)
	-o: Folder in which results will be stored. If not specified, the current working directory will be used.
		Note: if the specified folder is not empty, files might be overwritten!
	-l: Run script in loop-mode. The script will then add an additional output containing summary statistics of the clusters.
		This is useful for testing a series of cut-off values for different parameters and comparing the performance.
	-x: Run script in debugging-mode. In that case the script will print the code that is executed and will not delete temporary files.

		
EOF
}

# Function to check if file exists
check_file() 
{
    local file="$1"
	if [ ! -e "${file}" ]; then
        printf "***Error: File ${file} does not exist.***\n\n" 
		usage
        exit 1
	fi
}

# Function to check if input is integer
check_input() 
{
	if ! [[ "$1" =~ ^[0-9]+$ ]]; then
		printf "\n***Error: Expected an integer as input.***\n"
		usage
		exit 1
	fi
}


# Function to check if program is installed
check_program() 
{
    local program=$1
    if ! command -v "${program}" &>/dev/null; then
        printf "\n***Error: cannot find the program ${program}. Not installed or not found in PATH.**\n\n"
        exit 1
    fi
}

# Function to extract the first column from a file
get_first_column() {
    awk '{print $1}' "$1"
}


# Initialize required files:
input_file_multimappers=""
input_file_uniquemappers=""
chrom_file=""


# Default variables (if not specified otherwise by user)
min_pirnas=10
max_distance=5000
min_unique_pirnas=5
min_unique_pirna_positions=5
min_size_cluster=1000
min_pirna_density=10
out_directory="${workdir}"
loop="no"
debug="no"

# Parse options
while getopts ":i:u:g:a:d:c:p:s:m:o:lx" options; do
    case $options in
		i) 
			input_file_multimappers="$OPTARG"
			;;
		u) 
			input_file_uniquemappers="$OPTARG"
			;;
		g) 
			chrom_file="$OPTARG"
			;;
        a) 
			min_pirnas="$OPTARG"
			check_input "$min_pirnas"
			;;
        d) 
			max_distance="$OPTARG"
			check_input "$max_distance"
			;;
        c) 
			min_unique_pirnas="$OPTARG"
			check_input "$min_unique_pirnas"
			;;
        p) 
			min_unique_pirna_positions="$OPTARG"
			check_input "$min_unique_pirna_positions"
			;;
        s) 
			min_size_cluster="$OPTARG"
			check_input "$min_size_cluster"
			;;
        m) 
			min_pirna_density="$OPTARG"
			check_input "$min_pirna_density"
			;;
		o) 
			out_directory="$OPTARG"
			;;
		l) 
			loop="yes"
			;;
		x) 
			debug="yes"
			printf "Running in debugging mode. \n"
			;;
		:)
			printf "***Error: Option -${OPTARG} requires an argument.***\n\n"
			usage
			exit 1
			;;
        \?) 
			printf "Invalid option: -$OPTARG" >&2
            usage
            exit 1
			;;
    esac
done


# Check if the required files are provided

required_files=("${input_file_multimappers}" "${input_file_uniquemappers}" "${chrom_file}")
for file in "${required_files[@]}"; do
	check_file "${file}"
done


# Check: are required programs installed:
# List of required programs
required_programs=( "bedtools" "awk" )

# Check each required program
for prog in "${required_programs[@]}"; do
    check_program "${prog}"
done


##########################
### Cluster annotation ###
##########################

# print executed code if script is run in debuagging mode
if [ "${debug}" = "yes" ]; then
	set -x
fi


# Check if output directory exists, otherwise make
if [ ! -d "${out_directory}" ]; then
  mkdir "${out_directory}"
fi
	
# make temporary folder:
if [ ! -d "${out_directory}"/tmp ]; then
  mkdir "${out_directory}"/tmp
fi


# output files:
# storing the actual annotation
out_file=cluster_annotation.a"$min_pirnas".d"$max_distance".c"$min_unique_pirnas".p"$min_unique_pirna_positions"
touch "${out_directory}"/"${out_file}" 
# log file storing error messages and general info
out_log=log."${out_file}"
touch "${out_directory}"/"${out_log}" 
# summary file with statistics on clusters
if [ "${loop}" = "yes" ] && [ ! -f stats."${out_file}" ]; then
	out_stats="${out_directory}"/stats."${out_file}"
	touch "${out_directory}"/stats."${out_file}"
	printf "Min_piRNAs\tMax_dist\tMin_unique_piRNAs\tMin_unique_pos\tMin_size\tMin_density\tNo_clusters\tAv_length\tFrac_genome\tFrac_piRNAs_clusters\n" >> "${out_directory}"/stats."${out_file}"
fi



cat << EOF 2>&1 | tee -a "${out_directory}"/"${out_log}"
Cluster annotation was performed on $(date).

Clusters were annotated using the following files as input:
	Total piRNA reads: ${input_file_multimappers}
	Uniquely mapping piRNA reads: ${input_file_uniquemappers}
	Chromosome file: ${chrom_file}
The following parameters were used to annotate clusters:
	Minimal number of piRNAs per window: ${min_pirnas}.
	Maximal distance between windows to be joined: ${max_distance} bp.
	Minimal number of uniquely mapping piRNAs per cluster: ${min_unique_pirnas}.
	Number of uniquely mapping positions per cluster: ${min_unique_pirna_positions}.
	Minimal length of a cluster: ${min_size_cluster} bp.
	Minimal piRNA density per cluster: ${min_pirna_density} piRNAs/kb.
The output will be stored in the folder: ${out_directory}.
	
...Start cluster annotation....

EOF


# TO DO: remove piRNA-sized reads selection from script, but put this before?


# Make windows for cluster prediction
printf "Dividing genome into windows...\n"

bedtools makewindows -g "$chrom_file" -w 5000 |
	sort -k 1,1 -k 2,2n \
	> "${out_directory}"/tmp/file_tmp_genomic_windows
		
# Calculate normalization factor to scale number of piRNAs to million mapped piRNAs 
no_pirnas=`wc -l "${input_file_multimappers}" | awk '{ print $1 }'`
norm_fac=`echo "scale = 4; $no_pirnas / 1000000" | bc`


### 1. Determine windows with minimal piRNA coverage
#	
# 1.1. count piRNAs on windows, normalize to million piRNAs and
# 1.2. filter based on min piRNA cutoff
printf "Counting piRNAs on windows...\n"

bedtools coverage \
		-a "${out_directory}"/tmp/file_tmp_genomic_windows \
		-b "${input_file_multimappers}" \
		-F 0.51 \
		-counts |
	# normalize to million mapped piRNAs
	awk -v n="${norm_fac}" -v CONVFMT=%.3g 'BEGIN{ OFS = "\t" }{ $4 = $4/ n} 1' |
	awk -v a="${min_pirnas}" '( $4 >= a )' \
	> "${out_directory}"/tmp/file_tmp_filtered_windows

# 1.3. merge windows that reached threshold and are within the max distance to each other
printf "Merging windows...\n"
mergeBed \
		-i "${out_directory}"/tmp/file_tmp_filtered_windows \
		-d "${max_distance}" \
		-c 1 \
		-o count \
		> "${out_directory}"/tmp/file_tmp_merged_windows
				
# 1.4. count how many piRNAs are located on windows (total)
bedtools coverage \
		-a "${out_directory}"/tmp/file_tmp_merged_windows \
		-b "${input_file_multimappers}" \
		-F 0.51 \
		-counts > "${out_directory}"/tmp/file_tmp_windows_total_piRNAs


	
### 2. Filter clusters based on number of unique piRNAs and number of individual piRNA mapping positions
#
# 2.1. Count uniquely mapping piRNAs on filtered windows 
printf "Filtering clusters...\n"

bedtools coverage \
	-a "${out_directory}"/tmp/file_tmp_windows_total_piRNAs \
	-b "${input_file_uniquemappers}" \
	> "${out_directory}"/tmp/file_tmp_windows_uniq_piRNAs
	
# 2.2 normalize to million mapped piRNAs and
# 2.3 filter based on cutoff for minimal unique piRNAs and minimal unique piRNA mapping positions
awk -v n="$norm_fac" -v CONVFMT=%.3g 'BEGIN{ OFS = "\t" } 
		{ $5 = $5 / n; $6 = $6/ n} 1' "${out_directory}"/tmp/file_tmp_windows_uniq_piRNAs |
	# number of non-zero bases (in column 7) reported by bedtools coverage is
	# equal to the number of positions, as piRNA length is 1
	awk -v u="$min_unique_pirnas" -v p="$min_unique_pirna_positions" '( $6 >= u && $7 >= p )' |
	cut -d $'\t' -f 1-7 \
	> "${out_directory}"/tmp/file_tmp_merged_windows_filtered


### 3. Determine start/ end of cluster by first/last piRNA position on merged windows
#
# 3.1 Determine 5' most piRNA position
	# use only start position of merged windows
printf "Determining start and end positions of clusters...\n"
	
awk '{ OFS = "\t" }
		{ print $1, $2, $2 }' "${out_directory}"/tmp/file_tmp_merged_windows_filtered |
    closestBed \
		-a stdin \
		-b "${input_file_multimappers}" \
		-D ref \
		-iu \
		-t first |
    awk 'BEGIN { OFS = "\t" }
		{ $2 = $5; $3 = $6 } 1' |
	cut -d $'\t' -f 1-3 \
	> "${out_directory}"/tmp/file_tmp_clusters_start

# 3.2 Determine 3' most piRNA position
	# use only end position of merged windows
awk '{ OFS= "\t" }
		{ print $1, $3-1, $3 }' "${out_directory}"/tmp/file_tmp_merged_windows_filtered |
	closestBed \
		-a stdin \
		-b "${input_file_multimappers}" \
		-D ref \
		-fu \
		-t first |
	awk 'BEGIN { OFS = "\t" }
        { $2 = $5; $3 = $6 } 1' |
	cut -d $'\t' -f 1-3 \
	> "${out_directory}"/tmp/file_tmp_clusters_end


### 4. Combine information to get clusters (and info on mapping piRNAs)	

# 4.1 First a sanity check to make sure that files to be merged contain same clusters in each line:
mapfile -t col_file_clusters < <(get_first_column "${out_directory}"/tmp/file_tmp_merged_windows_filtered)
mapfile -t col_file_start < <(get_first_column "${out_directory}"/tmp/file_tmp_clusters_start)
mapfile -t col_file_end < <(get_first_column "${out_directory}"/tmp/file_tmp_clusters_end)

# 4.1.1 files with clusters, start and end positions need to have same number of lines and the dame first column:
lines_file_clusters=${#col_file_clusters[@]}
lines_file_start=${#col_file_start[@]}
lines_file_end=${#col_file_end[@]}

if [ "$lines_file_clusters" -ne "$lines_file_start" ] || [ "$lines_file_clusters" -ne "$lines_file_end" ] || [ "$lines_file_start" -ne "$lines_file_end" ]; then
    printf "Error: Files with clusters, start and end positions have different number of lines! Something went wrong! \n" 2>&1 | tee -a "${out_directory}"/"${out_log}"
    exit 1
fi

# 4.1.2 Check if first column is the same in all three files
if [ "${col_file_clusters[*]}" != "${col_file_start[*]}" ] || [ "${col_file_clusters[*]}" != "${col_file_end[*]}" ] || [ "${col_file_start[*]}" != "${col_file_end[*]}" ]; then
    printf "Error: First column of files with clusters, start and end positions is different! Something went wrong! \n" 2>&1 | tee -a "${out_directory}"/"${out_log}"
    exit 1
fi


# 4.2 Merge info on start and end positions into one file 
paste \
	<(cut -f 1,2 "${out_directory}"/tmp/file_tmp_clusters_start) \
	<(cut -f 3 "${out_directory}"/tmp/file_tmp_clusters_end) \
	<(cut -f 4-7 "${out_directory}"/tmp/file_tmp_merged_windows_filtered) \
	> "${out_directory}"/tmp/file_tmp_clusters
	
### 5. Filter clusters for minmal length and piRNA density
awk 'BEGIN{ OFS = "\t" }
		{ print $1, $2, $3, $4, $3-$2, "+", $5, $6, $7}' "${out_directory}"/tmp/file_tmp_clusters | # calculate cluster length, re-arrange column (and add strand for next step below)
	awk 'BEGIN { OFS = "\t" } 
		{ $10 = $7 / $5 *1000 } 1' | # piRNAs per kb
	awk -v s="$min_size_cluster" -v d="$min_pirna_density" '( $5 >= s && $10 >= d )' \
		> "${out_directory}"/tmp/file_tmp_clusters_filtered


### 6. Calculate piRNAs mapping on each strand
# 6.1 count how many piRNAs are on + strand, then calculate the - strand from that
printf "Calculating mapping statistics...\n"

bedtools coverage \
		-a "${out_directory}"/tmp/file_tmp_clusters_filtered \
		-b "${input_file_multimappers}" \
		-counts |
	bedtools coverage \
		-a stdin \
		-b "${input_file_multimappers}" \
		-counts \
		-s |
	awk 'BEGIN{ OFS="\t" }
		{ $13 = $11 - $12;
		  $14 = $12/ $11 * 100;
		  $15 = $13/ $11 *100 }
		{ print $1, $2, $3, $5, $11, $10, $7, $8, $9, $14, $15 }' |
	sort -nrk 4 \
	> "${out_directory}"/"${out_file}"

# Add a header to the out_file for easier interpretation
printf 'Chrom\tStart\tEnd\tCluster_length\tTotal_piRNAs\tpiRNA_density\tpiRNAs_per_million\tUnique_piRNAs_per_million\tUnique_piRNA_positions\tPerc_piRNAs_plus\tPerc_piRNAs_min\n' | 
	cat - "${out_directory}"/"${out_file}" \
	> "${out_directory}"/tmp/temp && mv "${out_directory}"/tmp/temp "${out_directory}"/"${out_file}"


### 7. calculate summary statistics for the annotation
genome_length=`awk '{ sum += $2 } END { print sum }' "${chrom_file}"`
total_length_clusters=`awk '{ sum += $4 } END { print sum }' "${out_directory}"/"${out_file}"`
no_cluster=`wc -l "${out_directory}"/"${out_file}" | cut -d $' ' -f 1`
no_pirnas_in_clusters=`awk '{ sum += $5 } END { print sum }' "${out_directory}"/"${out_file}"`
#no_pirnas_total_normalized=`echo "scale=4; ${no_pirnas} / ${norm_fac}" | bc`
fraction_pirnas_in_clusters=`echo "scale=4; ${no_pirnas_in_clusters} / ${no_pirnas} * 100" | bc`
average_length_clusters=`echo "scale=1; ${total_length_clusters} / ${no_cluster}" | bc`
genomic_space=`echo "scale=5; ${total_length_clusters} / ${genome_length} * 100 " | bc`  



# if run in loop-mode, print summary statistics in separate file for easy comparison
if [ "${loop}" = "yes" ]; then
	printf "${min_pirnas}\t${max_distance}\t${min_unique_pirnas}\t${min_unique_pirna_positions}\t${min_size_cluster}\t${min_pirna_density}\t$no_cluster\t$average_length_clusters\t$genomic_space\t$fraction_pirnas_in_clusters\n" >> "${out_directory}"/stats."${out_file}"
fi

# clean up temporary files
if [ ! "${debug}" = "yes" ]; then
	rm -r "${out_directory}"/tmp
fi



cat << EOF 2>&1 | tee -a "${out_directory}"/"${out_log}"

Summary of the results:

	Number of clusters annotated: ${no_cluster}
	Fraction of piRNAs in clusters: ${fraction_pirnas_in_clusters}% 
	Average length clusters: ${average_length_clusters} bp
	Fraction of genome covered by clusters: ${genomic_space}%
	
The output file ${out_directory}/${out_file} contains information on:

	Chrom: chromosome or scaffold name
	Start:  start coordinate of cluster (zero-based)
	End: end coordinate of cluster 
	Cluster_length: length of cluster [bp]
	Total_piRNAs: number of piRNAs (not normalized) mapping to a cluster
	piRNA_density: piRNAs per million piRNAs per kb mapping to a cluster
	piRNAs_per_million: number of piRNAs per million piRNAs mapping to a cluster
	Unique_piRNAs_per_million: number of unambigously mapping piRNAs per million piRNAs mapping to a cluster
	Unique_piRNA_positions: number of positions at which  piRNAs map unambigously 
	Perc_piRNAs_plus: fraction of piRNAs in cluster mapping on + strand [%]
	Perc_piRNAs_min: fraction of piRNAs in cluster mapping on - strand [%]
	
	
...Cluster annotation finished...

EOF

