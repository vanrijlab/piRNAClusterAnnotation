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

Summary: This script extracts the 5' ends of piRNA reads from a bam file. 

Usage: 
	bash $(basename "$0") [options] -i bamfile 
	
	-i: bam file
		
Options:
	-m: Minimal size of piRNA reads (default: 25)
	-M: maximum size of piRNA reads (default: 30)
	-b: store unfiltered bed file
	-o: Folder in which file will be stored. If not specified, the current working directory will be used.
	-x: Run script in debugging-mode. In that case the script will print the code that is executed.

Output:
	bed file containing 5' positions of piRNA-sized reads
	
EOF
}


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

# Initialize required files:
input_bam=""

# Default variables (if not specified otherwise by user)
min_size=25
max_size=30
out_directory="${workdir}"
debug="no"
keep_bed="no"

# Parse options
while getopts ":i:m:M:o:bx" options; do
    case $options in
		i) 
			input_bam="$OPTARG"
			;;
		m) 
			min_size="$OPTARG"
			check_input "$min_size"
			;;
		M) 
			max_size="$OPTARG"
			check_input "$max_size"
			;;
		o) 
			out_directory="$OPTARG"
			;;
		b) 
			keep_bed="yes"
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


# Check if the required bam file is provided
check_file "${input_bam}"



# Check: are required programs installed:
# List of required programs
required_programs=( "awk" "bedtools")

# Check each required program
for prog in "${required_programs[@]}"; do
    check_program "${prog}"
done


# print executed code if script is run in debuagging mode
if [ "${debug}" = "yes" ]; then
	set -x
fi


# Check if output directory exists, otherwise make
if [ ! -d "${out_directory}" ]; then
  mkdir "${out_directory}"
fi

# Extract filename of input bam file
file="$(basename -- $input_bam)"
filename="$(echo ${file%.*})"

# Convert bam file to bed file
bedtools bamtobed -i "${input_bam}" > "${out_directory}"/"${filename}".bed

# Extract piRNA-sized reads and determine 5' ends

awk -v n="${min_size}" -v m="${max_size}" 'BEGIN { OFS = "\t" }
        { $7 = $3 - $2 }
		{ if(( $7>=n ) && ( $7<=m )){ print }}' "${out_directory}"/"${filename}".bed | 
	awk 'BEGIN { OFS = "\t" } 
		{ if( $6=="+" ){ $3= $2 + 1 } else{ $2= $3 - 1 }} 1' |
	sort -k1,1 -k2,2n |
	awk 'NF > 0' > "${out_directory}"/"${filename}".piRNAs.bed
	

# delete bed file with all small RNA reads
if [ ! "${keep_bed}" = "yes" ]; then
	rm "${out_directory}"/"${filename}".bed
fi
