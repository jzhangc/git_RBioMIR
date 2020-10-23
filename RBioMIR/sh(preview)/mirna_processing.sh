#!/bin/bash
# Name: mir_processing.sh
# Version: 0.1.9
# Usage: ./mir_rocessing.sh <INPUTS> [OPTIONS]
# NOTE: In shell, 0 is true, and 1 is false - reverted from other languages like R and python.

# ------ Version history ------
# 0.1.9: Optional -o flag added for setting an output direcotry.
# 0.1.8: Quotations added for calling filenaem variables
# 0.1.7: Small cosmatic and editorial changes
# 0.1.6: Linux system support added
# 0.1.5: Config file support added
# 0.1.0-0.1.4: Initial versions

# ------ variables ------
# internal system variables
VERSION="0.1.9"
CURRENT_DAY=$(date +%d-%b-%Y)
PLATFORM="Unknown UNIX or UNIX-like system"
UNAMESTR=`uname`  # use `uname` variable to detect OS type
if [ $UNAMESTR == "Darwin" ]; then
	PLATFORM="macOS"
elif [ $UNAMESTR == "Linux" ]; then
	PLATFORM="Linux"
fi
HELP="\n
Format: ./mir_rocessing.sh <INPUTS> [OPTIONS] \n
Current version: $VERSION\n
\n
-h, --help: This help information. \n
--version: Display current version number. \n
\n
<INPUTS>: Mandatory \n
-i <dir>: Input raw file folder \n
-a <file>: Adapter file with full path\n
-n <file>: Negative reference with full path\n
-m <file>: Mature miRNA reference with full path\n
-s <file>: Stemloop/hairpin reference with full path\n
\n
[OPTIONS]: Optional \n
-o <dir>: Optional output directory. Default is the current folder. \n
-C <int>: Cut nucleotides from 5'\n
-c <int>: Cut nucleotides from 3'\n
-t: Change U to T for the reference files. Original file renamed with .bak \n
-p <int>: CPU core (thread) number\n"
CITE="Written by Jing Zhang PhD
Contact: jzha9@uwo.ca, jzhangcad@gmail.com
To cite in your research: TBA\n"

# colours
COLOUR_YELLOW="\033[1;33m"
COLOUR_ORANGE="\033[0;33m"
COLOUR_RED="\033[0;31m"
COLOUR_GREEN_L="\033[1;32m"
NO_COLOUR="\033[0;0m"

# flag variables with initial values
IFLAG=1 # initiate mandatory variable check variable. initial value 1 (false)
AFLAG=1 # initiate mandatory variable check variable. initial value 1 (false)
NFLAG=1 # initiate mandatory variable check variable. initial value 1 (false)
MFLAG=1 # initiate mandatory variable check variable. initial value 1 (false)
SFLAG=1 # initiate mandatory variable check variable. initial value 1 (false)
NTCUT=1  # initiate flag if to cut the nucleotides. initial value 1(false)
U2T=1  # initiate flag if to change U to T for reference files. initial value 1 (false)


# initial application variables from cammand flags
OUT_DIR=.
CORES=1  # initiate CPU core number as 1
NT5=0  # initiate variable for cutting nucleotides from 5'
NT3=0  # initiate variable for cutting nucleotides from 3'

# flag check and set/reset application variables from cammand flags
if [ $# -eq 0 ]; then # without any argument: display help file
	echo -e $HELP
	echo -e "\n"
  echo -e "=========================================================================="
	echo -e "${COLOUR_YELLOW}$CITE${NO_COLOUR}\n"
  exit 0  # exit 0: terminating without error. FYI exit 1 - exit with error, exit 2 - exit with message
else
	case "$1" in  # "one off" flags
		-h|--help)
			echo -e $HELP
			echo -e "\n"
			echo -e "=========================================================================="
			echo -e "${COLOUR_YELLOW}$CITE${NO_COLOUR}\n"
			exit 0
			;;
		--version)
			echo -e "Current version: $VERSION\n"
			exit 0
			;;
	esac

	while getopts ":ti:a:o:C:c:n:m:s:p:" opt; do
		case $opt in
			t)  # flag to change U to T for reference files
				U2T=0
				;;
			i)
				RAW_FILES=$OPTARG
				IFLAG=0
				;;
			a)
				ADAPTER=$OPTARG
				if ! [ -f "$ADAPTER" ]; then
					# >&2 means assign file descripter 2 (stderr). >&1 means assign to file descripter 1 (stdout)
					echo -e "${COLOUR_RED}ERROR: -a requires an adapter file in fasta or fa format, or file not found.${NO_COLOUR}\n" >&2
					exit 1  # exit 1: terminating with error
				fi
				AFLAG=0
				;;
			o)
				OUT_DIR=$OPTARG
				if ! [ -d "$OPTARG" ]; then
					echo -e "${COLOUR_YELLOW}\nWARNING: -o output direcotry not found. use the current directory instead.${NO_COLOUR}\n"
					OUT_DIR=.
				fi
				;;
			C)  # flag to decide if to cut the nucleotides
				NT5=$OPTARG
				NTCUT=0  # set to 0(true) with the presence of -C
				;;
			c)  # flag to decide if to cut the nucleotides
				NT3=-$OPTARG
				NTCUT=0  # set to 0(true) with the presence of -c
				;;
			n)  # flag to build bowtie index for negative ref
				NREF=$OPTARG
				if ! [ -f "$NREF" ]; then  # if not ...
				echo -e "${COLOUR_RED}ERROR: -n requires a negative reference file in fasta or fa format, or file not found.${NO_COLOUR}\n" >&2
				exit 1
				fi
				NFLAG=0
				# extract reference file names
				NREF_NAME=`basename "$NREF"`
				NREF_PATH=`dirname "$NREF"`
				NREF_NAME_WO_EXT="${NREF_NAME%%.*}"
				NREF_PATH_NAME_WO_EXT="$NREF_PATH/$NREF_NAME_WO_EXT"
				;;
			m)  # flag to build bowtie index for mature miRNA ref
				MREF=$OPTARG
				if ! [ -f "$MREF" ]; then  # if not ...
					echo -e "${COLOUR_RED}ERROR: -m requires a mature miRNA reference file in fasta or fa format, or file not found.${NO_COLOUR}\n" >&2
					exit 1
				fi
				MFLAG=0
				# extract reference file names
				MREF_NAME=`basename "$MREF"`
				MREF_PATH=`dirname "$MREF"`
				MREF_NAME_WO_EXT="${MREF_NAME%%.*}"
				MREF_PATH_NAME_WO_EXT="$MREF_PATH/$MREF_NAME_WO_EXT"
				;;
			s)  # flag to build bowtie index for stemloop/hairpin ref
				SREF=$OPTARG
				if ! [ -f "$SREF" ]; then  # if not ...
					echo -e "${COLOUR_RED}ERROR: -s requires a stemloop/hairpin file in fasta or fa format, or file not found.${NO_COLOUR}\n" >&2
					exit 1
				fi
				SFLAG=0
				# extract reference file names
				SREF_NAME=`basename "$SREF"`
				SREF_PATH=`dirname "$SREF"`
				SREF_NAME_WO_EXT="${SREF_NAME%%.*}"
				SREF_PATH_NAME_WO_EXT="$SREF_PATH/$SREF_NAME_WO_EXT"
				;;
			p)
				CORES=$OPTARG
				;;
			# \?) # only for invalid option warning. below ie. "*)" is a better option
				# echo "Invalid option: -$OPTARG" >&2
				# exit 1
				# ;;
			:)
				echo -e "${COLOUR_RED}ERROR: Option -$OPTARG requires an argument.${NO_COLOUR}\n" >&2
				exit 1
				;;
			*)  # if the input option not defined
				echo ""
				echo -e "${COLOUR_RED}ERROR: Invalid option: -$OPTARG${NO_COLOUR}\n" >&2
				echo -e $HELP
				echo -e "=========================================================================="
				echo -e "${COLOUR_YELLOW}$CITE${NO_COLOUR}\n"
				exit 1
				;;
			esac
		done
fi

if [[ $IFLAG -eq 1 || $AFLAG -eq 1 || $MFLAG -eq 1 ||$NFLAG -eq 1 || $SFLAG -eq 1 ]]; then
	echo -e "${COLOUR_RED}ERROR: -i, -a, -m, -n, -s flags are mandatory. Use -h or --help to see help info.${NO_COLOUR}\n" >&2
	exit 1
fi

# application variables from config file
# and their default settings if no config file
if [ -f ./mir_config ]; then  # variables read from the configeration file
  source mir_config
  ## below: -z tests if the variable has zero length. returns True if zero.
  if [[ -z $BASE_QUALITY || -z $MIN_LENGTH || -z $SEED || -z $SEED_MISMATCH ]]; then
    CONFIG_MSG="Config file detected. But one or more vairables missing.\n\tProceed with default settings."
		BASE_QUALITY=30  # cutadapt -q flag: minimum base quality threshold
		MIN_LENGTH=15  # cutadapt -m flag: minimum read length threshold
		SEED=20  # bowtie -l flag: seed length for -n
		SEED_MISMATCH=0  # bowtie -n flag: max mismatches in seed (can be 0-3)
  else
    CONFIG_MSG="Config file detected and loaded."
  fi
else
  CONFIG_MSG="Config file not detected. Proceed with default settings."
	BASE_QUALITY=30  # cutadapt -q flag: minimum base quality threshold
	MIN_LENGTH=15  # cutadapt -m flag: minimum read length threshold
	SEED=20  # bowtie -l flag: seed length for -n
	SEED_MISMATCH=0  # bowtie -n flag: max mismatches in seed (can be 0-3)
fi

# ------ functions ------
# functions to check if commands are installed
check_brew_app (){
	# Homebrew applications
	for i in bowtie fastqc samtools; do
		echo -en "	$i..."
		# sleep 1 # for testing only
		if hash $i 2>/dev/null; then
			echo -e "Pass!"
		else
			if [ $UNAMESTR == "Darwin" ]; then
				echo -e "Fail!"
				echo -e "\t-------------------------------------"
				echo -en "\t\tChecking Homebrew..."
				if hash brew 2>/dev/null; then
					echo -e "Pass!"
					echo -e "\t\tInstall $i through Homebrew..."
					brew install $i
					echo -e "\t-------------------------------------"
				else
					echo -e "Fail!\n"
					echo -e "${COLOUR_RED}ERROR: Homebrew not found for $i installation. Script terminated.${NO_COLOUR}\n" >&2
					echo -e "${COLOUR_RED}Please visit https://brew.sh to install Homebrew first.${NO_COLOUR}"
					exit 1
				fi
			else
				echo -e "Fail!\n"
				echo -e "${COLOUR_RED}Please install $i first.${NO_COLOUR}\n"
				exit 1
			fi
		fi
	done
}

check_other_app (){
	# conda/pip applications
	for i in cutadapt multiqc; do
		echo -en "	$i..."
		# sleep 1 # for testing only
		if hash $i 2>/dev/null; then
			echo -e "Pass!"
		else
			echo -e "Fail!"
			if hash conda 2>/dev/null; then
				echo -e "\t-------------------------------------"
				echo -e "\t\tPackage manager conda detected!"
				echo -e "\t\tInstall $i through conda..."
				conda install -c bioconda $i
				echo -e "\t-------------------------------------\n"
			elif hash pip 2>/dev/null; then
				echo -e "\t-------------------------------------"
				echo -e "\t\tPackage manager pip detected!"
				echo -e "\t\tInstall $i through pip..."
				pip install $i
				echo -e "\t-------------------------------------\n"
			else
				echo -e "Package manager not found for $i installation. Script terminated. \n" >&2
				echo -e "Please install either conda (https://conda.io) or pip first"
				exit 1
			fi
		fi
	done
}

# function to qc
qc (){
	# loop through fastq files for fastqc
	for i in "$1"/*.fastq.gz "$1"/*.fastq; do
		# below: || will run the following command if the test statement returns false
		[ -f "$i" ] || continue  # this line prevents from processing *.fastq or *.fastq.gz if file doesn't exist
		echo -en "\tRunning QC on `basename $i`..."
		fastqc "$i" --outdir="${OUT_DIR}"/QC_OUT/FASTQC_OUT 2>>"${OUT_DIR}"/LOG/qc.log 1>>"${OUT_DIR}"/LOG/qc.log # fastqc uses &2 for display
		echo -e "\n" >>"${OUT_DIR}"/LOG/qc.log
		echo -e "Done!"
	done
}

# function to trim adapters
# arguments: $1-input directory; $2-adapter file; $3-cores (i.e. $CORES)
adapter_trim (){
	# loop through fastq files for cutadapt
	for i in "$1"/*.fastq.gz "$1"/*.fastq; do
		[ -f "$i" ] || continue
		outname=`basename "$i"`  # NOTE: use ``, not "" or ''
		echo -en "\tTrimming adapters from $outname..."
		# below: cutadapt uses &2 for display system message, &1 for on-screen display messages
		cutadapt -b file:$2 -o "${OUT_DIR}"/TRIM_OUT/adapter_trimmed_$outname \
		"$i" \
		-m $MIN_LENGTH -q $BASE_QUALITY \
		--cores=$3 2>> "${OUT_DIR}"/LOG/adapter_trim_sys.log 1>> "${OUT_DIR}"/LOG/adapter_trim_results.log  # $3: third argument
		echo -e "Done!"

		if [ $NTCUT -eq 0 ]; then  # cut NTs if -c is present
			inname=adapter_trimmed_$outname
			echo -en "\tTrimming nucleotides from $inname..."
			cutadapt -u $NT5 -u $NT3 \
			-o "${OUT_DIR}"/TRIM_OUT/nt_trimmed_$outname \
			"${OUT_DIR}"/TRIM_OUT/$inname \
			-m $MIN_LENGTH --cores=$3 2>> "${OUT_DIR}"/LOG/adapter_trim_sys.log 1>> "${OUT_DIR}"/LOG/adapter_trim_results.log
			echo -e "Done!"
		fi
		echo -e "\n" >> "${OUT_DIR}"/LOG/adapter_trim_results.log
	done
}

# bowtie functions
idx_check (){
	# change U to T
	if [ $U2T -eq 0 ]; then
		echo -en "\tChanging U to T for mature miRNA reference..."
		sed -i '.bak' 's/U/T/g' $MREF
		echo -e "Done!"
		echo -en "\tChanging U to T for stemloop/hairpin reference..."
		sed -i '.bak' 's/U/T/g' $SREF
		echo -e "Done!"
		echo -en "\tChanging U to T for negative reference..."
		sed -i '.bak' 's/U/T/g' $NREF
		echo -e "Done!"
		echo -e "	-------------------------------------"
	fi

	# test if index exists, build if not
	echo -en "\tChecking bowtie index for negative reference..."
	if [ `ls "$NREF_PATH"/$NREF_NAME_WO_EXT.*.ebwt 2>/dev/null | wc -l` -gt 0 ]; then
		echo -e "Found!"
	else
		echo -e "Not found!"
		echo -en "\tBuilding bowtie index for negative reference..."
		bowtie-build $NREF $NREF_PATH_NAME_WO_EXT 1>>"${OUT_DIR}"/LOG/bowtie_idx.log
		echo -e "Done!"
	fi

	echo -en "\tChecking bowtie index for mature miRNA reference..."
	if [ `ls "$MREF_PATH"/$MREF_NAME_WO_EXT.*.ebwt 2>/dev/null | wc -l` -gt 0 ]; then
		echo -e "Found!"
	else
		echo -e "Not found!"
		echo -en "\tBuilding bowtie index for mature miRNA reference..."
		bowtie-build $MREF $MREF_PATH_NAME_WO_EXT 1>>"${OUT_DIR}"/LOG/bowtie_idx.log
		echo -e "Done!"
	fi

	echo -en "\tChecking bowtie index for stemloop/hairpin reference..."
	if [ `ls "$SREF_PATH"/$SREF_NAME_WO_EXT.*.ebwt 2>/dev/null | wc -l` -gt 0 ]; then
		echo -e "Found!"
	else
		echo -e "Not found!"
		echo -en "\tBuilding bowtie index for stemloop/hairpin reference..."
		bowtie-build $SREF $SREF_PATH_NAME_WO_EXT 1>>"${OUT_DIR}"/LOG/bowtie_idx.log
		echo -e "Done!"
	fi
}

neg_align (){
	if [ $NTCUT -eq 0 ]; then   # look for nt trimmed files first
		for i in "$1"/nt*.fastq.gz "$1"/nt*.fastq; do
			[ -f "$i" ] || continue

			neg_outname=`basename "$i"`
			neg_outname_wo_ext="${neg_outname##*_}"  # no extra crap with extension
			neg_outname_wo_ext="${neg_outname_wo_ext%%.*}"   # no extension

			echo -en "\tRemoving negative reads from $neg_outname..."
			echo -e "File: $neg_outname" >> "${OUT_DIR}"/LOG/bowtie_neg_results.log  # file name for log
			# below: bowtie command
			bowtie -q -p $CORES $NREF_PATH_NAME_WO_EXT \
			"$i" \
			--un "${OUT_DIR}"/BOWTIE_OUT/NEG_UNALIGN_OUT/neg_unaglin_$neg_outname_wo_ext.fastq \
			--al "${OUT_DIR}"/BOWTIE_OUT/NEG_ALIGN_OUT/neg_$neg_outname_wo_ext.fastq \
			1>/dev/null 2>>"${OUT_DIR}"/LOG/bowtie_neg_results.log  # discard flying message
			echo -e "\n" >> "${OUT_DIR}"/LOG/bowtie_neg_results.log
			echo -e "Done!"
		done
	else
		for i in "$1"/adapter*.fastq.gz "$1"/adapter*.fastq; do
			[ -f "$i" ] || continue

			neg_outname=`basename "$i"`
			neg_outname_wo_ext="${neg_outname##*_}"  # no extra crap with extension
			neg_outname_wo_ext="${neg_outname_wo_ext%%.*}" # no extension

			echo -en "\tRemoving negative reads from $neg_outname..."   # file name for log
			echo -e "File: $neg_outname" >> "${OUT_DIR}"/LOG/bowtie_neg_results.log
			# below: bowtie command
			bowtie -q -p $CORES $NREF_PATH_NAME_WO_EXT \
			"$i" \
			--un "${OUT_DIR}"/BOWTIE_OUT/NEG_UNALIGN_OUT/neg_unalign_$neg_outname_wo_ext.fastq \
			--al "${OUT_DIR}"/BOWTIE_OUT/NEG_ALIGN_OUT/neg_$neg_outname_wo_ext.fastq \
			1>/dev/null 2>>"${OUT_DIR}"/LOG/bowtie_neg_results.log  # discard flying message
			echo -e "\n" >> "${OUT_DIR}"/LOG/bowtie_neg_results.log
			echo -e "Done!"
		done
	fi
}

mature_align (){
	for i in "${OUT_DIR}"/BOWTIE_OUT/NEG_UNALIGN_OUT/*.fastq; do
		[ -f "$i" ] || continue

		mature_outname=`basename "$i"`
		mature_outname_wo_ext="${mature_outname##*_}"  # no extra crap with extension
		mature_outname_wo_ext="${mature_outname_wo_ext%%.*}"  # get file name without extension

		echo -en "\tAligning $mature_outname to mature miRNA reference..."
		echo -e "File: $mature_outname" >> "${OUT_DIR}"/LOG/bowtie_mature_results.log
		# below: the way to do "command 2>file | command2": { command | command 2; } 2>file
		{ bowtie --best --strata -a -q -S -p $CORES -l $SEED -n $SEED_MISMATCH \
			--un "${OUT_DIR}"/BOWTIE_OUT/MATURE_UNALIGN_OUT/mature_unalign_$mature_outname_wo_ext.fastq \
			$MREF_PATH_NAME_WO_EXT "$i" | samtools view -SbhF4 -o "${OUT_DIR}"/BOWTIE_OUT/MATURE_ALIGN_OUT/mature_$mature_outname_wo_ext.bam; } \
		1>/dev/null 2>>"${OUT_DIR}"/LOG/bowtie_mature_results.log
		echo -e "\n" >>"${OUT_DIR}"/LOG/bowtie_mature_results.log
		echo -e "Done!"
	done
}

hairpin_align (){
	for i in "${OUT_DIR}"/BOWTIE_OUT/MATURE_UNALIGN_OUT/*.fastq; do
		[ -f "$i" ] || continue

		hairpin_outname=`basename "$i"`
		hairpin_outname_wo_ext="${hairpin_outname##*_}"  # no extra crap with extension
		hairpin_outname_wo_ext="${hairpin_outname_wo_ext%%.*}"  # get file name without extension

		echo -en "\tAligning $hairpin_outname to stemloop/hairpin reference..."
		echo -e "File: $hairpin_outname" >> "${OUT_DIR}"/LOG/bowtie_hairpin_results.log
		# below: the way to do "command 2>file | command2": { command | command 2; } 2>file
		{ bowtie --best --strata -a -q -S -p $CORES -l $SEED -n $SEED_MISMATCH \
			--un "${OUT_DIR}"/BOWTIE_OUT/HAIRPIN_UNALIGN_OUT/haripin_unalign_$hairpin_outname_wo_ext.fastq \
			$SREF_PATH_NAME_WO_EXT "$i" | samtools view -SbhF4 -o "${OUT_DIR}"/BOWTIE_OUT/HAIRPIN_ALIGN_OUT/hairpin_$hairpin_outname_wo_ext.bam; } \
		1>/dev/null 2>>"${OUT_DIR}"/LOG/bowtie_hairpin_results.log
		echo -e "\n" >>"${OUT_DIR}"/LOG/bowtie_hairpin_results.log
		echo -e "Done!"
	done
}

read_count (){
	for i in "$1"/*.bam; do
		[ -f "$i" ] || continue

		count_outname=`basename "$i"`
		count_outname_wo_ext="${count_outname%%.*}"  # get file name without extension

		echo -en "\tCounting reads for $count_outname..."
		echo -e "File: $count_outname" >> "${OUT_DIR}"/LOG/count_sys.log
		# below: the way to do "command 2>file | command2": { command | command 2; } 2>file
		{ samtools sort -n "$i" \
			| samtools view \
			| awk '{print $3}' | sort | uniq -c \
			| sort -nr > "${OUT_DIR}"/READCOUNT_OUT/count_$count_outname_wo_ext.txt; } \
		2>>"${OUT_DIR}"/LOG/count_sys.log
		echo -e "\n" >>"${OUT_DIR}"/LOG/count_sys.log
		echo -e "Done!"
	done
}

# timing function
# from: https://www.shellscript.sh/tips/hms/
hms(){
  # Convert Seconds to Hours, Minutes, Seconds
  # Optional second argument of "long" makes it display
  # the longer format, otherwise short format.
  local SECONDS H M S MM H_TAG M_TAG S_TAG
  SECONDS=${1:-0}
  let S=${SECONDS}%60
  let MM=${SECONDS}/60 # Total number of minutes
  let M=${MM}%60
  let H=${MM}/60

  if [ "$2" == "long" ]; then
    # Display "1 hour, 2 minutes and 3 seconds" format
    # Using the x_TAG variables makes this easier to translate; simply appending
    # "s" to the word is not easy to translate into other languages.
    [ "$H" -eq "1" ] && H_TAG="hour" || H_TAG="hours"
    [ "$M" -eq "1" ] && M_TAG="minute" || M_TAG="minutes"
    [ "$S" -eq "1" ] && S_TAG="second" || S_TAG="seconds"
    [ "$H" -gt "0" ] && printf "%d %s " $H "${H_TAG},"
    [ "$SECONDS" -ge "60" ] && printf "%d %s " $M "${M_TAG} and"
    printf "%d %s\n" $S "${S_TAG}"
  else
    # Display "01h02m03s" format
    [ "$H" -gt "0" ] && printf "%02d%s" $H "h"
    [ "$M" -gt "0" ] && printf "%02d%s" $M "m"
    printf "%02d%s\n" $S "s"
  fi
}

# ------ script ------
# start time
start_t=`date +%s`

# display messages
echo -e "\nYou are running ${COLOUR_GREEN_L}mir_processing.sh${NO_COLOUR}"
echo -e "Version: $VERSION"
echo -e "Current OS: $PLATFORM"
echo -e "Output folder: $OUT_DIR"
echo -e "Today is: $CURRENT_DAY\n"
echo -e "${COLOUR_ORANGE}$CITE${NO_COLOUR}\n"

# loading config files
echo -e "\n"
echo -e "Configeration file (mir_config)"
echo -e "=========================================================================="
echo -e "\t$CONFIG_MSG"
echo -e "\t\tMIN_LENGTH=$MIN_LENGTH: cutadapt -m flag value"
echo -e "\t\tBASE_QUALITY=$BASE_QUALITY: cutadapt -q flag value"
echo -e "\t\tSEED=$SEED: bowtie -q flag value"
echo -e "\t\tSEED_MISMATCH=$SEED_MISMATCH: bowtie -n flag value"
echo -e "=========================================================================="

# command check
echo -e "\n"
echo -e "Checking all the necessary applications:"
echo -e "=========================================================================="
check_brew_app
check_other_app
echo -e "=========================================================================="

# create folders
echo -e "\n"
echo -en "Making output file folders..."
[ -d "${OUT_DIR}"/LOG ] || mkdir "${OUT_DIR}"/LOG  # [ -d /some/directory ]: test directory
[ -d "${OUT_DIR}"/QC_OUT ] || mkdir "${OUT_DIR}"/QC_OUT
[ -d "${OUT_DIR}"/QC_OUT/FASTQC_OUT ] || mkdir "${OUT_DIR}"/QC_OUT/FASTQC_OUT
[ -d "${OUT_DIR}"/TRIM_OUT ] || mkdir "${OUT_DIR}"/TRIM_OUT
[ -d "${OUT_DIR}"/BOWTIE_OUT ] || mkdir "${OUT_DIR}"/BOWTIE_OUT
[ -d "${OUT_DIR}"/BOWTIE_OUT/NEG_ALIGN_OUT ] || mkdir "${OUT_DIR}"/BOWTIE_OUT/NEG_ALIGN_OUT
[ -d "${OUT_DIR}"/BOWTIE_OUT/NEG_UNALIGN_OUT ] || mkdir "${OUT_DIR}"/BOWTIE_OUT/NEG_UNALIGN_OUT
[ -d "${OUT_DIR}"/BOWTIE_OUT/MATURE_ALIGN_OUT ] || mkdir "${OUT_DIR}"/BOWTIE_OUT/MATURE_ALIGN_OUT
[ -d "${OUT_DIR}"/BOWTIE_OUT/MATURE_UNALIGN_OUT ] || mkdir "${OUT_DIR}"/BOWTIE_OUT/MATURE_UNALIGN_OUT
[ -d "${OUT_DIR}"/BOWTIE_OUT/HAIRPIN_ALIGN_OUT ] || mkdir "${OUT_DIR}"/BOWTIE_OUT/HAIRPIN_ALIGN_OUT
[ -d "${OUT_DIR}"/BOWTIE_OUT/HAIRPIN_UNALIGN_OUT ] || mkdir "${OUT_DIR}"/BOWTIE_OUT/HAIRPIN_UNALIGN_OUT
[ -d "${OUT_DIR}"/READCOUNT_OUT ] || mkdir "${OUT_DIR}"/READCOUNT_OUT
echo -e "Done!"
echo -e "=========================================================================="
echo -e "\tFolders created and their usage:"
echo -e "\t\tQC_OUT: Quality check results"
echo -e "\t\tTRIM_OUT: Adapter-trimmed files"
echo -e "\t\tBOWTIE_OUT: Alignment results"
echo -e "\t\tREADCOUNT_OUT: Readcount results"
echo -e "\t\tLOG: Processing log files"
echo -e "=========================================================================="

# qc for raw reads
echo -e "\n"
echo -e "Quality checking raw files:"
echo -e "=========================================================================="
qc $RAW_FILES
echo -e "=========================================================================="

# cut and process adapters
echo -e "\n"
echo -e "Adapter trimming (speed depending on hardware configurations):"
echo -e "=========================================================================="
adapter_trim $RAW_FILES $ADAPTER $CORES
echo -e "=========================================================================="

# qc for the adapter trimmed files
echo -e "\n"
echo -e "Quality checking adapter-trimmed files:"
echo -e "=========================================================================="
qc "${OUT_DIR}"/TRIM_OUT
echo -e "\t-------------------------------------"
if hash multiqc 2>/dev/null; then
	echo -en "\tCompiling QC results into a single report..."
	echo -e "\n" >>"${OUT_DIR}"/LOG/qc.log
	multiqc "${OUT_DIR}"/QC_OUT/FASTQC_OUT/ -o "${OUT_DIR}"/QC_OUT -n multiqc_all.html 1>>"${OUT_DIR}"/LOG/qc.log 2>>"${OUT_DIR}"/LOG/qc.log
	echo -e "Done!"
else
	echo -e "\tMultiQC missing, QC results compiling step skipped."
fi
echo -e "=========================================================================="

# bowtie alignment
echo -e "\n"
echo -e "Bowtie alignment (speed depending on hardware configurations):"
echo -e "=========================================================================="
idx_check
echo -e "\t-------------------------------------"
neg_align "${OUT_DIR}"/TRIM_OUT
echo -e "\t-------------------------------------"
mature_align
echo -e "\t-------------------------------------"
hairpin_align
echo -e "=========================================================================="

# Counting reads
echo -e "\n"
echo -e "Read counting (speed depending on hardware configurations):"
echo -e "=========================================================================="
read_count "${OUT_DIR}"/BOWTIE_OUT/MATURE_ALIGN_OUT
read_count "${OUT_DIR}"/BOWTIE_OUT/HAIRPIN_ALIGN_OUT
echo -e "=========================================================================="

# end time and display
end_t=`date +%s`
tot=`hms $((end_t-start_t))`
echo -e "\n"
echo -e "Total run time: $tot"
echo -e "\n"
