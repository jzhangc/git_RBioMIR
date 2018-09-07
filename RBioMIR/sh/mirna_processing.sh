#!/bin/bash
# name: mirna.processing.sh
#
# ------ usage -------
# ./mirna_rocessing.sh -c -r [raw file directory] -a [adapter file name and full path] -p [CPU core number]
# -h: Help information
# -i <dir>: Input raw file path
# -a <file>: Adapter file with full path
# -C <int>: Cut nucleotides from 5'
# -c <int>: Cut nucleotides from 3'
# -n <file>: Negative reference with full path
# -m <file>: Mature miRNA reference with full path
# -s <file>: Stemloop/hairpin reference with full path
# -t: If to change U to T for the reference files
# -p <int>: Cpu core number
 
# ------ variables ------
# variables from command options
if [ $# -eq 0 ]; then # without any argument: display help file
	echo "TBC help messages from no options" >&2
    exit 0
else
	while getopts ":hti:a:C:c:n:m:p:" opt; do  # :hr: means :h and :r:
  		case $opt in
  			h)
  				echo -e "\n"
  				echo "TBC help messages from -h" >&2
  				echo -e "=========================================================================="
  				exit 1
  				;;
  			t)	# flag to change U to T for reference files
  				U2T=1
  				;;
    		i)
    			RAW_FILES=$OPTARG
    			;;
    		a)
    			ADAPTER=$OPTARG
    			if ! [ -f $ADAPTER ]; then  # if not ...
    			echo "-a requires an adapter file in fasta or fa format." >&2
    			exit 1
    			fi
    			;;
    		C)  # flag to decide if to cut the nucleotides
    			NTCUT=1  # set to 1(true) with the presence of -C
    			NT5=$OPTARG
    			;;
    		c)  # flag to decide if to cut the nucleotides
    			NTCUT=1  # set to 1(true) with the presence of -c
    			NT3=-$OPTARG
    			;;
    		n)  # flag to build bowtie index for negative ref
    			NEGREF=$OPTARG
    			if ! [ -f $NEGREF ]; then  # if not ...
    			echo "-n requires a negative reference file in fasta or fa format." >&2
    			exit 1
    			fi
    			    			
    			# extract reference file names
				NEGREF_PATH_NAME_WO_EXT="${NEGREF%%.*}"
    			;;
    		m)  # flag to build bowtie index for mature miRNA ref
    			MREF=$OPTARG
    			if ! [ -f $MREF ]; then  # if not ...
    			echo "-m requires a reference genome file in fasta or fa format." >&2
    			exit 1
    			fi
    			
    			# extract reference file names
    			MREF_PATH_NAME_WO_EXT="${MREF%%.*}"
    			;;
    		s)  # flag to build bowtie index for stemloop/hairpin ref
    			SREF=$OPTARG
    			if ! [ -f $SREF ]; then  # if not ...
    			echo "-m requires a reference genome file in fasta or fa format." >&2
    			exit 1
    			fi
    			
    			# extract reference file names
    			SREF_PATH_NAME_WO_EXT="${SREF%%.*}"
    			;;
    		p)
    			CORES=$OPTARG
    			;;
    		# \?) # only for invalid option warning. below ie. "*)" is a better option 
    	 		#echo "Invalid option: -$OPTARG" >&2
    	 		#exit 1
      	 		#;;
    		:)
    			echo "Option -$OPTARG requires an argument." >&2
      			exit 1
      			;;
    		*)  # if the input option not defined
        		echo "Invalid option: -$OPTARG" >&2
  				echo -e "\n"
  				echo "TBC help messages from -*" >&2
  				echo -e "=========================================================================="
  				exit 1
        		;;
  		esac
	done
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
			echo -e "Fail!"
			echo -e "	-------------------------------------"
			echo -en "		Check Homebrew..."
			if hash brew 2>/dev/null; then
				echo -e "Pass!"
				echo -e "		Install $i through Homebrew..."
				echo -e "\n"
				brew install $i	
				echo -e "	-------------------------------------"
			else
				echo -e "Fail!\n"
				echo -e "Package manager not found for $i installation. Script terminated."
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
				echo -e "	-------------------------------------"
				echo -e "		Package manager conda detected!"
				echo -e "		Install $i through conda..."
				conda install -c bioconda $i
				echo -e "	-------------------------------------\n"	
			elif hash pip 2>/dev/null; then
				echo -e "	-------------------------------------"
				echo -e "		Package manager pip detected!"
				echo -e "		Install $i through pip..."
				pip install $i
				echo -e "	-------------------------------------\n"
			else
				echo -e "Package manager not found for $i installation. Script terminated. \n"
				exit 1	
			fi
		fi
	done
}


# function to qc
qc (){	
	# loop through fastq files for fastqc
	for i in $1/*.fastq.gz $1/*.fastq; do
		# below: || will run the following command if the test statement returns false
		[ -f "$i" ] || break  # this line prevents from processing *.fastq or *.fastq.gz if file doesn't exist
		echo -en "	Runing QC on `basename $i`..."
		fastqc "$i" --outdir=./QC_OUT/ 2>>./LOG/qc.log 1>>./LOG/qc.log # fastqc uses &2 for display
		echo -e "\n" >>./LOG/qc.log
		echo -e "Done!"
	done
}

# function to trim adapters
adapter_trim (){
	# loop through fastq files for cutadapt
	for i in $1/*.fastq.gz $1/*.fastq; do
		[ -f "$i" ] || break
		outname=`basename "$i"`  # NOTE: use ``, not "" or ''
		echo -en "	Trimming adapters from $outname..."
		# below: cutadapt uses &2 for display system message, &1 for on-screen display messages
		cutadapt -b file:$2 -o ./TRIM_OUT/adapter_trimmed_$outname \
		"$i" \
		-m 15 -q 30 \
		--cores=$3 2>> ./LOG/adapter_trim_sys.log 1>> ./LOG/adapter_trim_results.log
		echo -e "Done!"
		
		if [ $NTCUT -eq 1 ]; then  # cut NTs if -c is present
			inname=adapter_trimmed_$outname
			echo -en "	Trimming nucleotides from $inname..."
			cutadapt -u $NT5 -u $NT3 \
			-o ./TRIM_OUT/nt_trimmed_$outname \
			./TRIM_OUT/$inname \
			-m 15 --cores=$3 2>> ./LOG/adapter_trim_sys.log 1>> ./LOG/adapter_trim_results.log
			echo -e "Done!"
		fi		
	done	
}

# bowtie functions
idx_check (){
	# change U to T
	if [ $U2T -eq 1 ]; then
		echo -en "	Change Us to Ts for mature miRNA reference..."
		sed -i '.bak' 's/U/T/g' $MREF
		echo -e "Done!"
		echo -en "	Change Us to Ts for stemloop/hairpin reference..."
		sed -i '.bak' 's/U/T/g' $SREF
		echo -e "Done!"
		echo -en "	Change Us to Ts for negative reference..."
		sed -i '.bak' 's/U/T/g' $NEGREF
		echo -e "Done!"
		echo -e "	-------------------------------------"
	fi
	
	# test if index exists, build if not
	neg_path=`dirname $NEGREF`  # extract negative ref path
	neg_name=`basename $NEGREF_PATH_NAME_WO_EXT`	
	echo -en "	Checking bowtie index for negative reference..."	
	if [ `ls $neg_path/$neg_name.*.ebwt 2>/dev/null | wc -l` -gt 0 ]; then
		echo -e "Found!"
	else
		echo -e "Not found!"
		echo -en "	Building bowtie index for negative reference..."
		bowtie-build $NEGREF $NEGREF_PATH_NAME_WO_EXT 1>>bowtie_idx.log
		echo -e "Done!"
	fi
	
	mature_path=`dirname $MREF`  # extract mature miRNA ref path
	mature_name=`basename $MREF_PATH_NAME_WO_EXT`
	echo -en "	Checking bowtie index for mature miRNA reference..."	
	if [ `ls $mature_path/$mature_name.*.ebwt 2>/dev/null | wc -l` -gt 0 ]; then
		echo -e "Found!"
	else
		echo -e "Not found!"
		echo -en "	Building bowtie index for mature miRNA reference..."
		bowtie-build $MREF $MREF_PATH_NAME_WO_EXT 1>>./LOG/bowtie_idx.log
		echo -e "Done!"
	fi
	
	haripin_path=`dirname $SREF`  # extract mature miRNA ref path
	hairpin_name=`basename $SREF_PATH_NAME_WO_EXT`
	echo -en "	Checking bowtie index for stemloop/hairpin reference..."	
	if [ `ls $haripin_path/$hairpin_name.*.ebwt 2>/dev/null | wc -l` -gt 0 ]; then
		echo -e "Found!"
	else
		echo -e "Not found!"
		echo -en "	Building bowtie index for stemloop/hairpin reference..."
		bowtie-build $SREF $SREF_PATH_NAME_WO_EXT 1>>./LOG/bowtie_idx.log
		echo -e "Done!"
	fi
}

neg_align (){	
	if [ $NTCUT -eq 1 ]; then   # look for nt trimmed files first
		for i in $1/nt*.fastq.gz $1/nt*.fastq; do	
			[ -f "$i" ] || break
			neg_outname=`basename "$i"`
			neg_outname_wo_ext="${neg_outname##*_}"  # no extra crap with extension
			neg_outname_wo_ext="${neg_outname_wo_ext%%.*}"   # no extension
			echo -en "	Removing negative reads from $neg_outname..."
			echo -e "File: $neg_outname" >> ./LOG/bowtie_neg_results.log  # file name for log
			# below: bowtie command
			bowtie -q -p $CORES $NEGREF_PATH_NAME_WO_EXT \
			"$i" \
			--un ./BOWTIE_OUT/NEG_UNALIGN_OUT/neg_unaglin_$neg_outname_wo_ext.fastq \
			--al ./BOWTIE_OUT/NEG_ALIGN_OUT/neg_$neg_outname_wo_ext.fastq \
			1>/dev/null 2>>./LOG/bowtie_neg_results.log  # discard flying message
			echo -e "\n" >> ./LOG/bowtie_neg_results.log
			echo -e "Done!"			
		done
	else 
		for i in $1/adapter*.fastq.gz $1/adapter*.fastq; do
			[ -f "$i" ] || break
			neg_outname=`basename "$i"`
			neg_outname_wo_ext="${neg_outname##*_}"  # no extra crap with extension
			neg_outname_wo_ext="${neg_outname_wo_ext%%.*}" # no extension
			echo -en "	Removing negative reads from $neg_outname..."   # file name for log
			echo -e "File: $neg_outname" >> ./LOG/bowtie_neg_results.log
			# below: bowtie command
			bowtie -q -p $CORES $NEGREF_PATH_NAME_WO_EXT \
			"$i" \
			--un ./BOWTIE_OUT/NEG_UNALIGN_OUT/neg_unalign_$neg_outname_wo_ext.fastq \
			--al ./BOWTIE_OUT/NEG_ALIGN_OUT/neg_$neg_outname_wo_ext.fastq \
			1>/dev/null 2>>./LOG/bowtie_neg_results.log  # discard flying message
			echo -e "\n" >> ./LOG/bowtie_neg_results.log
			echo -e "Done!"
		done			
	fi	
}

mature_align (){
	for i in ./BOWTIE_OUT/NEG_UNALIGN_OUT/*.fastq; do	
		[ -f "$i" ] || break
		mature_outname=`basename "$i"`
		mature_outname_wo_ext="${mature_outname##*_}"  # no extra crap with extension
		mature_outname_wo_ext="${mature_outname_wo_ext%%.*}"  # get file name without extension		
		echo -en "	Aligning $mature_outname to mature miRNA reference..."
		echo -e "File: $mature_outname" >> ./LOG/bowtie_mature_results.log
		# below: the way to do "command 2>file | command2": { command | command 2; } 2>file
		{ bowtie --best --strata -a -q -S -l 20 -n 0 -v 2 \
			--un ./BOWTIE_OUT/MATURE_UNALIGN_OUT/mature_unalign_$mature_outname_wo_ext.fastq -p $CORES \
			$MREF_PATH_NAME_WO_EXT "$i" | samtools view -SbhF4 -o ./BOWTIE_OUT/MATURE_ALIGN_OUT/mature_$mature_outname_wo_ext.bam; } \
		1>/dev/null 2>>./LOG/bowtie_mature_results.log
		echo -e "\n" >>./LOG/bowtie_mature_results.log
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
echo -e "\n"
echo -e "Script written by Jing Zhang PhD"
echo -e "Contact: jzha9@uwo.ca, jzhangcad@gmail.com"
echo -e "To cite in your research: TBA"

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
[ -d ./LOG ] || mkdir ./LOG  # [ -d /some/directory ]: test directory
[ -d ./QC_OUT ] || mkdir ./QC_OUT
[ -d ./QC_OUT/MULTIQC_OUT ] || mkdir ./QC_OUT/MULTIQC_OUT
[ -d ./TRIM_OUT ] || mkdir ./TRIM_OUT
[ -d ./BOWTIE_OUT ] || mkdir ./BOWTIE_OUT
[ -d ./BOWTIE_OUT/NEG_ALIGN_OUT ] || mkdir ./BOWTIE_OUT/NEG_ALIGN_OUT
[ -d ./BOWTIE_OUT/NEG_UNALIGN_OUT ] || mkdir ./BOWTIE_OUT/NEG_UNALIGN_OUT
[ -d ./BOWTIE_OUT/MATURE_ALIGN_OUT ] || mkdir ./BOWTIE_OUT/MATURE_ALIGN_OUT
[ -d ./BOWTIE_OUT/MATURE_UNALIGN_OUT ] || mkdir ./BOWTIE_OUT/MATURE_UNALIGN_OUT
[ -d ./BOWTIE_OUT/HAIRPIN_ALIGN_OUT ] || mkdir ./BOWTIE_OUT/HAIRPIN_ALIGN_OUT
[ -d ./BOWTIE_OUT/HAIRPIN_UNALIGN_OUT ] || mkdir ./BOWTIE_OUT/HAIRPIN_UNALIGN_OUT
[ -d ./READCOUNTS_OUT ] || mkdir ./READCOUNTS_OUT
[ -d ./READCOUNTS_OUT/SORTED_BAM ] || mkdir ./READCOUNTS_OUT/SORTED_BAM
echo -e "Done!"
echo -e "=========================================================================="
echo -e "	Folders created and their usage:"
echo -e "		QC_OUT: Quality check results"
echo -e "		TRIM_OUT: Adapter-trimmed files"
echo -e "		BOWTIE_OUT: Alignment results"
echo -e "		READCOUNTS_OUT: Readcount results"
echo -e "		LOG: Processing log files"
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
qc ./TRIM_OUT
echo -e "	-------------------------------------"
echo -en "	Compiling QC results into a single report..."
echo -e "\n" >>./LOG/qc.log
multiqc ./QC_OUT/ -o ./QC_OUT/MULTIQC_OUT -n multiqc_all.html 1>>./LOG/qc.log 2>>./LOG/qc.log
echo -e "Done!"
echo -e "=========================================================================="

# bowtie alignment
echo -e "\n"
echo -e "Bowtie alignment (speed depending on hardware configurations):"
echo -e "=========================================================================="
idx_check
echo -e "	-------------------------------------"
neg_align ./TRIM_OUT
echo -e "	-------------------------------------"
mature_align
echo -e "=========================================================================="

# end time and display
end_t=`date +%s`
tot=`hms $((end_t-start_t))`
echo -e "\n"
echo -e "Total run time: $tot"
echo -e "\n"