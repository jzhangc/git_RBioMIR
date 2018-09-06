#!/bin/bash
# name: mirna.processing.sh
#
# ------ usage -------
# ./mirna_rocessing.sh -c -r [raw file directory] -a [adapter file name and full path] -p [CPU core number]
# -h: show help information
# -C: cut nucleotides from 5'
# -c: cut nucleotides from 3'
# -r: raw file path
# -a: adapter file with full path
# -p: cpu core number
 
# ------ variables ------
# script internal variables
NTCUT=0  # initiate flag if to cut the nucleotides. initial value 0(false)
U2T=0  # initiate flag if to change U to T for reference files. initial value 0(false)
NT5=0  # initiate variable for cutting nucleotides from 5'
NT3=0  # initiate variable for cutting nucleotides from 3'
NEGIDX=0  # initiate flag to build bowtie negative reference index
POSIDX=0 # initiate flag to build bowtie positive reference index

# variables from command options
if [ $# -eq 0 ]; then # without any argument: display help file
	echo "TBC help messages from no options" >&2
    exit 0
else
	while getopts ":hti:a:C:c:n:r:p:" opt; do  # :hr: means :h and :r:
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
    			NEGIDX=1
    			NEGREF=$OPTARG
    			if ! [ -f $NEGREF ]; then  # if not ...
    			echo "-n requires a negative reference file in fasta or fa format." >&2
    			exit 1
    			fi
    			;;
    		r)  # flag to build bowtie index for negative ref
    			POSIDX=1
    			POSREF=$OPTARG
    			if ! [ -f $POSREF ]; then  # if not ...
    			echo "-r requires a reference genome file in fasta or fa format." >&2
    			exit 1
    			fi
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
# function to check if commands are installed
command_check (){	
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
	
	# conda/pip applications
	for j in cutadapt multiqc; do
		echo -en "	$j..."
		# sleep 1 # for testing only
		if hash $j 2>/dev/null; then
			echo -e "Pass!"
		else
			echo -e "Fail!"
			if hash conda 2>/dev/null; then
				echo -e "	-------------------------------------"
				echo -e "		Package manager conda detected!"
				echo -e "		Install $j through conda..."
				conda install -c bioconda $j
				echo -e "	-------------------------------------\n"	
			elif hash pip 2>/dev/null; then
				echo -e "	-------------------------------------"
				echo -e "		Package manager pip detected!"
				echo -e "		Install $j through pip..."
				pip install $j
				echo -e "	-------------------------------------\n"
			else
				echo -e "Package manager not found for multiqc installation. Script terminated. \n"
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
		echo -e "Done!"
	done
}

# function to trim adapters
adapter_trim (){
	# loop through fastq files for cutadapt
	for i in $1/*.fastq.gz $1/*.fastq; do
		[ -f "$i" ] || break
		outname=`basename $i`  # NOTE: use ``, not "" or ''
		echo -en "	Trimming adapters from `basename $i`..."
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

# bowtie function
align_foo (){
	# change U to T
	if [ $U2T -eq 1 ]; then
		echo -en "	Change Us to Ts for positive reference..."
		sed -i '.bak' 's/U/T/g' $POSREF
		echo -e "Done!"
		echo -en "	Change Us to Ts for positive reference..."
		sed -i '.bak' 's/U/T/g' $NEGREF
		echo -e "Done!"
		echo -e "	-------------------------------------"
	fi
	
	# extract reference file names
	neg_name="${NEGREF%.*}"
	pos_name="${POSREF%.*}"
	
	# build bowtie index if chosen
	if [ $NEGIDX -eq 1 ]; then
		echo -en "	Building bowtie index for negative reference..."
		bowtie-build $NEGREF $neg_name 1>>bowtie_idx.log
		echo -e "Done!"
	fi
	
	if [ $POSIDX -eq 1 ]; then
		echo -en "	Building bowtie index for positive reference..."
		bowtie-build $POSREF $pos_name 1>>./LOG/bowtie_idx.log
		echo -e "Done!"
	fi
	
	# bowtie negative align
	echo -e "	-------------------------------------"			
	if [ $NTCUT -eq 1 ]; then   # look for nt trimmed files first
		for i in $1/nt*.fastq.gz $1/nt*.fastq; do	
			[ -f "$i" ] || break
			neg_outname=`basename $i`
			echo -en "	Removing negative reads from $neg_outname..."
			echo -e "File: $neg_outname" >> ./LOG/bowtie_neg_results.log  # file name for log
			# below: bowtie command
			bowtie -p $CORES -q $neg_name \
			"$i" \
			--un ./BOWTIE_OUT/NEG_UNALIGN_OUT/neg_flt_$neg_outname \
			--al ./BOWTIE_OUT/NEG_ALIGN_OUT/neg_$neg_outname \
			1>/dev/null 2>>./LOG/bowtie_neg_results.log  # discard flying message
			echo -e "\n" >> ./LOG/bowtie_neg_results.log
			echo -e "Done!"			
		done
	else 
		for i in $1/adapter*.fastq.gz $1/adapter*.fastq; do
			[ -f "$i" ] || break
			neg_outname=`basename $i`
			echo -en "	Removing negative reads from $neg_outname..."   # file name for log
			echo -e "File: $neg_outname" >> ./LOG/bowtie_neg_results.log
			# below: bowtie command
			bowtie -p $CORES -q $neg_name \
			"$i" \
			--un ./BOWTIE_OUT/NEG_UNALIGN_OUT/neg_flt_$neg_outname \
			--al ./BOWTIE_OUT/NEG_ALIGN_OUT/neg_$neg_outname \
			1>/dev/null 2>>./LOG/bowtie_neg_results.log  # discard flying message
			echo -e "\n" >> ./LOG/bowtie_neg_results.log
			echo -e "Done!"
		done			
	fi
	
	# bowtie positive align
	echo -e "	-------------------------------------"
	for j in ./BOWTIE_OUT/NEG_UNALIGN_OUT/*.fastq.gz ./BOWTIE_OUT/NEG_UNALIGN_OUT/*.fastq; do	
		[ -f "$j" ] || break
		pos_outname=`basename $j`
		pos_outname=${pos_outname%.*.*}  # get file name without extension
		echo -en "	Aligning $pos_outname to positive reference..."
		echo -e "File: $pos_outname" >> ./LOG/bowtie_pos_results.log
		bowtie -p $CORES -q -S -l 20 -n 0 -v 2 -a --best --strata $pos_name "$j" --al ./BOWTIE_OUT/POS_ALIGN_OUT/pos_$pos_outname.sam --un ./BOWTIE_OUT/POS_UNALIGN_OUT/pos_unaligned_$pos_outname.sam 1>/dev/null 2>>./LOG/bowtie_pos_results.log # discard flying message
		echo -e "\n" >> ./LOG/bowtie_pos_results.log
		echo -e "Done!"		
	done	
}


# ------ script ------
echo -e "\n"
echo -e "Script written by Jing Zhang PhD"
echo -e "Contact: jzha9@uwo.ca, jzhangcad@gmail.com"
echo -e "To cite in your research: TBA"
# command check
echo -e "\n"
echo -e "Checking all the necessary applications:"
echo -e "=========================================================================="
command_check
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
[ -d ./BOWTIE_OUT/POS_ALIGN_OUT ] || mkdir ./BOWTIE_OUT/POS_ALIGN_OUT
[ -d ./BOWTIE_OUT/POS_UNALIGN_OUT ] || mkdir ./BOWTIE_OUT/POS_UNALIGN_OUT
[ -d ./BOWTIE_OUT/HAIRPIN_UNALIGN_OUT ] || mkdir ./BOWTIE_OUT/HAIRPIN_UNALIGN_OUT
[ -d ./READCOUNTS_OUT ] || mkdir ./READCOUNTS_OUT
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
multiqc ./QC_OUT/ -o ./QC_OUT/MULTIQC_OUT -n multiqc_all.html 1>>./LOG/qc.log 2>>./LOG/qc.log
echo -e "Done!"
echo -e "=========================================================================="

# bowtie alignment
echo -e "\n"
echo -e "Bowtie alignment (speed depending on hardware configurations):"
echo -e "=========================================================================="
align_foo ./TRIM_OUT
echo -e "=========================================================================="