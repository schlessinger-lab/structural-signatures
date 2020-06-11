#!/bin/bash

####################################
### Structural Signatures Pipeline #
### Written by: Rayees Rahman    ###
### For Mount Sinai DToXs        ###
### Version: 2 				     ###
### Dev: 6.24.19 				 ###
### Usage: See Below             ###
####################################

####################################
### Environment Variables:       ###
####################################
    
ref_cnt=""

if [ -z "$homed" ]; 
then 
	printf "\nThe variable: homed needs to be defined, please define the variable as such:\n\n\033[34mecho 'export homed="/path/to/structural-signatures"' >> ~/.bashrc \n\n\033[0m"; 
	exit  
else 
	printf "\nInstallation directory:\n\t\033[34m$homed\n\033[0m"
fi

db=$homed/database/structure_database.db 
if [[ -f "$db" ]]
then
	printf "Structure database located at:\n\t\033[33m$db\n\033[0m"
else 
	printf "\n\033[31mStructure database not found, please untar the structure database in:\n\t$homed/database/\n\033[0m\n"
	exit 
fi  

####################################
### Options:			         ###
####################################

u="Usage:\n\t -i input gene list \n\t -t domain enrichment (domain) scop enrichment (fold) or both (both)\n\t -n name type: gene name (gn) or uniprot id (uid)\n\t -o output file name \n)";  
while getopts ":i:o:t:n:d:b:f:g:T:p:lh:e:P:y:c:r:zh:" opt; do
    case $opt in
		i) 	input=$OPTARG ;;  #input
		t) 	type=$OPTARG ;; #domain or fold 
        n)  name=$OPTARG ;; #gene name (gn) or uid (uid) 
		o)  output=$OPTARG ;; #output file name 
		r) 	ref_db=$OPTARG ;; 
        \?)
        printf "$" 
        exit 1 ;;
    esac
done

####################################
### Checks   			         ###
####################################

if [ -z $input ] || [ -z $type ] || [ -z $name ] || [ -z $output ]
then
	printf "$u" ;
	exit  1
fi

if [[ -f "$input" ]]
then
	printf "" ; 
else 
	printf "\033[31mNo such file:\033[0m $input\n" ;
	exit 
fi 

if [ -z $ref_db ] 
then 
	printf "Using default human proteome background structure counts at:\n\033[32m\t$homed/bin/files/backgrounds/human_proteome/\n\033[0m\n"
	ref_cnt="$homed/bin/files/backgrounds/human_proteome/human_proteome"
else 
	printf "Using custom reference background structure counts at\n\033[32m\t$homed/bin/files/backgrounds/$ref_db/\n\033[0m\n"
	printf "Checking custom reference counts for correct formatting"	
	if [[ -f "$homed/bin/files/backgrounds/$ref_db/$ref_db.background.ipr.domain.cnt" ]] && [[ -f "$homed/bin/files/backgrounds/$ref_db/$ref_db.background.scop.family.cnt" ]]  && [[ -f "$homed/bin/files/backgrounds/$ref_db/$ref_db.background.scop.superfam.cnt" ]] && [[ -f "$homed/bin/files/backgrounds/$ref_db/$ref_db.background.scop.fold.cnt" ]] 
	then
		printf "."
		for f in $homed/bin/files/backgrounds/$ref_db/$ref_db.* 
		do 
			if [[ `(head -1 $f | tr "," "\n" | wc -l)` != 2 ]] 
			then
				printf "\n\033[31mIncorrect number of columns for:\033[0m $f\nEach background file should have two columns deliminated by a ',', the first column is the structure id, the second are the counts\n"
				exit 
			else 
				printf "."
				if [[ `(head -1 $f  | grep "^IPR\|^[abcdefg]\.[0-9]" | wc -l)` != 1 ]] 
				then 
					printf "\n\033[31mIncorrect columns for:\033[0m $f\nEach background file should have two columns deliminated by a ',', the first column is the structure id (starting with IPR for domains or a.2, etc. for SCOPe ), the second are the counts\n"
					exit 
				else 
					printf "."
				fi
			fi
		done 
		printf "\n\033[32mOK 👍\033[0m\n"
		ref_cnt="$homed/bin/files/backgrounds/$ref_db/$ref_db"
	else 
		printf "\n\033[31mDid not find background files, make sure the file names follow this format:\033[0m\n\t$ref_db.background.ipr.domain.cnt\n\t$ref_db.background.scop.family.cnt\n\t$ref_db.background.scop.superfam.cnt\n\t$ref_db.background.scop.fold.cnt"
		printf "\nAnd are located in this directory: \033[31m$homed/bin/files/backgrounds/$ref_db/\n\033[0m"
		exit 
	fi	
fi

####################################
### Main     			         ###
####################################

$homed/bin/scripts/get_struct_from_db.pl $homed/database/structure_database.db $type $name $homed/bin/files/ParentChildTreeFile.txt $input $output $homed/bin/files/available_uniprot_ids.csv yes #2> /dev/null
if [[ -f  "./$output.found.genes" ]]
then
	printf "Generating enrichments\n"
else 
	printf "\n\033[31mSomething went wrong, please check if the correct perl modules are available (DBI and parallel::forkmanager).\n\033[0m"
	exit	
fi 

numgenes=$( cat ./$output.found.genes ) 
if [ $type == "domain" ] 
then
    Rscript $homed/bin/scripts/compute_representation.R $output.domain.cnt $ref_cnt.background.ipr.domain.cnt $numgenes $output domain  #2> /dev/null
	# remove files ----------------------------------------
	rm $output.domain.cnt
	rm $output.found.genes
	rm $output.scop.structure.info.csv
	# Check if all files are generated ------------------------------
	if [[ -f $output.uid.converted ]] && [[ -f $output.ipr.info.csv ]] && [[ -f $output-domain-enrichments.csv ]]
	then 
		printf "\033[32mSuccess! 👏\n\033[0m"
	else
		printf "\033[31mSomething went wrong, not all of the output was generated\033[0m\n"
		printf "There should be the following files:\n\t$output.uid.converted\n$output.ipr.info.csv\n\t$output-domain-enrichments.csv"
		printf "Please check your inputs and try running again."
	fi 
elif [ $type == "fold" ]
then

	cat ./$output.scop.fold.cnt | cut -f1,2 -d"," > $output.tmp.fold 
	Rscript $homed/bin/scripts/compute_representation.R $output.tmp.fold $ref_cnt.background.scop.fold.cnt $numgenes $output fold #2> /dev/null
	cat ./$output.scop.superfam.cnt | cut -f1,2 -d"," > $output.tmp.superfam
	Rscript $homed/bin/scripts/compute_representation.R $output.tmp.superfam $ref_cnt.background.scop.superfam.cnt $numgenes $output superfam  #2> /dev/null
	cat ./$output.scop.family.cnt | cut -f1,2 -d"," > $output.tmp.fam
	Rscript $homed/bin/scripts/compute_representation.R $output.tmp.fam $ref_cnt.background.scop.family.cnt $numgenes $output family  #2> /dev/null
	# remove files ----------------------------------------
	rm ./$output.tmp.fold
	rm ./$output.scop.fold.cnt 
	rm ./$output.tmp.superfam 
	rm ./$output.scop.superfam.cnt
	rm ./$output.tmp.fam
	rm ./$output.scop.family.cnt
	rm $output.found.genes
	rm ./$output.scop.class.cnt
	rm $output.ipr.info.csv
	# Check if all files are generated ------------------------------
	if [[ -f $output.uid.converted ]] && [[ -f $output.scop.structure.info.csv ]] && [[ -f $output-fold-enrichments.csv ]] && [[ -f $output-family-enrichments.csv ]] && [[ -f $output-superfam-enrichments.csv ]]
	then 
		printf "\033[32mSuccess! 👏\n\033[0m"
	else
		printf "\033[31mSomething went wrong, not all of the output was generated\033[0m\n"
		printf "There should be the following files:\n\t$output.uid.converted\n\t$output.scop.structure.info.csv \n\t$output-fold-enrichments.csv\n\t$output-family-enrichments.csv\n\t$output-superfam-enrichments.csv"
		printf "\nPlease check your inputs and try running again\n"
	fi 
else 
	Rscript $homed/bin/scripts/compute_representation.R $output.domain.cnt $ref_cnt.background.ipr.domain.cnt $numgenes $output domain  #2> /dev/null
	cat ./$output.scop.fold.cnt | cut -f1,2 -d"," > $output.tmp.fold 
	Rscript $homed/bin/scripts/compute_representation.R $output.tmp.fold $ref_cnt.background.scop.fold.cnt $numgenes $output fold  #2> /dev/null
	cat ./$output.scop.superfam.cnt | cut -f1,2 -d"," > $output.tmp.superfam
	Rscript $homed/bin/scripts/compute_representation.R $output.tmp.superfam $ref_cnt.background.scop.superfam.cnt $numgenes $output superfam  #2> /dev/null
	cat ./$output.scop.family.cnt | cut -f1,2 -d"," > $output.tmp.fam
	Rscript $homed/bin/scripts/compute_representation.R $output.tmp.fam $ref_cnt.background.scop.family.cnt $numgenes $output family #2> /dev/null
	# remove files ----------------------------------------
	rm ./$output.domain.cnt
	rm ./$output.tmp.fold
	rm ./$output.scop.fold.cnt
	rm ./$output.tmp.superfam 
	rm ./$output.scop.superfam.cnt 
	rm ./$output.tmp.fam
	rm ./$output.scop.family.cnt
	rm ./$output.found.genes
	rm ./$output.scop.class.cnt
	# Check if all files are generated ------------------------------
	if [[ -f $output.uid.converted ]] && [[ -f $output.scop.structure.info.csv ]] && [[ -f $output-fold-enrichments.csv ]] && [[ -f $output-family-enrichments.csv ]] && [[ -f $output-superfam-enrichments.csv ]]  && [[ -f $output.ipr.info.csv ]] && [[ -f $output-domain-enrichments.csv ]]
	then 
		printf "\033[32mSuccess! 👏\n\033[0m"
	else
		printf "\033[31mSomething went wrong, not all of the output was generated\033[0m\n"
		printf "There should be the following files:\n\t$output.ipr.info.csv\n\t$output-domain-enrichments.csv\n\t$output.uid.converted\n\t$output.scop.structure.info.csv \n\t$output-fold-enrichments.csv\n\t$output-family-enrichments.csv\n\t$output-superfam-enrichments.csv"
		printf "\nPlease check your inputs and try running again\n"
	fi 
fi
exit ; 

