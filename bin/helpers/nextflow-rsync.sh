#!/usr/bin/env bash

# Author: 	Rachel Honeyghan-Williams, Lipovsek Lab, UCL Ear Institute
# Description:	Copy all files from the qc-pipeline safely to the RDSS
# Notes: 	1) 	Must be run in nextflow folder containing "work" and "results" directories
#		2)	Both gzip and rsync use checksums to confirm file integrity

##########################################
# using the script                       #
##########################################

# syntax
# bash nextflow-rsync.sh <RESULTS_DIR> <DESTINATION>

# example 
# (where <RDSS_DIRECTORY> is the full path beginning with /rdss/)
# nextflow-rsync.sh results <USERNAME>@rdp-ssh.arc.ucl.ac.uk:<RDSS_DIRECTORY>

# arguments
# DESTINATION	- 	destination directory for file transfer	(e.g. sjjgrww@rdp-ssh.arc.ucl.ac.uk:/rdss/rd01/ritd-ag-project-rd01le-mlipo86/atlas/qc/mouse/)

##########################################
# file compression explanation           #
##########################################

# find command
# find the absolute paths for each unique file in the work directories and the results directory | pass to tar
# exec 		- 	run a command on each result:
#		-	realpath	-	gives absolute real path, following symlinks
#		-	'{}'		-	current path
#		-	';'		-	end of exec expression
# sort		- 	sort by name
# uniq		- 	only keep unique file names

# command subsititution (brackets)
# supplies a list of all work directories used in the nextflow pipeline (nextflow log last) and the results directory

# tar command
# c		-	create a new tar archive		(a.k.a. tarball)
# z		-	compress files using gzip
# v		-	verbose output
# f		-	name of new tar archive			(qc-pipeline-<DATE>.tar.gz)
# T		-	take list of files from stdin		(i.e. results of the find command)

##########################################
# file transfer explanation              #
##########################################

# rsync command (use command 'rsync --help' for more details)
# a 		-	archive					(-rlptgoD)
# v 		-	verbose output
# h 		-	human readable file sizes
# P 		-	show progress bar
# L 		-	copy files content instead of symlink	(symlinks to absolute paths, e.g. '/home/rachel/.../work/' instead of './work', will be broken on transfer without this option, -a override)
# r 		-	copy all files within directories	
# no-p		- 	do NOT preserve permissions		(necessary to allow file sharing on RDSS, -a override)
# no-g 		- 	do NOT preserve groups			(necessary to allow file sharing on RDSS, -a override)
# $1		-	destination directory			(will be replaced with <DESTINATION> when you run the script)

##########################################
# script definition			 #
##########################################

# set logging variables
d=$(date +%Y%m%d-%H%M%S)
logFile="nextflow-rsync-$d.log"
outFile="qc-pipeline-$d.tar.gz"

# log start
printf "\n
###############################################################################################\n
\n
nextflow-rsync.sh\n
\n
###############################################################################################\n
run time 			: $d\n
source 				: $(pwd)\n
results dir         : $1\n
destination 			: $2\n
###############################################################################################\n
Beginning file compression...\n
###############################################################################################\n\n
" >> $logFile

# file compression
find $( nextflow log last ; echo $1 ) -exec realpath --relative-to=. '{}' ';' | sort | uniq | tar -czvf $outFile -T - >> $logFile

# log end of compression, start of file transfer
printf "\n
###############################################################################################\n
File compression complete.\n
Beginning file transfer...\n
###############################################################################################\n\n
" >> $logFile

# file transfer
rsync -avhPLr --no-p --no-g $outFile $2 >> $logFile

# log end of file transfer
printf "\n
###############################################################################################\n
File transfer complete.\n
###############################################################################################\n\n
" >> $logFile

# recommended follow-up to permanently delete compressed files:
# nextflow clean
# rm -r work
