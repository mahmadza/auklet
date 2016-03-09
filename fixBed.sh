#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
  echo >&2 "
$(basename $0) - fix BED file according to genome assembly

USAGE: cat *.bed | $(basename $0) -g <genome assembly>

"
  exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "g:" o
do
  case "$o" in
      g) genome="$OPTARG";;
     \?) exit 1;;
  esac
done

################################################################################
# Run program
################################################################################


awk -vFS="\t" -vOFS="\t" -vgenome_file=/GENOME_DIRECTORY/${genome}.chrom.sizes 'BEGIN{
    while((getline<genome_file)>0)
        chrsize[$1]=$2
}

($1 in chrsize){
    if($2<0 && $3>0) $2=0
    if($2<0 && $3<0) next                   #remove illegal lines

    if($3>chrsize[$1]) $3=chrsize[$1]       #correct ends > chromosome size

    if($2>chrsize[$1]) next                 #remove lines which start > chromosome size
    if($2>=$3) next                         #skip illegal lines

    print
}'
