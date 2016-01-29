#!/bin/bash
set -o errexit
set -o pipefail


################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
  echo >&2 "
$(basename $0) - expand the region around BED start and end column, and fix it

USAGE: cat *.bed | $(basename $0) -r <radius of expansion length> -g <genome assembly>

"
  exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "r:g:" o
do
  case "$o" in
      r) radius="$OPTARG";;
      g) genome="$OPTARG";;
     \?) exit 1;;
  esac
done

################################################################################
# Run program
################################################################################



awk -vflank=$radius -vOFS="\t" -vgenome_file=/groups/stark/genomes/chrom/${genome}.chrom.sizes 'BEGIN{
     while((getline<genome_file)>0)
          chrsize[$1]=$2  
}

($1 in chrsize){

     if($2>=$3) next     #skip illegal lines from input
     if($2<0 && $3<0) next

     $2=$2-flank
     $3=$3+flank

     #correct region if doesnt make sense
     if($2<0) $2=0
     if($3>chrsize[$1]) $3=chrsize[$1]

     print

}' | bedSort stdin stdout
