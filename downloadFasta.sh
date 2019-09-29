#!/bin/bash
set -o errexit
set -o pipefail

module load ncbi-blast+/2.4.0


################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
  echo >&2 "
$(basename $0) - download sequences using fastacmd

USAGE: cat *.bed | $(basename $0) -r <nucleotides around the fed coordinate to be downloaded> -g <genome assembly>

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

while read line; do

    chr=$(echo $line | awk '{print $1}')
    begin=$(echo $line | awk -vleft=$radius '{print $2+1-left}')        #input is bed file, realign for fastacmd
    finish=$(echo $line | awk -vright=$radius '{print $3+right}')
    strand=$(echo $line | awk '{if($6=="+"){print 1}else{print 2}}')    #fastacmd recognizes +strand as 1, -strand as 2
    id=$(echo $line | awk -vOFS="_" -vleft=$begin -vright=$finish '{if($6=="+") $NF="plus"; if($6=="-") $NF="minus";$2=left-1;$3=right;print $0}')
    #keep the id as bed

    #access database, then remove the line ">lcl|chr2L:132426-132526 No definition line found"
    #and put my own ID
    #$genome needs to point to the genome index. modify appropriately
    fastacmd -d $genome -s $chr -L$begin,$finish -S $strand | \
            awk -vchrmsm=$chr -vstart=$begin -vfin=$finish -vname=$id '(!/^>/){line=line $1}
                END{print ">"name; print line}'
                
done
