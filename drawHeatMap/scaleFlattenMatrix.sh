#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
  echo >&2 "
$(basename $0) - Trim boundary values and put symmteric min, max values as a fake row of a density heatmap matrix

USAGE: $(basename $0) -i <.txt heatmap matrix, densityHeatMap.sh output> -t <max value, int>
"
  exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################


while getopts "i:t:" o
do
  case "$o" in
    i) input="$OPTARG";;
    t) trim_value="$OPTARG";;
   \?) exit 1;;
  esac
done

if [ -z "$input" ]; then
     echo >&2 "ERROR: -i is required!"
     exit 1
fi

if [ ! -f "$input" ]; then
    echo >&2 "ERROR: file $input not found!"
    exit 1
fi

if [ -z "$trim_value" ]; then
     echo >&2 "ERROR: -t(trimming value) is required!"
     exit 1
fi

#check columns
uniq_col_num=$( cat $input | awk '{print NF}' | sort -k1,1n | uniq | wc -l )

if [ "$uniq_col_num" != "1" ]; then
    echo >&2 "ERROR: number of columns are not the same!"
    exit 1
fi

################################################################################
# Run program
################################################################################
#run check to make sure columns>3

cat $input | \
  awk -vboundary=$trim_value -vOFS="\t" '{

               for (i=1;i<=NF;i++)
               {
                    if($i!="NA" && $i>boundary) $i=boundary                         #if value more than +boundary
                         else
                            if ($i!="NA" && $i<(-boundary)) $i=-boundary          #if value less than -boundary
                              #else tmp=$i                                      #if value is within (-boundary,+boundary)
               }
               print $0
            }
            
            END{
            
                #put fake row
                line=boundary "\t" (-boundary)
                
                for(i=3;i<=NF;i++)
                     line=line"\t0"
                print line
            
            }'

