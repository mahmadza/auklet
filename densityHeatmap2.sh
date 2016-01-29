#!/bin/bash
set -o errexit
set -o pipefail


################################################################################
# Parse input and check for errors
################################################################################

genome=dm3                                                               #default assembly
radius=20000                                                          #window size around the peak summit
extend_length=0                                                           #extend length of the reads. default=0 for paired-end reads
down_sample_size=50                                                       #size of downsampling of the individual nucleotide positions
pseudo_count=1                                                            #pseudo-count to be added to input in cof_ratio.py
fixbed_sh=/groups/stark/zabidi/work/TCT_STARR/utils/fixBed.sh           #path to fixBed.sh (similar to bedClip, but correct lines > or < chromosomes, and removes $3>$2 lines)
transpose_awk=/groups/stark/zabidi/work/TCT_STARR/utils/transpose.awk     #path to transposase.awk
z_score=1.67

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
  echo >&2 "
$(basename $0) - Calculate confidence ratio of read density from STARR-Seq or ChIP for visualization in heatmap

USAGE: $(basename $0) -r <reference experiment bed file> -e <experiment bigbed file> -i <input bigbed file> -o <output matrix .txt file> [OPTIONS]

-g species assembly. default $species
-w bp window around the summit. default $window_size
-x bp extend length. default for paired-end: $extend_length (600 for single-end STARR-Seq data, 150 for ChIP-Seq, 120 for DHS)
-d bp down sampling size. default $down_sample_size
-p pseudo count to be used during confidence ratio. default $pseudo_count
"
  exit 1
fi

################################################################################
# Parse input and check for errors
################################################################################

while getopts "r:e:i:o:g:w:x:d:p" o
do
    case "$o" in
        r) reference="$OPTARG";;
        e) experiment="$OPTARG";;
        i) input="$OPTARG";;
        o) outfile="$OPTARG";;
        g) genome="$OPTARG";;
        w) radius="$OPTARG";;
        x) extend_length="$OPTARG";;
        d) down_sample_size="$OPTARG";;
        p) pseudo_count="$OPTARG";;
        \?) exit 1;;
    esac
done

if [ -z "$reference" ]; then
    echo >&2 "ERROR: -r is required!"
    exit 1
fi

if [ ! -f "$reference" ]; then
    echo >&2 "ERROR: file $ref_experiment not found"
    exit 1
fi

if [ ! -f "$experiment" ]; then
    echo >&2 "ERROR: file $experiment not found"
    exit 1
fi

if [ ! -f "$input" ]; then
    echo >&2 "ERROR: file $input not found"
    exit 1
fi

if [ -z "$outfile" ]; then
    echo >&2 "ERROR: -o is required"
    exit 1
fi

################################################################################
# Run program
################################################################################


#1)create window of $window_size around the summit of the ref_experiment, sorted according to the peak rank
window=$(mktemp)
cat $reference | \
    awk -vOFS="\t" '{print $1,$7,$8,$4}' | \
        /groups/stark/zabidi/work/TCT_STARR/utils/expandBed.sh -g $genome -r $radius > $window


##########################
#   extend reads
##########################

#2)extend reads of experiment as well as input to $extend_length
extended_reads=$(mktemp)
if [ $extend_length -ne "0" ]; then
    process_id=""
    bigBedToBed $experiment stdout | \
        awk -vOFS="\t" -vext=$extend_length '{ if($6=="+"){$3=$2+ext} else {$2=$3-ext}; print}' | \
            $fixbed_sh -g $genome | bedSort stdin ${extended_reads}_experiment &
    process_id="$! "
    bigBedToBed $input stdout | \
        awk -vOFS="\t" -vext=$extend_length '{ if($6=="+"){$3=$2+ext} else {$2=$3-ext}; print}' | \
                $fixbed_sh -g $genome | bedSort stdin ${extended_reads}_input &
    process_id="$process_id $! "
else
    process_id=""
    bigBedToBed $experiment stdout > ${extended_reads}_experiment &
    process_id="$! "
    bigBedToBed $input stdout > ${extended_reads}_input &
    process_id="$process_id $! "
fi
wait $process_id


##########################
#   grab coverage
##########################
#print to different files
#easier to correct for the peaks that are too close to chromosome edges
process_id=""
coverage=$(mktemp)
for signal in experiment input; do
    bedtools coverage -a ${extended_reads}_${signal} -b $window -d | \
        awk  -vhandle1=$coverage -vhandle2=$signal '{print $0 > handle1"_"handle2"_"$4 }' &
    process_id="$process_id $! "
done
wait $process_id


##########################
#   correct coverage
##########################
#correct the peaks that are too close to chromosome edges
#artificially add 0's
for signal in experiment input; do

    #correct peaks that are too close to chromosome starts
    cat $reference | \
        awk -vflank=$radius -vOFS="\t" -vgenome_file=/groups/stark/genomes/chrom/${genome}.chrom.sizes 'BEGIN{
             while((getline<genome_file)>0)
                  chrsize[$1]=$2  
            }
            {    
                $2=$2-flank
                if($2<0) print $4
            }' | \
                while read peak; do
                    all_lines=$( cat ${coverage}_experiment_${peak} | wc -l )
                    cat ${coverage}_experiment_${peak} | \
                        awk -vOFS="\t" -vx=$radius -vhave=$all_lines \
                            '(NR==1){
                                total=x*2+1                             #total nucleotides that shouldve been covered
                                need=total-have                         #nucleotides that are lacking
                                for(i=1;i<=need;i++)
                                    print $1,$2,$3,$4,i,0            #print filler lines
                                
                                $5=NR+i-1
                                print $0                                   #print the real line
                                next
                            }
                            {
                                $5=NR+i-1                               #now correct the bin numbering, even though wont be used
                                print $0
                            }' > ${coverage}_${signal}_${peak}_corr
                    mv ${coverage}_${signal}_${peak}_corr ${coverage}_${signal}_${peak}
                done
                
    #correct peaks that are too close to chromosome ends
    cat $reference | \
        awk -vflank=$radius -vOFS="\t" -vgenome_file=/groups/stark/genomes/chrom/${genome}.chrom.sizes 'BEGIN{
             while((getline<genome_file)>0)
                  chrsize[$1]=$2  
            }
            {    
                $3=$3+flank
                if($3>chrsize[$1]) print $4
            }' | \
                while read peak; do
                    cat ${coverage}_${signal}_${peak} | \
                        awk -vOFS="\t" -vx=$radius '{print}
                            END{
                                final=x*2+1
                                for(i=(NR+1);i<=final;i++)
                                    print $1,$2,$3,$4,i,0
                            }' > ${coverage}_${signal}_${peak}_corr
                        mv ${coverage}_${signal}_${peak}_corr ${coverage}_${signal}_${peak}
                done
done


##########################
#  downsample coverage
##########################
#grab peak names
#and sort them here
all_peaks=$( cat $reference | awk -vOFS="\t" '{split($4,x,"_"); print $4,x[2]}' | sort -k2,2n | cut -f1 )

for signal in experiment input; do
    for peak in $all_peaks; do
        cat ${coverage}_${signal}_${peak} | \
            awk -vOFS="\t" -vevery=$down_sample_size \
                '{
                    inside_bin++
                    total_bin+=$6
                    
                    if(inside_bin>=every)
                    {
                        print $1,$2,$3,$4,++bin_no,total_bin/every
                        total_bin=0
                        inside_bin=0
                    }
                }
                END{
                    print $1,$2,$3,$4,++bin_no,total_bin/inside_bin
                }'
    done > ${coverage}_${signal}_downsize &             #print to just two separate files
    process_id="$process_id $! "
done
wait $process_id




##########################
#  downcorrect matrices
##########################

conf_matrix=$(mktemp)
totalreads_experiment=$(bigBedInfo $experiment | awk '($1=="itemCount:"){gsub(",","",$2);print $2}')
totalreads_input=$(bigBedInfo $input | awk '($1=="itemCount:"){gsub(",","",$2);print $2}')
#determine which is bigger, $totalreads_experiment or $totalreads_input
#normalize the bigger one to the smaller one
#and then round both of the matrices
if [ $totalreads_experiment -gt $totalreads_input ]; then
    paste ${coverage}_experiment_downsize ${coverage}_input_downsize | \
        awk -vupper=$totalreads_input -vlower=$totalreads_experiment -vadd=$pseudo_count -vOFS="\t" \
            '{print int($6*(upper/lower)+0.5)+add,int($12+0.5)+add,$4}' > ${coverage}_normalized
elif [ $totalreads_input -gt $totalreads_experiment ]; then
    paste ${coverage}_experiment_downsize ${coverage}_input_downsize | \
        awk -vupper=$totalreads_experiment -vlower=$totalreads_input -vadd=$pseudo_count -vOFS="\t" \
            '{print int($6+0.5)+add,int($12*(upper/lower)+0.5)+add,$4}' > ${coverage}_normalized
else
    #both have exactly the same read count
    paste ${coverage}_experiment_downsize ${coverage}_input_downsize | \
        awk -vadd=$pseudo_count -vOFS="\t" '{print int($6+0.5)+add,int($12+0.5)+add,$4}' > ${coverage}_normalized       
fi


############################
# compute confidence ratio
#   and output
############################
#compute confidence ratio
#create marix and print
cat ${coverage}_normalized | \
    conf_ratio.py -z $z_score | \
                awk -vOFS="\t" \
                    '(NR==1){
                        peak_seen[$4]++
                        line=$1
                        next}
                    (NR>1 && !peak_seen[$4]++){
                        print line
                        line=$1
                        next}
                    {line=line"\t"$1}
                        END{print line}' > $outfile
                        

############################
#       cleanup
############################
rm $window ${coverage}_normalized
for signal in experiment input; do
    
    rm ${extended_reads}_${signal} ${coverage}_${signal}_downsize
    
    for peaks in $all_peaks; do
        rm ${coverage}_${signal}_${peaks}
    done

done

#exit success
exit 0
