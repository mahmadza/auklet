#convert fasta sequences to tab format
#fasta format:
#>SEQ_ID
#AATTCCGGTT
#tab format:
#SEQ_ID     AATTCCGGTT
#this format is handy to quickly many sequences in a tabular-like format

awk -v OFS="\t" '{

    #if sees ">", it means starting of a new sequence
    if($1 ~ /^>/)
    {
        sub(">","",$1)
        name[++count]=$1
    }

    else
        seq[count]=seq[count] $1
        #concatenate sequence
}

END{
    #print sequence
    for(i=1;i<=count;i++)
        print name[i],seq[i]
}'
