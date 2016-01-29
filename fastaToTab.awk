#convert fasta sequences to tab format
#fasta format:
#>SEQ_ID
#AATTCCGGTT
#tab format:
#SEQ_ID     AATTCCGGTT


awk -v OFS="\t" '{

    if($1 ~ /^>/)
    {
        sub(">","",$1)
        name[++count]=$1
    }

    else
        seq[count]=seq[count] $1

}

END{
    for(i=1;i<=count;i++)
        print name[i],seq[i]
}'
