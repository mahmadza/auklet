#design primers to amplify complete sequence
#assume the sequences are in column $1 in file $sequence_file

#grab total number of sequences
all_seq=$(cat $sequence_file | wc -l)

#temporary file
cand_list=$(mktemp)

for seq_now in $(seq 1 $all_seq); do

        #design forward primer
        cut -f1 $sequence_file | \
            awk -vOFS="\t" -vmaxlen=30 -vseek=$seq_now '(NR==seek){
        
                stretch=substr($1,1,maxlen)
        
                #print all Tm till maxlen
                for(i=1;i<=length(stretch);i++)
                {
                    if(substr(stretch,i,1)=="A" || substr(stretch,i,1)=="T")
                    Tm+=2
                    if(substr(stretch,i,1)=="C" || substr(stretch,i,1)=="G")
                    Tm+=4
                    print i,Tm,substr(stretch,1,i)
                }
                }' | \
                awk -vOFS="\t" '{
        
                    #print only primers which 3-ends are C or G
                    last_nt=substr($3,length($3),1)
                    if(last_nt=="G" || last_nt=="C")
                        print
                    }' | \
                    awk -vmin_len=15 '(length($3)>=min_len)' > ${cand_list}_${seq_now}_fwd
                    
        #design reverse primer
        cut -f1 $sequence_file | rev | tr ATGC TACG | \
            awk -vOFS="\t" -vmaxlen=30 -vseek=$seq_now '(NR==seek){
        
                stretch=substr($1,1,maxlen)
        
                #print all Tm till maxlen
                for(i=1;i<=length(stretch);i++)
                {
                    if(substr(stretch,i,1)=="A" || substr(stretch,i,1)=="T")
                    Tm+=2
                    if(substr(stretch,i,1)=="C" || substr(stretch,i,1)=="G")
                    Tm+=4
                    print i,Tm,substr(stretch,1,i)
                }
                }' | \
                awk -vOFS="\t" '{
                
                    #print only primers which 3-ends are C or G
                    last_nt=substr($3,length($3),1)
                    if(last_nt=="G" || last_nt=="C")
                        print
                    }' | \
                    awk -vmin_len=15 '(length($3)>=min_len)' > ${cand_list}_${seq_now}_rev

        #grab total number of forward primers
        number_of_primers=$(cat ${cand_list}_${seq_now}_fwd | wc -l)
    
        #for each forward primer, find the best reverse
        for now in $(seq 1 $number_of_primers); do
            cat ${cand_list}_${seq_now}_fwd | \
                awk -vOFS="\t" -vrev_file=${cand_list}_${seq_now}_rev -vp=$now -vTm_opt=62 \
                    'BEGIN{
                        while((getline<rev_file)>0)
                        {
                            rev_primer[++count]=$3
                            rev_Tm[count]=$2
                        }
                        }
                    (NR==p){
                        for(i=1;i<=count;i++)
                        {
                            Tm_diff=($2-rev_Tm[i]>0)?($2-rev_Tm[i]):-($2-rev_Tm[i])

                            diff_Tm_opt_fwd=($2-Tm_opt>=0)?($2-Tm_opt):-($2-Tm_opt)
                            diff_Tm_opt_rev=(rev_Tm[i]-Tm_opt>=0)?(rev_Tm[i]-Tm_opt):-(rev_Tm[i]-Tm_opt)

                            print $0,length(rev_primer[i]),rev_Tm[i],rev_primer[i],Tm_diff,diff_Tm_opt_fwd,diff_Tm_opt_rev
                        }
                    }' | sort -k7,7n -k8,8n -k9,9n | head -n1
        
        done | sort -k8,8n -k9,9n -k1,1nr -k4,4nr | head -n1
    
    done > $results file here
