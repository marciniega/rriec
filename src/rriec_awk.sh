#!/usr/bin/env bash
set -e

get_stats () { awk '(NR==1){for (i=2;i<=NF;i++){lab[i]=$i}}
                       (NR>1){c++;for (i=2;i<=NF;i++){sum[i]+=$i;sumsq[i]+=($i)^2}}
                        END{for (i=2;i<=NF;i++){
                             ave=sum[i]/(c);
                             sutwo=(sum[i]^2)/(c);
                             std=sqrt((sumsq[i]-sutwo)/(c));
                             printf "%-12s  %10.3f +/- %10.3f\n",lab[i],ave,std}}' $1 ;}

get_stats_res () { awk '/#total/{f++;r=1}
                        (f==1){k++;n[k]=$1}
                        {sum[r]+=$2;sumsq[r]+=($2)^2;r++}
                        END{for (i=1;i<r;i++){
                                ave=sum[i]/f;
                                sutwo=(sum[i]^2)/f;
                                std=sqrt((sumsq[i]-sutwo)/f);
                                printf "%-12s  %10.3f +/- %10.3f\n",n[i],ave,std}}' $1 ;}

get_value_from_keys () { awk -v t="$1" -v kt=2.49434 '(NR==FNR&&!(/\#/||/\&/||/\@/)){c++;a[c]=$NF;next}
                             {b[FNR]=a[$1];avg+=b[FNR]}
                             END{avg=avg/FNR;
                                  for (x=1;x<=FNR;x++){
                                  delta=b[x]-avg;
                                  val=delta/kt;
                                  r+=exp(val)};
                                  print t,kt*log(r/FNR)}' "$2" "$3";}

get_value () { awk -v t="$1" -v kt=2.49434 '(!(/\#/||/\&/||/\@/)){c++;b[c]=$NF;avg+=$NF}
                             END{avg=avg/c;
                                  for (x=1;x<=c;x++){
                                  delta=b[x]-avg;
                                  val=delta/kt;
                                  r+=exp(val)};
                                  print t,kt*log(r/c)}' "$2";}

get_bootstrap_from_keys () { awk '(NR==1){print $0;next}
                                  (NR==FNR&&!(/\#/||/\&/||/\@/)){c++;a[c]=$0;next}
                                  (NR>FNR){print a[$1]}' "$1" "$2";}

compute_entropy_res () { awk -v kt=2.49434 -v r=$2 -v t=$3 'BEGIN{tt=(t=="")?0:t}
                                                       (NR==1){f=1;v=1;c=0;next}
 						  	(/Frame/&&f!=$2){val[f]=v;f=$2;v=1;c=0;next}
   							{c++;v*=(c==r)?1:$2}
   							END{val[f]=v;for (x=1;x<=f;x++){sum+=val[x]};
                					ave=sum;
                					const=-kt*log(f);
                					print const+(kt*log(ave))-tt}' $1 ;}

compute_entropy_res_pool () { awk -v kt=2.49434 -v r=$2 -v t=$3 -v n=$4 'BEGIN{tt=(t=="")?0:t}
                                                       (NR==1){f=1;v=1;c=0;next}
 						  	(/Frame/&&f!=$2){val[f]=v;f=$2;v=1;c=0;next}
   							{c++;v*=(c==r)?1:$2}
   							END{val[f]=v;for (x=1;x<=f;x++){sum+=val[x]};
                					ave=sum;
                					const=-kt*log(f);
                					print n"   "const+(kt*log(ave))-tt}' $1 > cerp_"$2".txt ;}

compute_partial_res_sums () { awk -v kt=2.49434 '(NR==1){for (x=3;x<=NF;x++){n[x-1]=$x}}
                                                 (NR>1){for (x=2;x<=NF;x++){a[NR][x]=$x}}
	                                         END{for (x=2;x<=NF;x++){
       	                                                 for (z=2;z<=NR;z++){
	     	                            	             ave[x]+=a[z][x]};
	     	                                         ave[x]/=(NR-1)};
	     	                                     for (z=2;z<=NR;z++){
                                                          print "Frame "z-1;
	     	                            	          for (x=2;x<=NF;x++){
	     	                            	              diff=a[z][x]-ave[x];
	     	                            	              val=diff/kt;
	     	                            	              print n[x],exp(val)}}
                                                    }' "$1" ;} 

checkint() { if ! [[ "$1" =~ ^[0-9]+$ ]]; then
	           echo "Sorry integers only for $2 "; 
		   exit 1;
		   fi ;}

checkfile() { if [ ! -f $1 ]; then
                echo "File $1 was not found.\n  Nothing was done.";
                exit 1 ;
              fi
}

example() {
        echo -e "example: $script_name -f contrib_MM.dat -s 300 -b 30 -c 8"
}

show_help() {
        echo "-f --contrb_file          Residue contribution file from g_mmpbsa. "
        echo "-s --start               Frame (not time) at which the analysis will start. Default [1]. "
        echo "-b --bootstrap           Number of bootstrap runs to perform. Default [30]."
        echo "-c --cores               Number of cores to use in the computation. Default [ nproc/2 ]."
        example
}
