#!/bin/bash

file1=$2'_SeqsQ15.txt'
file2=$2'_QscoresQ15.txt'
file3=$2'_SeqsQ15.fas'

sed -n 'n;p' $1 > dummyfile2.txt
sed -n 'p;n' dummyfile2.txt | sed 's/ACAGGAGACCAAA/ /' | sed 's/CTGTAGAGACCAAG/ /' | cut -d' ' -f2 > dummyfile3.txt
awk '!s{match($0, /ACAGGAGACCAAAATG.*TAACTGTAGAGACCAAG/);s=1;next} {print substr($0,RSTART+13,RLENGTH-27);s=0}' dummyfile2.txt > dummyfile4.txt
rm dummyfile2.txt
paste dummyfile4.txt dummyfile3.txt |grep -v -e "[\!\"\#\$\%\&\'\(\)\*\+\,\.\/]" | grep -v -e "[\-]" | grep -e "^.\{620,\}$" > dummyfile5.txt
rm dummyfile3.txt
rm dummyfile4.txt
cut -f1  dummyfile5.txt > $file2
cut -f2  dummyfile5.txt > $file1
rm dummyfile5.txt
awk 'sub("^", ">\n")' $file1 > $file3