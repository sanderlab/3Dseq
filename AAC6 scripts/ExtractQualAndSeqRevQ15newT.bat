#!/bin/bash

file1=$2'_SeqsQ15rev.txt'
file2=$2'_QscoresQ15rev.txt'
file3=$2'_SeqsQ15rev.fas'

sed -n 'n;p' $1 > dummyfile2.txt
sed -n 'p;n' dummyfile2.txt | sed 's/CTTGGTCTCTACAG/ /' | sed 's/TTTGGTCTCCTGT/ /' | cut -d' ' -f2 > dummyfile3.txt
awk '!s{match($0, /CTTGGTCTCTACAGTTA.*CATTTTGGTCTCCTGT/);s=1;next} {print substr($0,RSTART+14,RLENGTH-27);s=0}' dummyfile2.txt > dummyfile4.txt
rm dummyfile2.txt
paste dummyfile4.txt dummyfile3.txt |grep -v -e "[\!\"\#\$\%\&\'\(\)\*\+\,\.\/]" | grep -v -e "[\-]" | grep -e "^.\{620,\}$" > dummyfile5.txt
rm dummyfile3.txt
rm dummyfile4.txt
cut -f1  dummyfile5.txt > $file2
cut -f2  dummyfile5.txt > $file1
rm dummyfile5.txt
awk 'sub("^", ">\n")' $file1 > $file3

