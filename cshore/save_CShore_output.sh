#!/bin/sh

cshorefiles="./O*"

for thisfile in $cshorefiles; do
   inname=$(basename "$thisfile" ".*")
   outname="master_$inname.txt"
   
   echo "$1: coast $2 profile $3" >> $outname
   echo >> $outname

   cat $thisfile >> $outname
   
   echo >> $outname
   echo "*******************************************************************************" >> $outname
done
 
