#!/bin/sh

#    This file is part of mpiMarkDup
#    
#    Copyright Institut Curie 2020
#    
#    This software is a computer program whose purpose is to sort SAM file and mark duplicates.
#    
#    You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
#    
#    The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND. Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data. 
#    
#    The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

# Author : Firmin Martin
# Last modification : 2018 Apr. 14
# Given a directory containing sam file, compare output of MarkDuplicates and mpiMarkDup
# SYNOPSIS 
#       ./runCmpOnDir.sh DIR  
# NOTE
#       User need modify *user variables* of cmpWithMD before using

LOGNAME=`date '+%Y-%m-%d_%H:%M:%S'`cmp.log

########################
######## SCRIPT ########
########################

FILES="$1/*.sam"

DIRSUM=0
DIREQ=0
DIRNE=0

for f in $FILES
do 
   echo "Processing `basename $f` ..."  
   . ./cmpWithMD.sh $f 2> /dev/null
   DIRSUM=$(( DIRSUM + SUM ))
   DIREQ=$(( DIREQ + EQ ))
   DIRNE=$(( DIRNE + NE ))
   echo -e "`basename $f`\t$SUM\t$EQ\t$NE\t$OPTICAL" >> "$LOGNAME"
done
PERC="`awk "BEGIN { print ($DIREQ/$DIRSUM)*100 }"`%"
STAT="\ntotal reads = $DIRSUM, equal = $DIREQ, not equal = $DIRNE, equal ratio = $PERC"
for (( i = 0; i < ${#STAT}; i++ )); do
    printf "-"  
done
echo -e $STAT
echo -e "\nSum\t$DIRSUM\t$DIREQ\t$DIRNE\t$DIROPTICAL" >> "$LOGNAME"
