#!/bin/bash

#Format a bed12 file to a gtf file using kent tools
# for format description, see UCSC format web page

# USAGE

usage() {
    echo "#" >&2
    echo -e "# USAGE: `basename $0` <file.bed12> \n OR\n    cat <my_file> | `basename $0`">&2
}


BINUCSC=~tderrien/progs/UCSC/exec/
infile=$1

##################################

if [ -p /dev/stdin ];then

    #from STDIN
    cat /dev/stdin | $BINUCSC/bedToGenePred stdin stdout | $BINUCSC/genePredToGtf file stdin stdout -utr -source=genePredToGtf -honorCdsStat


elif [ $# -eq  1 ];then

   $BINUCSC/bedToGenePred $infile stdout | $BINUCSC/genePredToGtf file stdin stdout -utr -source=genePredToGtf -honorCdsStat


else
    echo "#Error! no argument  file or file empty !"
    usage;

fi
