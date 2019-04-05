#./dream_yara_mapper /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa  -t 1 -b 1 -ft none -e 3 -s 0 -o empty.sam -vv

LOGNAME="$1"
INDEXPATH="$2"
#/srv/public/svnbngk/Data/hg38_N_index/
READSPATH="$3"
#/srv/public/svnbngk/Data/reads/illumina/illumina_1.fa
MTEMPPATH="$4"
#/srv/public/svnbngk/Data/hg38_N_index/
SPEC="$5"
#mappability or mappability100
BIN=1

echo "Start Benchmark" > $LOGNAME
echo "Start Benchmark"
echo "Test All Mapping with Mappability" >> $LOGNAME
echo "Test All Mapping with Mappability"

for e in `seq 1 3`;
do
	echo "Error: "$e
	echo "Error: "$e >> $LOGNAME
	MTEMP="$MTEMPPATH"
	MTEMP+="$SPEC"
	MTEMP2="$MTEMP"
	MTEMP2+="E"
	MTEMP+="H"
	MTEMP+="$e/"
	MTEMP2+="$e/"
	echo "$MTEMP"
	echo "$MTEMP2"
	echo "$MTEMP" >> $LOGNAME
	./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e $e -s $e -m $MTEMP -o empty.sam -vv -rb 100000 >> $LOGNAME
	echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
	echo "$MTEMP2" >> $LOGNAME
	./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e $e -s $e -m $MTEMP2 -o empty.sam -vv -rb 100000 >> $LOGNAME
	echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
done

#./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b 2 -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000 >> $LOGNAME
echo "End of Benchmark" >> $LOGNAME
echo "End of Benchmark"
lscpu >> $LOGNAME

#./dream_yara_mapper ../../Data/myhg5_N_index/ ../../Data/reads/myhg5_N/default_10k.fa -t 1 -b 2 -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000
#./dream_yara_mapper ../../Data/myhg5_N_index/ ../../Data/reads/myhg5_N/default_10k.fa -t 1 -b 2 -ft none -e 3 -s 3 -m ../../Data/myhg5_N_index/mappabilityH3/ -o empty.sam -vv -rb 100000
