#./dream_yara_mapper /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa  -t 1 -b 1 -ft none -e 3 -s 0 -o empty.sam -vv

LOGNAME=$1
INDEXPATH=$2
#/srv/public/svnbngk/Data/hg38_N_index/
READSPATH=$3
#/srv/public/svnbngk/Data/reads/illumina/illumina_1.fa
MAPPABILITY=$4
#/srv/public/svnbngk/Data/hg38_N_index/

BIN=1

./dream_yara_mapper /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa -t 1 -b 1 -ft none -e 3 -s 0 -o empty.sam -v -rb 500000

echo "Start Benchmark" > $LOGNAME
echo "Start Benchmark"
echo "Test All Mapping" >> $LOGNAME
echo "Test All Mapping"
for e in `seq 0 3`;
do
    echo "Error: "$e
    echo "Error: "$e >> $LOGNAME
    ./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e $e -s $e -o empty.sam -vv -rb 100000 -HD >> $LOGNAME
    echo "------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
done
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME

echo "End of Benchmark"

lscpu >> $LOGNAME

#./dream_yara_mapper ../../Data/myhg5_N_index/ ../../Data/reads/myhg5_N/default_10k.fa -t 1 -b 2 -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000
#./dream_yara_mapper ../../Data/myhg5_N_index/ ../../Data/reads/myhg5_N/default_10k.fa -t 1 -b 2 -ft none -e 3 -s 3 -m ../../Data/myhg5_N_index/mappabilityH3/ -o empty.sam -vv -rb 100000
