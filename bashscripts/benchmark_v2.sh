#./dream_yara_mapper /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa  -t 1 -b 1 -ft none -e 3 -s 0 -o empty.sam -vv

LOGNAME=$1
INDEXPATH=$2
#/srv/public/svnbngk/Data/hg38_N_index/
READSPATH=$3
#/srv/public/svnbngk/Data/reads/illumina/illumina_1.fa
MAPPABILITY=$4
#/srv/public/svnbngk/Data/hg38_N_index/

BIN=1
#Load Index and reference once before starting benchmark
./dream_yara_mapper /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa -t 1 -b 1 -ft none -e 3 -s 0 -o empty.sam -v -rb 500000

echo "Start Benchmark" > $LOGNAME
echo "Start Benchmark"
echo "Test All Mapping" >> $LOGNAME
echo "Test All Mapping"
for e in `seq 0 4`;
do
    echo "Error: "$e
    echo "Error: "$e >> $LOGNAME
    ./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e $e -s $e -o empty.sam -vv -rb 100000 >> $LOGNAME
    echo "------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
done
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME

echo "Test Stratified Mapping" >> $LOGNAME
echo "Test Stratified Mapping"

for s in `seq 0 2`;
do
    echo "Strata: "$s
    echo "Strata: "$s >> $LOGNAME
    ./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e 3 -s $s -o empty.sam -vv -rb 100000 >> $LOGNAME
    echo "------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
done
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
echo "End of Benchmark"
echo "Testing Influence of Options"

echo "No Suffix Filter" >> $LOGNAME
echo "No Suffix Filter"
./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000 -nS >> $LOGNAME
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
echo "No Delay" >> $LOGNAME
echo "No Delay"
./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000 -nD >> $LOGNAME
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
echo "No ITV" >> $LOGNAME
echo "No ITV"
./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000 -nD -nI >> $LOGNAME
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
echo "No ITV & No Suffix Filter" >> $LOGNAME
echo "No ITV & No Suffix Filter"
./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000 -nD -nI -nS >> $LOGNAME
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
echo "Hamming Distance" >> $LOGNAME
echo "Hamming Distance"
./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000 -HD >> $LOGNAME
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
echo "Hamming Distance & No Suffix Filter" >> $LOGNAME
echo "Hamming Distance & No Suffix Filter"
./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000 -HD -nS >> $LOGNAME
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME

echo "Start Yara Benchmark"
echo "Test All Mapping" >> $LOGNAME
echo "Test All Mapping"
for e in `seq 0 4`;
do
    echo "Error: "$e
    echo "Error: "$e >> $LOGNAME
    ./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e $e -s $e -o empty.sam -vv -rb 100000 -y full -of >> $LOGNAME
    echo "------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
done
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME

echo "Test Stratified Mapping" >> $LOGNAME
echo "Test Stratified Mapping"

for s in `seq 0 2`;
do
    echo "Strata: "$s
    echo "Strata: "$s >> $LOGNAME
    ./dream_yara_mapper $INDEXPATH $READSPATH -t 1 -b $BIN -ft none -e 3 -s $s -o empty.sam -vv -rb 100000 -y full -of >> $LOGNAME
    echo "------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
done
echo "--------------------------------------------------------------------------------------------------------------------------" >> $LOGNAME
echo "End of Benchmark"

lscpu >> $LOGNAME


#./dream_yara_mapper ../../Data/myhg5_N_index/ ../../Data/reads/myhg5_N/default_10k.fa -t 1 -b 2 -ft none -e 3 -s 3 -o empty.sam -vv -rb 100000
#./dream_yara_mapper ../../Data/myhg5_N_index/ ../../Data/reads/myhg5_N/default_10k.fa -t 1 -b 2 -ft none -e 3 -s 3 -m ../../Data/myhg5_N_index/mappabilityH3/ -o empty.sam -vv -rb 100000
