Misc: Various programms for read mapping with optimum search schemes
===================================================================================

Installation
-------------------------

The following instructions assume Linux or OS X. For more information, including Windows instructions, refer to the `SeqAn getting started tutorial <http://trac.seqan.de/wiki/Tutorial/GettingStarted>`_.

Software requirements
~~~~~~~~~~~~~~~~~~~~~

**A modern C++11 compiler with OpenMP 3.0 extensions is required to build Yara. If unsure, use GNU G++ 4.9 or newer.**

* Git.
* CMake 3.2 or newer.
* G++ 4.9 or newer.

Download
~~~~~~~~

Misc sources downloaded by executing:

::

  $git clone https://github.com/svnbgnk/misc.git


Configuration
~~~~~~~~~~~~~

Create a build project by executing CMake as follows:

::

  $ mkdir misc-build
  $ cd misc-build
  $ cmake ../misc

Build
~~~~~

Invoke make as follows:

::

  $ make

Usage
-----


testOSSAlignments
~~~~~~~~~~~~~~~~~~~

Count alternative Alignments of a single Occurrence when search with up 3 errors with my Version of the OSS

::

  $ testOSSAlignments -l 100 -mv -e 3




MapMap Short Instructions
~~~~~~~~~~~~~~~~~~~

MapMap setup (requires 6-7 GB of RAM)
::

 $ git clone --recurse-submodules https://github.com/svnbgnk/dream_yara.git
 $ cd dream_yara/include/seqan/
 $ git remote add upstream https://github.com/svnbgnk/seqan.git
 $ git checkout upstream/mappa 
 $ cd ../../..
 $ mkdir mapmap-build
 $ cd mapmap-build
 $ cmake ../dream_yara
 $ make all



Acquire hg38.fa
::

 $ /srv/public/svnbngk/Data/reference
 $ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
 $ gunzip hg38.fa.gz
 $ mkdir bin
 $ cp hg38.fa bin/0.fa

Build Index (requires 200GB of secondary memory)
::
 $ ./dream_yara_indexer --threads 8 --output-prefix /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reference/bin/*.fa -td /srv/public/svnbngk/tmp/

Computing sequence mappability and bit vectors with TH = 10
::
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 10 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability10E3

All mapping with up to 3 errors
::
 $ ./dream_yara_mapper /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbgnk/Data/reads/illumina/illumina_1.fa -t 1 -b 1 -ft none -e 3 -s 3 -o result.sam -vv 

Stratified all-mapping with strata 2 and 3 errors
::
 $ ./dream_yara_mapper /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbgnk/Data/reads/illumina/illumina_1.fa -t 1 -b 1 -ft none -e 3 -s 3 -o result.sam -vv

Mapping with sequence mappability up to 3 errors
::
 $ ./dream_yara_mapper /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbgnk/Data/reads/illumina/illumina_1.fa -t 1 -b 1 -ft none -e 3 -s 3 -m /srv/public/svnbngk/Data/hg38_N_index/mappability10E3/ -o result.sam -vv


Complete replication of results
-------------------------------

testOSSAlignments
~~~~~~~~~~~~~~~~~~~

Default OSS
::
 $ ./testOSSAlignments -l 100 -e 1
 $ ./testOSSAlignments -l 100 -e 2
 $ ./testOSSAlignments -l 100 -e 3
 $ ./testOSSAlignments -l 100 -e 4
With 1 read error
 $ ./testOSSAlignments -l 100 -e 4 -em -m 1

Simulating on MapMap OSS
 $ ./testOSSAlignments -l 100 -mv -e 1
 $ ./testOSSAlignments -l 100 -mv -e 2
 $ ./testOSSAlignments -l 100 -mv -e 3
 $ ./testOSSAlignments -l 100 -mv -e 4
With 1 read error
 $ ./testOSSAlignments -l 100 -mv -e 3 -em -m 1

MapMap
~~~~~~~~~~~~~~~~~~~

Computation of Sequence Mappability
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 5 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability5E3 
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 5 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability5E2 
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 1 -T 5 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability5E1 
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 5 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability5H3 
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 5 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability5H2 
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 1 -T 5 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability5H1

 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 10 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability10E3
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 10 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability10E2
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 1 -T 10 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability10E1
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 10 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability10H3
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 10 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability10H2
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 1 -T 10 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability10H1

 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 100 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability100E3
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 100 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability100E2
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 100 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability100E1
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 100 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability100H3
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 100 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability100H2
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 1 -T 100 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability100H1

 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 1000 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability1000E3
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 1000 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability1000E2
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 1 -T 1000 -s 0 -t 20 -o 35 -v -i -O /srv/public/svnbngk/Data/hg38_N_index/mappability1000E1
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 3 -T 1000 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability1000H3
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 2 -T 1000 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability1000H2
 $ ./dream_yara_mappability /srv/public/svnbngk/Data/hg38_N_index/ -b 1 -K 100 -E 1 -T 1000 -s 0 -t 20 -o 35 -v -O /srv/public/svnbngk/Data/hg38_N_index/mappability1000H1

Using bashscripts in ./bashscripts

 $ ./benchmark_v2.sh master.log /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa
 $ ./benchmark_hamming.sh masterHamming.log /srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa
 $ ./benchmark_f2.sh map5.log srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa /srv/public/svnbngk/Data/hg38_N_index/ mappability5
 $ ./benchmark_f2.sh map10.log srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa /srv/public/svnbngk/Data/hg38_N_index/ mappability10
 $ ./benchmark_f2.sh map100.log srv/public/svnbngk/Data/hg38_N_index/ /srv/public/svnbngk/Data/reads/illumina/illumina_1.fa /srv/public/svnbngk/Data/hg38_N_index/ mappability100



