#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>

#include<iostream>
#include<fstream>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace seqan;
using namespace std;





int main(int argc, char const * argv[])
{
    ArgumentParser parser("Search");


    addOption(parser, ArgParseOption("r1", "reads1", "Path to the reads file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "reads1");

    addOption(parser, ArgParseOption("r2", "reads2", "Path to the reads file", ArgParseArgument::INPUT_FILE, "IN"));
//     setRequired(parser, "reads2");

    addOption(parser, ArgParseOption("o", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("b", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString referencePath, readsPath1, readsPath2, outputPath;
    int batchSize1 = 100000;


//     getOptionValue(referencePath, parser, "reference");
    getOptionValue(readsPath1, parser, "reads1");
    getOptionValue(readsPath2, parser, "reads2");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(batchSize1, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");
    srand (time(NULL));



    int count1 = 0;
    string line;
    ifstream file(toCString(readsPath1));

    while (getline(file, line))
        count1++;

    int count2 = 0;
    ifstream file2(toCString(readsPath2));

    while (getline(file2, line))
        count2++;

    file.close();
    file2.close();


    StringSet<CharString> idsReads1;
    StringSet<Dna5String> reads1;
    SeqFileIn seqFileInReads1(toCString(readsPath1));

    StringSet<CharString> idsReads2;
    StringSet<Dna5String> reads2;
    SeqFileIn seqFileInReads2(toCString(readsPath2));

    SeqFileOut seqFileout(toCString(outputPath));


    int batchSize2 = (static_cast<float>(count2) / count1) * batchSize1;



    while(true){
        reserve(idsReads1, batchSize1);
        reserve(reads1, batchSize1);

        reserve(idsReads2, batchSize2);
        reserve(reads2, batchSize2);

        readRecords(idsReads1, reads1, seqFileInReads1, batchSize1);
        readRecords(idsReads2, reads2, seqFileInReads2, batchSize2);
        if(verbose)
            std::cout << "Batchsizes: " << length(reads1) << "\t" << length(reads2) << "\n";

        StringSet<Dna5String> mergeReads;
        StringSet<CharString> mergeIdsReads;
        reserve(mergeReads, batchSize1 + batchSize2);

        int i = 0, j = 0;

        while(i < length(reads1) || j < length(reads2)){
            float samplingProb1 = static_cast<float>(length(reads1) - i) / (length(reads1) - i + length(reads2) - j);
            if(verbose)
                std::cout << "Prob for taking read from first bin: " << samplingProb1 << "\n";
            //sampling
            if(i < length(reads1) && j < length(reads2)){

                float prob = (float) rand()/RAND_MAX;
                if(samplingProb1 > prob){
                    appendValue(mergeReads, reads1[i]);
                    appendValue(mergeIdsReads, idsReads1[i]);
                    ++i;
                }else{
                    appendValue(mergeReads, reads2[j]);
                    appendValue(mergeIdsReads, idsReads2[j]);
                    ++j;
                }
            }else if(i < length(reads1)){
                appendValue(mergeReads, reads1[i]);
                appendValue(mergeIdsReads, idsReads1[i]);
                ++i;
            }else if(j < length(reads2)){
                appendValue(mergeReads, reads2[j]);
                appendValue(mergeIdsReads, idsReads2[j]);
                ++j;
            }else{
                std::cout << "What?\n";
            }

            if(verbose){
                std::cout << mergeIdsReads[length(mergeIdsReads) - 1] << "\n";
                std::cout << mergeReads[length(mergeReads) - 1] << "\n\n";
            }
        }

        for(uint32_t k = 0; k < length(mergeReads); ++k){
            writeRecord(seqFileout, mergeIdsReads[k], mergeReads[k]);
        }

        if(length(reads1) < batchSize1 && length(reads2) < batchSize2)
            break;


        clear(idsReads1);
        clear(reads1);
        clear(idsReads2);
        clear(reads2);
    }


    close(seqFileout);

    std::cout << "Finished!\n";

    return 0;
}



