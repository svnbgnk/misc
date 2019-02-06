#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>
#include <typeinfo>
#include <set>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace seqan;
using namespace std;





int main(int argc, char const * argv[])
{
    ArgumentParser parser("Search");
/*
    addOption(parser, ArgParseOption("R", "reference", "Path to the reference", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "reference");*/

    addOption(parser, ArgParseOption("r", "reads", "Path to the reads", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "reads");

    addOption(parser, ArgParseOption("o", "output", "Path to output file prefix", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("sp", "splitPoint", "Length of first read in matepair", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("b", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("p", "samplingProb", "Probability of which a reads accepted for the ouput", ArgParseArgument::DOUBLE, "FLOAT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString referencePath, readsPath, outputPath;
    int splitPoint = 0;
    int batchSize = 100000;
    bool sampling = true;
    float samplingProb = 1;

//     getOptionValue(referencePath, parser, "reference");
    getOptionValue(readsPath, parser, "reads");
    getOptionValue(splitPoint, parser, "splitPoint");
    getOptionValue(samplingProb, parser, "samplingProb");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(batchSize, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");

    if(samplingProb == 1)
        sampling = false;

    bool splitHalf;
    if(splitPoint == 0)
        splitHalf = true;

    srand (time(NULL));

    StringSet<CharString> idsReads;
    StringSet<Dna5String> reads;
    SeqFileIn seqFileInReads(toCString(readsPath));

    CharString outputPath1 = outputPath;
    CharString outputPath2 = outputPath;
    outputPath1 += + "_1.fa";
    outputPath2 += + "_2.fa";

    SeqFileOut seqFileout1(toCString(outputPath1));
    SeqFileOut seqFileout2(toCString(outputPath2));



    while(true){
        reserve(idsReads, batchSize);
        reserve(reads, batchSize);

        readRecords(idsReads, reads, seqFileInReads, batchSize);
        std::cout << "batch: " << length(reads) << "\n";

        StringSet<Dna5String> reads1;
        StringSet<Dna5String> reads2;
        reserve(reads1, length(reads));
        reserve(reads2, length(reads));

        for(int i = 0; i < length(reads); ++i){
            //sampling
            if(sampling){
                float prob = (float) rand()/RAND_MAX;
                if(samplingProb < prob)
                    continue;
            }
            if(splitHalf){
                int readlength = length(reads[i]);
                splitPoint = readlength/2;
            }
            appendValue(reads1, prefix(reads[i], splitPoint));
            appendValue(reads2, suffix(reads[i], splitPoint));
            if(verbose){
                std::cout << idsReads[i] << "\n";
                std::cout << prefix(reads[i], splitPoint) << "\n" << suffix(reads[i], splitPoint) << "\n\n";
            }
        }

        for(uint32_t i = 0; i < length(reads1); ++i){
            writeRecord(seqFileout1, idsReads[i], reads1[i]);
        }


        for(uint32_t i = 0; i < length(reads2); ++i){
            writeRecord(seqFileout2, idsReads[i], reads2[i]);
        }

        if(length(reads) < batchSize)
            break;

        clear(idsReads);
        clear(reads);
    }


    close(seqFileout1);
    close(seqFileout2);

/*

//     String<long unsigned int, seqan::Alloc<> >
    auto limits = stringSetLimits(reference);
    typedef Pair <uint64_t, uint64_t>       TPos;

    //start loop
    for(uint32_t i = 0; i < length(reads); ++i)
    {

        uint32_t rint = rand() % total_length;
        TPos result;
        posLocalize(result, rint, limits);
        if(verbose)
            cout << "selected pos: " << result << "\n";

        uint32_t result_end = getSeqOffset(result) + length(reads[i]);
        if(length(reference[getSeqNo(result)]) > result_end)
        {
            uint32_t k = 0;
            for(uint32_t j = getSeqOffset(result); j < result_end; ++j){
                reference[getSeqNo(result)][j] = reads[i][k];
                ++k;
                if(verbose)
                    std::cout << reference[getSeqNo(result)][j];
            }
            if(verbose)
                std::cout << endl;
        }else{
            --i;
            if(verbose)
                std::cout << "not valid position try again" << "\n"; //very unlikly for longer sequences
        }
    }


    SeqFileOut seqFileout(toCString(outputPath));
    for(uint32_t i = 0; i < length(reference); ++i){
        writeRecord(seqFileout, ids[i], reference[i]);
    }*/
//     writeRecord(seqFileout, idsReads, sampledReads);

//     close(seqFileout);
    std::cout << "Finished!\n";

    return 0;
}



