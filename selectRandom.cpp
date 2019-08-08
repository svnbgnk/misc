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

    addOption(parser, ArgParseOption("r2", "reads2", "Path to the reads paired reads if the exist", ArgParseArgument::INPUT_FILE, "IN"));

    addOption(parser, ArgParseOption("o", "output", "Path to output incase paired _1.fasta and _2 fasta is appended", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("b", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("p", "samplingProb", "Probability of which a reads accepted for the ouput", ArgParseArgument::DOUBLE, "FLOAT"));
//     setRequired(parser, "samplingProb");

    addOption(parser, ArgParseOption("s", "stopEarly", "Maximum Number of reads that should be read", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString referencePath, readsPath, readsPath2 = "", outputPath;
    int batchSize = 100000;
    uint64_t stopEarly = 0;
    bool sampling = true;
    float samplingProb = 1;

//     getOptionValue(referencePath, parser, "reference");
    getOptionValue(readsPath, parser, "reads");
    getOptionValue(readsPath2, parser, "reads2");

    getOptionValue(samplingProb, parser, "samplingProb");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(batchSize, parser, "batchSize");
    getOptionValue(stopEarly, parser, "stopEarly");
    bool verbose = isSet(parser, "verbose");

    if(samplingProb == 1)
        sampling = false;

    srand (time(NULL));

    if(length(readsPath2) == 0){
        StringSet<CharString> idsReads;
        StringSet<Dna5String> reads;
        SeqFileIn seqFileInReads(toCString(readsPath));
        SeqFileOut seqFileout(toCString(outputPath));
        uint64_t k = 0;
        while(true){
            reserve(idsReads, batchSize);
            reserve(reads, batchSize);
            readRecords(idsReads, reads, seqFileInReads, batchSize);
            if(verbose)
                std::cout << "batch: " << length(reads) << "\n";

            std::vector<uint32_t> sampledPos;

            for(int i = 0; i < length(reads); ++i){
                //sampling
                if(sampling){
                    float prob = (float) rand()/RAND_MAX;
                    if(samplingProb > prob)
                        sampledPos.push_back(i);
                }
                else
                {
                    sampledPos.push_back(i);
                }
            }

            if(verbose)
                std::cout << "Number of sampled reads: " << sampledPos.size() << "\n";

            for(uint32_t i = 0; i < length(sampledPos); ++i){
                writeRecord(seqFileout, idsReads[sampledPos[i]], reads[sampledPos[i]]);
            }

            k+= batchSize;
            if(length(reads) < batchSize || (stopEarly != 0 && stopEarly <= k))
                break;

            clear(idsReads);
            clear(reads);
            if(verbose)
                std::cout << "Load next batch\n";

        }
        close(seqFileout);
    }
    else
    {
        StringSet<CharString> idsReads;
        StringSet<Dna5String> reads;
        StringSet<CharString> idsReads2;
        StringSet<Dna5String> reads2;



        CharString outputPath1 = outputPath;
        outputPath1 += "_1.fasta";
        CharString outputPath2 = outputPath;
        outputPath2 += "_2.fasta";

        SeqFileIn seqFileInReads(toCString(readsPath));
        SeqFileOut seqFileout(toCString(outputPath1));

        SeqFileIn seqFileInReads2(toCString(readsPath2));
        SeqFileOut seqFileout2(toCString(outputPath2));

        uint64_t k = 0;
        while(true){
            reserve(idsReads, batchSize);
            reserve(reads, batchSize);
            reserve(idsReads2, batchSize);
            reserve(reads2, batchSize);
            readRecords(idsReads, reads, seqFileInReads, batchSize);
            readRecords(idsReads2, reads2, seqFileInReads2, batchSize);
            if(verbose)
                std::cout << "batch: " << length(reads) << "\n";

            std::vector<uint32_t> sampledPos;

            if(length(reads) != length(reads2)){
                std::cout << length(reads) << "\t" << length(reads2) << "\n";
                std::cerr << "Unequal amount of left paired and right paired reads\n";
                exit(0);
            }

            for(int i = 0; i < length(reads); ++i){
                //sampling
                if(sampling){
                    float prob = (float) rand()/RAND_MAX;
                    if(samplingProb > prob)
                        sampledPos.push_back(i);
                }
                else
                {
                    sampledPos.push_back(i);
                }
            }

            if(verbose)
                std::cout << "Number of sampled reads: " << sampledPos.size() << "\n";


            for(uint32_t i = 0; i < sampledPos.size(); ++i){
                writeRecord(seqFileout, idsReads[sampledPos[i]], reads[sampledPos[i]]);
                writeRecord(seqFileout2, idsReads2[sampledPos[i]], reads2[sampledPos[i]]);
            }

            k+= batchSize;
            if(length(reads) < batchSize || (stopEarly != 0 && stopEarly <= k))
                break;
            clear(idsReads);
            clear(reads);
            clear(idsReads2);
            clear(reads2);
            if(verbose)
                std::cout << "Load next batch\n";
        }
        close(seqFileout);
        close(seqFileout2);
    }

    std::cout << "Finished!\n";

    return 0;
}



