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

    addOption(parser, ArgParseOption("o", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("b", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("p", "samplingProb", "Probability of which a reads accepted for the ouput", ArgParseArgument::DOUBLE, "FLOAT"));
    setRequired(parser, "samplingProb");

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString referencePath, readsPath, outputPath;
    int batchSize = 100000;
    bool sampling = true;
    float samplingProb = 1;

//     getOptionValue(referencePath, parser, "reference");
    getOptionValue(readsPath, parser, "reads");
    getOptionValue(samplingProb, parser, "samplingProb");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(batchSize, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");

    if(samplingProb == 1)
        sampling = false;

    srand (time(NULL));

    StringSet<CharString> idsReads;
    StringSet<Dna5String> reads;
    SeqFileIn seqFileInReads(toCString(readsPath));

    SeqFileOut seqFileout(toCString(outputPath));


    while(true){
        reserve(idsReads, batchSize);
        reserve(reads, batchSize);
        readRecords(idsReads, reads, seqFileInReads, batchSize);
        std::cout << "batch: " << length(reads) << "\n";

        StringSet<Dna5String> reads1;
        reserve(reads1, length(reads));

        for(int i = 0; i < length(reads); ++i){
            //sampling
            if(sampling){
                float prob = (float) rand()/RAND_MAX;
                if(samplingProb < prob)
                    continue;
            }

            appendValue(reads1, reads[i]);
            if(verbose){
                std::cout << idsReads[i] << "\n";
                std::cout << reads[i] << "\n";
            }
        }

        for(uint32_t i = 0; i < length(reads1); ++i){
            writeRecord(seqFileout, idsReads[i], reads1[i]);
        }

        if(length(reads) < batchSize)
            break;



        clear(idsReads);
        clear(reads);
        if(verbose)
            std::cout << "Load next batch\n";
    }


    close(seqFileout);

    std::cout << "Finished!\n";

    return 0;
}



