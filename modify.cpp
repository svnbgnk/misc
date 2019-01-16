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

    addOption(parser, ArgParseOption("R", "reference", "Path to the reference", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "reference");

    addOption(parser, ArgParseOption("r", "reads", "Path to the reads", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "reads");

    addOption(parser, ArgParseOption("O", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString referencePath, readsPath, outputPath;

    getOptionValue(referencePath, parser, "reference");
    getOptionValue(readsPath, parser, "reads");
    getOptionValue(outputPath, parser, "output");
    bool verbose = isSet(parser, "verbose");


    srand (time(NULL));

    StringSet<CharString> idsReads;
    StringSet<DnaString> reads;
    SeqFileIn seqFileInReads(toCString(readsPath));

    readRecords(idsReads, reads, seqFileInReads);

    StringSet<CharString> ids;
    StringSet<DnaString> reference;
    SeqFileIn seqFileIn(toCString(referencePath));

    readRecords(ids, reference, seqFileIn);
    uint32_t total_length = lengthSum(reference);



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
    }
//     writeRecord(seqFileout, idsReads, sampledReads);

    close(seqFileout);
    std::cout << "Finished!\n";

    return 0;
}



