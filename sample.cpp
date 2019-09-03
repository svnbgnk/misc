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
#include <algorithm>

using namespace seqan;
using namespace std;





int main(int argc, char const * argv[])
{
    ArgumentParser parser("Search");

    addOption(parser, ArgParseOption("R", "reference", "Path to the reference", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "reference");

    addOption(parser, ArgParseOption("O", "output", "Path to output directory", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("n", "number", "number of read Samples", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("l", "length", "length of read Samples", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("e", "errors", "errors in modified reads", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("nm", "nmod", "number of modified reads from Sample reads", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("pi", "probInsertion", "probability Insertion", ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, ArgParseOption("pd", "probDeletion", "probability Deletion", ArgParseArgument::DOUBLE, "DOUBLE"));

    addOption(parser, ArgParseOption("nn", "noN", "Do not allow n in the sampled reads"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString referencePath, outputPath;
    int n = 10, l = 100, e = 3, nmod = 10;
    float probInsertion = 0.0, probDeletion = 0.0;


    getOptionValue(referencePath, parser, "reference");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(n, parser, "number");
    getOptionValue(nmod, parser, "nmod");
    getOptionValue(l, parser, "length");
    getOptionValue(e, parser, "errors");
    getOptionValue(probInsertion, parser, "probInsertion");
    getOptionValue(probDeletion, parser, "probDeletion");
    bool noN = isSet(parser, "noN");
    bool verbose = isSet(parser, "verbose");


    srand (time(NULL));
    StringSet<CharString> idsReads;
    StringSet<Dna5String> sampledReads;

    SeqFileIn seqFileIn(toCString(referencePath));
    StringSet<CharString> ids;
    StringSet<Dna5String> reference;

    readRecords(ids, reference, seqFileIn);

    uint32_t total_length = lengthSum(reference);

//     String<long unsigned int, seqan::Alloc<> >
    auto limits = stringSetLimits(reference);
    typedef Pair <uint64_t, uint64_t>       TPos;

    //start loop
    uint32_t k = 0;
    for(uint32_t i = 0; length(sampledReads) < n; ++i)
    {

        uint32_t rint = rand() % total_length;
        TPos result;
        posLocalize(result, rint, limits);
        bool valid = false;
        if(verbose)
            cout << "result: " << result << "\n";
        if((length(reference[getSeqNo(result)]) > getSeqOffset(result) + l))
            valid = true;
        Dna5String tmpread = infix(reference[getSeqNo(result)], getSeqOffset(result), getSeqOffset(result) + l);
        if(valid && noN){
            for(uint32_t p = 0; p < length(tmpread); ++p){
                if (tmpread[p] == 'N'){
                    valid = false;
                    if (verbose)
                        std::cout << "Contains N: " << tmpread << "\n";
                    break;
                }
            }
        }
        if(valid)
        {
            if(verbose)
                std::cout << "valid" << "\n";

            if (verbose)
                std::cout << "sampled read: " << tmpread << "\n";
            appendValue(sampledReads, tmpread);
            appendValue(idsReads, to_string(k));
            ++k;
        }else{
            if(verbose)
                std::cout << "not valid try again" << "\n";
        }
    }

    //save sampled Reads
    auto outputPathSampled = outputPath;
    outputPathSampled += ".fa";
    SeqFileOut seqFileout(toCString(outputPathSampled));
    for(uint32_t i = 0; i < length(sampledReads); ++i){
        writeRecord(seqFileout, idsReads[i], sampledReads[i]);
    }
//     writeRecord(seqFileout, idsReads, sampledReads);
    close(seqFileout);

    if (nmod == 0)
        return 0;

    StringSet<DnaString> modifiedReads;


    for(uint32_t k = 0; k < length(sampledReads); ++k)
    {
	if(verbose)
             std::cout << "At sampled read: " << k << "\n";
        for(uint32_t j = 0; j < nmod; ++j)
        {
            DnaString tmpread;
            if(j <= nmod/2){
                tmpread = sampledReads[k];
            }else{
                if(verbose)
                    cout << "reverseComplement\n";
                Dna5StringReverseComplement myModifier(sampledReads[k]);
                tmpread = myModifier;
            }

            std::vector<int> v;
            v.reserve(e);
            for(int i = 0; i < e; ++i){
                uint32_t rpos = rand() % length(tmpread);
                while(std::find(v.begin(), v.end(), rpos) != v.end()){
                    rpos = rand() % length(tmpread);
                }
                v.push_back(rpos);
                if(verbose)
                    std::cout << rpos << "\n";


                float prob = (float) rand()/RAND_MAX;
                if(verbose)
                    std::cout << "prob: " << prob << "\n";

                //Mismatch
                if(probInsertion + probDeletion < prob){
                    int rValue = rand() % 3;
                    int cValue = tmpread[rpos].value;
                    if(rValue >=  cValue)
                        ++rValue;
                    Dna rBase(rValue);
                    tmpread[rpos] = rBase;
                    if(verbose){
                        std::cout << "Mismatch: \n";
                        std::cout << "rBase " << rBase << "\n";
                    }
                }else if(probInsertion > prob){
                    //Insertion

                    int rValue = rand() % 4;
                    Dna rBase(rValue);
                    for(uint32_t h = length(tmpread) - 1; h > rpos; --h){
                        tmpread[h] = tmpread[h - 1];
                    }
                    tmpread[rpos] = rBase;
                    if(verbose){
                        std::cout << "Insertion: \n";
                        std::cout << "rBase " << rBase << "\n";
                    }
                }else{
                    //Deletion
                    if(verbose)
                        std::cout << "Deletion: \n";
                    for(uint32_t h = rpos; h < length(tmpread) - 1; ++h){
                        tmpread[h] = tmpread[h + 1];
                    }
                    tmpread[length(tmpread) - 1] = 'A';
                }
            }
            if(verbose){
                std::cout << "sampledReads: " << "\n" << sampledReads[k] << "\n";
                std::cout << tmpread << "\n";
            }
            appendValue(modifiedReads, tmpread);
        }
    }

    //save sampled Reads
    auto outputPathModified = outputPath;
    outputPathModified += "mod.fa";
    SeqFileOut seqFileout2(toCString(outputPathModified));
    for(uint32_t i = 0; i < length(modifiedReads); ++i){
        writeRecord(seqFileout2, to_string(i), modifiedReads[i]);
    }
    close(seqFileout2);

    return 0;
}



