#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include <numeric>

#include <iostream>
#include <seqan/align.h>
#include <typeinfo>
#include <set>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

using namespace seqan;
using namespace std;


template <typename TChar, typename TConfig>
void generateText(String<TChar, TConfig> & text, unsigned const length)
{
    unsigned alphabetSize = ValueSize<TChar>::VALUE;
    unsigned minChar = MinValue<TChar>::VALUE;

    resize(text, length);

    for (unsigned i = 0; i < length; ++i)
    {
        text[i] = std::rand() % alphabetSize - minChar;
    }
}

template <typename TChar, typename TConfig>
void mutate(String<TChar, TConfig> & read, int const errors, bool const editMuations){
    bool verbose = false;
    std::vector<int> v;
    float probInsertion = 0.15;
    float probDeletion = 0.15;
    v.reserve(errors);
    if(verbose)
        std::cout << read << "\n";
    for(int i = 0; i < errors; ++i){
        uint32_t rpos = rand() % length(read);
        while(std::find(v.begin(), v.end(), rpos) != v.end()){
            rpos = rand() % length(read);
        }
        v.push_back(rpos);
        if(verbose)
            std::cout << rpos << "\n";
        float prob = (float) rand()/RAND_MAX;

        //Mismatch
        if(!editMuations || probInsertion + probDeletion < prob){
            int rValue = rand() % 3;
            int cValue = read[rpos].value;
            if(rValue >=  cValue)
                ++rValue;
            Dna rBase(rValue);
            read[rpos] = rBase;
            if(verbose){
                std::cout << "Mismatch: \n";
                std::cout << "rBase " << rBase << "\n";
            }
        }

        else if(probInsertion > prob)
        {
            //Insertion
            int rValue = rand() % 4;
            Dna rBase(rValue);
            for(uint32_t h = length(read) - 1; h > rpos; --h){
                read[h] = read[h - 1];
            }
            read[rpos] = rBase;
            if(verbose){
                std::cout << "Insertion: \n";
                std::cout << "rBase " << rBase << "\n";
            }
        }
        else
        {
            //Deletion
            if(verbose)
                std::cout << "Deletion: \n";
                for(uint32_t h = rpos; h < length(read) - 1; ++h){
                    read[h] = read[h + 1];
                }
                read[length(read) - 1] = 'A';
        }
    }
    if(verbose)
        std::cout << read << "\n";
}

template <typename TDelegate,
          typename TText, typename TIndexSpec,
          typename TNeedle,
          typename TDistanceTag>
inline void find(const int minErrors,
     const int maxErrors,
     TDelegate & delegate,
     Index<TText, BidirectionalIndex<TIndexSpec> > & index,
     TNeedle const & needle,
     TDistanceTag const & )
{
    switch (maxErrors)
    {
        case 1: find<0, 1>(delegate, index, needle, TDistanceTag());
                break;
        case 2: find<0, 2>(delegate, index, needle, TDistanceTag());
                break;
        case 3: find<0, 3>(delegate, index, needle, TDistanceTag());
                break;
        case 4: find<0, 4>(delegate, index, needle, TDistanceTag());
                break;
        default: std::cerr << "E = " << maxErrors << " not yet supported.\n";
                exit(1);
    }
}

template<typename TVector>
void printv(TVector & a){
    for(int i = 0; i < a.size(); ++i){
        std::cout << a[i] << ", ";
    }
    std::cout << "\n";
}

template<typename TPair>
void printPair(TPair & a){
    std::cout << "<" << a.i1 << ", " << a.i2 << ">";
}

struct SARange
{
    Pair<uint32_t, uint32_t> range;
    uint32_t repLength;
    uint8_t errors;
    uint32_t shift = 0;
};

struct Hit
{
    uint32_t occ;
    uint32_t occEnd;
    uint8_t errors;
    uint32_t needleId;
};

void print(Hit & a)
{
    std::cout << a.needleId << " @ <" << a.occ << ", " << a.occEnd << ">" << "\t + " << (int)a.errors << "\n";
}

int main(int argc, char const * argv[])
{
    ArgumentParser parser("Test");

    addOption(parser, ArgParseOption("l", "length", "length of the read", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("e", "errors", "Number of allowed Errors", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("m", "mutations", "Number of mutations done on the read", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("i", "iterations", "Number of iterations testing search", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("HD", "hammingDistance", ""));

    addOption(parser, ArgParseOption("em", "editMuations", "Allow Insertions and Deletions as mutations"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString referencePath, readsPath, outputPath;
    int length = 100, errors = 2, mutations = 0, iterations = 500;

//     getOptionValue(referencePath, parser, "reference");
    getOptionValue(length, parser, "length");
    getOptionValue(errors, parser, "errors");
    getOptionValue(mutations, parser, "mutations");
    getOptionValue(iterations, parser, "iterations");
    bool hammingDistance = isSet(parser, "hammingDistance");
    bool editMuations = isSet(parser, "editMuations");
    bool verbose = isSet(parser, "verbose");

    time_t seed = std::time(nullptr);
    std::srand(seed);

    typedef DnaString TText;
    typedef FastFMIndexConfig<void, uint32_t, 2, 1> TMyFastConfig;
    typedef Index<TText, BidirectionalIndex<FMIndex<void, TMyFastConfig> > > TIndex;
    typedef Iter<TIndex, VSTree<TopDown<> > > TIter;

    std::vector<std::vector<int> > hitCounters;

    hitCounters.resize(errors + 1);

    for(int i = 0; i < iterations; ++i){

        TText text;
        generateText(text, length + 2 * errors);
        TText read = infix(text, errors, length + errors);
        if(mutations > 0)
            mutate(read, mutations, editMuations);

        if(verbose)
            std::cout << text << "\nRead: " << read << "\n";
        TIndex index(text);
        TIter it(index);
        std::vector<int> hitCounter(errors + 1, 0);
        std::vector<Hit> hits;
        std:vector<SARange> ranges;
        /*
        auto delegate = [&hitCounter, &ranges](auto & iter, DnaString const & needle, uint8_t const cerrors)
        {
//             Pair<int8_t, int8_t> range(iter.fwdIter.vDesc.range.i1, iter.fwdIter.vDesc.range.i2);
            SARange range;
            range.errors = cerrors;
            range.range = Pair<int8_t, int8_t> (iter.fwdIter.vDesc.range.i1, iter.fwdIter.vDesc.range.i2);
            ranges.push_back(range);
            for (auto occ : getOccurrences(iter)){
                ++hitCounter[cerrors];
            }
        };*/


        auto delegate = [&hitCounter, &hits](auto & index, auto & saRange, uint32_t const & needleId)
        {
//             int shift = saRange.shift;
            for (uint32_t i = saRange.range.i1; i < saRange.range.i2; ++i){
                Hit hit;
                hit.occ = index.fwd.sa[i]; // + shift;
                hit.occEnd = hit.occ + saRange.repLength;
                hit.errors = saRange.errors;
                hit.needleId = needleId;
                hits.push_back(hit);
                ++hitCounter[saRange.errors];
            }
        };

        auto delegateRange = [&ranges](auto & iter, DnaString const & needle, uint8_t const cerrors)
        {
//             Pair<int8_t, int8_t> range(iter.fwdIter.vDesc.range.i1, iter.fwdIter.vDesc.range.i2);
            SARange range;
            range.repLength = repLength(iter);
            range.errors = cerrors;
            range.range = Pair<int8_t, int8_t> (iter.fwdIter.vDesc.range.i1, iter.fwdIter.vDesc.range.i2);
            ranges.push_back(range);
        };



        if(hammingDistance)
            find(0, errors, delegateRange, index, read, HammingDistance());
        else
            find(0, errors, delegateRange, index, read, EditDistance());
        if(verbose)
            printv(hitCounter);
        std::sort(ranges.begin(), ranges.end(), [](SARange & x, SARange & y){
            if(x.range.i1 == y.range.i1)
                return(x.errors < y.errors);
            else
                return(x.range.i1 < y.range.i1);
        });

        if(verbose){
            for(int l = 0; l < ranges.size(); ++l){
                auto cRange = ranges[l];
                printPair(cRange.range);
                std::cout << "\t" << (int)cRange.errors << "\n";
            }
        }

        //locate //call global delegate
        std::vector <uint32_t> lastEnds(errors + 1, 0);
        for(int l = 0; l < ranges.size(); ++l){
            auto cRange = ranges[l];
            uint8_t cE = cRange.errors;
            if(lastEnds[cE] < cRange.range.i1)
            {
                delegate(index, cRange, 0);
                lastEnds[cE] = cRange.range.i2;
                // <2,5> 2 <3,8> 1 <6, 7> 2 miss reporting last Interval with 2 errors
                //but maybe we are sure a range with a smaller error has to smaller than range with larger error
//                 for(int e = cE; e < errors + 1; ++e)
//                     lastEnds[e] = cRange.range.i2;
            }
            else if(lastEnds[cE] < cRange.range.i2)
            {
                cRange.range.i1 = lastEnds[cE];
                delegate(index, cRange, 0);
                lastEnds[cE] = cRange.range.i2;
                //is this one the same case??
//                 for(int e = cE; e < errors + 1; ++e)
//                     lastEnds[e] = cRange.range.i2;
            }
        }


        if(verbose){
            for(int j = 0; j < hits.size(); ++j)
                print(hits[j]);
        }
        for(int j = 0; j < hitCounter.size(); ++j){
            hitCounters[j].push_back(hitCounter[j]);
        }

    }


    cout << "Number of Iterations: " << iterations << "\n";
    for(int i = 0; i < errors + 1; ++i){
        cout << "Number of Alignments  with " << i << " errors:\n";
        double average = accumulate( hitCounters[i].begin(), hitCounters[i].end(), 0.0)/hitCounters[i].size();
        cout << "Average \t\t\t\t\t" << average << endl;

        std::vector<double> diff(hitCounters[i].size());
        std::transform(hitCounters[i].begin(), hitCounters[i].end(), diff.begin(), [average](double x) { return x - average; });
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / hitCounters[i].size());
        cout << "Standart diviation \t\t\t\t" << stdev << endl;
    }


    std::cout << "Finished!\n";

    return 0;
}



