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



template<typename TString, typename TNeedle, typename TError>
inline bool inTextVerificationE(StringSet<TString> const & reference, StringSet<TNeedle> const & needles, TError maxErrors)
{/*
    typedef MapperTraits<TSpec, TConfig>                        TTraits;
    typedef typename TConfig::TContigsLen                       TContigsLen;
    typedef typename TConfig::TContigsSize                      TContigsSize;

    typedef typename TTraits::TContigSeqs                       TContigSeqs;
    typedef typename StringSetPosition<TContigSeqs const>::Type TContigPos;
    typedef typename InfixOnValue<TContigSeqs const>::Type      TContigSeqsInfix;
    typedef typename InfixOnValue<TNeedle const>::Type          TNeedleInfix;
    typedef ModifiedString<TNeedle, ModReverse>                 TNeedleInfixRev;
    typedef ModifiedString<TContigSeqsInfix, ModReverse>        TStringInfixRev;

    TContigSeqs & contigSeqs = me.contigs.seqs;

    TContigsSize seqNo = getMember(match, ContigId());
    TContigsLen seqOffset = getMember(match, ContigBegin());
    TContigsLen seqOffsetEnd = getMember(match, ContigEnd());
    TContigSeqsInfix text = infix(contigSeqs[seqNo], seqOffset, seqOffsetEnd);*/

    typedef uint32_t                       TContigsLen;
    TString text = reference[0];
    TNeedle needle = needles[0];

/*
    typedef Finder<TString>                        TFinder;
//     typedef Finder<TStringInfixRev>                       TFinder2;
    typedef AlignTextBanded<FindInfix,
                            NMatchesNone_,
                            NMatchesNone_>                     TMyersSpecInfix;
    typedef Myers<TMyersSpecInfix, True, void>                 TAlgorithmInfix;
    typedef PatternState_<TNeedle, TAlgorithmInfix>   TPatternInfix;
    TPatternInfix patternInfix;
    */

    typedef Finder<TString>                                     TMyersFinder;
    typedef PatternState_<TNeedle, Myers<AlignTextBanded<FindInfix, NMatchesNone_, NMatchesNone_>, True, void> > TPatternState;
    TPatternState patternState;


    typedef PatternState_<TNeedle, Myers<AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>, True, void> > TPatternStatePrefix;
    TPatternStatePrefix patternPrefixState;

    typedef Finder<ModifiedString<TString, ModReverse>>                                 TMyersRevFinder;
    typedef PatternState_<ModifiedString<TNeedle, ModReverse>, Myers<AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>, True, void> > TPatternStatePrefixRev;
    TPatternStatePrefixRev patternStatePrefixRev;

    uint8_t minErrors = maxErrors + 1;
    TMyersFinder finderInfix(text);

    std::cout << "ref: " << text << "\n";
    std::cout << "ned: " << needle << "\n";

    while (find(finderInfix, needle, patternState, -static_cast<int>(length(needle))))
    {
        auto currentErrors = -getScore(patternState);
        std::cout << "Loop E: " << currentErrors << "\n";
//         std::cout << currentErrors << "\n";
        if(minErrors > currentErrors)
            minErrors = currentErrors;
    }


    std::cout << "minErrors " << (int)minErrors << "\n";



    if(length(text) <= length(needle)){
        Gaps<TString> textGaps(text);
        Gaps<TNeedle> needleGaps(needle);
        Score<int, Simple> scoringScheme(0, -999, -1000);
        int score = globalAlignment(textGaps, needleGaps, scoringScheme, AlignConfig<true, false, false, true>()) / -999;
        clipSemiGlobal(textGaps, needleGaps);
        std::cout << textGaps << "\n" << needleGaps << "\n";
        std::cout << beginPosition(textGaps) << "\n" << endPosition(textGaps) << "\n";
        std::cout << beginPosition(needleGaps) << "\n" << endPosition(needleGaps) << "\n";
        std::cout << "Score " << (int)score << "\n";
        return score < maxErrors;
    }

//     if(minErrors <= maxErrors)
//         match.errors = minErrors;

    TMyersFinder finder(text);
    uint8_t mErrors = maxErrors * 4;
    TContigsLen endPos = 0;
    while (find(finder, needle, patternPrefixState, -static_cast<int>(maxErrors * 4)))
    {
        int currentEnd = position(finder) + 1;
        uint16_t currentErrors = -getScore(patternPrefixState);
//         if (getValue(text, currentEnd) != back(needle))
//             ++currentErrors;
        std::cout << "Loop: " << currentErrors << "\t" << currentEnd << "\n";
        if (currentErrors <= mErrors)
        {
            mErrors = currentErrors;
            endPos = currentEnd;
        }
    }
    std::cout << "BEst: " << (int)mErrors << "\t" << endPos << "\n";

//     TString infixPrefix = infix(text, 0, endPos);
//     std::cout << "Prefix: " << infixPrefix << "\n";

    ModifiedString<TString, ModReverse> infixRev(text/*infixPrefix*/);
    ModifiedString<TNeedle, ModReverse> needleRev(needle);

    std::cout << "revPrefices: \n" << infixRev << "\n" << needleRev << "\n";


    TMyersRevFinder finder2(infixRev);

    mErrors = maxErrors * 4;
    TContigsLen startPos = endPos;

    while (find(finder2, needleRev, patternStatePrefixRev, -static_cast<int>(maxErrors * 4)))
    {
        int currentEnd = position(finder2) + 1;
        uint16_t currentErrors = -getScore(patternStatePrefixRev);
        std::cout << "Loop rev: " << currentErrors << "\t" << currentEnd << "\n";
        if (currentErrors <= mErrors)
        {
            mErrors = currentErrors;
            startPos = currentEnd;
        }
    }
    if(mErrors >= maxErrors * 4)
    {
        std::cout << "no alignment between refPrefixes smaller maxErrors * 4 was found \n";
    }

    std::cout << "BEst: " << (int)mErrors << "\t" << endPos << "\n";


    std::cout << "start: " << (length(infixRev) - startPos) << "\t endPos: " << endPos << "\t Errors: " << (int)minErrors << "\n";

//     std::cout << "revSuffix: " << (endPos - startPos) << "\t endPos: " << endPos << "\t Errors: " << (int)minErrors << "\n";

    std::cout << infix(text, (length(infixRev) - startPos), endPos) << "\n" << needle << "\n";


}

int main(int argc, char const * argv[])
{


        StringSet<Dna5String> reference;
        StringSet<Dna5String> needles;

// myers fails third alignment (can be fixed by appending random bases to reference at the beginning)
    Dna5String needle = "CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACGCTAACCCTAACCCTAACCCTAACCCCAACCCTA";
    appendValue(needles, needle);
    appendValue(reference, "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC");


// myers fails in third alignment (can be fixed by appending random bases to reference at the beginning)
/*
    Dna5String needle = "TTTGTGTGTCTATCTGTATTTTTGGTCTAAATATTTGGCCCAACAGAGGTTAAAGGCTTTGATGTTCTCAGCAAAAGCCTTGTGAGATCTCTAGTTTATT";
    appendValue(needles, needle);
    appendValue(reference, "TTTGTGTGTCTATCTGTATTTTTGGTCTAAAATATTTGGCCCAACAGAGGTTAAAGGCTTTGATGGTCTCAGCAAAAGCCTTGTGAGATCTCTAGTTTAT");*/


/*
// myers does not fail
    Dna5String needle = "CATAAGTTGCATTACTTCAGCGTCCCAACTGCACCCTTACCACGAAGACAGGTTTGTCCATTCCCATACTGCGGCGTTGGCAGGGGGTTCGCATGTCCCA";
    appendValue(needles, needle);
    appendValue(reference, "ACCATAAGTTGCATTACTTCAGCGTCCCAACTGCACCCTTACCACGAAGACAGGTTTGTCCATTCCCATACTGCGGCGTTGGCAGGGGGTTCGCATGTCCCATA");
*/



    std::cout << "Call: \n";
    inTextVerificationE(reference, needles, 3);

    return 0;
}



