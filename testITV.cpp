#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>

using namespace seqan;
using namespace std;

template<typename TString, typename TNeedle, typename TError>
inline bool verification(TString & text, TNeedle & needle, TError maxErrors)
{
    typedef Finder<TString>                                     TMyersFinder;
    typedef PatternState_<TNeedle, Myers<AlignTextBanded<FindInfix, NMatchesNone_, NMatchesNone_>, True, void> > TPatternState;
    TPatternState patternState;


    typedef PatternState_<TNeedle, Myers<AlignTextBanded<FindPrefix, NMatchesNone_, NMatchesNone_>, True, void> > TPatternStatePrefix;
    TPatternStatePrefix patternPrefixState;

    uint8_t minErrors = maxErrors + 1;

    std::cout << "Seqs: \n" << "reference: " << text << "\nneedle:    " << needle << "\n";

    TMyersFinder finder(text);
    uint8_t mErrors = maxErrors * 4;
    int startPos = 99;
    while (find(finder, needle, patternPrefixState, -static_cast<int>(maxErrors * 4)))
    {
        int currentEnd = position(finder) + 1;
        int currentErrors = -getScore(patternPrefixState);
        std::cout << "Loop errors: " << currentErrors << "\tselected end: " << currentEnd << "\n";
        if (currentErrors <= mErrors)
        {
            mErrors = currentErrors;
            startPos = currentEnd;
        }
    }
    std::cout << "BEST error: " << (int)mErrors << "\n";

    int score = globalAlignmentScore(needle, text, MyersBitVector());
    std::cout << "score global Alignment function :" << score << "\n";

}

int main(int argc, char const * argv[])
{

    // The bug only occurres when length(reference) == length(needle);

//     Dna5String reference = "TAAATATTATT";
//     Dna5String needle =    "TAAAATATTAT";

     Dna5String reference = "TAAATATTAT";
     Dna5String needle =    "TAAATAATTAT";

//     Dna5String reference = "TAAATAATTAT";
//     Dna5String needle =    "TAAATATTAT";


    verification(reference, needle, 3);

// This works
//     Dna5String needle =      "TATTATAAAAT";
//     Dna5String reference =   "TTATTATAAAT";

    return 0;
}

