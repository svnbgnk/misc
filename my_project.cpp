#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/seq_io.h>
#include "common.h"
#include "extension.h"

#include <iostream>
#include <seqan/align.h>
#include <typeinfo>
#include <set>

using namespace seqan;
using namespace std;

myGlobalParameters params;
int global;

/*
class myBamAlignmentRecord
{
public:
    CharString qName;               // QNAME
    uint16_t flag;                  // FLAG
    int32_t rID;                    // REF
    int32_t beginPos;               // POS
    uint8_t mapQ;                   // MAPQ mapping quality, 255 for invalid
    uint16_t bin;                   // bin for indexing
    String<CigarElement<> > cigar;  // CIGAR string
    int32_t rNextId;                // RNEXT (0-based)
    int32_t pNext;                  // PNEXT (0-based)
    int32_t tLen;                   // TLEN
    CharString seq;                 // SEQ, as in SAM/BAM file.
    CharString qual;                // Quality string as in SAM (Phred).
    CharString tags;                // Tags, raw as in BAM.

    // Constants for marking pos, reference id and length members invalid (== 0).
    static int32_t const INVALID_POS = -1;
    static int32_t const INVALID_REFID = -1;
    static int32_t const INVALID_LEN = 0;
};
*/

/*
template<typename T>
T * ptr(T * obj) { return obj; } //obj is already pointer, return it!

vector<int> & getElement(vector<vector<int>> & v, int e){
    return v[e];
}


template<typename TValue1, typename TValue2>
TValue1 & getFirst(std::vector<std::pair<TValue1, TValue2> *>& b, int e)
{
//     std::cout << b[e]->first << "\n";
    return b[e]->first;
}

template<typename TValue1, typename TValue2>
TValue1 & getFirst(std::vector<std::pair<TValue1, TValue2> >& b, int e)
{
    return b[e].first;
}

template<typename Tbit>
void my_f(std::vector<Tbit> & b){
    auto b1 = getFirst(b, 1);
    std::cout << "b1: " << b1 << "\n";
}

struct smallHit{
    Pair <uint32_t, uint32_t> occ;
//     uint8_t errors;
//     DnaString read;
};*/




int main(int argc, char const ** argv)
{

    typedef String<Dna5> TSeq;

    TSeq text = "AAAAAAAANNAAGC";
    TSeq read = "AAAAAANNAAAAGA";

    StringSet<TSeq> contigs;

    appendValue(contigs, text);
    appendValue(contigs, read);

    typedef typename InfixOnValue<StringSet<TSeq> const>::Type   TContigSeqsInfix;

    typedef typename InfixOnValue<TSeq const>::Type   TSeqsInfix;


    TContigSeqsInfix seq1 = infix(contigs[0], 1, length(contigs[0]) - 2);
    TContigSeqsInfix seq2 = infix(contigs[1], 1, length(contigs[1]) - 2);

    TSeqsInfix seq3 = infix(text, 1, length(text) - 2);
    TSeqsInfix seq4 = infix(read, 1, length(read) - 2);

    auto seq5 = infix(text, 1, length(text) - 2);
    auto seq6 = infix(read, 1, length(read) - 2);

    std::cout << seq1 << "\t" << seq2 << "\n";

//     int score = globalAlignmentScore(contigs[0], contigs[1], MyersBitVector());
//     int score = globalAlignmentScore(seq1, seq2, MyersBitVector());
//     int score = globalAlignmentScore(seq3, seq4, MyersBitVector());
    int score = globalAlignmentScore(seq5, seq6, MyersBitVector());
//     int score = globalAlignmentScore(contigs[0], contigs[1], MyersBitVector());
    cout << "score: " << score << endl;

/*
    StringSet<Dna5String> contig;
    Dna5String a = "ACGT";
    Dna5String c = "ACGTACGT";
    Dna5String b = "ACGTNACGTacgtACGT";

    std::cout << length(b) << "\t" << b << "\n";

    appendValue(contig, a);
//     appendValue(contig, b);
//     appendValue(contig, c);

    std::vector<uint32_t> sl;
    sl.push_back(0);
    for(uint32_t i = 0; i < countSequences(contig); ++i){
        sl.push_back(seqan::length(contig[i]));
    }

    std::cout << "SL: \t";
    for(int i = 0; i < sl.size(); ++i)
        std::cout << sl[i] << "\n";

    auto mylimits = stringSetLimits(contig);
    std::cout << "\nlim: \n";
    for(int i = 0; i < length(mylimits); ++i){
        std::cout <<  mylimits[i] << "\t";
    }*/
/*
    std::pair<int, int> a (5,7);
//     std::pair<int, int> * pa = &a;
//     std::cout << pa->first << "\n";
//     std::cout << a.first << "\n";
    std::vector<std::pair<int, int> > b;
    std::vector<std::pair<int, int> *> p;
    b.push_back(a);
    b.push_back(a);
    p.push_back(&a);
    p.push_back(&a);

    my_f(b);
    my_f(p);*/


/*
    vector<int> a(500000000, 1);
    std::cout << "finished vector0" << "\n";
    vector<vector<int>> v;
    v.push_back(a);
    v.push_back(a);
    std::cout << "finished vector" << "\n";

    vector<int> & b = getElement(v, 1);
    std::cout << "accessing vector" << "\n";
    std::cout << "f " << b[0] << "\n";

    std::cout << "f " << getElement(v, 1)[1] << "\n";
//     std::cout << b << "\n";*/



/*
    text = "AGTCGGATCTACTG";
    read =   "TCGCAACTAC";
    score = globalAlignmentScore(text, read, MyersBitVector()) + 4;
    cout << "score: " << score << endl;

    text = "AGTCGGATCTACTG";
    read =   "TCGCATCTGC";
    score = globalAlignmentScore(text, read, MyersBitVector()) + 4;
    cout << "score: " << score << endl;*/
/*
    cout << "Test" << endl;
    int max_e = 2; //max number of Insertion + Deletions
    for(int e = 0; e <= max_e; ++e){
        cout << "E: " << e << endl;;
        for(int i = 0; i <= e; ++i){
            //i is number of insertions
            int d = e - i //number of deletions
            cout << "Number of insertions: " << i << endl;
            for(unsigned pos = 0; pos < 2; ++pos){
                if(pos == 0){
                    auto const & tmp = infix(ex_needle, 2 - i, length(needle) + 2 - d);
                }else{
                    auto const & tmp = infix(ex_needle, 2 + d, length(needle) + 2 + i);
                }
                uint8_t tmpscore = 0 - globalAlignmentScore(tmp, needle, MyersBitVector());
                if(tmpscore <= max_e)
                    delegateDirect(sa_info, needle, score);

            }
        }
        cout << endl;
    }*/

//             for(unsigned positions = 0; positions <= ndistributions; ++positions){
//                 bitset<4> bitset1{types};
//                 bitset<4> bitset2{positions};
//                 cout << bitset1 << endl;
//                 cout << bitset2 << endl;
//                 cout << endl;
//             }
//
    unsigned short1 = 4;
    bitset<16> bitset1{short1};   // the bitset representation of 4
    cout << bitset1 << endl;
    short1 = (short1 << 1);
    bitset<16>bitset2{short1};
    cout << bitset2 << endl;
/*
    String<int> locations;
    for (unsigned i = 0; i < length(text) - length(pattern); ++i)
    {
        // Compute the MyersBitVector in current window of text.
        TSequence tmp = infix(text, i, i + length(pattern));

        // Report hits with at most 2 errors.
        if (globalAlignmentScore(tmp, pattern, MyersBitVector()) >= -2)
        {
            appendValue(locations, i);
        }
    }
*/



    /*
    {
    typedef String<char> TSequence;                             // sequence type
    typedef StringSet<TSequence> TStringSet;                    // container for strings
    typedef StringSet<TSequence, Dependent<> > TDepStringSet;   // dependent string set
    typedef Graph<Alignment<TDepStringSet> > TAlignGraph;       // alignment graph

    TSequence seq1 = "blablubalu";
    TSequence seq2 = "abba";

    TStringSet sequences;
    appendValue(sequences, seq1);
    appendValue(sequences, seq2);

    TAlignGraph alignG(sequences);

    int score = globalAlignment(alignG, Score<int, Simple>(1, -1, -1), AlignConfig<true, true, true, true>(), LinearGaps());
    std::cout << "Score: " << score << std::endl;
    std::cout << alignG << std::endl;
    }*/

/*
    typedef String<Dna> TSequence;                             // sequence type                  // container for strings
    typedef Align<TSequence, ArrayGaps> TAlign;

    TSequence seq1 = "AAATGACGGATTG";
    TSequence seq2 = "AGTCGGATCTACTG";

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    int score = globalAlignment(align, Score<int, Simple>(4, -2, -2, -4), AffineGaps());
    std::cout << "Score: " << score << std::endl;
    std::cout << align << std::endl;
*/

/*
    typedef String<Dna> TSequence;                 // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign;     // align type

    TSequence seq1 = "ACCGGTTGCAAGT";
    TSequence seq2 = "ACCGCCGCATTGT";

//     TAlign align;
//     resize(rows(align), 2);
//     assignSource(row(align, 0), seq1);
//     assignSource(row(align, 1), seq2);
    int score = globalAlignmentScore(seq1, seq2, MyersBitVector());
    std::cout << "Score: " << score << std::endl;
//     std::cout << align << std::endl;
*/
    /*

    CharString samFileOutName = CharString("/home/sven/devel/my_project-build/ex.sam");

    myGlobalParameters myobject;
    params.flipdensity = 0.9;
    myobject.flipdensity = 0.8;
    cout << params.flipdensity << endl;
    cout << "Use Function: " << endl;
//     my::printAllP(params);
    params.print();
    cout << "Try function" << endl;
    cout << my::calcsomething(0.5, 0.5) << endl;
    global = 10;
    seqan::testglobal();
*/
//     BamFileOut samFileOut(std::cout, Sam());


//
//     typedef typename BamHeaderRecord::TTag    TTag;
/*
@HD	VN:1.3	SO:coordinate
@SQ	SN:ref	LN:45
@SQ	SN:ref2	LN:40
*/
 /*
    BamFileOut samFileOut(toCString(samFileOutName));
    BamHeader header;
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, BamHeaderRecord::TTag("VN", "3.0"));
    appendValue(firstRecord.tags, BamHeaderRecord::TTag("SO", "coordinate"));
    appendValue(header, firstRecord);

    BamHeaderRecord seqRecord;
    seqRecord.type = BAM_HEADER_REFERENCE;
    appendValue(seqRecord.tags, BamHeaderRecord::TTag("SN", "ref"));
    appendValue(seqRecord.tags, BamHeaderRecord::TTag("LN", "999"));
    appendValue(header, seqRecord);


    BamHeaderRecord thirdRecord;
    thirdRecord.type = BAM_HEADER_REFERENCE;
    appendValue(thirdRecord.tags, BamHeaderRecord::TTag("SN", "ref2"));
    appendValue(thirdRecord.tags, BamHeaderRecord::TTag("LN", "999"));
    appendValue(header, thirdRecord);


    writeHeader(samFileOut, header);





    vector<Pair<DnaString, Pair<uint32_t, uint32_t> > > occ;
    ;
    Pair<uint32_t, uint32_t> cord = Pair <uint32_t, uint32_t>(0, 32);
    occ.push_back(Pair < DnaString, Pair <uint32_t, uint32_t> > (DnaString("AAACATGCCCTCCCTTCCCCTGTGAGATTACTCGTGGGCAGGGCTTGTCCCTAGGGGGGCCCTGTCATTCAGTGGTTAACCCCGTGGTCCCAGGGTTTAC"), cord));
    cord = Pair <uint32_t, uint32_t>(0, 12);
    occ.push_back(Pair < DnaString, Pair <uint32_t, uint32_t> > (DnaString("TGTGCTAAGAATGTATTACTCTCTTTAGCAAATAACAGCCATAAATTTTGTTTCCAAAAATCTTTATCTTACTTGACAAGTGGGTGGAGACTGCTGAGCA"), cord));
    cord = Pair <uint32_t, uint32_t>(1, 62);
    occ.push_back(Pair < DnaString, Pair <uint32_t, uint32_t> > (DnaString("CCATAGTAGGTGCTCGGTAAATGTCAAAGATTTCCCAAATATTTGTATTGCTGCCAAGCTGAAGCATTTAGAAATGCATCAGTTTCATCTCCTCTTTGTG"), cord));
    cord = Pair <uint32_t, uint32_t>(1, 92);
    occ.push_back(Pair < DnaString, Pair <uint32_t, uint32_t> > (DnaString("AAACATGCCCTCCCTTCCCCTGTGAGATTACTCGTGGGCAGGGCTTGTCCCTAGGGGGGCCCTGTCATTCAGTGGTTAACCCCGTGGTCCCAGGGTTTAC"), cord));

    StringSet<CharString> contigNameStore;
    CharString id = "ref";
    appendValue(contigNameStore, id);
    id = "ref2";
    appendValue(contigNameStore, id);
    NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);



// SeqFileIn genomeFileIn(toCString(genomePath));
// while (!atEnd(genomeFileIn))
//     {
//         CharString id;
//         Dna5String genome;
//         readRecord(id, genome, genomeFileIn);
//         appendValue(contigNameStore, id);
//         // throw away seq
// }

    for(int i = 0; i < 4; ++i){

        BamAlignmentRecord record;


        record.qName = CharString("name" + to_string(i));
        record.rID = 0; //Segementaion fault //occ[i].i2.i1
        record.beginPos = occ[i].i2.i2;
        record.seq = occ[i].i1;

        record.mapQ = 255;
        record.rNextId = BamAlignmentRecord::INVALID_REFID;
        record.pNext = BamAlignmentRecord::INVALID_POS;
        record.tLen = BamAlignmentRecord::INVALID_LEN;
        record.qual = "IIIIIIIIIIIIIIIIIIIIIII";


        BamTagsDict tagsDict(record.tags);
        //getCigarString(...)

        setTagValue(tagsDict, "NM", 2);
        // => tags: "NM:i:2"
        setTagValue(tagsDict, "NH", 1);
        // => tags: "NM:i:2 NH:i:1"
        setTagValue(tagsDict, "NM", 3);
        // => tags: "NM:i:3 NH:i:1"



        try{
        writeRecord(samFileOut, record);
        }
        catch (ParseError const & e)
        {
            std::cerr << "ERROR: input header is badly formatted. " << e.what() << "\n";
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
        }
    }
     */

/*

    CharString bamFileInName = getAbsolutePath("demos/tutorial/sam_and_bam_io/example.sam");

    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamFileInName)))
    {
        std::cerr << "ERROR: could not open input file " << bamFileInName << ".\n";
        return 1;
    }

//     BamFileOut samFileOut(context(bamFileIn), std::cout, Sam());

    BamHeader header;
    readHeader(header, bamFileIn);
    BamAlignmentRecord record;
    uint32_t count = 0;
    while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            cout << record.seq << endl;
            count += static_cast<int>(hasFlagUnmapped(record));
//             cout << std::string(record.cigar) << endl;

//             writeRecord(samFileOut, record);
        }
    cout << "Count: " << count << endl;
  */


/*
    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

    TBamContext const & bamContext = context(bamFileIn);

    for(int i = 0; i< seqan::length(contigNames(bamContext)); ++i){
        cout << contigNames(bamContext)[i] << "\t" << contigLengths(bamContext)[i] << "\n";
    }*/



    /*
    try
    {
        BamHeader header;
        BamAlignmentRecord record;
        readHeader(header, bamFileIn);
        writeHeader(samFileOut, header);
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            writeRecord(samFileOut, record);
        }
    }
    catch (ParseError const & e)
    {
        std::cerr << "ERROR: input header is badly formatted. " << e.what() << "\n";
    }
    catch (IOError const & e)
    {
        std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
    }
    */


    /*
    CharString bamFileInName = getAbsolutePath("demos/tutorial/file_io_overview/example.bam");
    CharString samFileOutName = getAbsolutePath("demos/tutorial/file_io_overview/example.sam");

    // Open input BAM file.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamFileInName)))
    {
        std::cerr << "ERROR: could not open input file " << bamFileInName << ".\n";
        return 1;
    }

    // Open output SAM file.
//     BamFileOut samFileOut(context(bamFileIn), toCString(samFileOutName));
    BamFileOut samFileOut(context(bamFileIn), std::cout, Sam());

    // Copy header.
    BamHeader header;
    try
    {
        readHeader(header, bamFileIn);
        writeHeader(samFileOut, header);
    }
    catch (ParseError const & e)
    {
        std::cerr << "ERROR: input header is badly formatted. " << e.what() << "\n";
    }
    catch (IOError const & e)
    {
        std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
    }

    // Copy all records.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        try
        {
            readRecord(record, bamFileIn);
            writeRecord(samFileOut, record);
        }
        catch (ParseError const & e)
        {
            std::cerr << "ERROR: input record is badly formatted. " << e.what() << "\n";
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }*/

    return 0;
}

/*
Errors: 1   < AAACATGCCCTCCCTTCCCCTGTGAGATTACTCGTGGGCAGGGCTTGTCCCTAGGGGGGCCCTGTCATTCAGTGGTTAACCCCGTGGTCGCAGGGTTTAC , < 0 , 97437165 > >
AAACATGCCCTCCCTTCCCCTGTGAGATTACTCGTGGGCAGGGCTTGTCCCTAGGGGGGCCCTGTCATTCAGTGGTTAACCCCGTGGTCCCAGGGTTTAC
Errors: 1   < TGTGCTAAGAATGTATTAATCTCTTTAGCAAATAACAGCCATAAATTTTGTTTCCAAAAATCTTTATCTTACTTGACAAGTGGGTGGAGACTGCTGAGCA , < 0 , 97482311 > >
TGTGCTAAGAATGTATTACTCTCTTTAGCAAATAACAGCCATAAATTTTGTTTCCAAAAATCTTTATCTTACTTGACAAGTGGGTGGAGACTGCTGAGCA
Errors: 0   < CCATAGTAGGTGCTCGGTAAATGTCAAAGATTTCCCAAATATTTGTATTGCTGCCAAGCTGAAGCATTTAGAAATGCATCAGTTTCATCTCCTCTTTGTG , < 0 , 97547742 > >
CCATAGTAGGTGCTCGGTAAATGTCAAAGATTTCCCAAATATTTGTATTGCTGCCAAGCTGAAGCATTTAGAAATGCATCAGTTTCATCTCCTCTTTGTG
Errors: 2   < GTGCACATGCTCACACACATCCTCACACACATCCTTACACACCCTCACCCACATGCACTCACACACATGCACACACACTCCCTCACTCATGCACACATAC , < 0 , 97770616 > >
GTGCACATGCTCACACACATCCTCACACACTTCCTCACACACCCTCACCCACATGCACTCACACACATGCACACACACTCCCTCACTCATGCACACATAC*/

    /*

    */


