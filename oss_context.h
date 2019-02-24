#ifndef OSSBESTX_H_
#define OSSBESTX_H_

#include "common.h"

using namespace seqan;
struct modusParameters{
public:
    bool nomappability;
    bool directsearch;
    bool compmappable;
    bool suspectunidirectional;

    bool testflipdensity;
    uint32_t step;
    uint32_t distancetoblockend;
    uint32_t directsearch_th;
    uint32_t directsearchblockoffset;
    float filter_th;
    float invflipdensity;
    uint32_t intervalsize;

    modusParameters(){
        setdefault();
    }

    void setdefault(){
        nomappability = true;
        directsearch = true;
        compmappable = true;
        suspectunidirectional = true;

        testflipdensity = true;
        //binaryNumber //has to be 2^x - 1 for fast modulo calculation
        step = 0b11;
        distancetoblockend = 2;

        directsearchblockoffset = 0;
        directsearch_th = 2;
        filter_th = 0.5;

        invflipdensity = 0.5;

        intervalsize = 3;
    }

    void print(){
        std::cout << "Cases Enabled: " << "\n";
        std::cout << nomappability << " " << directsearch << " " << compmappable << " " << suspectunidirectional << "\n";
        std::cout << "Params: " << "\n";

        std::cout << "step: " << step << "\n";
        std::cout << "distancetoblockend: " << distancetoblockend << "\n";
        std::cout << "directsearchblockoffset: " << directsearchblockoffset << "\n";
        std::cout << "directsearch_th: " << directsearch_th << "\n";
        std::cout << "filter_th: " << filter_th << "\n";
        std::cout << "invflipdensity: " << invflipdensity << "\n";
        std::cout << "intervalsize: " << intervalsize << "\n";
    }
};

// template <typename TTraits>
class OSSContext
{
public:
    //Parameters
    modusParameters normal;
    modusParameters comp;
    modusParameters uni;


    bool bestXMapper = false; //still needed multiple searches
    bool oneSSBestXMapper = false;
/*
    // Shared-memory read-write data.
    ReadsContext<void, void>       ctx;
    StringSet<DnaString> & reads;
    std::vector<hit> & hits;
    std::vector<hit> & dhits;
    std::vector<uint32_t>  readOccCount;*/

    // Shared-memory read-only data.
//     TContigSeqs const & contigSeqs;
    bool filterDelegate = true;
    bool trackReadCount = false;
    bool itv = true;
    uint32_t itvOccThreshold = 10;

//     uint32_t readCount = length(reads); //if readCount is not set than all reads are assumend to be on one strand
    uint8_t maxError = 0;
    uint8_t strata = 99;
    uint32_t readLength;
    uint32_t numberOfSequences;
//     std::vector<std::pair<int, bool> > bitvectorMeta;

/*
    OSSContext(StringSet<DnaString> inreads,
               std::vector<hit> & inhits,
               std::vector<hit> & indhits) :
        reads(inreads),
        hits(inhits),
        dhits(indhits)
//         readOccCount(inreadOccCount)
    {}*/

    void setdefault(){
        normal.setdefault();
        comp.setdefault();
        uni.setdefault();
    }

    void print(){
        std::cout << "Normal: ";
        normal.print();
        std::cout << "Comp: ";
        comp.print();
        std::cout << "Uni: ";
        uni.print();
    }

    void loadInputParameters(uint8_t inMaxError, uint8_t inStrata, uint32_t inReadLength, uint32_t inNumberOfSequences){
        maxError = inMaxError;
        strata = inStrata;
        readLength = inReadLength;
        numberOfSequences = inNumberOfSequences;
    }

    template <size_t nbrBlocks, typename TSALength>
    bool itvCondition(OptimalSearch<nbrBlocks> const & s,
                      uint8_t const blockIndex,
                      TSALength ivalOne)
    {
        return(itv && ivalOne < itvOccThreshold/*(static_cast<int>(s.pi.size()) - blockIndex - 1 + normal.directsearchblockoffset) * normal.directsearch_th*/);
    }


    template <typename TIter,
              size_t nbrBlocks>
    bool itvConditionComp(TIter iter,
                          uint32_t const needleLeftPos,
                          uint32_t const needleRightPos,
                          uint8_t const errors,
                          OptimalSearch<nbrBlocks> const & s,
                          uint8_t const blockIndex)
    {
        return(itv && iter.fwdIter.vDesc.range.i2 - iter.fwdIter.vDesc.range.i1 < itvOccThreshold/*< (s.pi.size() - blockIndex - 1 + comp.directsearchblockoffset) * comp.directsearch_th*/);
    }

    template <typename TSALength>
    bool itvConditionUni(uint8_t const blockSize,
                         uint8_t const blockIndex,
                         TSALength ivalOne)
    {
        return(itv && ivalOne < itvOccThreshold/*(static_cast<int>(blockSize) - blockIndex - 1 + uni.directsearchblockoffset) * uni.directsearch_th*/);
    }

    template<size_t nbrBlocks>
    bool inBlockCheckMappabilityCondition(uint32_t needleLeftPos,
                                          uint32_t needleRightPos,
                                           OptimalSearch<nbrBlocks> const & s,
                                          uint8_t blockIndex)
    {
        uint32_t step = (needleRightPos - needleLeftPos - 1);
        if(normal.distancetoblockend > step)
            return false;
        uint32_t prevBlocklength = (blockIndex > 0) ? s.blocklength[blockIndex - 1] : 0;
        uint32_t nextBlocklength = s.blocklength[blockIndex];

        bool enoughDistanceToBlockEnds = step + normal.distancetoblockend < nextBlocklength && step - normal.distancetoblockend > prevBlocklength;
        return(((step & normal.step) == 0) && enoughDistanceToBlockEnds);
    }
};

#endif
