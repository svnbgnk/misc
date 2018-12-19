#ifndef COMMON_H_
#define COMMON_H_

#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace seqan;


struct myGlobalParameters{
public:
    float flipdensity = 0.5;
    int intervalsize = 3;
    
    void print(){
        cout << flipdensity << "\n" << intervalsize << "\n";
    }
};

extern myGlobalParameters params;
extern int global;



namespace my{
void printAllP(myGlobalParameters m);

float calcsomething(float a, float b);
}

namespace seqan{

void testglobal();
    
using TMyFastConfig = seqan::FastFMIndexConfig<void, uint32_t, 2, 1>;
using TIndexConfig = seqan::BidirectionalIndex<seqan::FMIndex<void, TMyFastConfig> >;

}
#endif
