#include "common.h"



namespace seqan{
void testglobal(){
    cout << global << endl;
}

}

namespace my{
    
void printAllP(myGlobalParameters m)
{
    cout << m.flipdensity << "\n" << m.intervalsize << "\n" << endl;
}

float calcsomething(float a, float b){
    float calc = (a + b) * params.flipdensity;
    return(calc);
}
}