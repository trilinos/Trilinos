#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <iomanip>
using namespace std;

#include "TSF_RDP_MultiVector.h"

// Constuctor takes two vectors.  Adder adds them together, result in V1
class WhizBangVectorAdder {

public:

  WhizBangVectorAdder(TSF_RDP_MultiVector * V1, TSF_RDP_MultiVector * V2 );
  
    ~WhizBangVectorAdder();
  
    void Add();  // Adds V1 and V2
  
private:
  
  TSF_RDP_MultiVector * V1_;  
  TSF_RDP_MultiVector * V2_;

};
