
#include "WhizBangVectorAdder.h"

// Constuctor takes two vectors.  Adder adds them together, result in V1

WhizBangVectorAdder::WhizBangVectorAdder(TSF_RDP_MultiVector * V1, TSF_RDP_MultiVector * V2 ):
  V1_(V1),
  V2_(V2)
{}
  
WhizBangVectorAdder::~WhizBangVectorAdder() {}
  
void WhizBangVectorAdder::Add() { 

  V1_->Update(1.0, *V2_, 1.0);
  return;
}
