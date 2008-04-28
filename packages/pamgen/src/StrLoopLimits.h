// $Id$

#ifndef STRLOOPLIMITS_H
#define STRLOOPLIMITS_H

// LoopLimits
namespace PAMGEN_NEVADA {
struct LoopLimits {
  int is;  // Start of i-loop
  int ie;  // End   of i-loop
  int il;  // Last i-iterate for WorkSet's
  int js;
  int je;
  int jl;
  int ks;
  int ke;
  int kl;
  int jstride;  // j-stride size
  int kstride;  // k-stride size for the array
  int total;   // Total array count for strip-mining

  LoopLimits():
  is(0),
  ie(0),
  il(0),
  js(0),
  je(0),
  jl(0),
  ks(0),
  ke(0),
  kl(0),
  jstride(0),
  kstride(0),
  total(0)
  {}

};
} // end namespace PAMGEN_NEVADA 

#endif // end of STRLOOPLIMITS_H
