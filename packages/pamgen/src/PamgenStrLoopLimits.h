// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef PAMSTRLOOPLIMITS_H
#define PAMSTRLOOPLIMITS_H

// LoopLimits
namespace PAMGEN_NEVADA {
struct LoopLimits {
  long long is;  // Start of i-loop
  long long ie;  // End   of i-loop
  long long il;  // Last i-iterate for WorkSet's
  long long js;
  long long je;
  long long jl;
  long long ks;
  long long ke;
  long long kl;
  long long jstride;  // j-stride size
  long long kstride;  // k-stride size for the array
  long long total;   // Total array count for strip-mining

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
