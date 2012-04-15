/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _CoefAccess_h_
#define _CoefAccess_h_


#include <cstdlib>

class CoefAccess {
 public:
  CoefAccess() : patternID_(-1), numRowIDs_(0), rowIDs_(NULL),
    numColIDsPerRow_(0), colIDs_(NULL), numRowCoefs_(0), numColCoefs_(0),
    coefs_(NULL) {}

  CoefAccess(const CoefAccess& src)
    {
      *this = src;
    }

  CoefAccess& operator=(const CoefAccess& src)
    {
      patternID_ = src.patternID_;

      numRowIDs_ = src.numRowIDs_;
      numColIDsPerRow_ = src.numColIDsPerRow_;
      numRowCoefs_ = src.numRowCoefs_;
      numColCoefs_ = src.numColCoefs_;

      if (numRowIDs_ > 0) {
	rowIDs_ = new GlobalID[numRowIDs_];
	for(int i=0; i<numRowIDs_; i++) rowIDs_[i] = src.rowIDs_[i];
      }

      if (numColIDsPerRow_ > 0 && numRowIDs_ > 0) {
	int len = numRowIDs_*numColIDsPerRow_;
	colIDs_ = new GlobalID[len];
	for(int i=0; i<len; i++) colIDs_[i] = src.colIDs_[i];
      }

      if (numRowCoefs_ > 0 && numColCoefs_ > 0) {
	int len = numRowCoefs_*numColCoefs_;
	coefs_ = new double[len];
	for(int i=0; i<len; i++) coefs_[i] = src.coefs_[i];
      }

      return(*this);
    }

  ~CoefAccess()
    {
      delete [] rowIDs_; delete [] colIDs_; delete [] coefs_;
      numRowIDs_ = 0; numColIDsPerRow_ = 0; numRowCoefs_ = 0; numColCoefs_ = 0;
    }

  int patternID_;

  int numRowIDs_;
  GlobalID* rowIDs_;

  int numColIDsPerRow_;
  GlobalID* colIDs_;

  int numRowCoefs_;
  int numColCoefs_;

  double* coefs_;
};

#endif // _CoefAccess_h_
