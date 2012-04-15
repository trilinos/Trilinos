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

#ifndef _AccessPattern_h_
#define _AccessPattern_h_


#include <cstdlib>

class AccessPattern {
 public:
  AccessPattern() : ID_(-1), numRowIDs_(0), numFieldsPerRow_(NULL),
    rowFieldIDs_(NULL), numColIDsPerRow_(0), numFieldsPerCol_(NULL),
    colFieldIDs_(NULL), interleaveStrategy_(0) {}

  ~AccessPattern()
    {
      int i;
      for(i=0; i<numRowIDs_; i++) delete [] rowFieldIDs_[i];
      for(i=0; i<numColIDsPerRow_; i++) delete [] colFieldIDs_[i];

      delete [] rowFieldIDs_;
      delete [] colFieldIDs_;
      delete [] numFieldsPerRow_;
      delete [] numFieldsPerCol_;
      numRowIDs_ = 0;
      numColIDsPerRow_ = 0;
    }

  int ID_;
  int numRowIDs_;
  int* numFieldsPerRow_;
  int** rowFieldIDs_;
  int numColIDsPerRow_;
  int* numFieldsPerCol_;
  int** colFieldIDs_;
  int interleaveStrategy_;
};

#endif // _AccessPattern_h_
