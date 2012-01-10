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


#include <fei_macros.hpp>
#include <fei_defs.h>
#include <test_utils/ElemBlock.hpp>
#include <cstdlib>

//==============================================================================
ElemBlock::ElemBlock()
 : blockID_(0),
   numElements_(0),
   numNodesPerElement_(0),
   numFieldsPerNode_(NULL),
   nodalFieldIDs_(NULL),
   elemIDs_(NULL),
   elemConn_(NULL),
   numStiffRows_(0),
   elemFormat_(0),
   elemStiff_(NULL),
   elemLoad_(NULL),
   numElemDOF_(0),
   elemDOFFieldIDs_(NULL),
   interleaveStrategy_(0),
   lumpingStrategy_(0)
{
}

//==============================================================================
ElemBlock::~ElemBlock() {
   deleteMemory();
}

//==============================================================================
void ElemBlock::deleteMemory() {
   for(int i=0; i<numElements_; i++) {
      for(int j=0; j<numStiffRows_; j++) {
         delete [] elemStiff_[i][j];
      }
      delete [] elemStiff_[i];
      delete [] elemLoad_[i];
      delete [] elemConn_[i];
   }

   for(int j=0; j<numNodesPerElement_; j++) {
      delete [] nodalFieldIDs_[j];
   }
   delete [] nodalFieldIDs_;
   delete [] numFieldsPerNode_;

   delete [] elemStiff_;
   delete [] elemLoad_;
   delete [] elemConn_;
   delete [] elemIDs_;

   if (numElemDOF_ > 0) {
      delete [] elemDOFFieldIDs_;
      numElemDOF_ = 0;
   }

   numElements_ = 0;
   numNodesPerElement_ = 0;
}

