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


#ifndef _fei_VectorTraits_CSVec_hpp_
#define _fei_VectorTraits_CSVec_hpp_

#include <fei_VectorTraits.hpp>
#include <fei_CSVec.hpp>

namespace fei {
  template<>
  struct VectorTraits<CSVec> {
    static const char* typeName()
      { return("CSVec"); }
    
    static int setValues(CSVec* vec, int firstLocalOffset,
                         double scalar, bool isSolnVector=false)
      {
        set_values(*vec, scalar);
        return(0);
      }

    static int putValuesIn(CSVec* vec,
                     int firstLocalOffset,
                     int numValues, const int* indices, const double* values,
                     bool sum_into,
                     bool isSolnVector=false,
                     int vectorIndex=0)
      {
        if (sum_into) {
          for(int i=0; i<numValues; ++i) {
            if (indices[i] < 0) continue;
            add_entry(*vec, indices[i], values[i]);
          }
        }
        else {
          for(int i=0; i<numValues; ++i) {
            if (indices[i] < 0) continue;
            put_entry(*vec, indices[i], values[i]);
          }
        }

        return( 0 );
      }

    static int copyOut(CSVec* vec,
                       int firstLocalOffset,
                       int numValues, const int* indices, double* values,
                       bool isSolnVector=false,
                       int vectorIndex=0)
      {
        for(int i=0; i<numValues; ++i) {
          try {
            values[i] = get_entry(*vec, indices[i]);
          }
          catch(...) {}
        }
        return(0);
      }

    static int update(CSVec* vec,
                      double a,
                      const CSVec* x,
                      double b)
    { return(-1); }

  };
}//namespace fei

#endif // _fei_VectorTraits_CSVec_hpp_

