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


#ifndef _fei_VectorTraits_Aztec_h_
#define _fei_VectorTraits_Aztec_h_

#ifdef HAVE_FEI_AZTECOO

//
//IMPORTANT NOTE: Make sure that wherever this file is included from, it
//appears BEFORE any include of fei_base.hpp or fei_Vector.hpp !!!
//
#include <fei_VectorTraits.hpp>
#include <fei_Include_Trilinos.hpp>

namespace fei {
  /** Declare an fei::Aztec_LSVector specialization of the
      fei::VectorTraits struct.

      This allows fei::Aztec_LSVector to be used as the template parameter of the
      fei::Vector class.
  */
  template<>
  struct VectorTraits<Aztec_LSVector> {
    static const char* typeName()
      { return("fei::Aztec_LSVector"); }

    static int setValues(Aztec_LSVector* vec, int firstLocalOffset,
                         double scalar, bool isSolnVector=false)
      {
        return( vec->put(scalar) );
      }

    //note that incoming indices are point-entry indices, not block-indices.
    static int putValuesIn(Aztec_LSVector* vec,
                           int firstLocalOffset,
                           int numValues,
                           const int* indices,
                           const double* values,
                           bool sum_into,
                           bool isSolnVector=false,
                           int /*vectorIndex=0*/)
      {
        double* localVecValues = vec->startPointer();
        if (sum_into) {
          for(int i=0; i<numValues; ++i) {
            localVecValues[indices[i]-firstLocalOffset] += values[i];
          }
        }
        else {
          for(int i=0; i<numValues; ++i) {
            localVecValues[indices[i]-firstLocalOffset] = values[i];
          }
        }
        return(0);
      }

    //note that incoming indices are point-entry indices, not block-indices.
    static int copyOut(Aztec_LSVector* vec,
                       int firstLocalOffset,
                       int numValues, const int* indices, double* values,
                       bool isSolnVector=false,
                       int vectorIndex=0)
      {
        double* localVecValues = vec->startPointer();
        for(int i=0; i<numValues; ++i) {
          values[i] = localVecValues[indices[i]-firstLocalOffset];
        }

        return(0);
      }

    static double* getLocalCoefsPtr(Aztec_LSVector* vec,
                                    bool isSolnVector=false,
                                    int vectorIndex=0)
      {
        return(vec->startPointer());
      }

    static int update(Aztec_LSVector* vec,
                      double a,
                      const Aztec_LSVector* x,
                      double b)
    {
      vec->scale(b);
      vec->addVec(a, x);
      return(0);
    }

  };//struct VectorTraits<Aztec_LSVector>
}//namespace fei

#endif //HAVE_FEI_AZTECOO

#endif // _fei_VectorTraits_Aztec_hpp_
