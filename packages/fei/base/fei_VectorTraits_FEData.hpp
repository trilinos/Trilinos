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


#ifndef _fei_VectorTraits_FEData_hpp_
#define _fei_VectorTraits_FEData_hpp_

//This file defines vector traits for FiniteElementData vectors
//(well, "vector-views" to be more precise).
//

#include <fei_VectorTraits.hpp>
#include <fei_FiniteElementData.hpp>

namespace fei {

  /** specialization for FiniteElementData */
  template<>
  struct VectorTraits<FiniteElementData>  {

    /** name of VectorTraits type */
    static const char* typeName()
      { return("FiniteElementData"); }

    /** set all vector values to specified scalar */
    static int setValues(FiniteElementData* vec, int firstLocalOffset,
			 double scalar, bool isSolnVector=false)
      {
	return(-1);
      }

    /** sum-into operation for vector data */
    static int putValuesIn(FiniteElementData* vec,
		     int firstLocalOffset,
		     int numValues, const int* indices, const double* values,
                     bool sum_into,
		     bool isSolnVector=false,
		     int vectorIndex=0)
      {
	return(-1);
      }

    /** copy out vector data */
    static int copyOut(FiniteElementData* vec,
		       int firstLocalOffset,
		       int numValues, const int* indices, double* values,
		       bool isSolnVector=false,
		       int vectorIndex=0)
      {
	return(-1);
      }

    /** vec = b*vec + a*x */
    static int update(FiniteElementData* vec,
		      double a,
		      const FiniteElementData* x,
		      double b)
    { return(-1); }

  };//struct VectorTraits
}//namespace fei

#endif // _fei_VectorTraits_FEData_hpp_
