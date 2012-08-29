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


#ifndef _fei_VectorTraits_LinSysCore_hpp_
#define _fei_VectorTraits_LinSysCore_hpp_


#include <fei_VectorTraits.hpp>
#include <fei_LinearSystemCore.hpp>

namespace fei {

/** This struct specialization defines vector traits for LinearSystemCore
  vectors (well, "vector-views" to be more precise).
*/
  template<>
  struct VectorTraits<LinearSystemCore>  {

   /** Return a string type-name for the vector. */
    static const char* typeName()
      { return("LinearSystemCore"); }

   /** Set a specified scalar value throughout the vector.
       */
    static int setValues(LinearSystemCore* vec, int firstLocalOffset,
			 double scalar, bool isSolnVector=false)
      {
	if (isSolnVector) {
	  //LinearSystemCore doesn't have a 'resetSolnVector()'.
	  return(-1);
	}
	int err = vec->resetRHSVector(scalar);
	return(err);
      }

   /** Sum values into the vector, adding to any
          that may already exist at the specified indices.
      */
    static int putValuesIn(LinearSystemCore* vec,
		     int firstLocalOffset,
		     int numValues, const int* indices, const double* values,
                     bool sum_into,
		     bool isSolnVector=false,
		     int vectorIndex=0)
    {
      int err = 0;
      if (isSolnVector) {
        if (sum_into) {
          return(-97);//LinearSystemCore allows 'put' (overwrite) operations on
          //the soln-vector, but not 'sumInto'.
        }
        else {
          err = vec->putInitialGuess(indices, values, numValues);
        }
      }
      else {
        if (sum_into) {
          err = vec->sumIntoRHSVector(numValues, values, indices);
        }
        else {
          err = vec->putIntoRHSVector(numValues, values, indices);
        }
      }
      return(err);
    }

   /** Copy values from the specified indices out into the user-allocated
          array 'values'.
      */
    static int copyOut(LinearSystemCore* vec,
		       int firstLocalOffset,
		       int numValues, const int* indices, double* values,
		       bool isSolnVector=false,
		       int vectorIndex=0)
      {
	int err = 0;
	if (isSolnVector) {
	  for(int i=0; i<numValues; ++i) {
	    if (vec->getSolnEntry(indices[i], values[i]) != 0) return(-1);
	  }
	}
	else {
	  err = vec->getFromRHSVector(numValues, values, indices);
	}

	return(err);
      }

    /** Update 'vec' = b*'vec' + a*x
       */
    static int update(LinearSystemCore* vec,
		      double a,
		      const LinearSystemCore* x,
		      double b)
    { return(-1); }

  };//struct VectorTraits
}//namespace fei

#endif // _fei_VectorTraits_LinSysCore_hpp_
