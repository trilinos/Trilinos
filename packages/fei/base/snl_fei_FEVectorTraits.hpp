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


#ifndef _snl_fei_FEVectorTraits_hpp_
#define _snl_fei_FEVectorTraits_hpp_

#include <fei_macros.hpp>

namespace snl_fei {

  /** Internal implementation vector traits. Define a "template" for accessing
      vector data.
      Provide function stubs for default type "T", which will catch the
      use of any vector type for which specialized traits have not been defined.
  */
  template<class T>
  struct FEVectorTraits {

    /** Return a string type-name for the vector. */
    static const char* typeName()
      { return("unsupported"); }

    /** Reset (zero) the vector.
     */
    static int reset(T* /*vec*/)
      { return(-1); }

    /** Sum an element-vector contribution into the vector object */
    static int sumInElemVector(T* /*vec*/,
			       int /*elemBlockID*/,
			       int /*elemID*/,
			       int /*numNodes*/,
			       const int* /*nodeNumbers*/,
			       const int* /*dofPerNode*/,
             const int* /*dof_ids*/,
			       const double* /*coefs*/)
      { return(-1); }

    /** Copy data out of the vector object */
    static int copyOut(T* /*vec*/,
		       int /*nodeNumber*/,
		       int /*dofOffset*/,
		       double& /*value*/)
      { return( -1 ); }

  };//struct FEVectorTraits
}//namespace snl_fei

#endif // _snl_fei_FEVectorTraits_hpp_
