// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER


#ifndef PANZER_PTR_FROM_STL_VECTOR_HPP
#define PANZER_PTR_FROM_STL_VECTOR_HPP

#include <vector>

namespace panzer {

  /** Function to eliminate runtime errors triggered by compiling with
      checked_stl enabled builds.  We have cases where raw pointers
      must be passed to external libraries like Epetra.  If the raw
      data is stored in a stl vector object, it is legal to access the
      vector using &v[0].  In an effort to be efficient (for
      performance essential kernels) we do not check the size of the
      stl vector before grabbing the pointer. We just grab the pointer
      and pass the size to the epetra function.  However, if the
      vector size is zero, then the checked stl will throw an error.
      This is a case of checked stl forcing us to write less efficient
      code (by having to check the size).  To preserve efficiency, we
      use a define based on if checked stl is enabled.
   */
  template <typename T>
  inline T* ptrFromStlVector(std::vector<T>& v)
  {
#ifdef _GLIBCXX_DEBUG  
    if (v.size() > 0)
      return &v[0];
    
    return NULL;
#else
    return &v[0];
#endif
  }

}

#endif
