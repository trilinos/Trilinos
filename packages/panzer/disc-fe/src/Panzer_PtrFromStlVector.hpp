// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
