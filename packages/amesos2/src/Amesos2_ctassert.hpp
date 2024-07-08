// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
 * \file   Amesos2_ctassert.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Tue Jul 26 15:31:53 2011
 * 
 * \brief  Simple compile-time assertion class
 */

/* ETB: I think the Teuchos CompileTimeAssert class is backwards
 * (i.e. if the test evaluates to true, a compile-time error is
 * thrown), so I'm defining my own.
 *
 * This definition is based off a Dr Dobb's article on compile-time
 * assertions.
 */

#ifndef AMESOS2_CTASSERT_HPP
#define AMESOS2_CTASSERT_HPP

namespace Amesos2{
  
  template <bool assertion>
  struct ctassert {
    enum { N = 1 - 2*int(!assertion) };
    //  1 if assertion is true
    // -1 if assertion is false
    static char A[N];
  };

  template <bool assertion>
  char ctassert<assertion>::A[N];

}

#endif	// AMESOS2_CTASSERT_HPP
