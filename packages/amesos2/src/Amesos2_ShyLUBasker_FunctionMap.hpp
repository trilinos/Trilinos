// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_Baker_FunctionMap.hpp
   \author Joshua Dennis Booth <jdbooth@sandia.gov>
           Siva Rajamanickam <srajama@sandia.gov>

   \brief  Provides a mechanism to map function calls to the correct Solver
           function based on the scalar type of Matrices and MultiVectors
*/

#ifndef AMESOS2_SHYLUBASKER_FUNCTIONMAP_HPP
#define AMESOS2_SHYLUBASKER_FUNCTIONMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_ShyLUBasker_TypeMap.hpp"

//#include "shylubasker_decl.hpp"
//#include "shylubasker_def.hpp"
//#include "shylubasker_trilinos_decl.hpp"

namespace Amesos2 {

  /* ==================== Specializations ====================
   *
   * \cond ShyLUBasker_function_specializations
   */

  /**
   * \brief Pass function calls to ShyLUBasker based on data type.

   */
#ifdef HAVE_TEUCHOS_COMPLEX
  template <>
  struct FunctionMap<ShyLUBasker,Kokkos::complex<double>>
  {
    static std::complex<double> * convert_scalar(Kokkos::complex<double> * pData) {
      return reinterpret_cast<std::complex<double> *>(pData);
    }
  };

#endif // HAVE_TEUCHOS_COMPLEX

  // if not specialized, then assume generic conversion is fine
  template <typename scalar_t>
  struct FunctionMap<ShyLUBasker,scalar_t>
  {
    static scalar_t * convert_scalar(scalar_t * pData) {
      return pData; // no conversion necessary
    }
  };


  /* \endcond ShyLUBasker_function_specializations */


} // end namespace Amesos2

#endif  // AMESOS2_SHYLUBASKER_FUNCTIONMAP_HPP
