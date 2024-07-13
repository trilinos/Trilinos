// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
   \file   Amesos2_CssMKL_TypeMap.hpp
   \author John Doe <jd@sandia.gov>
   \date

   \brief Provides definition of CssMKL types as well as
          conversions and type traits.  For the purpose of
          demonstration, we assume that CssMKL has defined its own
          complex data-types called `complex' and `doublecomplex'.
*/

#ifndef AMESOS2_CLUSTERSPARSEMKL_TYPEMAP_HPP
#define AMESOS2_CLUSTERSPARSEMKL_TYPEMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <mkl_types.h>
#include <mkl_dss.h>

#include <Teuchos_as.hpp>
#ifdef HAVE_TEUCHOS_COMPLEX
#include <Teuchos_SerializationTraits.hpp>
#endif

#include "Amesos2_TypeMap.hpp"
#include "Amesos2_PardisoMKL_TypeMap.hpp"

namespace Amesos2 {

  // forward declaration due to circular reference
  template <class, class> class CssMKL;

  /* Specialize the Amesos::TypeMap struct for CssMKL types.
   *
   * Additional nested types may be added without harm.  For an example, look at
   * Amesos2_Superlu_TypeMap.hpp
   */

  template <>
  struct TypeMap<CssMKL,float>
  {
    typedef PMKL::_REAL_t type;
    typedef PMKL::_REAL_t magnitude_type;
  };


  template <>
  struct TypeMap<CssMKL,double>
  {
    typedef PMKL::_DOUBLE_PRECISION_t type;
    typedef PMKL::_DOUBLE_PRECISION_t magnitude_type;
  };

#ifdef HAVE_TEUCHOS_COMPLEX

  /*
   * We map the std complex types to the appropriate CssMKL complex
   * types.
   */

  template <>
  struct TypeMap<CssMKL,std::complex<float> >
  {
    typedef PMKL::_MKL_Complex8 type;
    typedef PMKL::_REAL_t magnitude_type;
  };


  template <>
  struct TypeMap<CssMKL,std::complex<double> >
  {
    typedef PMKL::_DOUBLE_COMPLEX_t type;
    typedef PMKL::_DOUBLE_PRECISION_t magnitude_type;
  };


  template <>
  struct TypeMap<CssMKL,PMKL::_MKL_Complex8>
  {
    typedef PMKL::_MKL_Complex8 type;
    typedef PMKL::_REAL_t magnitude_type;
  };


  template <>
  struct TypeMap<CssMKL,PMKL::_DOUBLE_COMPLEX_t>
  {
    typedef PMKL::_DOUBLE_COMPLEX_t type;
    typedef PMKL::_DOUBLE_PRECISION_t magnitude_type;
  };
#endif  // HAVE_TEUCHOS_COMPLEX

  template <>
  struct TypeMap<CssMKL,int>
  {
    typedef PMKL::_INTEGER_t type;
    //typedef int   type;
  };

  template <>
  struct TypeMap<CssMKL,long long int>
  {
    typedef long long int type;
  };

  /*
   * We check whether the size of long int is bigger than an int.  If
   * it is, then long int should be the same size as a long long int,
   * so we can safely promote.  Otherwise, long int will probably be
   * the same size as int, and we can safely treat it as such.
   */
  template <>
  struct TypeMap<CssMKL,long int>
  {
    typedef std::conditional_t<
      sizeof(int) < sizeof(long int),
      TypeMap<CssMKL,long long int>::type,
      TypeMap<CssMKL,int>::type > type;
  };

} // end namespace Amesos

#endif  // AMESOS2_CLUSTERSPARSEMKL_TYPEMAP_HPP
