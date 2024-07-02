// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef IFPACK2_DETAILS_LAPACKSUPPORTSSCALAR_HPP
#define IFPACK2_DETAILS_LAPACKSUPPORTSSCALAR_HPP

#include "Teuchos_ConfigDefs.hpp"

#ifdef HAVE_TEUCHOS_COMPLEX
#  include <complex>
#endif // HAVE_TEUCHOS_COMPLEX

namespace Ifpack2 {
namespace Details {

  /// \class LapackSupportsScalar
  /// \brief Type traits class that says whether Teuchos::LAPACK
  ///   has a valid implementation for the given ScalarType.
  template<class ScalarType>
  class LapackSupportsScalar {
  public:
    const static bool value = false;
  };

  template<>
  class LapackSupportsScalar<float> {
  public:
    const static bool value = true;
  };

  template<>
  class LapackSupportsScalar<double> {
  public:
    const static bool value = true;
  };

#ifdef HAVE_TEUCHOS_COMPLEX
  template<>
  class LapackSupportsScalar<std::complex<float> > {
  public:
    const static bool value = true;
  };

  template<>
  class LapackSupportsScalar<std::complex<double> > {
  public:
    const static bool value = true;
  };
#endif // HAVE_TEUCHOS_COMPLEX

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LAPACKSUPPORTSSCALAR_HPP
