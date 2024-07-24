// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_Combine.hpp
/// \brief Interface to TSQR's six computational kernels.

#ifndef TSQR_COMBINEFACTORY_HPP
#define TSQR_COMBINEFACTORY_HPP

#include "Tsqr_CombineDefault.hpp"
#include "Tsqr_CombineNative.hpp"
#include "Teuchos_TestForException.hpp"
#include <memory>
#include <string>

namespace TSQR {
  /// \class CombineFactory
  /// \brief Factory for creating Combine instances.
  /// \author Mark Hoemmen
  template<class Ordinal, class Scalar>
  class CombineFactory {
  public:
    /// \brief Given the maximum number of columns in either the
    ///   matrix to factor, or the matrix to which to apply a Q factor
    ///   or compute an explicit Q factor, return an appropriate
    ///   Combine implementation.
    static std::unique_ptr<Combine<Ordinal, Scalar>>
    create (const Ordinal maxNumCols)
    {
      // FIXME (mfh 19 Dec 2019) This _should_ depend on the BLAS
      // implementation.
      constexpr Ordinal blas_3_threshold = 32;
      if (maxNumCols >= blas_3_threshold) {
        using impl_type = CombineDefault<Ordinal, Scalar>;
        // NOTE (mfh 19 Dec 2019) We can't use std::make_unique yet,
        // because it requires C++14.
        return std::unique_ptr<impl_type> (new impl_type);
      }
      else {
        using impl_type = CombineNative<Ordinal, Scalar>;
        return std::unique_ptr<impl_type> (new impl_type);
      }
    }

    static std::unique_ptr<Combine<Ordinal, Scalar>>
    create (const std::string& combineType)
    {
      if (combineType == "CombineNative" ||
          combineType == "Native") {
        using impl_type = CombineNative<Ordinal, Scalar>;
        return std::unique_ptr<impl_type> (new impl_type);
      }        
      else if (combineType == "CombineDefault" ||
               combineType == "Default") {
        using impl_type = CombineDefault<Ordinal, Scalar>;
        return std::unique_ptr<impl_type> (new impl_type);
      }        
      else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument, "TSQR::CombineFactory: "
           "Invalid Combine subclass name \"" << combineType <<
           "\".");
      }
    }
  };

} // namespace TSQR

#endif // TSQR_COMBINEFACTORY_HPP
