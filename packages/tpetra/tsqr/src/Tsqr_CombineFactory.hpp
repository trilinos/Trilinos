//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

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
