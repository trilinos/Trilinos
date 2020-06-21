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

#ifndef TSQR_COMBINEUSER_HPP
#define TSQR_COMBINEUSER_HPP

#include "Tsqr_CombineFactory.hpp"

namespace TSQR {
namespace Impl {

/// \class CombineUser
/// \brief Private base class for TSQR classes that use Combine.
///
/// Classes that use Combine should inherit privately from this class,
/// in order to reuse getCombine.
template<class LocalOrdinal, class Scalar>
class CombineUser {
public:
  /// \brief Given the maximum number of columns that the caller
  ///   intends to give to Combine functions, return the best choice
  ///   of Combine implementation.
  Combine<LocalOrdinal, Scalar>&
  getCombine (const LocalOrdinal maxNumCols) const {
    if (combine_.get () == nullptr) {
      using factory_type = CombineFactory<LocalOrdinal, Scalar>;
      combine_ = factory_type::create (maxNumCols);
    }
    return *combine_;
  }

  //! Return a specific Combine implementation.
  Combine<LocalOrdinal, Scalar>&
  getCombine (const std::string& combineType) const {
    if (combine_.get () == nullptr) {
      using factory_type = CombineFactory<LocalOrdinal, Scalar>;
      combine_ = factory_type::create (combineType);
    }
    return *combine_;
  }

private:
  mutable std::unique_ptr<Combine<LocalOrdinal, Scalar>> combine_;
};

} // namespace Impl
} // namespace TSQR

#endif // TSQR_COMBINEUSER_HPP
