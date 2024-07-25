// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
