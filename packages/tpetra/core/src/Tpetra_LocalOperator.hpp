// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_LOCALOPERATOR_HPP
#define TPETRA_LOCALOPERATOR_HPP

#include "Tpetra_LocalOperator_fwd.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <type_traits>

namespace Tpetra {

  /// \class LocalOperator
  /// \brief Abstract interface for local operators (e.g., matrices
  ///   and preconditioners).
  ///
  /// \tparam Scalar The type of the entries of the input and output
  ///   (multi)vectors.
  /// \tparam Device The Kokkos Device type; must be a specialization
  ///   of Kokkos::Device.
  template<class Scalar, class Device>
  class LocalOperator {
  public:
    using scalar_type = typename Kokkos::ArithTraits<Scalar>::val_type;
    using array_layout = Kokkos::LayoutLeft;
    using device_type =
      Kokkos::Device<typename Device::execution_space,
                     typename Device::memory_space>;

    virtual ~LocalOperator () = default;

    virtual void
    apply (Kokkos::View<const scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
           Kokkos::View<scalar_type**, array_layout,
             device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
           const Teuchos::ETransp mode,
           const scalar_type alpha,
           const scalar_type beta) const = 0;

    virtual bool hasTransposeApply () const { return false; }
  };

} // namespace Tpetra

#endif // TPETRA_LOCALOPERATOR_HPP
