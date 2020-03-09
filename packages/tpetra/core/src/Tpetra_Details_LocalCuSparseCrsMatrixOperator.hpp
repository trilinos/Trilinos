// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_DETAILS_LOCALCUSPARSECRSMATRIXOPERATOR_HPP
#define TPETRA_DETAILS_LOCALCUSPARSECRSMATRIXOPERATOR_HPP

/// \file Tpetra_Details_LocalCuSparseCrsMatrixOperator.hpp
/// \brief Declaration of
///   Tpetra::Details::LocalCuSparseCrsMatrixOperator

#include "Tpetra_Details_LocalCuSparseCrsMatrixOperator_fwd.hpp"
#include "Tpetra_Details_LocalCrsMatrixOperatorWithSetup.hpp"
#ifdef HAVE_TPETRACORE_CUSPARSE
#include "Tpetra_Details_CuSparseHandle_fwd.hpp"
#include "Tpetra_Details_CuSparseMatrix_fwd.hpp"
#endif // HAVE_TPETRACORE_CUSPARSE

namespace Tpetra {
namespace Details {

/// \brief Subclass of LocalCrsMatrixOperatorWithSetup that may call
///   cuSPARSE, if that library is enabled in Trilinos, and if all the
///   template arguments match cuSPARSE's requirements.
template<class MultiVectorScalar, class MatrixScalar, class Device>
class LocalCuSparseCrsMatrixOperator :
  public LocalCrsMatrixOperatorWithSetup<
    MultiVectorScalar, MatrixScalar, Device>
{
private:
  using base_type =
    LocalCrsMatrixOperator<MultiVectorScalar, MatrixScalar, Device>;

public:
  using execution_space = typename Device::execution_space;
  using local_matrix_type = typename base_type::local_matrix_type;

  LocalCuSparseCrsMatrixOperator(
    const execution_space& /* execSpace */,
    const std::shared_ptr<local_matrix_type>& A)
    : base_type(A)
  {}
  ~LocalCuSparseCrsMatrixOperator() override = default;
};

#ifdef HAVE_TPETRACORE_CUSPARSE

/// \brief Common base class of the specializations of
///   LocalCuSparseCrsMatrixOperator that can actually use cuSPARSE.
///
/// LocalCuSparseCrsMatrixOperator can only use cuSPARSE with the
/// following template arguments:
///
/// float, float, Device<Cuda, Cuda::memory_space>, or
///
/// double, double, Device<Cuda, Cuda::memory_space>.
///
/// If you instantiate LocalCuSparseCrsMatrixOperator with any other
/// template arguments, it will just do whatever
/// LocalCrsMatrixOperatorWithSetup does.  That class may or may not
/// (currently does not) use cuSPARSE.
///
/// The reason for only using cuSPARSE for Scalar=float and
/// Scalar=double, is that Kokkos::complex and std::complex have less
/// stringent alignment requirements than CUDA's corresponding native
/// complex types.
///
/// The point of this common base class is to avoid code duplication
/// for the full specializations below.
template<class Scalar>
class LocalCuSparseCrsMatrixOperatorRealBase :
  public LocalCrsMatrixOperatorWithSetup<
    Scalar, Scalar,
    Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space> >
{
private:
  using device_type =
    Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>;
  using base_type =
    LocalCrsMatrixOperatorWithSetup<Scalar, Scalar, device_type>;
  using array_layout =
    typename ::Tpetra::LocalOperator<Scalar, device_type>::
      array_layout;
  using LO = ::Tpetra::Details::DefaultTypes::local_ordinal_type;

public:
  using execution_space = Kokkos::Cuda;
  using local_matrix_type = typename base_type::local_matrix_type;

  LocalCuSparseCrsMatrixOperatorRealBase(
    const execution_space& execSpace,
    const std::shared_ptr<local_matrix_type>& A);
  // This can't be = default, because else unique_ptr<CuSparseMatrix>
  // would demand a definition of CuSparseMatrix.
  ~LocalCuSparseCrsMatrixOperatorRealBase() override;

  void resumeFill() override;
  void
  setMinMaxNumberOfEntriesPerRow(const LO minNumEntPerRow,
                                 const LO maxNumEntPerRow) override;
  void fillComplete() override;

  void
  apply(Kokkos::View<const Scalar**, array_layout,
          device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > X,
        Kokkos::View<Scalar**, array_layout,
          device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Y,
        const Teuchos::ETransp mode,
        const Scalar alpha,
        const Scalar beta) const override;

private:
  LO minNumEntPerRow_ = 0;
  LO maxNumEntPerRow_ = 0;
  std::shared_ptr<CuSparseHandle> handle_;
  std::unique_ptr<CuSparseMatrix, decltype(&Impl::deleteCuSparseMatrix)> matrix_;
  Kokkos::View<LO*, device_type> ptr_;
};

//! Full specialization that uses cuSPARSE for Scalar=float.
template<>
class LocalCuSparseCrsMatrixOperator<float, float,
  Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>> :
  public LocalCuSparseCrsMatrixOperatorRealBase<float>
{
private:
  using device_type =
    Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>;
  using base_type = LocalCuSparseCrsMatrixOperatorRealBase<float>;

public:
  using execution_space = Kokkos::Cuda;
  using local_matrix_type = typename base_type::local_matrix_type;

  LocalCuSparseCrsMatrixOperator(
    const execution_space& execSpace,
    const std::shared_ptr<local_matrix_type>& A);
  ~LocalCuSparseCrsMatrixOperator() override = default;
};

//! Full specialization that uses cuSPARSE for Scalar=double.
template<>
class LocalCuSparseCrsMatrixOperator<double, double,
  Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>> :
  public LocalCuSparseCrsMatrixOperatorRealBase<double>
{
private:
  using device_type =
    Kokkos::Device<Kokkos::Cuda, Kokkos::Cuda::memory_space>;
  using base_type = LocalCuSparseCrsMatrixOperatorRealBase<double>;

public:
  using execution_space = Kokkos::Cuda;
  using local_matrix_type = typename base_type::local_matrix_type;

  LocalCuSparseCrsMatrixOperator(
    const execution_space& execSpace,
    const std::shared_ptr<local_matrix_type>& A);
  ~LocalCuSparseCrsMatrixOperator() override = default;
};

#endif // HAVE_TPETRACORE_CUSPARSE

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_LOCALCUSPARSECRSMATRIXOPERATOR_HPP
