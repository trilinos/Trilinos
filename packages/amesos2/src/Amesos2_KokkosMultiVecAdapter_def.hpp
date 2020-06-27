// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_KokkosMultiVecAdapter_def.hpp
  \author
  \date

  \brief  Amesos2::MultiVecAdapter specialization for the
          Kokkos::View class.
*/

#ifndef AMESOS2_KOKKOS_MULTIVEC_ADAPTER_DEF_HPP
#define AMESOS2_KOKKOS_MULTIVEC_ADAPTER_DEF_HPP

#include <type_traits>
#include "Amesos2_KokkosMultiVecAdapter_decl.hpp"
#include "Amesos2_Kokkos_View_Copy_Assign.hpp"


namespace Amesos2 {

  template <typename Scalar, typename ExecutionSpace >
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::MultiVecAdapter( const Teuchos::RCP<multivec_t>& m )
  : mv_(m)
  {}

  template <typename Scalar, typename ExecutionSpace >
  Scalar *
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::getMVPointer_impl() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "getMVPointer_impl not implemented.");
  }

  template <typename Scalar, typename ExecutionSpace >
  void
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::get1dCopy(const Teuchos::ArrayView<scalar_t>& av,
    size_t lda,
    Teuchos::Ptr<
      const Tpetra::Map<local_ordinal_t, global_ordinal_t,
                      MultiVecAdapter<Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace>>::node_t>> distribution_map,
                      EDistribution distribution) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "get1dCopy for kokkos not implemented.");
  }

  template <typename Scalar, typename ExecutionSpace >
  Teuchos::ArrayRCP<Scalar>
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::get1dViewNonConst (bool local)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Amesos2::MultiVecAdapter::get1dViewNonConst: "
      "Not implemented.");
  }

  template <typename Scalar, typename ExecutionSpace>
  void
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::put1dData(
      const Teuchos::ArrayView<const scalar_t>& new_data,
      size_t lda,
      Teuchos::Ptr<
        const Tpetra::Map<local_ordinal_t, global_ordinal_t,
                      MultiVecAdapter<Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace>>::node_t> > source_map,
                      EDistribution /* distribution */) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Amesos2::MultiVecAdapter::put1dData: "
      "Not implemented.");
  }

  template <typename Scalar, typename ExecutionSpace >
  std::string
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::description() const
  {
    std::ostringstream oss;
    oss << "Amesos2 adapter wrapping: ";
    oss << mv_->description();
    return oss.str();
  }


  template <typename Scalar, typename ExecutionSpace >
  void
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::describe (Teuchos::FancyOStream& os,
                                   const Teuchos::EVerbosityLevel verbLevel) const
  {
    mv_->describe (os, verbLevel);
  }

} // end namespace Amesos2

#endif // AMESOS2_KOKKOS_MULTIVEC_ADAPTER_DEF_HPP
