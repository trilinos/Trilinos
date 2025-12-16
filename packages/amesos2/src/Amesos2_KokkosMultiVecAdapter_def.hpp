// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  Teuchos::RCP< Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >
   MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::clone() const
  {
    using MV = Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace>;
    MV Y("clonedY", mv_->extent(0), mv_->extent(1));
    return Teuchos::rcp( &Y );
  }

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

  template <typename Scalar, typename ExecutionSpace>
  template<typename KV, typename host_ordinal_type_array>
  int
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::gather (
      KV& kokkos_new_view,
      host_ordinal_type_array &perm_g2l,
      host_ordinal_type_array &recvCountRows,
      host_ordinal_type_array &recvDisplRows,
      EDistribution distribution ) const
  {
    return -1;
  }

  template <typename Scalar, typename ExecutionSpace>
  template<typename KV, typename host_ordinal_type_array>
  int
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::scatter (
      KV& kokkos_new_view,
      host_ordinal_type_array &perm_g2l,
      host_ordinal_type_array &recvCountRows,
      host_ordinal_type_array &recvDisplRows,
      EDistribution distribution ) const
  {
    return -1;
  }

  template <typename Scalar, typename ExecutionSpace >
  std::string
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::description() const
  {
    std::ostringstream oss;
    size_t m = mv_->extent(0);
    size_t n = mv_->extent(1);
    oss << "Amesos2 adapter wrapping: ";
    oss << " Kokkos::View(" << std::to_string(m) << " x " << std::to_string(n) << ")";
    oss << " of type " << std::string(typeid(Scalar).name());
    return oss.str();
  }


  template <typename Scalar, typename ExecutionSpace >
  void
  MultiVecAdapter<
    Kokkos::View<Scalar**, Kokkos::LayoutLeft, ExecutionSpace> >::describe (Teuchos::FancyOStream& os,
                                   const Teuchos::EVerbosityLevel verbLevel) const
  {
    size_t m = mv_->extent(0);
    size_t n = mv_->extent(1);
    os << " Kokkos::View(" << std::to_string(m) << " x " << std::to_string(n) << ")";
    os << " of type " << std::string(typeid(Scalar).name()) << std::endl;
  }

} // end namespace Amesos2

#endif // AMESOS2_KOKKOS_MULTIVEC_ADAPTER_DEF_HPP
