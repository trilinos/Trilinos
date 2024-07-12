// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_OPERATORWITHDIAGONAL_HPP
#define TPETRA_OPERATORWITHDIAGONAL_HPP

#include <Teuchos_RCP.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_RowMatrix.hpp>
#include "Teuchos_VerbosityLevel.hpp"

namespace Tpetra {

template <class Scalar        = Tpetra::RowMatrix<>::scalar_type,
          class LocalOrdinal  = typename Tpetra::RowMatrix<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename Tpetra::RowMatrix<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class OperatorWithDiagonal : public Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using mv_type  = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using map_type = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

  using row_matrix_type = RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  using impl_scalar_type = typename row_matrix_type::impl_scalar_type;
  using mag_type         = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

  using local_inds_device_view_type =
      typename row_matrix_type::local_inds_device_view_type;
  using local_inds_host_view_type =
      typename row_matrix_type::local_inds_host_view_type;
  using nonconst_local_inds_host_view_type =
      typename row_matrix_type::nonconst_local_inds_host_view_type;

  using global_inds_device_view_type =
      typename row_matrix_type::global_inds_device_view_type;
  using global_inds_host_view_type =
      typename row_matrix_type::global_inds_host_view_type;
  using nonconst_global_inds_host_view_type =
      typename row_matrix_type::nonconst_global_inds_host_view_type;

  using values_device_view_type =
      typename row_matrix_type::values_device_view_type;
  using values_host_view_type =
      typename row_matrix_type::values_host_view_type;
  using nonconst_values_host_view_type =
      typename row_matrix_type::nonconst_values_host_view_type;

  //! @name Constructor/Destructor
  //@{

  virtual ~OperatorWithDiagonal() = default;

  size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
    throw std::runtime_error("Not implemented.");
  }

  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    throw std::runtime_error("Not implemented.");
  }

  size_t getGlobalMaxNumRowEntries() const {
    throw std::runtime_error("Not implemented.");
  }

  LocalOrdinal getBlockSize() const {
    throw std::runtime_error("Not implemented.");
  }

  size_t getLocalMaxNumRowEntries() const {
    throw std::runtime_error("Not implemented.");
  }

  bool hasColMap() const {
    return false;
  }

  bool isLocallyIndexed() const {
    return true;
  }

  bool isGloballyIndexed() const {
    return true;
  }

  bool isFillComplete() const {
    return true;
  }

  bool supportsRowViews() const {
    return false;
  }

  void
  getGlobalRowCopy(GlobalOrdinal GlobalRow,
                   nonconst_global_inds_host_view_type& Indices,
                   nonconst_values_host_view_type& Values,
                   size_t& NumEntries) const {
    throw std::runtime_error("Not implemented.");
  }

  void
  getLocalRowCopy(LocalOrdinal LocalRow,
                  nonconst_local_inds_host_view_type& Indices,
                  nonconst_values_host_view_type& Values,
                  size_t& NumEntries) const {
    throw std::runtime_error("Not implemented.");
  }

  void
  getGlobalRowView(GlobalOrdinal GlobalRow,
                   global_inds_host_view_type& indices,
                   values_host_view_type& values) const {
    throw std::runtime_error("Not implemented.");
  }

  void
  getLocalRowView(LocalOrdinal LocalRow,
                  local_inds_host_view_type& indices,
                  values_host_view_type& values) const {
    throw std::runtime_error("Not implemented.");
  }

  void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
    throw std::runtime_error("Not implemented.");
  }

  void rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
    throw std::runtime_error("Not implemented.");
  }

  mag_type getFrobeniusNorm() const {
    return 0.;
  }
};
}  // namespace Tpetra

#endif  // TPETRA_OPERATORWITHDIAGONAL_HPP
