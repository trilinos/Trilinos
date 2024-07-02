// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef AMESOS2_MATRIXTRAITS_HPP
#define AMESOS2_MATRIXTRAITS_HPP

#include "Amesos2_config.h"

#include <Tpetra_CrsMatrix.hpp>


#ifdef HAVE_AMESOS2_EPETRA
#  include <Epetra_RowMatrix.h>
#  include <Epetra_CrsMatrix.h>
// #  include <Epetra_MsrMatrix.h>
#  include <Epetra_VbrMatrix.h>
// and perhaps some others later...
#endif

#include "Amesos2_Util.hpp"

namespace Amesos2 {

  // The declaration
  template <class Matrix>
  struct MatrixTraits {};

  /*******************
   * Specializations *
   *******************/

  template < typename Scalar,
             typename LocalOrdinal,
             typename GlobalOrdinal,
             typename Node >
  struct MatrixTraits<
    Tpetra::RowMatrix<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node> > {
    typedef Scalar scalar_t;
    typedef LocalOrdinal local_ordinal_t;
    typedef GlobalOrdinal global_ordinal_t;
    typedef Node node_t;

    typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>  matrix_type;
    typedef typename matrix_type::impl_scalar_type                    impl_scalar_type;

    typedef typename matrix_type::nonconst_global_inds_host_view_type global_host_idx_type;
    typedef typename matrix_type::nonconst_values_host_view_type      global_host_val_type;

    typedef row_access major_access;
  };

  template < typename Scalar,
             typename LocalOrdinal,
             typename GlobalOrdinal,
             typename Node >
  struct MatrixTraits<
    Tpetra::CrsMatrix<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node> > {
    typedef Scalar scalar_t;
    typedef LocalOrdinal local_ordinal_t;
    typedef GlobalOrdinal global_ordinal_t;
    typedef Node node_t;

    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>  matrix_type;
    typedef typename matrix_type::impl_scalar_type                    impl_scalar_type;

    typedef typename matrix_type::nonconst_global_inds_host_view_type global_host_idx_type;
    typedef typename matrix_type::nonconst_values_host_view_type      global_host_val_type;

    typedef row_access major_access;
  };

  template < typename Scalar,
             typename LocalOrdinal,
             typename DeviceType >
  struct MatrixTraits<
    KokkosSparse::CrsMatrix<Scalar,
                      LocalOrdinal,
                      DeviceType> > {
    typedef Scalar scalar_t;
    typedef LocalOrdinal local_ordinal_t;
    typedef Tpetra::Map<>::global_ordinal_type global_ordinal_t;
    typedef LocalOrdinal global_size_t;

    typedef KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, DeviceType>  matrix_type;
    typedef Scalar impl_scalar_type;
    typedef Tpetra::Map<>::node_type node_t;

    typedef typename matrix_type::HostMirror::index_type  global_host_idx_type;
    typedef typename matrix_type::HostMirror::values_type global_host_val_type;

    typedef row_access major_access;
  };

#ifdef HAVE_AMESOS2_EPETRA

  template <>
  struct MatrixTraits<Epetra_RowMatrix> {
    typedef double scalar_t;
    typedef double impl_scalar_type;
    typedef int local_ordinal_t;
    typedef Tpetra::Map<>::global_ordinal_type global_ordinal_t;
    typedef Tpetra::Map<>::node_type node_t;

    typedef Epetra_RowMatrix matrix_type;
    typedef matrix_type local_matrix_t;
    typedef int* sparse_ptr_type;
    typedef int* sparse_idx_type;
    typedef double* sparse_values_type;

    typedef Kokkos::View<global_ordinal_t*, Kokkos::HostSpace> global_host_idx_type;
    typedef Kokkos::View<scalar_t*,         Kokkos::HostSpace> global_host_val_type;

    typedef row_access major_access;
  };

  template <>
  struct MatrixTraits<Epetra_CrsMatrix> {
    typedef double scalar_t;
    typedef double impl_scalar_type;
    typedef int local_ordinal_t;
    typedef Tpetra::Map<>::global_ordinal_type global_ordinal_t;
    typedef Tpetra::Map<>::node_type node_t;

    typedef Epetra_CrsMatrix matrix_type;
    typedef matrix_type local_matrix_t;
    typedef int* sparse_ptr_type;
    typedef int* sparse_idx_type;
    typedef double* sparse_values_type;

    typedef Kokkos::View<global_ordinal_t*, Kokkos::HostSpace> global_host_idx_type;
    typedef Kokkos::View<scalar_t*,         Kokkos::HostSpace> global_host_val_type;

    typedef row_access major_access;
  };

  // template <>
  // struct MatrixTraits<Epetra_MsrMatrix> {
  //   typedef double scalar_t;
  //   typedef int local_ordinal_t;
  //   typedef Tpetra::Map<>::global_ordinal_type global_ordinal_t;
  //   typedef Tpetra::Map<>::node_type node_t;

  //   typedef row_access major_access;
  // };

  template <>
  struct MatrixTraits<Epetra_VbrMatrix> {
    typedef double scalar_t;
    typedef double impl_scalar_type;
    typedef int local_ordinal_t;
    typedef Tpetra::Map<>::global_ordinal_type global_ordinal_t;
    typedef Tpetra::Map<>::node_type node_t;

    typedef Epetra_VbrMatrix matrix_type;
    typedef matrix_type local_matrix_t;
    typedef int* sparse_ptr_type;
    typedef int* sparse_idx_type;
    typedef double* sparse_values_type;

    typedef row_access major_access;
  };

#endif

}

#endif  // AMESOS2_MATRIXTRAITS_HPP
