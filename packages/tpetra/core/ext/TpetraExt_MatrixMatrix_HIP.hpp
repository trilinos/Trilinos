// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_HIP_DEF_HPP
#define TPETRA_MATRIXMATRIX_HIP_DEF_HPP

#include "Tpetra_Details_IntRowPtrHelper.hpp"

#ifdef HAVE_TPETRA_INST_HIP
namespace Tpetra {
namespace MMdetails {

template <>
struct KokkosKernelsSPGEMMBackend<Tpetra::KokkosCompat::KokkosHIPWrapperNode> {
  static std::string parameter_prefix() { return "hip"; }
  static std::string algorithm_label() { return "HIP"; }

  template <class MatrixType>
  static void pre_spgemm(MatrixType&) {}
};

/*********************************************************************************************************/
// MMM KernelWrappers for Partial Specialization to HIP
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
struct KernelWrappers<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode, LocalOrdinalViewType> {
  using Node = Tpetra::KokkosCompat::KokkosHIPWrapperNode;

  static inline void mult_A_B_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                       const LocalOrdinalViewType& Acol2Brow,
                                                       const LocalOrdinalViewType& Acol2Irow,
                                                       const LocalOrdinalViewType& Bcol2Ccol,
                                                       const LocalOrdinalViewType& Icol2Ccol,
                                                       CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                       Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > Cimport,
                                                       const std::string& label                           = std::string(),
                                                       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
    Tpetra::MMdetails::kokkos_kernels_mult_A_B_newmatrix(
        Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
  }

  static inline void mult_A_B_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                   const LocalOrdinalViewType& Acol2Brow,
                                                   const LocalOrdinalViewType& Acol2Irow,
                                                   const LocalOrdinalViewType& Bcol2Ccol,
                                                   const LocalOrdinalViewType& Icol2Ccol,
                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                   Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > Cimport,
                                                   const std::string& label                           = std::string(),
                                                   const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
    Tpetra::MMdetails::host_mult_A_B_reuse(
        Aview, Bview, Acol2Brow, Acol2Irow,
        Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
  }
};

// Jacobi KernelWrappers for Partial Specialization to HIP
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalOrdinalViewType> {
  using Node = Tpetra::KokkosCompat::KokkosHIPWrapperNode;

  static inline void jacobi_A_B_newmatrix_kernel_wrapper(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                         const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Dinv,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                         const LocalOrdinalViewType& Acol2Brow,
                                                         const LocalOrdinalViewType& Acol2Irow,
                                                         const LocalOrdinalViewType& Bcol2Ccol,
                                                         const LocalOrdinalViewType& Icol2Ccol,
                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                         Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > Cimport,
                                                         const std::string& label                           = std::string(),
                                                         const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
    // Node-specific code
    using Teuchos::RCP;

    // Options
    // int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer // unreferenced
    std::string myalg("KK");
    if (!params.is_null()) {
      if (params->isParameter("hip: jacobi algorithm"))
        myalg = params->get("hip: jacobi algorithm", myalg);
    }

    if (myalg == "MSAK") {
      ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_newmatrix_MultiplyScaleAddKernel(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else if (myalg == "KK") {
      kokkos_kernels_jacobi_A_B_newmatrix(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else {
      throw std::runtime_error("Tpetra::MatrixMatrix::Jacobi newmatrix unknown kernel");
    }

    Tpetra::Details::ProfilingRegion("TpetraExt: Jacobi: Newmatrix HIPESFC");

    // Final Fillcomplete
    RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
    labelList->set("Timer Label", label);
    if (!params.is_null()) labelList->set("compute global constants", params->get("compute global constants", true));

    // NOTE: MSAK already fillCompletes, so we have to check here
    if (!C.isFillComplete()) {
      RCP<const Export<LocalOrdinal, GlobalOrdinal, Node> > dummyExport;
      C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport, dummyExport, labelList);
    }
  }

  static inline void jacobi_A_B_reuse_kernel_wrapper(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                     const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Dinv,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                     const LocalOrdinalViewType& Acol2Brow,
                                                     const LocalOrdinalViewType& Acol2Irow,
                                                     const LocalOrdinalViewType& Bcol2Ccol,
                                                     const LocalOrdinalViewType& Icol2Ccol,
                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                     Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > Cimport,
                                                     const std::string& label                           = std::string(),
                                                     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
    host_jacobi_A_B_reuse(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow,
                          Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
  }

}  // namespace MMdetails
}  // namespace Tpetra

#endif  // HIP

#endif
