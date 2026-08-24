// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_SYCL_DEF_HPP
#define TPETRA_MATRIXMATRIX_SYCL_DEF_HPP

#include "Tpetra_Details_IntRowPtrHelper.hpp"

#ifdef HAVE_TPETRA_INST_SYCL
namespace Tpetra {
namespace MMdetails {

template <>
struct KokkosKernelsSPGEMMBackend<Tpetra::KokkosCompat::KokkosSYCLWrapperNode> {
  static std::string parameter_prefix() { return "sycl"; }
  static std::string algorithm_label() { return "SYCL"; }

  template <class MatrixType>
  static void pre_spgemm(MatrixType&) {}
};

/*********************************************************************************************************/
// MMM KernelWrappers for Partial Specialization to SYCL
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
struct KernelWrappers<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSYCLWrapperNode, LocalOrdinalViewType> {
  using Node = Tpetra::KokkosCompat::KokkosSYCLWrapperNode;
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
        Aview, Bview, targetMapToOrigRow_dev, targetMapToImportRow_dev,
        Bcol2Ccol_dev, Icol2Ccol_dev, C, Cimport, label, params);
  }
};

// Jacobi KernelWrappers for Partial Specialization to SYCL
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers2<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSYCLWrapperNode, LocalOrdinalViewType> {
  using Node = Tpetra::KokkosCompat::KokkosSYCLWrapperNode;
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
#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
    using Teuchos::TimeMonitor;
    Teuchos::RCP<TimeMonitor> MM;
#endif

    // Node-specific code
    using Teuchos::RCP;

    // Options
    // int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer // unreferenced
    std::string myalg("KK");
    if (!params.is_null()) {
      if (params->isParameter("sycl: jacobi algorithm"))
        myalg = params->get("sycl: jacobi algorithm", myalg);
    }

    if (myalg == "MSAK") {
      ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_newmatrix_MultiplyScaleAddKernel(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else if (myalg == "KK") {
      kokkos_kernels_jacobi_A_B_newmatrix(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else {
      throw std::runtime_error("Tpetra::MatrixMatrix::Jacobi newmatrix unknown kernel");
    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix SYCLESFC"))));
#endif

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
    host_jacobi_A_B_reuse(omega, Dinv, Aview, Bview, targetMapToOrigRow_dev, targetMapToImportRow_dev, Bcol2Ccol_dev, Icol2Ccol_dev, C, Cimport, label, params);
  }
};

}  // namespace MMdetails
}  // namespace Tpetra

#endif  // SYCL

#endif
