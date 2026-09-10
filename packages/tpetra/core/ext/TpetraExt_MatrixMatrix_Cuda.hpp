// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_CUDA_DEF_HPP
#define TPETRA_MATRIXMATRIX_CUDA_DEF_HPP

#include "Tpetra_Details_IntRowPtrHelper.hpp"

#ifdef HAVE_TPETRA_INST_CUDA
namespace Tpetra {
namespace MMdetails {

template <>
struct KokkosKernelsSPGEMMBackend<Tpetra::KokkosCompat::KokkosCudaWrapperNode> {
  static std::string parameter_prefix() { return "cuda"; }
  static std::string algorithm_label() { return "Cuda"; }

  template <class MatrixType>
  static void pre_spgemm(MatrixType& Bmerged) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && ((CUDA_VERSION < 11000) || (CUDA_VERSION >= 11040))
    using device_t = typename MatrixType::device_type;
    if constexpr (std::is_same_v<typename device_t::execution_space, Kokkos::Cuda>) {
      if (!KokkosSparse::isCrsGraphSorted(Bmerged.graph.row_map, Bmerged.graph.entries)) {
        Import_Util::sortCrsMatrix(Bmerged);
      }
    }
#else
    (void)Bmerged;
#endif
  }
};

/*********************************************************************************************************/
// MMM KernelWrappers for Partial Specialization to CUDA
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
struct KernelWrappers<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosCudaWrapperNode, LocalOrdinalViewType> {
  using Node = Tpetra::KokkosCompat::KokkosCudaWrapperNode;

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

// Jacobi KernelWrappers for Partial Specialization to Cuda
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers2<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosCudaWrapperNode, LocalOrdinalViewType> {
  using Node = Tpetra::KokkosCompat::KokkosCudaWrapperNode;

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
    using Teuchos::rcp;
    RCP<Tpetra::Details::ProfilingRegion> MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: MMM: Jacobi CudaWrapper"));

    // Options
    // int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer // unreferenced
    std::string myalg("KK");
    if (!params.is_null()) {
      if (params->isParameter("cuda: jacobi algorithm"))
        myalg = params->get("cuda: jacobi algorithm", myalg);
    }

    if (myalg == "MSAK") {
      ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_newmatrix_MultiplyScaleAddKernel(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else if (myalg == "KK") {
      kokkos_kernels_jacobi_A_B_newmatrix(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else {
      throw std::runtime_error("Tpetra::MatrixMatrix::Jacobi newmatrix unknown kernel");
    }

    MM = Teuchos::null;
    MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: Jacobi: Newmatrix CudaESFC"));

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
    host_jacobi_A_B_reuse(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
  }
};

}  // namespace MMdetails
}  // namespace Tpetra

#endif  // CUDA

#endif
