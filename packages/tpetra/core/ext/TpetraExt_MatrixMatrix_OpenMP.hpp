// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_OPENMP_DEF_HPP
#define TPETRA_MATRIXMATRIX_OPENMP_DEF_HPP

#ifdef HAVE_TPETRA_INST_OPENMP
namespace Tpetra {
namespace MMdetails {

template <>
struct KokkosKernelsSPGEMMBackend<Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> {
  static std::string parameter_prefix() { return "openmp"; }
  static std::string algorithm_label() { return "OpenMP"; }

  template <class MatrixType>
  static void pre_spgemm(MatrixType&) {}
};

/*********************************************************************************************************/
// MMM KernelWrappers for Partial Specialization to OpenMP
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode, LocalOrdinalViewType> {
  static inline void mult_A_B_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                       const LocalOrdinalViewType& Acol2Brow,
                                                       const LocalOrdinalViewType& Acol2Irow,
                                                       const LocalOrdinalViewType& Bcol2Ccol,
                                                       const LocalOrdinalViewType& Icol2Ccol,
                                                       CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                       Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                       const std::string& label                           = std::string(),
                                                       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
    std::string myalg("SPGEMM_DEFAULT");

    if (!params.is_null()) {
      if (params->isParameter("openmp: algorithm"))
        myalg = params->get("openmp: algorithm", myalg);
    }

    if (myalg == "LTG") {
      // Use the LTG kernel if requested
      ::Tpetra::MatrixMatrix::ExtraKernels::mult_A_B_newmatrix_LowThreadGustavsonKernel(Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else {
      Tpetra::MMdetails::kokkos_kernels_mult_A_B_newmatrix(
          Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    }

#if 0
  {
    Teuchos::ArrayRCP< const size_t > Crowptr;
    Teuchos::ArrayRCP< const LocalOrdinal > Ccolind;
    Teuchos::ArrayRCP< const Scalar > Cvalues;
    C.getAllValues(Crowptr,Ccolind,Cvalues);

    // DEBUG
    int MyPID = C->getComm()->getRank();
    printf("[%d] Crowptr = ",MyPID);
    for(size_t i=0; i<(size_t) Crowptr.size(); i++) {
      printf("%3d ",(int)Crowptr.getConst()[i]);
    }
    printf("\n");
    printf("[%d] Ccolind = ",MyPID);
    for(size_t i=0; i<(size_t)Ccolind.size(); i++) {
      printf("%3d ",(int)Ccolind.getConst()[i]);
    }
    printf("\n");
    fflush(stdout);
    // END DEBUG
  }
#endif
  }

  static inline void mult_A_B_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                   const LocalOrdinalViewType& Acol2Brow,
                                                   const LocalOrdinalViewType& Acol2Irow,
                                                   const LocalOrdinalViewType& Bcol2Ccol,
                                                   const LocalOrdinalViewType& Icol2Ccol,
                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                   Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                   const std::string& label                           = std::string(),
                                                   const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
    using Teuchos::TimeMonitor;
    Teuchos::RCP<TimeMonitor> MM;
#endif

    // Lots and lots of typedefs
    using Teuchos::RCP;

    // Options
    int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer
    std::string myalg("LTG");
    if (!params.is_null()) {
      if (params->isParameter("openmp: algorithm"))
        myalg = params->get("openmp: algorithm", myalg);
      if (params->isParameter("openmp: team work size"))
        team_work_size = params->get("openmp: team work size", team_work_size);
    }

    if (myalg == "LTG") {
      // Use the LTG kernel if requested
      ::Tpetra::MatrixMatrix::ExtraKernels::mult_A_B_reuse_LowThreadGustavsonKernel(Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else {
      throw std::runtime_error("Tpetra::MatrixMatrix::MMM reuse unknown kernel");
    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Reuse OpenMPESFC"))));
#endif
    C.fillComplete(C.getDomainMap(), C.getRangeMap());
  }
};

// Jacobi KernelWrappers for Partial Specialization to OpenMP
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers2<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode, LocalOrdinalViewType> {
  static inline void jacobi_A_B_newmatrix_kernel_wrapper(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                         const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Dinv,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                         const LocalOrdinalViewType& Acol2Brow,
                                                         const LocalOrdinalViewType& Acol2Irow,
                                                         const LocalOrdinalViewType& Bcol2Ccol,
                                                         const LocalOrdinalViewType& Icol2Ccol,
                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                         Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
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
    int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer
    std::string myalg("LTG");
    if (!params.is_null()) {
      if (params->isParameter("openmp: jacobi algorithm"))
        myalg = params->get("openmp: jacobi algorithm", myalg);
      if (params->isParameter("openmp: team work size"))
        team_work_size = params->get("openmp: team work size", team_work_size);
    }

    if (myalg == "LTG") {
      // Use the LTG kernel if requested
      ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_newmatrix_LowThreadGustavsonKernel(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else if (myalg == "MSAK") {
      ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_newmatrix_MultiplyScaleAddKernel(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else if (myalg == "KK") {
      kokkos_kernels_jacobi_A_B_newmatrix(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else {
      throw std::runtime_error("Tpetra::MatrixMatrix::Jacobi newmatrix unknown kernel");
    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix OpenMPESFC"))));
#endif

    // Final Fillcomplete
    RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
    labelList->set("Timer Label", label);
    if (!params.is_null()) labelList->set("compute global constants", params->get("compute global constants", true));

    // NOTE: MSAK already fillCompletes, so we have to check here
    if (!C.isFillComplete()) {
      RCP<const Export<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > dummyExport;
      C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport, dummyExport, labelList);
    }
  }

  static inline void jacobi_A_B_reuse_kernel_wrapper(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                     const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Dinv,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                     const LocalOrdinalViewType& Acol2Brow,
                                                     const LocalOrdinalViewType& Acol2Irow,
                                                     const LocalOrdinalViewType& Bcol2Ccol,
                                                     const LocalOrdinalViewType& Icol2Ccol,
                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                     Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                     const std::string& label                           = std::string(),
                                                     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
    using Teuchos::TimeMonitor;
    Teuchos::RCP<TimeMonitor> MM;
#endif

    // Lots and lots of typedefs
    using Teuchos::RCP;

    // Options
    int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer
    std::string myalg("LTG");
    if (!params.is_null()) {
      if (params->isParameter("openmp: jacobi algorithm"))
        myalg = params->get("openmp: jacobi algorithm", myalg);
      if (params->isParameter("openmp: team work size"))
        team_work_size = params->get("openmp: team work size", team_work_size);
    }

    if (myalg == "LTG") {
      // Use the LTG kernel if requested
      ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_reuse_LowThreadGustavsonKernel(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
    } else {
      throw std::runtime_error("Tpetra::MatrixMatrix::Jacobi reuse unknown kernel");
    }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Reuse OpenMPESFC"))));
#endif
    C.fillComplete(C.getDomainMap(), C.getRangeMap());
  }
};

// Triple-Product KernelWrappers for Partial Specialization to OpenMP
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers3<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode, LocalOrdinalViewType> {
  static inline void mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Rview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                         const LocalOrdinalViewType& Acol2Prow,
                                                         const LocalOrdinalViewType& Acol2PIrow,
                                                         const LocalOrdinalViewType& Pcol2Ccol,
                                                         const LocalOrdinalViewType& PIcol2Ccol,
                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                         Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                         const std::string& label                           = std::string(),
                                                         const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static inline void mult_R_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Rview,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                     const LocalOrdinalViewType& Acol2Prow,
                                                     const LocalOrdinalViewType& Acol2PIrow,
                                                     const LocalOrdinalViewType& Pcol2Ccol,
                                                     const LocalOrdinalViewType& PIcol2Ccol,
                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                     Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                     const std::string& label                           = std::string(),
                                                     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static inline void mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                          CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                          const LocalOrdinalViewType& Acol2Prow,
                                                          const LocalOrdinalViewType& Acol2PIrow,
                                                          const LocalOrdinalViewType& Pcol2Ccol,
                                                          const LocalOrdinalViewType& PIcol2Ccol,
                                                          CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                          Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                          const std::string& label                           = std::string(),
                                                          const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static inline void mult_PT_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                      CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                      const LocalOrdinalViewType& Acol2Prow,
                                                      const LocalOrdinalViewType& Acol2PIrow,
                                                      const LocalOrdinalViewType& Pcol2Ccol,
                                                      const LocalOrdinalViewType& PIcol2Ccol,
                                                      CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                      Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                      const std::string& label                           = std::string(),
                                                      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
};

/*********************************************************************************************************/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers3<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode, LocalOrdinalViewType>::mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Rview,
                                                                                                                                                                    CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                                                                                                                                    CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                                                                                                                                    const LocalOrdinalViewType& Acol2Prow,
                                                                                                                                                                    const LocalOrdinalViewType& Acol2PIrow,
                                                                                                                                                                    const LocalOrdinalViewType& Pcol2Accol,
                                                                                                                                                                    const LocalOrdinalViewType& PIcol2Accol,
                                                                                                                                                                    CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                                                                                                                                    Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                                                                                                                                    const std::string& label,
                                                                                                                                                                    const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM;
#endif

  // Node-specific code
  std::string nodename("OpenMP");

  // Options
  std::string myalg("LTG");

  if (!params.is_null()) {
    if (params->isParameter("openmp: rap algorithm"))
      myalg = params->get("openmp: rap algorithm", myalg);
  }

  if (myalg == "LTG") {
    // Use the LTG kernel if requested
    ::Tpetra::MatrixMatrix::ExtraKernels::mult_R_A_P_newmatrix_LowThreadGustavsonKernel(Rview, Aview, Pview, Acol2Prow, Acol2PIrow, Pcol2Accol, PIcol2Accol, Ac, Acimport, label, params);
  } else {
    throw std::runtime_error("Tpetra::MatrixMatrix::R_A_P newmatrix unknown kernel");
  }
}

/*********************************************************************************************************/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers3<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode, LocalOrdinalViewType>::mult_R_A_P_reuse_kernel_wrapper(
    CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Rview,
    CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
    CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,

    const LocalOrdinalViewType& Acol2Prow,
    const LocalOrdinalViewType& Acol2Irow,
    const LocalOrdinalViewType& Pcol2Ccol,
    const LocalOrdinalViewType& Icol2Ccol,
    CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
    Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
    const std::string& label,
    const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM;
#endif

  // Lots and lots of typedefs
  using Teuchos::RCP;

  // Options
  std::string myalg("LTG");
  if (!params.is_null()) {
    if (params->isParameter("openmp: rap algorithm"))
      myalg = params->get("openmp: rap algorithm", myalg);
  }

  if (myalg == "LTG") {
    // Use the LTG kernel if requested
    ::Tpetra::MatrixMatrix::ExtraKernels::mult_R_A_P_reuse_LowThreadGustavsonKernel(Rview, Aview, Pview, Acol2Prow, Acol2Irow, Pcol2Ccol, Icol2Ccol, C, Cimport, label, params);
  } else {
    throw std::runtime_error("Tpetra::MatrixMatrix::R_A_P newmatrix unknown kernel");
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Reuse OpenMPESFC"))));
#endif
  C.fillComplete(C.getDomainMap(), C.getRangeMap());
}

/*********************************************************************************************************/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers3<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode, LocalOrdinalViewType>::mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,

                                                                                                                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                                                                                                                                     const LocalOrdinalViewType& Acol2Prow,
                                                                                                                                                                     const LocalOrdinalViewType& Acol2PIrow,
                                                                                                                                                                     const LocalOrdinalViewType& Pcol2Accol,
                                                                                                                                                                     const LocalOrdinalViewType& PIcol2Accol,
                                                                                                                                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                                                                                                                                     Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                                                                                                                                     const std::string& label,
                                                                                                                                                                     const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM;
#endif

  // Node-specific code
  std::string nodename("OpenMP");

  // Options
  std::string myalg("LTG");

  if (!params.is_null()) {
    if (params->isParameter("openmp: ptap algorithm"))
      myalg = params->get("openmp: ptap algorithm", myalg);
  }

  if (myalg == "LTG") {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose"))));
#endif

    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using SC = Scalar;

    // We don't need a kernel-level PTAP, we just transpose here
    using transposer_type =
        RowMatrixTransposer<SC, LO, GO, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
    transposer_type transposer(Pview.origMatrix, label + "XP: ");
    RCP<ParameterList> transposeParams(new ParameterList);
    if (!params.is_null()) {
      transposeParams->set("compute global constants",
                           params->get("compute global constants: temporaries",
                                       false));
    }
    transposeParams->set("sort", false);
    RCP<CrsMatrix<SC, LO, GO, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Ptrans =
        transposer.createTransposeLocal(transposeParams);
    CrsMatrixStruct<SC, LO, GO, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> Rview;
    Rview.origMatrix = Ptrans;

    using ::Tpetra::MatrixMatrix::ExtraKernels::mult_R_A_P_newmatrix_LowThreadGustavsonKernel;
    mult_R_A_P_newmatrix_LowThreadGustavsonKernel(Rview, Aview, Pview, Acol2Prow, Acol2PIrow, Pcol2Accol,
                                                  PIcol2Accol, Ac, Acimport, label, params);
  } else {
    throw std::runtime_error("Tpetra::MatrixMatrix::PT_A_P newmatrix unknown kernel");
  }
}

/*********************************************************************************************************/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers3<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode, LocalOrdinalViewType>::mult_PT_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,

                                                                                                                                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                                                                                                                                 const LocalOrdinalViewType& Acol2Prow,
                                                                                                                                                                 const LocalOrdinalViewType& Acol2PIrow,
                                                                                                                                                                 const LocalOrdinalViewType& Pcol2Accol,
                                                                                                                                                                 const LocalOrdinalViewType& PIcol2Accol,
                                                                                                                                                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                                                                                                                                 Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                                                                                                                                 const std::string& label,
                                                                                                                                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM;
#endif

  // Node-specific code
  std::string nodename("OpenMP");

  // Options
  std::string myalg("LTG");

  if (!params.is_null()) {
    if (params->isParameter("openmp: ptap algorithm"))
      myalg = params->get("openmp: ptap algorithm", myalg);
  }

  if (myalg == "LTG") {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose"))));
#endif

    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using SC = Scalar;

    // We don't need a kernel-level PTAP, we just transpose here
    using transposer_type =
        RowMatrixTransposer<SC, LO, GO, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>;
    transposer_type transposer(Pview.origMatrix, label + "XP: ");
    RCP<ParameterList> transposeParams(new ParameterList);
    if (!params.is_null()) {
      transposeParams->set("compute global constants",
                           params->get("compute global constants: temporaries",
                                       false));
    }
    transposeParams->set("sort", false);
    RCP<CrsMatrix<SC, LO, GO, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Ptrans =
        transposer.createTransposeLocal(transposeParams);
    CrsMatrixStruct<SC, LO, GO, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> Rview;
    Rview.origMatrix = Ptrans;

    using ::Tpetra::MatrixMatrix::ExtraKernels::mult_R_A_P_reuse_LowThreadGustavsonKernel;
    mult_R_A_P_reuse_LowThreadGustavsonKernel(Rview, Aview, Pview, Acol2Prow, Acol2PIrow, Pcol2Accol,
                                              PIcol2Accol, Ac, Acimport, label, params);
  } else {
    throw std::runtime_error("Tpetra::MatrixMatrix::PT_A_P reuse unknown kernel");
  }
  Ac.fillComplete(Ac.getDomainMap(), Ac.getRangeMap());
}

}  // namespace MMdetails
}  // namespace Tpetra

#endif  // OpenMP

#endif
