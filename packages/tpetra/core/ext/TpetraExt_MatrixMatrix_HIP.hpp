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
  static inline void mult_A_B_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                       const LocalOrdinalViewType& Acol2Brow,
                                                       const LocalOrdinalViewType& Acol2Irow,
                                                       const LocalOrdinalViewType& Bcol2Ccol,
                                                       const LocalOrdinalViewType& Icol2Ccol,
                                                       CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                       Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                       const std::string& label                           = std::string(),
                                                       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static inline void mult_A_B_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                   const LocalOrdinalViewType& Acol2Brow,
                                                   const LocalOrdinalViewType& Acol2Irow,
                                                   const LocalOrdinalViewType& Bcol2Ccol,
                                                   const LocalOrdinalViewType& Icol2Ccol,
                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                   Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                   const std::string& label                           = std::string(),
                                                   const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
};

// Jacobi KernelWrappers for Partial Specialization to HIP
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers2<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode, LocalOrdinalViewType> {
  static inline void jacobi_A_B_newmatrix_kernel_wrapper(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                         const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Dinv,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                         const LocalOrdinalViewType& Acol2Brow,
                                                         const LocalOrdinalViewType& Acol2Irow,
                                                         const LocalOrdinalViewType& Bcol2Ccol,
                                                         const LocalOrdinalViewType& Icol2Ccol,
                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                         Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                         const std::string& label                           = std::string(),
                                                         const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static inline void jacobi_A_B_reuse_kernel_wrapper(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                     const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Dinv,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                     const LocalOrdinalViewType& Acol2Brow,
                                                     const LocalOrdinalViewType& Acol2Irow,
                                                     const LocalOrdinalViewType& Bcol2Ccol,
                                                     const LocalOrdinalViewType& Icol2Ccol,
                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                     Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                     const std::string& label                           = std::string(),
                                                     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static inline void jacobi_A_B_newmatrix_KokkosKernels(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                        const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Dinv,
                                                        CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                        CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                        const LocalOrdinalViewType& Acol2Brow,
                                                        const LocalOrdinalViewType& Acol2Irow,
                                                        const LocalOrdinalViewType& Bcol2Ccol,
                                                        const LocalOrdinalViewType& Icol2Ccol,
                                                        CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                        Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                        const std::string& label                           = std::string(),
                                                        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
};

/*********************************************************************************************************/
// AB NewMatrix Kernel wrappers (KokkosKernels/HIP Version)
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode, LocalOrdinalViewType>::mult_A_B_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                                                                                                                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                                                                                                                              const LocalOrdinalViewType& Acol2Brow,
                                                                                                                                                              const LocalOrdinalViewType& Acol2Irow,
                                                                                                                                                              const LocalOrdinalViewType& Bcol2Ccol,
                                                                                                                                                              const LocalOrdinalViewType& Icol2Ccol,
                                                                                                                                                              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                                                                                                                              Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                                                                                                                              const std::string& label,
                                                                                                                                                              const Teuchos::RCP<Teuchos::ParameterList>& params) {
  Tpetra::MMdetails::kokkos_kernels_mult_A_B_newmatrix(
      Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
}

/*********************************************************************************************************/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode, LocalOrdinalViewType>::mult_A_B_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                                                                                                                          CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                                                                                                                          const LocalOrdinalViewType& targetMapToOrigRow_dev,
                                                                                                                                                          const LocalOrdinalViewType& targetMapToImportRow_dev,
                                                                                                                                                          const LocalOrdinalViewType& Bcol2Ccol_dev,
                                                                                                                                                          const LocalOrdinalViewType& Icol2Ccol_dev,
                                                                                                                                                          CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                                                                                                                          Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                                                                                                                          const std::string& label,
                                                                                                                                                          const Teuchos::RCP<Teuchos::ParameterList>& params) {
  Tpetra::MMdetails::host_mult_A_B_reuse(
      Aview, Bview, targetMapToOrigRow_dev, targetMapToImportRow_dev,
      Bcol2Ccol_dev, Icol2Ccol_dev, C, Cimport, label, params);
}

/*********************************************************************************************************/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers2<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode, LocalOrdinalViewType>::jacobi_A_B_newmatrix_kernel_wrapper(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                                                                                                                                 const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Dinv,
                                                                                                                                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                                                                                                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                                                                                                                                 const LocalOrdinalViewType& Acol2Brow,
                                                                                                                                                                 const LocalOrdinalViewType& Acol2Irow,
                                                                                                                                                                 const LocalOrdinalViewType& Bcol2Ccol,
                                                                                                                                                                 const LocalOrdinalViewType& Icol2Ccol,
                                                                                                                                                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                                                                                                                                 Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                                                                                                                                 const std::string& label,
                                                                                                                                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {
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
    jacobi_A_B_newmatrix_KokkosKernels(omega, Dinv, Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C, Cimport, label, params);
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
    RCP<const Export<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > dummyExport;
    C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport, dummyExport, labelList);
  }
}

/*********************************************************************************************************/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers2<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode, LocalOrdinalViewType>::jacobi_A_B_reuse_kernel_wrapper(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                                                                                                                             const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Dinv,
                                                                                                                                                             CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                                                                                                                             CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                                                                                                                             const LocalOrdinalViewType& targetMapToOrigRow_dev,
                                                                                                                                                             const LocalOrdinalViewType& targetMapToImportRow_dev,
                                                                                                                                                             const LocalOrdinalViewType& Bcol2Ccol_dev,
                                                                                                                                                             const LocalOrdinalViewType& Icol2Ccol_dev,
                                                                                                                                                             CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                                                                                                                             Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                                                                                                                             const std::string& label,
                                                                                                                                                             const Teuchos::RCP<Teuchos::ParameterList>& params) {
  host_jacobi_A_B_reuse(omega, Dinv, Aview, Bview, targetMapToOrigRow_dev, targetMapToImportRow_dev, Bcol2Ccol_dev, Icol2Ccol_dev, C, Cimport, label, params);
}

/*********************************************************************************************************/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class LocalOrdinalViewType>
void KernelWrappers2<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode, LocalOrdinalViewType>::jacobi_A_B_newmatrix_KokkosKernels(typename Teuchos::ScalarTraits<Scalar>::magnitudeType omega,
                                                                                                                                                                const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Dinv,
                                                                                                                                                                CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Aview,
                                                                                                                                                                CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& Bview,
                                                                                                                                                                const LocalOrdinalViewType& Acol2Brow,
                                                                                                                                                                const LocalOrdinalViewType& Acol2Irow,
                                                                                                                                                                const LocalOrdinalViewType& Bcol2Ccol,
                                                                                                                                                                const LocalOrdinalViewType& Icol2Ccol,
                                                                                                                                                                CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>& C,
                                                                                                                                                                Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > Cimport,
                                                                                                                                                                const std::string& label,
                                                                                                                                                                const Teuchos::RCP<Teuchos::ParameterList>& params) {
  // Check if the diagonal entries exist in debug mode
  const bool debug = Tpetra::Details::Behavior::debug();
  if (debug) {
    auto rowMap = Aview.origMatrix->getRowMap();
    Tpetra::Vector<Scalar> diags(rowMap);
    Aview.origMatrix->getLocalDiagCopy(diags);
    size_t diagLength = rowMap->getLocalNumElements();
    Teuchos::Array<Scalar> diagonal(diagLength);
    diags.get1dCopy(diagonal());

    for (size_t i = 0; i < diagLength; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(diagonal[i] == Teuchos::ScalarTraits<Scalar>::zero(),
                                 std::runtime_error,
                                 "Matrix A has a zero/missing diagonal: " << diagonal[i] << std::endl
                                                                          << "KokkosKernels Jacobi-fused SpGEMM requires nonzero diagonal entries in A" << std::endl);
    }
  }

  using Teuchos::RCP;
  using Teuchos::rcp;
  RCP<Tpetra::Details::ProfilingRegion> MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: Jacobi: Newmatrix KokkosKernels"));

  // Usings
  using device_t       = typename Tpetra::KokkosCompat::KokkosHIPWrapperNode::device_type;
  using matrix_t       = typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode>::local_matrix_device_type;
  using graph_t        = typename matrix_t::StaticCrsGraphType;
  using lno_view_t     = typename graph_t::row_map_type::non_const_type;
  using int_view_t     = Kokkos::View<int*,
                                  typename lno_view_t::array_layout,
                                  typename lno_view_t::memory_space,
                                  typename lno_view_t::memory_traits>;
  using c_lno_view_t   = typename graph_t::row_map_type::const_type;
  using lno_nnz_view_t = typename graph_t::entries_type::non_const_type;
  using scalar_view_t  = typename matrix_t::values_type::non_const_type;

  // KokkosKernels handle
  using handle_t = typename KokkosKernels::Experimental::KokkosKernelsHandle<
      typename lno_view_t::const_value_type, typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type,
      typename device_t::execution_space, typename device_t::memory_space, typename device_t::memory_space>;
  using int_handle_t = typename KokkosKernels::Experimental::KokkosKernelsHandle<
      typename int_view_t::const_value_type, typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type,
      typename device_t::execution_space, typename device_t::memory_space, typename device_t::memory_space>;

  // Merge the B and Bimport matrices
  const matrix_t Bmerged = Tpetra::MMdetails::merge_matrices(Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C.getColMap()->getLocalNumElements());

  // Get the properties and arrays of input matrices
  const matrix_t Amat = Aview.origMatrix->getLocalMatrixDevice();
  const matrix_t Bmat = Bview.origMatrix->getLocalMatrixDevice();

  typename handle_t::nnz_lno_t AnumRows = Amat.numRows();
  typename handle_t::nnz_lno_t BnumRows = Bmerged.numRows();
  typename handle_t::nnz_lno_t BnumCols = Bmerged.numCols();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmerged.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmerged.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmerged.values;

  // Arrays of the output matrix
  lno_view_t row_mapC(Kokkos::ViewAllocateWithoutInitializing("row_mapC"), AnumRows + 1);
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  // Options
  int team_work_size = 16;
  std::string myalg("SPGEMM_DEFAULT");
  if (!params.is_null()) {
    if (params->isParameter("hip: algorithm"))
      myalg = params->get("hip: algorithm", myalg);
    if (params->isParameter("hip: team work size"))
      team_work_size = params->get("hip: team work size", team_work_size);
  }

  // Get the algorithm mode
  std::string nodename("HIP");
  std::string alg = nodename + std::string(" algorithm");
  if (!params.is_null() && params->isParameter(alg)) myalg = params->get(alg, myalg);
  KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(myalg);

  Tpetra::Details::IntRowPtrHelper<decltype(Bmerged)> irph(Bmerged.nnz(), Bmerged.graph.row_map);
  const bool useIntRowptrs =
      irph.shouldUseIntRowptrs() &&
      Aview.origMatrix->getApplyHelper()->shouldUseIntRowptrs();

  const Scalar jacobiOmega = omega * Teuchos::ScalarTraits<Scalar>::one();

  if (useIntRowptrs) {
    int_handle_t kh;
    kh.create_spgemm_handle(alg_enum);
    kh.set_team_work_size(team_work_size);

    int_view_t int_row_mapC(Kokkos::ViewAllocateWithoutInitializing("int_row_mapC"), AnumRows + 1);

    auto Aint = Aview.origMatrix->getApplyHelper()->getIntRowptrMatrix(Amat);
    auto Bint = irph.getIntRowptrMatrix(Bmerged);
    {
      Tpetra::Details::ProfilingRegion MM2("TpetraExt: Jacobi: Newmatrix KokkosKernels symbolic int");

      KokkosSparse::spgemm_symbolic(&kh, AnumRows, BnumRows, BnumCols,
                                    Aint.graph.row_map, Aint.graph.entries, false,
                                    Bint.graph.row_map, Bint.graph.entries, false,
                                    int_row_mapC);
    }

    Tpetra::Details::ProfilingRegion MM2("TpetraExt: Jacobi: Newmatrix KokkosKernels numeric int");

    size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size) {
      entriesC = lno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC  = scalar_view_t(Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }

    if (c_nnz_size) {
      // even though there is no TPL for this, we have to use the same handle that was used in the symbolic phase,
      // so need to have a special int-typed call for this as well.
      KokkosSparse::Experimental::spgemm_jacobi(&kh, AnumRows, BnumRows, BnumCols,
                                                Aint.graph.row_map, Aint.graph.entries, Amat.values, false,
                                                Bint.graph.row_map, Bint.graph.entries, Bint.values, false,
                                                int_row_mapC, entriesC, valuesC,
                                                jacobiOmega, Dinv.getLocalViewDevice(Access::ReadOnly));
    }
    // transfer the integer rowptrs back to the correct rowptr type
    Kokkos::parallel_for(
        int_row_mapC.size(), KOKKOS_LAMBDA(int i) { row_mapC(i) = int_row_mapC(i); });
    kh.destroy_spgemm_handle();
  } else {
    handle_t kh;
    kh.create_spgemm_handle(alg_enum);
    kh.set_team_work_size(team_work_size);

    {
      Tpetra::Details::ProfilingRegion MM2("TpetraExt: Jacobi: Newmatrix KokkosKernels symbolic non-int");

      KokkosSparse::spgemm_symbolic(&kh, AnumRows, BnumRows, BnumCols,
                                    Amat.graph.row_map, Amat.graph.entries, false,
                                    Bmerged.graph.row_map, Bmerged.graph.entries, false,
                                    row_mapC);
    }

    Tpetra::Details::ProfilingRegion MM2("TpetraExt: Jacobi: Newmatrix KokkosKernels numeric non-int");

    size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size) {
      entriesC = lno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC  = scalar_view_t(Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
      KokkosSparse::Experimental::spgemm_jacobi(&kh, AnumRows, BnumRows, BnumCols,
                                                Amat.graph.row_map, Amat.graph.entries, Amat.values, false,
                                                Bmerged.graph.row_map, Bmerged.graph.entries, Bmerged.values, false,
                                                row_mapC, entriesC, valuesC,
                                                jacobiOmega, Dinv.getLocalViewDevice(Access::ReadOnly));
    }
    kh.destroy_spgemm_handle();
  }

  MM = Teuchos::null;
  MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: Jacobi: Newmatrix HIPSort"));

  // Sort & set values
  if (params.is_null() || params->get("sort entries", true))
    Import_Util::sortCrsEntries(row_mapC, entriesC, valuesC);
  C.setAllValues(row_mapC, entriesC, valuesC);

  MM = Teuchos::null;
  MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: Jacobi: Newmatrix HIPESFC"));

  // Final Fillcomplete
  Teuchos::RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
  labelList->set("Timer Label", label);
  if (!params.is_null()) labelList->set("compute global constants", params->get("compute global constants", true));
  Teuchos::RCP<const Export<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> > dummyExport;
  C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport, dummyExport, labelList);
}

}  // namespace MMdetails
}  // namespace Tpetra

#endif  // HIP

#endif
