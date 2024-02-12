// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <algorithm>  // shuffle
#include <vector>     // vector
#include <random>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_IO.hpp>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
// MueLu
#include "MueLu.hpp"
#include "MueLu_TestHelpers.hpp"

#if defined(HAVE_TPETRA_INST_OPENMP)
#include "Xpetra_TpetraVector.hpp"
#include "TpetraExt_MatrixMatrix_def.hpp"
#endif

// =========================================================================
// Support Routines
// =========================================================================
template <class View1, class View2>
inline void copy_view_n(const int n, const View1 x1, View2 x2) {
  Kokkos::parallel_for(
      n, KOKKOS_LAMBDA(const size_t i) {
        x2[i] = x1[i];
      });
}

template <class View1, class View2>
inline void copy_view(const View1 x1, View2 x2) {
  copy_view_n(x1.extent(0), x1, x2);
}

template <class V1, class V2>
void print_crs_graph(std::string name, const V1 rowptr, const V2 colind) {
  printf("%s rowptr[%d] = ", name.c_str(), rowptr.extent(0));
  for (size_t i = 0; i < rowptr.extent(0); i++)
    printf(" %d", (int)rowptr[i]);
  printf("\n%s colind[%d] = ", name.c_str(), colind.extent(0));
  for (size_t i = 0; i < colind.extent(0); i++)
    printf(" %d", (int)colind[i]);
  printf("\n");
}

// =========================================================================
// MKL Testing
// =========================================================================
#ifdef HAVE_MUELU_MKL
#include "mkl.h"

std::string mkl_error(sparse_status_t code) {
  switch (code) {
    case SPARSE_STATUS_SUCCESS:
      return std::string("Success");
    case SPARSE_STATUS_NOT_INITIALIZED:
      return std::string("Empty handle or matrix array");
    case SPARSE_STATUS_ALLOC_FAILED:
      return std::string("Memory allocation failed");
    case SPARSE_STATUS_INVALID_VALUE:
      return std::string("Input contains an invalid value");
    case SPARSE_STATUS_EXECUTION_FAILED:
      return std::string("Execution failed");
    case SPARSE_STATUS_INTERNAL_ERROR:
      return std::string("Internal error");
    case SPARSE_STATUS_NOT_SUPPORTED:
      return std::string("Operation not supported");
  };
}

// mkl_sparse_spmm
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MM2_MKL(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B1, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B2, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C, Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > &Ccolmap, std::string algorithm_name) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  // NOTE: MKL uses int as its columnn index type and the either double or float for its Scalar type

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
  RCP<TimeMonitor> tm;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef typename crs_matrix_type::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;
  typedef Tpetra::Map<LO, GO, NO> map_type;
  typedef typename map_type::local_map_type local_map_type;

  typedef typename Kokkos::View<MKL_INT *, typename lno_nnz_view_t::array_layout, typename lno_nnz_view_t::device_type> mkl_int_type;

  RCP<const crs_matrix_type> Au  = Utilities::Op2TpetraCrs(rcp(&A, false));
  RCP<const crs_matrix_type> B1u = Utilities::Op2TpetraCrs(rcp(&B1, false));
  RCP<const crs_matrix_type> B2u = Utilities::Op2TpetraCrs(rcp(&B2, false));
  RCP<const crs_matrix_type> Cu  = Utilities::Op2TpetraCrs(rcp(&C, false));
  RCP<crs_matrix_type> Cnc       = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);

  const KCRS &Amat  = Au->getLocalMatrixDevice();
  const KCRS &B1mat = B1u->getLocalMatrixDevice();
  const KCRS &B2mat = B2u->getLocalMatrixDevice();

  if (A.getLocalNumRows() != C.getLocalNumRows()) throw std::runtime_error("C is not sized correctly");

  c_lno_view_t Arowptr = Amat.graph.row_map, B1rowptr = B1mat.graph.row_map, B2rowptr = B2mat.graph.row_map;
  lno_view_t Crowptr("Crowptr", C.getLocalNumRows() + 1);
  c_lno_nnz_view_t Acolind = Amat.graph.entries, B1colind = B1mat.graph.entries, B2colind = B2mat.graph.entries;
  lno_nnz_view_t Ccolind;
  const scalar_view_t Avals = Amat.values, B1vals = B1mat.values, B2vals = B2mat.values;
  scalar_view_t Cvals;
  RCP<const Tpetra::Map<LO, GO, Node> > Ccolmap_t = Xpetra::toTpetra(Ccolmap);
  local_map_type Bcolmap_local                    = B1u->getColMap()->getLocalMap();
  local_map_type Icolmap_local                    = B2u->getColMap()->getLocalMap();
  local_map_type Ccolmap_local                    = Ccolmap_t->getLocalMap();

  sparse_matrix_t AMKL;
  sparse_matrix_t B1MKL;
  sparse_matrix_t B2MKL;
  sparse_matrix_t CMKL;
  sparse_matrix_t Temp1MKL, Temp2MKL;

  sparse_status_t result;
  // **********************************
  // Copy in the matrix data for MKL
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM2 MKL " + algorithm_name + ": Matrix CopyIn")));

  mkl_int_type ArowptrMKL("Arowptr", Arowptr.extent(0));
  mkl_int_type AcolindMKL("Acolind", Acolind.extent(0));
  mkl_int_type B1colindMKL("B1colind", B1colind.extent(0));
  mkl_int_type B2colindMKL("B2colind", B2colind.extent(0));

  // Generate new rowptrs that have all the rows for both B and I (as MKL requires)
  size_t Nb = B1rowptr.extent(0) - 1, Ni = B2rowptr.extent(0) - 1;
  mkl_int_type B1rowptrMKL("B1rowptr", Nb + Ni + 1);
  mkl_int_type B2rowptrMKL("B2rowptr", Nb + Ni + 1);
  Kokkos::parallel_for(
      Nb + 1, KOKKOS_LAMBDA(const size_t i) {
        B1rowptrMKL[i] = B1rowptr[i];
        if (i != Nb) B2rowptrMKL[i] = 0;
      });
  Kokkos::parallel_for(
      Ni + 1, KOKKOS_LAMBDA(const size_t i) {
        if (i != 0) B1rowptrMKL[Nb + i] = B1rowptr[Nb - 1];
        B2rowptrMKL[Nb + i] = B2rowptr[i];
      });

  // Reindex the column indices into C's column map
  Kokkos::parallel_for(
      B1colind.extent(0), KOKKOS_LAMBDA(const size_t i) {
        B1colindMKL[i] = Ccolmap_local.getLocalElement(Bcolmap_local.getGlobalElement(B1colind(i)));
      });

  Kokkos::parallel_for(
      B2colind.extent(0), KOKKOS_LAMBDA(const size_t i) {
        B2colindMKL[i] = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(B2colind(i)));
      });

  // Copy in A
  copy_view(Arowptr, ArowptrMKL);
  copy_view(Acolind, AcolindMKL);

  if (std::is_same<Scalar, double>::value) {
    mkl_sparse_d_create_csr(&AMKL, SPARSE_INDEX_BASE_ZERO, Au->getLocalNumRows(), Au->getLocalNumCols(), ArowptrMKL.data(), ArowptrMKL.data() + 1, AcolindMKL.data(), (double *)Avals.data());
    mkl_sparse_d_create_csr(&B1MKL, SPARSE_INDEX_BASE_ZERO, Nb + Ni, Ccolmap->getLocalNumElements(), B1rowptrMKL.data(), B1rowptrMKL.data() + 1, B1colindMKL.data(), (double *)B1vals.data());
    mkl_sparse_d_create_csr(&B2MKL, SPARSE_INDEX_BASE_ZERO, Nb + Ni, Ccolmap->getLocalNumElements(), B2rowptrMKL.data(), B2rowptrMKL.data() + 1, B2colindMKL.data(), (double *)B2vals.data());
  } else
    throw std::runtime_error("MKL Type Mismatch");
  tm = Teuchos::null;
  Au->getComm()->barrier();

  tm = Teuchos::null;
  Au->getComm()->barrier();

  if (algorithm_name == "MULT_ADD") {
    // **********************************
    // Multiply #1 (A*B1)
    {
      TimeMonitor tm1(*TimeMonitor::getNewTimer("MM2 MKL MULT_ADD: Multiply 1"));
      result = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, AMKL, B1MKL, &Temp1MKL);
      typename KCRS::execution_space().fence();
      if (result != SPARSE_STATUS_SUCCESS) throw std::runtime_error("MKL Multiply 1 failed: " + mkl_error(result));
    }

    // **********************************
    // Multiply #2 (A*B2)
    {
      TimeMonitor tm2(*TimeMonitor::getNewTimer("MM2 MKL MULT_ADD: Multiply 2"));
      result = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, AMKL, B1MKL, &Temp2MKL);
      typename KCRS::execution_space().fence();
      if (result != SPARSE_STATUS_SUCCESS) throw std::runtime_error("MKL Multiply 2 failed: " + mkl_error(result));
    }

    // **********************************
    // Add (A*B1) + (A*B2)
    {
      TimeMonitor tm3(*TimeMonitor::getNewTimer("MM2 MKL MULT_ADD: Add"));
      result = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, Temp1MKL, 1.0, Temp2MKL, &CMKL);
      typename KCRS::execution_space().fence();
      if (result != SPARSE_STATUS_SUCCESS) throw std::runtime_error("MKL Add failed: " + mkl_error(result));
    }
  } else if (algorithm_name == "ADD_MULT") {
    // **********************************
    // Add B1 + B2
    {
      TimeMonitor tm4(*TimeMonitor::getNewTimer("MM2 MKL ADD_MULT: Add"));
      result = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, B1MKL, 1.0, B2MKL, &Temp1MKL);
      typename KCRS::execution_space().fence();
      if (result != SPARSE_STATUS_SUCCESS) throw std::runtime_error("MKL Add failed: " + mkl_error(result));
    }

    // **********************************
    // Multiply A*(B1+B2)
    {
      TimeMonitor tm5(*TimeMonitor::getNewTimer("MM2 MKL ADD_MULT: Multiply"));
      result = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, AMKL, Temp1MKL, &CMKL);
      typename KCRS::execution_space().fence();
      if (result != SPARSE_STATUS_SUCCESS) throw std::runtime_error("MKL Multiply failed: " + mkl_error(result));
    }
  } else
    throw std::runtime_error("Invalid MKL algorithm");

  // **********************************
  // Copy out the data for MKL
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM2 MKL " + algorithm_name + ": Copy Out")));
  sparse_index_base_t c_indexing;
  MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
  double *values;
  mkl_sparse_d_export_csr(CMKL, &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values);
  size_t cnnz = rows_end[c_rows - 1];
  Kokkos::resize(Ccolind, cnnz);
  Kokkos::resize(Cvals, cnnz);
  if (c_rows != A.getLocalNumRows() || c_rows + 1 != Crowptr.extent(0)) throw std::runtime_error("C row size mismatch");
  copy_view_n(c_rows, rows_start, Crowptr);
  Crowptr(c_rows) = rows_end[c_rows - 1];
  copy_view_n(cnnz, columns, Ccolind);
  copy_view_n(cnnz, values, Cvals);

  Cnc->replaceColMap(Ccolmap_t);
  Cnc->setAllValues(Crowptr,
                    Ccolind,
                    Cvals);

  mkl_sparse_destroy(AMKL);
  mkl_sparse_destroy(B1MKL);
  mkl_sparse_destroy(B2MKL);
  mkl_sparse_destroy(CMKL);
  mkl_sparse_destroy(Temp1MKL);
  if (algorithm_name == "MULT_ADD") mkl_sparse_destroy(Temp2MKL);

  tm = Teuchos::null;
  Au->getComm()->barrier();
}

#endif

// =========================================================================
// Tpetra Kernel Testing
// =========================================================================
#include "Tpetra_Import_Util2.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_decl.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_def.hpp"

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MM2_Wrapper(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B1, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B2, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C, Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > &Ccolmap, std::string algorithm_name, int team_work_size) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
  RCP<TimeMonitor> tm;

  std::string name = algorithm_name;

  if (lib == Xpetra::UseTpetra) {
#if defined(HAVE_TPETRA_INST_OPENMP)
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
    typedef Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> import_type;
    typedef typename crs_matrix_type::local_matrix_device_type KCRS;
    typedef typename KCRS::device_type device_t;
    typedef typename KCRS::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename KCRS::values_type::non_const_type scalar_view_t;
    typedef Kokkos::View<LO *, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;
    typedef Tpetra::Map<LO, GO, NO> map_type;
    typedef typename map_type::local_map_type local_map_type;
    typedef typename Node::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
    LocalOrdinal LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
    RCP<const import_type> Cimport;
    RCP<const crs_matrix_type> Au  = Utilities::Op2TpetraCrs(rcp(&A, false));
    RCP<const crs_matrix_type> B1u = Utilities::Op2TpetraCrs(rcp(&B1, false));
    RCP<crs_matrix_type> B2u       = Teuchos::rcp_const_cast<crs_matrix_type>(Utilities::Op2TpetraCrs(rcp(&B2, false)));
    RCP<const crs_matrix_type> Cu  = Utilities::Op2TpetraCrs(rcp(&C, false));
    RCP<crs_matrix_type> Cnc       = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);

    // **********************************
    // Copy in the data for Kernel Wrapper
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM2 " + name + ": CopyIn")));

    Tpetra::CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node> Aview, Bview;
    Aview.origMatrix   = Au;
    Aview.origRowMap   = Au->getRowMap();
    Aview.rowMap       = Au->getRowMap();
    Aview.colMap       = Au->getColMap();
    Aview.domainMap    = Au->getDomainMap();
    Aview.importColMap = Teuchos::null;

    Bview.origMatrix   = B1u;
    Bview.origRowMap   = B1u->getRowMap();
    Bview.rowMap       = B1u->getRowMap();
    Bview.colMap       = B1u->getColMap();
    Bview.domainMap    = B1u->getDomainMap();
    Bview.importColMap = Teuchos::null;

    // A bunch of this stuff isn't getting set "right", but we fake it for testing purposes here...
    Bview.importMatrix = B2u;
    Bview.importColMap = B2u->getColMap();

    RCP<const Tpetra::Map<LO, GO, Node> > Ccolmap_t = Xpetra::toTpetra(Ccolmap);
    Cnc->replaceColMap(Ccolmap_t);

    local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
    local_map_type Irowmap_local = Bview.importMatrix->getRowMap()->getLocalMap();

    local_map_type Acolmap_local = Aview.colMap->getLocalMap();
    local_map_type Bcolmap_local = Bview.origMatrix->getColMap()->getLocalMap();
    local_map_type Icolmap_local = Bview.importMatrix->getColMap()->getLocalMap();
    local_map_type Ccolmap_local = Ccolmap_t->getLocalMap();

    // Bcol2Ccol
    lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"), Bview.origMatrix->getColMap()->getLocalNumElements());
    lo_view_t Icol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Icol2Ccol"), Bview.importMatrix->getColMap()->getLocalNumElements());
    Kokkos::parallel_for(
        "Tpetra::mult_A_B_newmatrix::Bcol2Ccol_getGlobalElement", range_type(0, Bview.origMatrix->getColMap()->getLocalNumElements()), KOKKOS_LAMBDA(const LO i) {
          Bcol2Ccol(i) = Ccolmap_local.getLocalElement(Bcolmap_local.getGlobalElement(i));
        });
    Kokkos::parallel_for(
        "Tpetra::mult_A_B_newmatrix::Icol2Ccol_getGlobalElement", range_type(0, Bview.importMatrix->getColMap()->getLocalNumElements()), KOKKOS_LAMBDA(const LO i) {
          Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
        });

    // Acol2Brow
    lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"), Aview.colMap->getLocalNumElements());
    lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"), Aview.colMap->getLocalNumElements());
    Kokkos::parallel_for(
        "Tpetra::mult_A_B_newmatrix::construct_tables", range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex() + 1), KOKKOS_LAMBDA(const LO i) {
          GO aidx  = Acolmap_local.getGlobalElement(i);
          LO B_LID = Browmap_local.getLocalElement(aidx);
          if (B_LID != LO_INVALID) {
            targetMapToOrigRow(i)   = B_LID;
            targetMapToImportRow(i) = LO_INVALID;
          } else {
            LO I_LID                = Irowmap_local.getLocalElement(aidx);
            targetMapToOrigRow(i)   = LO_INVALID;
            targetMapToImportRow(i) = I_LID;
          }
        });
    tm = Teuchos::null;
    Au->getComm()->barrier();

    // **********************************
    // Multiply
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("openmp: algorithm", algorithm_name);
    params->set("openmp: team work size", team_work_size);

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("MM2 ") + name + ": Kernel")));
    Tpetra::MMdetails::KernelWrappers<SC, LO, GO, Node, lno_nnz_view_t>::mult_A_B_newmatrix_kernel_wrapper(Aview, Bview, targetMapToOrigRow, targetMapToImportRow, Bcol2Ccol, Icol2Ccol, *Cnc, Cimport, name, params);

    tm = Teuchos::null;
    Au->getComm()->barrier();

#endif
  }
}

// =========================================================================
// =========================================================================
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using std::endl;
  using Teuchos::RCP;  // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    //    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> pOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream &out      = *pOut;
    out.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    //    GO nx = 50, ny = 50, nz = 50;
    // Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace3D"); // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);  // manage parameters of Xpetra

    bool printTimings = true;
    clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
    int nrepeat = 100;
    clp.setOption("nrepeat", &nrepeat, "repeat the experiment N times");
    // the kernels
    bool do_mkl_ma                = true;
    bool do_mkl_am                = true;
    bool do_kk_mem                = true;
    bool do_kk_dense              = true;
    bool do_kk_default            = true;
    bool do_ltg                   = true;
    double percent_rows_to_remove = 10.0;
    clp.setOption("remove-rows", &percent_rows_to_remove, "% of rows to remove from B to make B1 & B2");

#ifndef HAVE_MUELU_MKL
    do_mkl_ma = false;
    do_mkl_am = false;
#endif
    clp.setOption("mkl_ma", "nomkl_ma", &do_mkl_ma, "Evaluate MKL Mult-Add");
    clp.setOption("mkl_am", "nomkl_am", &do_mkl_am, "Evaluate MKL Add-Mult");
    clp.setOption("kk_mem", "nokk_mem", &do_kk_mem, "Evaluate KK Mem");
    clp.setOption("kk_dense", "nokk_dense", &do_kk_dense, "Evaluate KK Dense");
    clp.setOption("kk_default", "nokk_default", &do_kk_default, "Evaluate KK Default");
    clp.setOption("ltg", "noltg", &do_ltg, "Evaluate LTG");
    int kk_team_work_size = 16;

    std::string matrixFileNameA = "A.mm";
    clp.setOption("matrixfileA", &matrixFileNameA, "matrix market file containing matrix");
    std::string matrixFileNameB = "B.mm";
    clp.setOption("matrixfileB", &matrixFileNameB, "matrix market file containing matrix");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    if (percent_rows_to_remove < 0.0) percent_rows_to_remove = 0.0;
    if (percent_rows_to_remove > 100.0) percent_rows_to_remove = 100.0;

#ifndef HAVE_MUELU_MKL
    if (do_mkl_ma || do_mkl_am) {
      out << "MKL was requested, but this kernel is not available. Disabling..." << endl;
      do_mkl_ma = false;
      do_mkl_am = false;
    }
#endif

    // simple hack to randomize order of experiments
    enum class Experiments { MKL_MULT_ADD = 0,
                             MKL_ADD_MULT,
                             KK_MEM,
                             KK_DENSE,
                             KK_DEFAULT,
                             LTG };
    std::vector<Experiments> my_experiments;
    // add the experiments we will run

#ifdef HAVE_MUELU_MKL
    if (do_mkl_ma) my_experiments.push_back(Experiments::MKL_MULT_ADD);  // MKL Mult-Add
    if (do_mkl_am) my_experiments.push_back(Experiments::MKL_ADD_MULT);  // MKL Add-Mult
#endif

    // assume these are available
    if (do_kk_mem) my_experiments.push_back(Experiments::KK_MEM);          // KK Mem
    if (do_kk_dense) my_experiments.push_back(Experiments::KK_DENSE);      // KK Dense
    if (do_kk_default) my_experiments.push_back(Experiments::KK_DEFAULT);  // KK Default
    if (do_ltg) my_experiments.push_back(Experiments::LTG);                // LTG

    out << "========================================================" << endl
        << xpetraParameters
        // << matrixParameters
        << "========================================================" << endl
        << "Template Types:" << endl
        << "  Scalar:        " << Teuchos::demangleName(typeid(SC).name()) << endl
        << "  LocalOrdinal:  " << Teuchos::demangleName(typeid(LO).name()) << endl
        << "  GlobalOrdinal: " << Teuchos::demangleName(typeid(GO).name()) << endl
        << "  Node:          " << Teuchos::demangleName(typeid(NO).name()) << endl
        << "Sizes:" << endl
        << "  Scalar:        " << sizeof(SC) << endl
        << "  LocalOrdinal:  " << sizeof(LO) << endl
        << "  GlobalOrdinal: " << sizeof(GO) << endl
        << "========================================================" << endl
        << "Matrix:        " << Teuchos::demangleName(typeid(Matrix).name()) << endl
        << "Vector:        " << Teuchos::demangleName(typeid(Vector).name()) << endl
        << "Hierarchy:     " << Teuchos::demangleName(typeid(Hierarchy).name()) << endl
        << "========================================================" << endl;

#if defined(HAVE_TPETRA_INST_OPENMP)
    out << "Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency() = " << Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency() << endl
        << "========================================================" << endl;
#endif

    // At the moment, this test only runs on one MPI rank
    if (comm->getSize() != 1) exit(1);

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time")));

    comm->barrier();

    RCP<Matrix> A = Xpetra::IO<SC, LO, GO, Node>::Read(std::string(matrixFileNameA), lib, comm);
    RCP<Matrix> B = Xpetra::IO<SC, LO, GO, Node>::Read(std::string(matrixFileNameB), lib, comm);
    RCP<Matrix> B1, B2, C;

    RCP<const Map> Browmap = B->getRowMap();
    RCP<const Map> Bcolmap = B->getColMap();

    // Do some surgery to generate B1 and B2
    const GO GO_INVALID = Teuchos::OrdinalTraits<GO>::invalid();
    size_t brows        = Browmap->getLocalNumElements();
    size_t b2rows       = std::max((size_t)0, std::min((size_t)(brows * percent_rows_to_remove / 100.0), brows));
    size_t b1rows       = brows - b2rows;
    {
      Teuchos::Array<GO> b1list(b1rows), b2list(b2rows);
      for (size_t i = 0, i1 = 0, i2 = 0; i < brows; i++) {
        if (i < b1rows) {
          b1list[i1] = Browmap->getGlobalElement(i);
          i1++;
        } else {
          b2list[i2] = Browmap->getGlobalElement(i);
          i2++;
        }
      }

      RCP<const Map> B1rowmap = Xpetra::MapFactory<LO, GO, Node>::Build(lib, GO_INVALID, b1list(), 0, comm);
      RCP<const Map> B2rowmap = Xpetra::MapFactory<LO, GO, Node>::Build(lib, GO_INVALID, b2list(), 0, comm);

      RCP<const Import> B1import = Xpetra::ImportFactory<LO, GO, Node>::Build(Browmap, B1rowmap);
      RCP<const Import> B2import = Xpetra::ImportFactory<LO, GO, Node>::Build(Browmap, B2rowmap);
      RCP<const Map> dummy;
      B1 = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(B, *B1import, dummy, B1rowmap);
      B2 = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(B, *B2import, dummy, B2rowmap);
    }

    globalTimeMonitor = Teuchos::null;
    comm->barrier();

    out << "Matrix Read complete." << endl
        << "Matrix A:" << endl
        << *A
        << "========================================================" << endl
        << "Matrix B1: rows " << b1rows << endl
        << *B1
        << "========================================================" << endl
        << "Matrix B2: rows " << b2rows << endl
        << *B2
        << "========================================================" << endl;

    // random source
    std::random_device rd;
    std::mt19937 random_source(rd());

    // no need for a barrier, because the randomization process uses a collective.
    if (!my_experiments.empty()) {
      for (int i = 0; i < nrepeat; i++) {
        // randomize the experiments
        if (comm->getRank() == 0) {
          std::shuffle(my_experiments.begin(), my_experiments.end(), random_source);
        }

        // Broadcast this ordering to the other processes
        comm->broadcast(0,
                        static_cast<int>(sizeof(Experiments::LTG) * my_experiments.size()),
                        reinterpret_cast<char *>(my_experiments.data()));

        // loop over the randomized experiments
        for (const auto &experiment_id : my_experiments) {
          switch (experiment_id) {
            // MKL MULT_ADD
            case Experiments::MKL_MULT_ADD:
#ifdef HAVE_MUELU_MKL
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM2 MKL MULT_ADD: Total"));
                MM2_MKL(*A, *B1, *B2, *C, Bcolmap, std::string("MULT_ADD"));
              }
#endif
              break;
            // MKL ADD_MULT
            case Experiments::MKL_ADD_MULT:
#ifdef HAVE_MUELU_MKL
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM2 MKL ADD_MULT: Total"));
                MM2_MKL(*A, *B1, *B2, *C, Bcolmap, std::string("ADD_MULT"));
              }
#endif
              break;
            // KK Algorithms (KK Memory)
            case Experiments::KK_MEM:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM2 SPGEMM_KK_MEMORY: Total"));
                MM2_Wrapper(*A, *B1, *B2, *C, Bcolmap, std::string("SPGEMM_KK_MEMORY"), kk_team_work_size);
              }
              break;
            // KK Algorithms (KK Dense)
            case Experiments::KK_DENSE:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM2 SPGEMM_KK_DENSE: Total"));
                MM2_Wrapper(*A, *B1, *B2, *C, Bcolmap, std::string("SPGEMM_KK_DENSE"), kk_team_work_size);
              }
              break;
            // KK Algorithms (KK Default)
            case Experiments::KK_DEFAULT:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM2 SPGEMM_KK: Total"));
                MM2_Wrapper(*A, *B1, *B2, *C, Bcolmap, std::string("SPGEMM_KK"), kk_team_work_size);
              }
              break;
            // LTG
            case Experiments::LTG:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM2 LTG: Total"));
                MM2_Wrapper(*A, *B1, *B2, *C, Bcolmap, std::string("LTG"), kk_team_work_size);
              }
              break;

            default:
              std::cerr << "Unknown experiment ID encountered: " << (int)experiment_id << std::endl;
          }
          comm->barrier();
        }  // end random exp loop
      }    // end repeat
    }      // end ! my_experiments.empty()

    if (printTimings) {
      TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
