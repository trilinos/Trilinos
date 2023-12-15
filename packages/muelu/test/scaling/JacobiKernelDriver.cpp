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

// =========================================================================
// MKL Testing
// =========================================================================
#ifdef HAVE_MUELU_MKL
#include "mkl.h"

// mkl_sparse_spmm
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Jacobi_MKL_SPMM(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C, Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &D, Scalar omega) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  // NOTE: MKL uses int as its columnn index type and the either double or float for its Scalar type

  // Jacobi Does ( I - omega D A) B

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
  RCP<TimeMonitor> tm;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;
  typedef typename crs_matrix_type::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  typedef typename vector_type::device_type device_type;
  typedef typename Kokkos::View<MKL_INT *, typename lno_nnz_view_t::array_layout, typename lno_nnz_view_t::device_type> mkl_int_type;

  RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A, false));
  RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B, false));
  RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C, false));
  RCP<crs_matrix_type> Cnc      = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);
  RCP<const vector_type> Du     = Xpetra::toTpetra(D);

  const KCRS &Amat = Au->getLocalMatrixDevice();
  const KCRS &Bmat = Bu->getLocalMatrixDevice();

  if (A.getLocalNumRows() != C.getLocalNumRows()) throw std::runtime_error("C is not sized correctly");

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  lno_view_t Crowptr("Crowptr", C.getLocalNumRows() + 1);
  c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  lno_nnz_view_t Ccolind;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  scalar_view_t Cvals;

  auto Dvals = Du->getLocalViewDevice(Tpetra::Access::ReadOnly);

  sparse_matrix_t AMKL;
  sparse_matrix_t BMKL;
  sparse_matrix_t CMKL;
  sparse_matrix_t DMKL;

  sparse_matrix_t XTempMKL, YTempMKL;

  sparse_status_t result;
  // **********************************
  // Copy in the matrix data for MKL
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Jacobi MKL: Matrix CopyIn")));

  mkl_int_type ArowptrMKL("Arowptr", Arowptr.extent(0));
  mkl_int_type BrowptrMKL("Browptr", Browptr.extent(0));

  mkl_int_type AcolindMKL("Acolind", Acolind.extent(0));
  mkl_int_type BcolindMKL("Bcolind", Bcolind.extent(0));

  copy_view(Arowptr, ArowptrMKL);
  copy_view(Browptr, BrowptrMKL);
  copy_view(Acolind, AcolindMKL);
  copy_view(Bcolind, BcolindMKL);

  if (std::is_same<Scalar, double>::value) {
    mkl_sparse_d_create_csr(&AMKL, SPARSE_INDEX_BASE_ZERO, Au->getLocalNumRows(), Au->getLocalNumCols(), ArowptrMKL.data(), ArowptrMKL.data() + 1, AcolindMKL.data(), (double *)Avals.data());
    mkl_sparse_d_create_csr(&BMKL, SPARSE_INDEX_BASE_ZERO, Bu->getLocalNumRows(), Bu->getLocalNumCols(), BrowptrMKL.data(), BrowptrMKL.data() + 1, BcolindMKL.data(), (double *)Bvals.data());
  } else
    throw std::runtime_error("MKL Type Mismatch");
  tm = Teuchos::null;
  Au->getComm()->barrier();

  // **********************************
  // Copy in the vector-as-matrix data for MKL
  tm       = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Jacobi MKL: Scaled D-Vector CopyIn")));
  size_t n = Dvals.extent(0);
  mkl_int_type DrowptrMKL("Drowptr", n + 1);
  mkl_int_type DcolindMKL("Dcolind", n);
  scalar_view_t DvalsMKL("Dvals", n);
  Kokkos::parallel_for(
      n + 1, KOKKOS_LAMBDA(const size_t i) {
        DrowptrMKL[i] = i;
        if (i < n) {
          DcolindMKL[i] = i;
          DvalsMKL[i]   = -omega * Dvals(i, 0);
        }
      });

  // NOTE: No sanity checks here.  We did that above
  mkl_sparse_d_create_csr(&DMKL, SPARSE_INDEX_BASE_ZERO, Dvals.extent(0), Dvals.extent(0), DrowptrMKL.data(), DrowptrMKL.data() + 1, DcolindMKL.data(), DvalsMKL.data());

  tm = Teuchos::null;
  Au->getComm()->barrier();

  // **********************************
  // Multiply (A*B)
  tm     = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Jacobi MKL: Multiply")));
  result = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, AMKL, BMKL, &XTempMKL);
  typename KCRS::execution_space().fence();
  if (result != SPARSE_STATUS_SUCCESS) throw std::runtime_error("MKL Multiply failed");

  // **********************************
  // Scale (-omegaD) * AB)
  tm     = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Jacobi MKL: Scale-Via-Multiply")));
  result = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, DMKL, XTempMKL, &YTempMKL);
  typename KCRS::execution_space().fence();
  if (result != SPARSE_STATUS_SUCCESS) throw std::runtime_error("MKL Scale failed");

  // **********************************
  // Add B - ((-omegaD) * AB))
  tm     = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Jacobi MKL: Add")));
  result = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, BMKL, 1.0, YTempMKL, &CMKL);
  typename KCRS::execution_space().fence();
  if (result != SPARSE_STATUS_SUCCESS) throw std::runtime_error("MKL Add failed");

  // **********************************
  // Copy out the data for MKL
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Jacobi MKL: Copy Out")));
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

  Cnc->replaceColMap(Bu->getColMap());
  Cnc->setAllValues(Crowptr,
                    Ccolind,
                    Cvals);

  mkl_sparse_destroy(AMKL);
  mkl_sparse_destroy(BMKL);
  mkl_sparse_destroy(CMKL);
  mkl_sparse_destroy(XTempMKL);
  mkl_sparse_destroy(YTempMKL);

  tm = Teuchos::null;
  Au->getComm()->barrier();
}

#endif

// =========================================================================
// LTG Testing
// =========================================================================
#include "Tpetra_Import_Util2.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_decl.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_def.hpp"

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Jacobi_Wrapper(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C, const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &D, Scalar omega, std::string jacobi_algorithm_name, std::string spgemm_algorithm_name, int team_work_size) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
  RCP<TimeMonitor> tm;

  std::string name = jacobi_algorithm_name + "/" + spgemm_algorithm_name;

  if (lib == Xpetra::UseTpetra) {
#if defined(HAVE_TPETRA_INST_OPENMP)
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
    typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;
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
    RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A, false));
    RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B, false));
    RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C, false));
    RCP<crs_matrix_type> Cnc      = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);

    RCP<const vector_type> Du = Xpetra::toTpetra(D);

    // **********************************
    // Copy in the data for Jacobi Kernel Wrapper
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Jacobi " + name + ": CopyIn")));

    Tpetra::CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node> Aview, Bview;
    Aview.origMatrix   = Au;
    Aview.origRowMap   = Au->getRowMap();
    Aview.rowMap       = Au->getRowMap();
    Aview.colMap       = Au->getColMap();
    Aview.domainMap    = Au->getDomainMap();
    Aview.importColMap = Teuchos::null;

    Bview.origMatrix   = Bu;
    Bview.origRowMap   = Bu->getRowMap();
    Bview.rowMap       = Bu->getRowMap();
    Bview.colMap       = Bu->getColMap();
    Bview.domainMap    = Bu->getDomainMap();
    Bview.importColMap = Teuchos::null;

    // Because we're in serial...
    Cnc->replaceColMap(Bu->getColMap());

    // Bcol2Ccol is trivial
    lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"), Bview.colMap->getLocalNumElements()), Icol2Ccol;
    const LO colMapSize = static_cast<LO>(Bview.colMap->getLocalNumElements());
    Kokkos::parallel_for(
        "Tpetra::mult_A_B_newmatrix::Bcol2Ccol_fill",
        Kokkos::RangePolicy<execution_space, LO>(0, colMapSize),
        KOKKOS_LAMBDA(const LO i) {
          Bcol2Ccol(i) = i;
        });

    // Acol2Brow
    local_map_type Acolmap_local = Aview.colMap->getLocalMap();
    local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
    lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"), Aview.colMap->getLocalNumElements());
    lo_view_t targetMapToImportRow;
    Kokkos::parallel_for(
        "Tpetra::mult_A_B_newmatrix::construct_tables", range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex() + 1), KOKKOS_LAMBDA(const LO i) {
          GO aidx  = Acolmap_local.getGlobalElement(i);
          LO B_LID = Browmap_local.getLocalElement(aidx);
          if (B_LID != LO_INVALID) {
            targetMapToOrigRow(i) = B_LID;
            //        targetMapToImportRow(i) = LO_INVALID;
          } else {
            // This shouldn't happen here
          }
        });
    tm = Teuchos::null;
    Au->getComm()->barrier();

    // **********************************
    // Multiply
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set("openmp: jacobi algorithm", jacobi_algorithm_name);
    params->set("openmp: algorithm", spgemm_algorithm_name);
    params->set("openmp: team work size", team_work_size);

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(std::string("Jacobi ") + name + ": Kernel")));
    Tpetra::MMdetails::KernelWrappers2<SC, LO, GO, Node, lno_nnz_view_t>::jacobi_A_B_newmatrix_kernel_wrapper(omega, *Du, Aview, Bview, targetMapToOrigRow, targetMapToImportRow, Bcol2Ccol, Icol2Ccol, *Cnc, Cimport, name, params);

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
    bool do_mkl        = true;
    bool do_kk_mem     = true;
    bool do_kk_dense   = true;
    bool do_kk_default = true;
    bool do_ltg_msak   = true;
    bool do_ltg_jac    = true;

#ifndef HAVE_MUELU_MKL
    do_mkl = false;
#endif
    clp.setOption("mkl", "nomkl", &do_mkl, "Evaluate MKL SpMM");
    clp.setOption("kk_mem", "nokk_mem", &do_kk_mem, "Evaluate KK Mem");
    clp.setOption("kk_dense", "nokk_dense", &do_kk_dense, "Evaluate KK Dense");
    clp.setOption("kk_default", "nokk_default", &do_kk_default, "Evaluate KK Default");
    clp.setOption("ltg_msak", "noltg_msak", &do_ltg_msak, "Evaluate LTG (MSAK)");
    clp.setOption("ltg_jac", "noltg_jac", &do_ltg_jac, "Evaluate LTG (Jacobi)");

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

#ifndef HAVE_MUELU_MKL
    if (do_mkl == true) {
      out << "MKL was requested, but this kernel is not available. Disabling..." << endl;
      do_mkl = false;
    }
#endif

    // simple hack to randomize order of experiments
    enum class Experiments { MKL_SPMM = 0,
                             KK_MEM,
                             KK_DENSE,
                             KK_DEFAULT,
                             LTG_MSAK,
                             LTG_JAC };
    std::vector<Experiments> my_experiments;
    // add the experiments we will run

#ifdef HAVE_MUELU_MKL
    if (do_mkl) my_experiments.push_back(Experiments::MKL_SPMM);  // MKL SPMM
#endif

    // assume these are available
    if (do_kk_mem) my_experiments.push_back(Experiments::KK_MEM);          // KK Mem
    if (do_kk_dense) my_experiments.push_back(Experiments::KK_DENSE);      // KK Dense
    if (do_kk_default) my_experiments.push_back(Experiments::KK_DEFAULT);  // KK Default
    if (do_ltg_msak) my_experiments.push_back(Experiments::LTG_MSAK);      // LTG MSAK
    if (do_ltg_jac) my_experiments.push_back(Experiments::LTG_JAC);        // LTG Jacobi

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
    RCP<Matrix> C;

    // Jacobi-specific
    RCP<Vector> Diagonal = Xpetra::VectorFactory<SC, LO, GO, Node>::Build(A->getRowMap(), true);
    A->getLocalDiagCopy(*Diagonal);
    RCP<Vector> D = Xpetra::VectorFactory<SC, LO, GO, Node>::Build(A->getRowMap(), true);
    D->reciprocal(*Diagonal);
    SC omega = (4.0 / 3.0) * Teuchos::ScalarTraits<Scalar>::one();

    globalTimeMonitor = Teuchos::null;
    comm->barrier();

    out << "Matrix Read complete." << endl
        << "Matrix A:" << endl
        << *A
        << "========================================================" << endl
        << "Matrix B:" << endl
        << *B
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
                        static_cast<int>(sizeof(Experiments::LTG_MSAK) * my_experiments.size()),
                        reinterpret_cast<char *>(my_experiments.data()));

        // loop over the randomized experiments
        for (const auto &experiment_id : my_experiments) {
          switch (experiment_id) {
            // MKL_SPMM
            case Experiments::MKL_SPMM:
#ifdef HAVE_MUELU_MKL
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("JAC MKL: Total"));
                Jacobi_MKL_SPMM(*A, *B, *C, *D, omega);
              }
#endif
              break;
            // KK Algorithms (KK Memory)
            case Experiments::KK_MEM:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("JAC MSAK/SPGEMM_KK_MEMORY: Total"));
                Jacobi_Wrapper(*A, *B, *C, *D, omega, std::string("MSAK"), std::string("SPGEMM_KK_MEMORY"), kk_team_work_size);
              }
              break;
            // KK Algorithms (KK Dense)
            case Experiments::KK_DENSE:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("JAC MSAK/SPGEMM_KK_DENSE: Total"));
                Jacobi_Wrapper(*A, *B, *C, *D, omega, std::string("MSAK"), std::string("SPGEMM_KK_DENSE"), kk_team_work_size);
              }
              break;
            // KK Algorithms (KK Default)
            case Experiments::KK_DEFAULT:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("JAC SPGEMM_KK: Total"));
                Jacobi_Wrapper(*A, *B, *C, *D, omega, std::string("MSAK"), std::string("SPGEMM_KK"), kk_team_work_size);
              }
              break;
            // LTG (MSAK)
            case Experiments::LTG_MSAK:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("JAC MSAK/LTG: Total"));
                Jacobi_Wrapper(*A, *B, *C, *D, omega, std::string("MSAK"), std::string("LTG"), 0);
              }
              break;
            // LTG (Jacobi)
            case Experiments::LTG_JAC:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("JAC JLTG: Total"));
                Jacobi_Wrapper(*A, *B, *C, *D, omega, std::string("LTG"), std::string("--"), 0);
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
