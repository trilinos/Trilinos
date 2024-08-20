// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stdlib.h>                         // for exit, EXIT_FA...
#include <Teuchos_StandardCatchMacros.hpp>  // for TEUCHOS_STAND...
#include <iostream>                         // for endl, ostring...
#include <limits>                           // for numeric_limits
#include <random>                           // for mt19937, rand...
#include <stdexcept>                        // for runtime_error
#include <string>                           // for string
#include <fstream>                          // write CSV
#include "MueLu_ConfigDefs.hpp"
#include "KokkosKernels_Handle.hpp"          // for KokkosKernels...
#include "KokkosSparse_spgemm_handle.hpp"    // for StringToSPGEM...
#include "KokkosSparse_spgemm_numeric.hpp"   // for spgemm_numeric
#include "KokkosSparse_spgemm_symbolic.hpp"  // for spgemm_symbolic
#include "Teuchos_CommandLineProcessor.hpp"  // for CommandLinePr...
#include "Teuchos_ENull.hpp"                 // for ENull::null
#include "Teuchos_FancyOStream.hpp"          // for FancyOStream
#include "Teuchos_OrdinalTraits.hpp"         // for OrdinalTraits
#include "Teuchos_RCP.hpp"                   // for rcp_const_cast
#include "Teuchos_RCPDecl.hpp"               // for RCP, rcp
#include "Teuchos_ScalarTraitsDecl.hpp"      // for ScalarTraits
#include "Tpetra_CrsMatrix_decl.hpp"         // for CrsMatrix
#include "Tpetra_Import_decl.hpp"            // for Import
#include <Tpetra_Import_Util2.hpp>           // Import_Util::sortCrsEntries
#include "Tpetra_Map_decl.hpp"               // for Map
#include <Xpetra_IO.hpp>
#include "Xpetra_BlockedMap.hpp"    // for LO, GO, NO
#include "Xpetra_Map.hpp"           // for UnderlyingLib
#include "Xpetra_Parameters.hpp"    // for Parameters
#include "Xpetra_MatrixMatrix.hpp"  //Xpetra SpGEMM
namespace Teuchos {
class ParameterList;
}  // namespace Teuchos
namespace Teuchos {
class TimeMonitor;
}  // namespace Teuchos
namespace Teuchos {
template <typename Ordinal>
class Comm;
}  // namespace Teuchos
namespace Xpetra {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class Matrix;
}  // namespace Xpetra

#ifdef HAVE_VTUNE
#include "ittnotify.h"
__itt_domain *vtune_domain_modLTG;
__itt_string_handle *vtune_task_modLTG_setup;
__itt_string_handle *vtune_task_modLTG_mult;
__itt_string_handle *vtune_task_modLTG_copy;
#endif

// MueLu
#include "MueLu.hpp"
#include "MueLu_TestHelpers.hpp"

#include <sys/types.h>
#include <unistd.h>

int getRssHWM() {
  static const std::string mem_file = "/proc/self/status";
  std::ifstream infile(mem_file);

  std::string line;
  while (std::getline(infile, line)) {
    if (line.rfind("VmHWM", 0) == 0) {
      // VmHWM:      1764 kB
      std::stringstream ss(line);
      std::string key;
      int value;
      std::string unit;
      ss >> key >> value >> unit;
      infile.close();
      return (value);
    }
  }
  infile.close();
  return (-1);
}

#if USE_DESCRIPTIVE_STATS == 1
#include "descriptive_stats/descriptive_stats.hpp"
#include <iomanip>

// this will generate a pedantic warning
#define DescriptiveTime(x) ({                                            \
  int rssStart = -1;                                                     \
  int rssEnd   = -1;                                                     \
  if (trackRssHWM) {                                                     \
    rssStart = getRssHWM();                                              \
  }                                                                      \
  auto t0 = Time::now();                                                 \
  x;                                                                     \
  auto t1 = Time::now();                                                 \
  if (trackRssHWM) {                                                     \
    rssEnd = getRssHWM();                                                \
    mem_vector.push_back(static_cast<double>(rssEnd - rssStart));        \
  }                                                                      \
  double_secs ds = t1 - t0;                                              \
  if (descStats)                                                         \
    timing_vector.push_back(std::chrono::duration_cast<ns>(ds).count()); \
})

std::map<std::string, std::vector<double> > my_experiment_timings;

std::chrono::steady_clock::time_point ltg_core_start;
std::chrono::steady_clock::time_point ltg_sort_start;
std::chrono::steady_clock::time_point ltg_copy_start;
std::chrono::steady_clock::time_point ltg_mult_start;

std::chrono::steady_clock::time_point ltg_core_stop;
std::chrono::steady_clock::time_point ltg_sort_stop;
std::chrono::steady_clock::time_point ltg_copy_stop;
std::chrono::steady_clock::time_point ltg_mult_stop;

// MKL
std::chrono::steady_clock::time_point mkl_mult_start;
std::chrono::steady_clock::time_point mkl_copyin_start;
std::chrono::steady_clock::time_point mkl_copyout_start;
std::chrono::steady_clock::time_point mkl_sort_start;

std::chrono::steady_clock::time_point mkl_mult_stop;
std::chrono::steady_clock::time_point mkl_copyin_stop;
std::chrono::steady_clock::time_point mkl_copyout_stop;
std::chrono::steady_clock::time_point mkl_sort_stop;

// ViennaCL
std::chrono::steady_clock::time_point vcl_mult_start;
std::chrono::steady_clock::time_point vcl_mult_stop;

// KK D
std::chrono::steady_clock::time_point kk_dense_mult_start;
std::chrono::steady_clock::time_point kk_dense_copyout_start;

std::chrono::steady_clock::time_point kk_dense_mult_stop;
std::chrono::steady_clock::time_point kk_dense_copyout_stop;

// KK M
std::chrono::steady_clock::time_point kk_mem_mult_start;
std::chrono::steady_clock::time_point kk_mem_copyout_start;

std::chrono::steady_clock::time_point kk_mem_mult_stop;
std::chrono::steady_clock::time_point kk_mem_copyout_stop;

// KK D
std::chrono::steady_clock::time_point kk_def_mult_start;
std::chrono::steady_clock::time_point kk_def_copyout_start;

std::chrono::steady_clock::time_point kk_def_mult_stop;
std::chrono::steady_clock::time_point kk_def_copyout_stop;

// Xpetra
std::chrono::steady_clock::time_point xpetra_mult_start;
std::chrono::steady_clock::time_point xpetra_mult_stop;

#undef remove_pedantic_macro_warning_type
#else
#define DescriptiveTime(x) x
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

enum class Experiments { ViennaCL = 0,
                         MKL_SPMM,
                         KK_MEM,
                         KK_DENSE,
                         KK_DEFAULT,
                         LTG,
                         SERIAL,
                         XPETRA };

static inline std::string to_string(const Experiments &experiment_id) {
  switch (experiment_id) {
    case Experiments::ViennaCL:
      return ("ViennaCL");
    case Experiments::MKL_SPMM:
      return ("MKL");
    case Experiments::KK_MEM:
      return ("KK_MEM");
    case Experiments::KK_DENSE:
      return ("KK_DENSE");
    case Experiments::KK_DEFAULT:
      return ("KK_DEFAULT");
    case Experiments::LTG:
      return ("LTG");
    case Experiments::SERIAL:
      return ("SERIAL");
    case Experiments::XPETRA:
      return ("XPETRA");
    default:
      return ("UNDEFINED");
  }
}

// =========================================================================
// ViennaCL Testing
// =========================================================================
#ifdef HAVE_MUELU_VIENNACL
#define VIENNACL_WITH_OPENMP
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/host_based/common.hpp"

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Multiply_ViennaCL(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
  RCP<TimeMonitor> tm;

  // NOTE: ViennaCL matrices use "unsigned int" for rowptr and colind and are templated on scalar type (yay); which means (for us) POD only.

  if (lib == Xpetra::UseTpetra) {
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
    typedef typename crs_matrix_type::local_matrix_type KCRS;
    typedef typename KCRS::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::row_map_type::const_type c_lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
    typedef typename KCRS::values_type::non_const_type scalar_view_t;
    typedef typename Kokkos::View<unsigned int *, typename lno_nnz_view_t::array_layout, typename lno_nnz_view_t::device_type> vcl_size_t_type;

    RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A, false));
    RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B, false));
    RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C, false));
    RCP<crs_matrix_type> Cnc      = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);

    const KCRS &Amat = Au->getLocalMatrixDevice();
    const KCRS &Bmat = Bu->getLocalMatrixDevice();

    using no_init_view = Kokkos::ViewAllocateWithoutInitializing;

    c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
    lno_view_t Crowptr(no_init_view("Crowptr"), C.getLocalNumRows() + 1);
    c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
    lno_nnz_view_t Ccolind;
    const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
    scalar_view_t Cvals;

    // **********************************
    // Copy in the data for ViennaCL

    vcl_size_t_type ArowptrVCL(no_init_view("Arowptr"), Arowptr.extent(0));
    vcl_size_t_type BrowptrVCL(no_init_view("Browptr"), Browptr.extent(0));

    vcl_size_t_type AcolindVCL(no_init_view("Acolind"), Acolind.extent(0));
    vcl_size_t_type BcolindVCL(no_init_view("Bcolind"), Bcolind.extent(0));

    viennacl::compressed_matrix<Scalar> AVCL, BVCL;
    {
      Teuchos::TimeMonitor tm_copyin(*TimeMonitor::getNewTimer("MM ViennaCL: CopyIn"));
      copy_view(Arowptr, ArowptrVCL);
      copy_view(Browptr, BrowptrVCL);
      copy_view(Acolind, AcolindVCL);
      copy_view(Bcolind, BcolindVCL);

      AVCL.set(ArowptrVCL.data(), AcolindVCL.data(), Avals.data(),
               Au->getLocalNumRows(), Au->getLocalNumCols(), Au->getLocalNumEntries());
      BVCL.set(BrowptrVCL.data(), BcolindVCL.data(), Bvals.data(),
               Bu->getLocalNumRows(), Bu->getLocalNumCols(), Bu->getLocalNumEntries());
    }
    Au->getComm()->barrier();

    viennacl::compressed_matrix<Scalar> CVCL;
    // **********************************
    // Multiply
    {
      Teuchos::TimeMonitor tm_mult(*TimeMonitor::getNewTimer("MM ViennaCL: Multiply"));
#ifdef USE_DESCRIPTIVE_STATS
      vcl_mult_start = std::chrono::steady_clock::now();
#endif

      CVCL = viennacl::linalg::prod(AVCL, BVCL);

#ifdef USE_DESCRIPTIVE_STATS
      vcl_mult_stop = std::chrono::steady_clock::now();
#endif
    }
    KCRS::execution_space::fence();

    // **********************************
    // Copy out the data for ViennaCL
    {
      Teuchos::TimeMonitor tm_copy(*TimeMonitor::getNewTimer("MM ViennaCL: Copy Out"));
      size_t cnnz = (size_t)CVCL.nnz();
      Ccolind     = lno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("Ccolind"), cnnz);
      Cvals       = scalar_view_t(Kokkos::ViewAllocateWithoutInitializing("Cvals"), cnnz);
#ifdef VIENNACL_WITH_CUDA
      const unsigned int *CrowptrVCL = viennacl::cuda_arg<unsigned int>(CVCL.handle1());
      const unsigned int *CcolindVCL = viennacl::cuda_arg<unsigned int>(CVCL.handle2());
      const Scalar *CvalsVCL         = viennacl::cuda_arg<Scalar>(CVCL.handle());
#else
      const unsigned int *CrowptrVCL = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(CVCL.handle1());
      const unsigned int *CcolindVCL = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(CVCL.handle2());
      const Scalar *CvalsVCL         = viennacl::linalg::host_based::detail::extract_raw_pointer<Scalar>(CVCL.handle());
#endif

      copy_view_n(Crowptr.extent(0), CrowptrVCL, Crowptr);
      copy_view_n(cnnz, CcolindVCL, Ccolind);
      copy_view_n(cnnz, CvalsVCL, Cvals);

      // set values
      {
        // Because we're in serial...
        Cnc->replaceColMap(Bu->getColMap());

        Cnc->setAllValues(Crowptr,
                          Ccolind,
                          Cvals);
      }
    }
#ifdef USE_DESCRIPTIVE_STATS
    {
      std::chrono::duration<double> ds;

      ds = vcl_mult_stop - vcl_mult_start;
      my_experiment_timings["ViennaCL: Call"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());
    }
#endif

    Au->getComm()->barrier();
  }
}

#endif

// =========================================================================
// MKL Testing
// =========================================================================
#ifdef HAVE_MUELU_MKL
#include "mkl_spblas.h"  // for mkl_sparse_de...
#include "mkl_types.h"   // for MKL_INT
//#include "mkl.h"

#define MMKD_MKL_ERROR_CHECK(rc)                                                \
  {                                                                             \
    if (mkl_rc != SPARSE_STATUS_SUCCESS) {                                      \
      std::stringstream ss;                                                     \
      switch (mkl_rc) {                                                         \
        case SPARSE_STATUS_NOT_INITIALIZED:                                     \
          ss << "[ERROR] mkl::mkl_sparse_spmm: SPARSE_STATUS_NOT_INITIALIZED";  \
          break;                                                                \
        case SPARSE_STATUS_ALLOC_FAILED:                                        \
          ss << "[ERROR] mkl::mkl_sparse_spmm: SPARSE_STATUS_ALLOC_FAILED";     \
          break;                                                                \
        case SPARSE_STATUS_INVALID_VALUE:                                       \
          ss << "[ERROR] mkl::mkl_sparse_spmm: SPARSE_STATUS_INVALID_VALUE";    \
          break;                                                                \
        case SPARSE_STATUS_EXECUTION_FAILED:                                    \
          ss << "[ERROR] mkl::mkl_sparse_spmm: SPARSE_STATUS_EXECUTION_FAILED"; \
          break;                                                                \
        case SPARSE_STATUS_INTERNAL_ERROR:                                      \
          ss << "[ERROR] mkl::mkl_sparse_spmm: SPARSE_STATUS_INTERNAL_ERROR";   \
          break;                                                                \
        case SPARSE_STATUS_NOT_SUPPORTED:                                       \
          ss << "[ERROR] mkl::mkl_sparse_spmm: SPARSE_STATUS_NOT_SUPPORTED";    \
          break;                                                                \
        default:                                                                \
          break;                                                                \
      }                                                                         \
      std::cerr << ss.str() << std::endl;                                       \
      return;                                                                   \
    }                                                                           \
  }

// mkl_sparse_spmm
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Multiply_MKL_SPMM(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  // NOTE: MKL uses int as its columnn index type and the either double or float for its Scalar type
  static constexpr bool sortCSR = false;

  RCP<TimeMonitor> tm;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef typename crs_matrix_type::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  typedef typename Kokkos::View<MKL_INT *, typename lno_nnz_view_t::array_layout, typename lno_nnz_view_t::device_type> mkl_int_type;

  using no_init_view = Kokkos::ViewAllocateWithoutInitializing;

  RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A, false));
  RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B, false));
  RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C, false));
  RCP<crs_matrix_type> Cnc      = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);

  const KCRS &Amat = Au->getLocalMatrixDevice();
  const KCRS &Bmat = Bu->getLocalMatrixDevice();
  if (A.getLocalNumRows() != C.getLocalNumRows()) throw std::runtime_error("C is not sized correctly");

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  lno_view_t Crowptr(no_init_view("Crowptr"), C.getLocalNumRows() + 1);
  c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;

  sparse_matrix_t AMKL;
  sparse_matrix_t BMKL;
  sparse_matrix_t CMKL;
  sparse_status_t mkl_rc;

  mkl_int_type ArowptrMKL(no_init_view("Arowptr"), Arowptr.extent(0));
  mkl_int_type BrowptrMKL(no_init_view("Browptr"), Browptr.extent(0));

  mkl_int_type AcolindMKL(no_init_view("Acolind"), Acolind.extent(0));
  mkl_int_type BcolindMKL(no_init_view("Bcolind"), Bcolind.extent(0));

  {
    // **********************************
    // Copy in the data for MKL
    Teuchos::TimeMonitor tm_copyin(*TimeMonitor::getNewTimer("MM MKL: CopyIn"));
#ifdef USE_DESCRIPTIVE_STATS
    mkl_copyin_start = std::chrono::steady_clock::now();
#endif
    copy_view(Arowptr, ArowptrMKL);
    copy_view(Browptr, BrowptrMKL);
    copy_view(Acolind, AcolindMKL);
    copy_view(Bcolind, BcolindMKL);

    if (std::is_same<Scalar, double>::value) {
      mkl_rc = mkl_sparse_d_create_csr(&AMKL, SPARSE_INDEX_BASE_ZERO,
                                       Au->getLocalNumRows(),
                                       Au->getLocalNumCols(),
                                       ArowptrMKL.data(),
                                       ArowptrMKL.data() + 1,
                                       AcolindMKL.data(),
                                       (double *)Avals.data());
      MMKD_MKL_ERROR_CHECK(mkl_rc);
      mkl_rc = mkl_sparse_d_create_csr(&BMKL, SPARSE_INDEX_BASE_ZERO,
                                       Bu->getLocalNumRows(),
                                       Bu->getLocalNumCols(),
                                       BrowptrMKL.data(),
                                       BrowptrMKL.data() + 1,
                                       BcolindMKL.data(),
                                       (double *)Bvals.data());
      MMKD_MKL_ERROR_CHECK(mkl_rc);
    } else
      throw std::runtime_error("MKL Type Mismatch");

#ifdef USE_DESCRIPTIVE_STATS
    mkl_copyin_stop = std::chrono::steady_clock::now();
#endif
  }
  Au->getComm()->barrier();

  // **********************************
  // Multiply
  {
    Teuchos::TimeMonitor tm_mult(*TimeMonitor::getNewTimer("MM MKL: Multiply"));
#ifdef USE_DESCRIPTIVE_STATS
    mkl_mult_start = std::chrono::steady_clock::now();
#endif
    mkl_rc = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, AMKL, BMKL, &CMKL);
    MMKD_MKL_ERROR_CHECK(mkl_rc);
#ifdef USE_DESCRIPTIVE_STATS
    mkl_mult_stop = std::chrono::steady_clock::now();
#endif
    // KCRS::execution_space::fence();
  }
  // **********************************

  // Copy out the data for MKL
  Teuchos::TimeMonitor tm_copy(*TimeMonitor::getNewTimer("MM MKL: Copy Out"));
#ifdef USE_DESCRIPTIVE_STATS
  mkl_copyout_start = std::chrono::steady_clock::now();
#endif

  sparse_index_base_t c_indexing;
  MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
  double *values;
  mkl_rc = mkl_sparse_d_export_csr(CMKL, &c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values);
  MMKD_MKL_ERROR_CHECK(mkl_rc);
  size_t cnnz            = rows_end[c_rows - 1];
  lno_nnz_view_t Ccolind = lno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("Ccolind"), cnnz);
  scalar_view_t Cvals    = scalar_view_t(Kokkos::ViewAllocateWithoutInitializing("Cvals"), cnnz);
  Kokkos::resize(Ccolind, cnnz);
  Kokkos::resize(Cvals, cnnz);
  if ((size_t)c_rows != A.getLocalNumRows() || (size_t)c_rows + 1 != Crowptr.extent(0)) throw std::runtime_error("C row size mismatch");
  copy_view_n(c_rows, rows_start, Crowptr);
  Crowptr(c_rows) = rows_end[c_rows - 1];
  copy_view_n(cnnz, columns, Ccolind);
  copy_view_n(cnnz, values, Cvals);

#ifdef USE_DESCRIPTIVE_STATS
  mkl_copyout_stop = std::chrono::steady_clock::now();
#endif

  if (sortCSR == true) {
#ifdef USE_DESCRIPTIVE_STATS
    mkl_sort_start = std::chrono::steady_clock::now();
#endif

    Tpetra::Import_Util::sortCrsEntries(Crowptr,
                                        Ccolind,
                                        Cvals);

#ifdef USE_DESCRIPTIVE_STATS
    mkl_sort_stop = std::chrono::steady_clock::now();
#endif
  }

  // set values
  {
    // Because we're in serial...
    Cnc->replaceColMap(Bu->getColMap());

    Cnc->setAllValues(Crowptr,
                      Ccolind,
                      Cvals);
  }

  mkl_rc = mkl_sparse_destroy(AMKL);
  MMKD_MKL_ERROR_CHECK(mkl_rc);
  mkl_rc = mkl_sparse_destroy(BMKL);
  MMKD_MKL_ERROR_CHECK(mkl_rc);
  mkl_rc = mkl_sparse_destroy(CMKL);
  MMKD_MKL_ERROR_CHECK(mkl_rc);

#ifdef USE_DESCRIPTIVE_STATS
  {
    std::chrono::duration<double> ds;

    ds = mkl_mult_stop - mkl_mult_start;
    my_experiment_timings["MKL: Call"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());

    ds = mkl_copyin_stop - mkl_copyin_start;
    my_experiment_timings["MKL: CopyIn"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());

    ds = mkl_copyout_stop - mkl_copyout_start;
    my_experiment_timings["MKL: CopyOut"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());

    ds = mkl_sort_stop - mkl_sort_start;
    my_experiment_timings["MKL: Sort"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());
  }
#endif
}

#undef MMKD_MKL_ERROR_CHECK

#endif

// =========================================================================
// Kokkos Kernels Testing
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Multiply_KokkosKernels(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B, Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C, std::string algorithm_name, int team_work_size) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();

  std::string prefix = std::string("MM KokkosKernels ") + algorithm_name + std::string(": ");

  if (lib == Xpetra::UseTpetra) {
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
    typedef typename crs_matrix_type::local_matrix_device_type KCRS;
    typedef typename KCRS::StaticCrsGraphType graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::row_map_type::const_type c_lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
    typedef typename KCRS::values_type::non_const_type scalar_view_t;
    typedef typename KCRS::device_type device_t;

    RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A, false));
    RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B, false));
    RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C, false));
    RCP<crs_matrix_type> Cnc      = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);

    const KCRS &Amat = Au->getLocalMatrixDevice();
    const KCRS &Bmat = Bu->getLocalMatrixDevice();

    c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
    lno_view_t Crowptr("Crowptr", A.getLocalNumRows() + 1);
    c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
    const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
    lno_nnz_view_t Ccolind;
    scalar_view_t Cvals;

    // KokkosKernelsHandle
    typedef KokkosKernels::Experimental::KokkosKernelsHandle<
        typename lno_view_t::const_value_type, typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type,
        typename device_t::execution_space, typename device_t::memory_space, typename device_t::memory_space>
        KernelHandle;
    KokkosSparse::SPGEMMAlgorithm alg_enum    = KokkosSparse::StringToSPGEMMAlgorithm(algorithm_name);
    typename KernelHandle::nnz_lno_t AnumRows = Au->getLocalNumRows();
    typename KernelHandle::nnz_lno_t BnumRows = Bu->getLocalNumRows();
    typename KernelHandle::nnz_lno_t BnumCols = Bu->getLocalNumCols();

    // **********************************
    // Multiply
    {
      Teuchos::TimeMonitor tm_mult(*TimeMonitor::getNewTimer(prefix + std::string("Multiply")));
#ifdef USE_DESCRIPTIVE_STATS
      switch (alg_enum) {
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK_DENSE:
          kk_dense_mult_start = std::chrono::steady_clock::now();
          break;
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK_MEMORY:
          kk_mem_mult_start = std::chrono::steady_clock::now();
          break;
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK:
          kk_def_mult_start = std::chrono::steady_clock::now();
          break;
        default:
          break;
      }
#endif
      KernelHandle kh;
      kh.create_spgemm_handle(alg_enum);
      kh.set_team_work_size(team_work_size);

      KokkosSparse::Experimental::spgemm_symbolic(&kh, AnumRows, BnumRows, BnumCols,
                                                  Arowptr, Acolind, false,
                                                  Browptr, Bcolind, false,
                                                  Crowptr);

      size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
      if (c_nnz_size) {
        using no_init_view = Kokkos::ViewAllocateWithoutInitializing;
        Ccolind            = lno_nnz_view_t(no_init_view("entriesC"), c_nnz_size);
        Cvals              = scalar_view_t(no_init_view("valuesC"), c_nnz_size);
      }
      KokkosSparse::Experimental::spgemm_numeric(&kh, AnumRows, BnumRows, BnumCols,
                                                 Arowptr, Acolind, Avals, false,
                                                 Browptr, Bcolind, Bvals, false,
                                                 Crowptr, Ccolind, Cvals);
      kh.destroy_spgemm_handle();
      typename KCRS::execution_space().fence();

#ifdef USE_DESCRIPTIVE_STATS
      switch (alg_enum) {
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK_DENSE:
          kk_dense_mult_stop = std::chrono::steady_clock::now();
          break;
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK_MEMORY:
          kk_mem_mult_stop = std::chrono::steady_clock::now();
          break;
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK:
          kk_def_mult_stop = std::chrono::steady_clock::now();
          break;
        default:
          break;
      }
#endif
    }  // scope time monitor mult
    Au->getComm()->barrier();

    // **********************************
    // Copy out the data for KokkosKernels
    {
      Teuchos::TimeMonitor tm_copy(*TimeMonitor::getNewTimer(prefix + std::string("Copy Out")));
#ifdef USE_DESCRIPTIVE_STATS
      switch (alg_enum) {
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK_DENSE:
          kk_dense_copyout_start = std::chrono::steady_clock::now();
          break;
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK_MEMORY:
          kk_mem_copyout_start = std::chrono::steady_clock::now();
          break;
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK:
          kk_def_copyout_start = std::chrono::steady_clock::now();
          break;
        default:
          break;
      }
#endif
      // Because we're in serial...
      Cnc->replaceColMap(Bu->getColMap());
      Cnc->setAllValues(Crowptr,
                        Ccolind,
                        Cvals);

#ifdef USE_DESCRIPTIVE_STATS
      std::chrono::duration<double> ds;
      switch (alg_enum) {
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK_DENSE:
          kk_dense_copyout_stop = std::chrono::steady_clock::now();

          ds = kk_dense_mult_stop - kk_dense_mult_start;
          my_experiment_timings["KK DENSE: Call"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());

          ds = kk_dense_copyout_stop - kk_dense_copyout_start;
          my_experiment_timings["KK DENSE: CopyOut"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());
          break;
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK_MEMORY:
          kk_mem_copyout_stop = std::chrono::steady_clock::now();

          ds = kk_mem_mult_stop - kk_mem_mult_start;
          my_experiment_timings["KK MEM: Call"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());

          ds = kk_mem_copyout_stop - kk_mem_copyout_start;
          my_experiment_timings["KK MEM: CopyOut"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());
          break;
        case KokkosSparse::SPGEMMAlgorithm::SPGEMM_KK:
          kk_def_copyout_stop = std::chrono::steady_clock::now();

          ds = kk_def_mult_stop - kk_def_mult_start;
          my_experiment_timings["KK DEFAULT: Call"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());

          ds = kk_def_copyout_stop - kk_def_copyout_start;
          my_experiment_timings["KK DEFAULT: CopyOut"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());
          break;
        default:
          break;
      }
#endif
    }  // scope for Teuchos time monitor
  }    // lib
}

// =========================================================================
// LTG Testing
// =========================================================================
#include "Tpetra_Import_Util2.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_decl.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_def.hpp"

// The LTG kernel is only defined for the Kokkos OpenMP node, so
// its test must only be enabled for OpenMP
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct LTG_Tests {
  static void Multiply_LTG(
      const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &,
      const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &,
      Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &) {}
};

#ifdef HAVE_TPETRA_INST_OPENMP

template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
struct LTG_Tests<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> {
  typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;
  static void Multiply_LTG(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A,
                           const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B,
                           Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C) {
#include <MueLu_UseShortNames.hpp>
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::TimeMonitor;

    Xpetra::UnderlyingLib lib = A.getRowMap()->lib();

    if (lib == Xpetra::UseTpetra) {
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
      RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A, false));
      RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B, false));
      RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C, false));
      RCP<crs_matrix_type> Cnc      = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);

      using no_init_view = Kokkos::ViewAllocateWithoutInitializing;

      // **********************************
      // Copy in the data for LTG

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

      lo_view_t Bcol2Ccol(no_init_view("Bcol2Ccol"), Bview.colMap->getLocalNumElements()), Icol2Ccol;
      const LO colMapSize = static_cast<LO>(Bview.colMap->getLocalNumElements());

      local_map_type Acolmap_local = Aview.colMap->getLocalMap();
      local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
      lo_view_t targetMapToOrigRow(no_init_view("targetMapToOrigRow"), Aview.colMap->getLocalNumElements());
      lo_view_t targetMapToImportRow;

      // Copy in work
      {
        Teuchos::TimeMonitor tm_copyin(*TimeMonitor::getNewTimer("MM LTG: CopyIn"));

        // Because we're in serial...
        Cnc->replaceColMap(Bu->getColMap());

        // Bcol2Ccol is trivial
        Kokkos::parallel_for(
            "Tpetra::mult_A_B_newmatrix::Bcol2Ccol_fill",
            Kokkos::RangePolicy<execution_space, LO>(0, colMapSize),
            KOKKOS_LAMBDA(const LO i) {
              Bcol2Ccol(i) = i;
            });

        // Acol2Brow
        Kokkos::parallel_for(
            "Tpetra::mult_A_B_newmatrix::construct_tables",
            range_type(Aview.colMap->getMinLocalIndex(),
                       Aview.colMap->getMaxLocalIndex() + 1),
            KOKKOS_LAMBDA(const LO i) {
              GO aidx  = Acolmap_local.getGlobalElement(i);
              LO B_LID = Browmap_local.getLocalElement(aidx);
              if (B_LID != LO_INVALID) {
                targetMapToOrigRow(i) = B_LID;
                //        targetMapToImportRow(i) = LO_INVALID;
              } else {
                // This shouldn't happen here
              }
            });

      }  // copyin
      Au->getComm()->barrier();

      // **********************************
      // Multiply
      Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
      params->set("sort entries", false);
      {
        Teuchos::TimeMonitor tm(*TimeMonitor::getNewTimer("MM LTG: Multiply"));
#ifdef USE_DESCRIPTIVE_STATS
        ltg_mult_start = std::chrono::steady_clock::now();
#endif
#ifdef HAVE_VTUNE
        __itt_resume();
#endif
        Tpetra::MatrixMatrix::ExtraKernels::mult_A_B_newmatrix_LowThreadGustavsonKernel(
            Aview, Bview,
            targetMapToOrigRow, targetMapToImportRow,
            Bcol2Ccol, Icol2Ccol,
            *Cnc, Cimport,
            std::string("LTG"), params);
#ifdef HAVE_VTUNE
        __itt_pause();
#endif
#ifdef USE_DESCRIPTIVE_STATS
        ltg_mult_stop = std::chrono::steady_clock::now();
        {
          std::chrono::duration<double> ds;
          ds = ltg_mult_stop - ltg_mult_start;
          my_experiment_timings["LTG: Call"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());
        }
#endif
      }
      Au->getComm()->barrier();
    }
  }
};
#endif  // HAVE OPENMP

// =========================================================================
// Xpetra Testing
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Multiply_Xpetra(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B, Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &C, Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::null) {
#ifdef USE_DESCRIPTIVE_STATS
  xpetra_mult_start = std::chrono::steady_clock::now();
#endif
  Teuchos::FancyOStream fos(Teuchos::rcpFromRef(std::cout));
  Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(A, false, B, false, C, fos, true, true, std::string(""), params);
#ifdef USE_DESCRIPTIVE_STATS
  xpetra_mult_stop = std::chrono::steady_clock::now();
  {
    std::chrono::duration<double> ds;
    ds = ltg_mult_stop - ltg_mult_start;
    my_experiment_timings["Xpetra: Call"].push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(ds).count());
  }
#endif
}

// =========================================================================
// Utilities
// =========================================================================

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void dump_matrix(const std::string &prefix,
                 const std::string &kernel_name,
                 const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> &C) {
  std::ostringstream ss;
  ss << prefix << "_" << kernel_name << ".mtx";
  Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(ss.str(), C);
}

inline bool file_exists(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

static int prior_pos;
void simple_progress_bar(const int iter, std::ostream &out, const int iter_max) {
  double progress = (double)iter / (double)iter_max;

  int barWidth = 70;

  if (iter == 0) prior_pos = 0;

  int new_pos = barWidth * progress;

  if (new_pos > prior_pos) {
    out << "[";
    for (int i = 0; i < barWidth; ++i) {
      if (i < new_pos)
        out << "=";
      else if (i == new_pos)
        out << ">";
      else
        out << " ";
    }
    out << "] " << int(progress * 100.0) << " %\r";
    out.flush();
  }
  prior_pos = new_pos;
}

void print_desc_stats(const std::string &csvFile, const bool writeCSV, std::ostream &out) {
#if USE_DESCRIPTIVE_STATS == 1
  typedef DescriptiveStats::ns ns;

  std::ostringstream oss;
  std::ostringstream oss_header;

  // it is safe to reuse a 'stats' data structure. It will be cleared before resuse
  DescriptiveStats::descriptive_stat_map_type stats;
  const double timer_ratio = double(ns::period::num) / double(ns::period::den);

  DescriptiveStats::print_descriptive_stats(oss_header, stats, timer_ratio, "header", true, true);
  DescriptiveStats::profile_timer(oss);

  // loop through and dump the stats to the stringstream
  for (auto &dataset : my_experiment_timings) {
    auto &measurements           = dataset.second;
    const bool compute_median_CI = true;
    DescriptiveStats::get_descriptive_stats(measurements, measurements.size(), stats, compute_median_CI);
    DescriptiveStats::print_descriptive_stats(oss, stats, timer_ratio, dataset.first, true);
  }

  // write to 'out'
  out << oss_header.str() << oss.str();

  // if rank == 0...
  if (writeCSV) {
    std::ofstream csv_fptr;
    const bool haveExistingCSV = file_exists(csvFile);
    csv_fptr.open(csvFile, std::ios_base::app);

    out << " Writing CSV [" << csvFile << "]...";
    if (!haveExistingCSV) {
      out << "adding header...";
      csv_fptr << oss_header.str();
    }
    csv_fptr << oss.str();
    csv_fptr.close();
    out << "Wrote data... done" << std::endl;
  }
#endif  // if have descriptive stats

}  // end print desc stats

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

#if USE_DESCRIPTIVE_STATS == 1
  typedef DescriptiveStats::Time Time;
  typedef DescriptiveStats::ns ns;
  typedef DescriptiveStats::double_secs double_secs;
  typedef DescriptiveStats::descriptive_stat_map_type descriptive_stat_map_type;
#endif

#ifdef HAVE_VTUNE
#pragma message "Pausing Vtune"
  __itt_pause();
#endif

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    // disable this annoyance
    Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> pOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream &out      = *pOut;
    out.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    clp.throwExceptions(false);
    Xpetra::Parameters xpetraParameters(clp);  // manage parameters of Xpetra

    bool printTimings = true;
    clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
    bool descStats = true;
    clp.setOption("descStats", "nodescStats", &descStats, "track all timings and compute descriptive statistics (uses more memory+time)");
    int nrepeat = 100;
    clp.setOption("nrepeat", &nrepeat, "repeat the experiment N times");
    bool printFNorm = true;
    clp.setOption("fnorm", "nofnorm", &printFNorm, "Compute the Frobenius norm after each multiplication");

    bool printTypes = true;
    clp.setOption("types", "noTypes", &printTypes, "Print the types being used");
    bool printMats = true;
    clp.setOption("printMats", "noPrintMats", &printMats, "Describe matrices");
    bool progressBar = true;
    clp.setOption("progressBar", "noProgressBar", &progressBar, "Show an ascii progress bar while executing.");
    bool writeCSV = false;
    clp.setOption("writeCSV", "noWriteCSV", &writeCSV, "Write a CSV file if descriptive stats is enabled. See --csvFile=NAME");
    std::string csvFile = "desc_stats.csv";
    clp.setOption("csvFile", &csvFile, "CSV file to store desc stats in.");

    bool trackRssHWM = false;
    clp.setOption("trackRssHWM", "noTrackRssHWM", &trackRssHWM, "Track the RSS high water mark (assumes no HP usage)");

    // the kernels
    bool is_serial     = comm->getSize() == 1;
    bool do_viennaCL   = is_serial;
    bool do_mkl        = is_serial;
    bool do_kk_mem     = is_serial;
    bool do_kk_dense   = is_serial;
    bool do_kk_default = is_serial;
    bool do_ltg        = is_serial;
    bool do_xpetra     = true;

#ifndef HAVE_MUELU_VIENNACL
    do_viennaCL = false;
#endif

#ifndef HAVE_MUELU_MKL
    do_mkl = false;
#endif
    clp.setOption("viennaCL", "noviennaCL", &do_viennaCL, "Evaluate ViennaCL");
    clp.setOption("mkl", "nomkl", &do_mkl, "Evaluate MKL SpMM");
    clp.setOption("kk_mem", "nokk_mem", &do_kk_mem, "Evaluate KK Mem");
    clp.setOption("kk_dense", "nokk_dense", &do_kk_dense, "Evaluate KK Dense");
    clp.setOption("kk_default", "nokk_default", &do_kk_default, "Evaluate KK Default");
    clp.setOption("ltg", "noltg", &do_ltg, "Evaluate LTG");
    clp.setOption("xpetra", "noxpetra", &do_ltg, "Evaluate Xpetra");

    std::string dump_result = "";
    clp.setOption("dump_result", &dump_result, "Dump the resulting matrices (once)");

    int kk_team_work_size = 16;

    std::string matrixFileNameA = "A.mm";
    clp.setOption("matrixfileA", &matrixFileNameA, "matrix market file containing matrix");
    std::string matrixFileNameB = "B.mm";
    clp.setOption("matrixfileB", &matrixFileNameB, "matrix market file containing matrix");

    clp.recogniseAllOptions(false);
    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

#ifndef HAVE_MUELU_VIENNACL
    if (do_viennaCL == true) {
      out << "ViennaCL was requested, but this kernel is not available. Disabling..." << endl;
      do_viennaCL = false;
    }
#endif

#ifndef HAVE_MUELU_MKL
    if (do_mkl == true) {
      out << "MKL was requested, but this kernel is not available. Disabling..." << endl;
      do_mkl = false;
    }
#endif

    // simple hack to randomize order of experiments
    std::vector<Experiments> my_experiments;
// add the experiments we will run
#ifdef HAVE_MUELU_VIENNACL
    if (do_viennaCL) my_experiments.push_back(Experiments::ViennaCL);  // ViennaCL
#endif

#ifdef HAVE_MUELU_MKL
    if (do_mkl) my_experiments.push_back(Experiments::MKL_SPMM);  // MKL SPMM
#endif

    // assume these are available
    if (do_kk_mem) my_experiments.push_back(Experiments::KK_MEM);          // KK Mem
    if (do_kk_dense) my_experiments.push_back(Experiments::KK_DENSE);      // KK Dense
    if (do_kk_default) my_experiments.push_back(Experiments::KK_DEFAULT);  // KK Default
    if (do_ltg) my_experiments.push_back(Experiments::LTG);                // LTG
    if (do_xpetra) my_experiments.push_back(Experiments::XPETRA);          // Xpetra
    if (printTypes) {
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
    }

    //    bool dump_matrices = false;
    //    if (dump_result != "") dump_matrices = true;

    // Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);

    std::stringstream ss;

// this stuff was designed to report nested parallelism it isn't super helpful
// but at some point I was working with OMP_NESTED=true, which is very different
// from Kokkos' team based parallelism with OpenMP
#ifdef HAVE_MUELU_OPENMP
#pragma omp parallel
    {
      int num_threads = omp_get_num_threads();

#pragma omp for ordered
      for (int n = 0; n < num_threads; n++) {
#pragma omp critical
        {
          ss << "Level: " << omp_get_level()
             << ", Active Level: " << omp_get_active_level()
             << ", Team #: " << omp_get_team_num() + 1 << " / " << omp_get_num_teams()
             << ", Thread: " << omp_get_thread_num() + 1 << " / " << omp_get_num_threads()
             << ", cpu: " << sched_getcpu()
             << std::endl;
        }

#pragma omp parallel
        {
          if (omp_get_active_level() == omp_get_level()) {
#pragma omp critical
            {
              ss << "Level: " << omp_get_level()
                 << ", Active Level: " << omp_get_active_level()
                 << ", Parent TID: " << omp_get_ancestor_thread_num(1) + 1
                 << ", Team #: " << omp_get_team_num() + 1 << " / " << omp_get_num_teams()
                 << ", Team Size: " << omp_get_team_size(2)
                 << ", Thread: " << omp_get_thread_num() + 1 << " / " << omp_get_num_threads()
                 << ", cpu: " << sched_getcpu()
                 << std::endl;
            }
          } else {
#pragma omp critical
            {
              ss << "Level: " << omp_get_level()
                 << ", Active Level: " << omp_get_active_level()
                 << ", Parent TID: " << omp_get_ancestor_thread_num(1) + 1
                 << ", Team #: " << omp_get_team_num() + 1 << " / " << omp_get_num_teams()
                 << ", Team Size: " << omp_get_team_size(2)
                 << ", Thread: " << omp_get_thread_num() + 1 << " / " << omp_get_num_threads()
                 << ", cpu: " << sched_getcpu()
                 << std::endl;
            }
          }
        }
      }
    }

    out << ss.str();
    out.flush();
#endif  // have openmp

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time")));

    comm->barrier();

    RCP<Matrix> A = Xpetra::IO<SC, LO, GO, Node>::Read(std::string(matrixFileNameA), lib, comm);
    RCP<Matrix> B;

    if (matrixFileNameB == "dupe") {
      B = A;

      if (verbose) out << "duplicating A" << endl;
    } else if (matrixFileNameB == "alias") {
      B = A;
      if (verbose) out << "Aliasing A" << endl;
    } else {
      B = Xpetra::IO<SC, LO, GO, Node>::Read(std::string(matrixFileNameB), lib, comm);
    }

    globalTimeMonitor = Teuchos::null;
    comm->barrier();

    if (verbose) {
      out << "Matrix Read complete." << endl;
    }

    if (printMats) {
      out << "Matrix A:" << endl
          << *A
          << "========================================================" << endl
          << "Matrix B:" << endl
          << *B
          << "========================================================" << endl;
      //      if (do_serial)
      //      out << "Matrix A:" << endl
      //          << *Aserial
      //          << "========================================================" << endl
      //          << "Matrix B:" << endl
      //          << *Bserial
      //          << "========================================================" << endl;
    }

    // random source
    std::random_device rd;
    std::mt19937 random_source(rd());

    // no need for a barrier, because the randomization process uses a collective.
    if (!my_experiments.empty()) {
      std::vector<double> mem_dummy;
      double prior_fnorm = std::numeric_limits<double>::max();
      int fnorm_failures = 0;

      globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Kernels Total Time")));
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
          RCP<Matrix> C;
          std::string kernel_name = to_string(experiment_id);
#if USE_DESCRIPTIVE_STATS == 1
          auto &timing_vector = my_experiment_timings[kernel_name + ": Total"];
          auto &mem_vector    = (trackRssHWM) ? my_experiment_timings[kernel_name + ": RssHWM"] : mem_dummy;
#endif

          switch (experiment_id) {
            // ViennaCL
            case Experiments::ViennaCL:
#ifdef HAVE_MUELU_VIENNACL
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM ViennaCL: Total"));
                DescriptiveTime(Multiply_ViennaCL(*A, *B, *C));
              }
#endif
              break;
            // MKL_SPMM
            case Experiments::MKL_SPMM:
#ifdef HAVE_MUELU_MKL
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM MKL: Total"));
                DescriptiveTime(Multiply_MKL_SPMM(*A, *B, *C));
              }
#endif
              break;
            // KK Algorithms (KK Memory)
            case Experiments::KK_MEM:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM SPGEMM_KK_MEMORY: Total"));
                DescriptiveTime(Multiply_KokkosKernels(*A, *B, *C, std::string("SPGEMM_KK_MEMORY"), kk_team_work_size));
              }
              break;
            // KK Algorithms (KK Dense)
            case Experiments::KK_DENSE:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM SPGEMM_KK_DENSE: Total"));
                DescriptiveTime(Multiply_KokkosKernels(*A, *B, *C, std::string("SPGEMM_KK_DENSE"), kk_team_work_size));
              }
              break;
            // KK Algorithms (KK Default)
            case Experiments::KK_DEFAULT:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM SPGEMM_KK: Total"));
                DescriptiveTime(Multiply_KokkosKernels(*A, *B, *C, std::string("SPGEMM_KK"), kk_team_work_size));
              }
              break;
            // LTG
            case Experiments::LTG:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM LTG: Total"));
                using ltg_tester = LTG_Tests<SC, LO, GO, Node>;
                DescriptiveTime(ltg_tester::Multiply_LTG(*A, *B, *C));
              }
              break;
            // Xpetra
            case Experiments::XPETRA:
              C = Xpetra::MatrixFactory<SC, LO, GO, Node>::Build(A->getRowMap(), 0);
              {
                TimeMonitor t(*TimeMonitor::getNewTimer("MM XPETRA: Total"));
                DescriptiveTime(Multiply_Xpetra(*A, *B, C));
              }
              break;

            default:
              std::cerr << "Unknown experiment ID encountered: " << (int)experiment_id << std::endl;
          }
          comm->barrier();

          if (printFNorm) {
            const auto current_fnorm = C->getFrobeniusNorm();
            const bool fnorm_error   = (std::abs(current_fnorm - prior_fnorm) / prior_fnorm) > 1.e-14;
            if (i == 0) {
              out << kernel_name << " "
                  << "Current Fnorm: " << current_fnorm
                  << std::endl;
            } else if (fnorm_error) {
              out << kernel_name << " "
                  << "Current Fnorm: " << current_fnorm << ", Difference from Prior: "
                  << std::abs(prior_fnorm - current_fnorm)
                  << std::endl;
              fnorm_failures++;
            }
            prior_fnorm = current_fnorm;
          }
          C = Teuchos::null;
        }  // end random exp loop

        if (progressBar) simple_progress_bar(i, out, nrepeat);
      }  // end repeat
    }    // end ! my_experiments.empty()

    globalTimeMonitor = Teuchos::null;

    if (printTimings) {
      TimeMonitor::summarize(comm.ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
#if USE_DESCRIPTIVE_STATS == 1

      if (descStats) print_desc_stats(csvFile, writeCSV, out);
#endif
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
