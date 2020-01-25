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
#include <Teuchos_StackedTimer.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
// MueLu
#include "MueLu.hpp"
#include "MueLu_TestHelpers.hpp"
#include <MatrixLoad.hpp>

#if defined(HAVE_MUELU_TPETRA)
#include "Xpetra_TpetraMultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "KokkosSparse_spmv.hpp"
#endif

#if defined(HAVE_MUELU_CUSPARSE)
#include "cusparse.h"
#endif

// =========================================================================
// Support Routines
// =========================================================================
template<class View1, class View2>
inline void copy_view_n(const int n, const View1 x1, View2 x2) {
  Kokkos::parallel_for(n,KOKKOS_LAMBDA(const size_t i) {
      x2[i] = x1[i];
    });
}

template<class View1, class View2>
inline void copy_view(const View1 x1, View2 x2) {
  copy_view_n(x1.extent(0),x1,x2);
}



template <class V1, class V2>
void print_crs_graph(std::string name, const V1 rowptr, const V2 colind) {
  printf("%s rowptr[%d] = ",name.c_str(),rowptr.extent(0));
  for(size_t i=0; i<rowptr.extent(0); i++)
    printf(" %d",(int)rowptr[i]);
  printf("\n%s colind[%d] = ",name.c_str(),colind.extent(0));
  for(size_t i=0; i<colind.extent(0); i++)
    printf(" %d",(int)colind[i]);
  printf("\n");
}



// =========================================================================
// MAGMASparse Testing
// =========================================================================
#if defined(HAVE_MUELU_MAGMASPARSE) && defined(HAVE_MUELU_TPETRA)
#include <magma_v2.h>
#include <magmasparse.h>
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class MagmaSparse_SpmV_Pack {
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> vector_type;

public:
  MagmaSparse_SpmV_Pack (const crs_matrix_type& A,
                      const vector_type& X,
                            vector_type& Y)
  {}

 ~MagmaSparse_SpmV_Pack() {};

  bool spmv(const Scalar alpha, const Scalar beta) { return (true); }
};

template<typename LocalOrdinal, typename GlobalOrdinal>
class MagmaSparse_SpmV_Pack<double,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosCudaWrapperNode> {
  // typedefs shared among other TPLs
  typedef double Scalar;
  typedef typename Kokkos::Compat::KokkosCudaWrapperNode Node;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> vector_type;
  typedef typename crs_matrix_type::local_matrix_type    KCRS;
  typedef typename KCRS::StaticCrsGraphType              graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type     c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::const_type     c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type     scalar_view_t;
  typedef typename Node::device_type device_type;

  typedef typename Kokkos::View<int*,
                                typename lno_nnz_view_t::array_layout,
                                typename lno_nnz_view_t::device_type> int_array_type;

public:
  MagmaSparse_SpmV_Pack (const crs_matrix_type& A,
                         const vector_type& X,
                               vector_type& Y)
  {
    // data access common to other TPLs
    const KCRS & Amat = A.getLocalMatrix();
    c_lno_view_t Arowptr = Amat.graph.row_map;
    c_lno_nnz_view_t Acolind = Amat.graph.entries;
    const scalar_view_t Avals = Amat.values;

    // type conversion
    Arowptr_int = int_array_type("Arowptr", Arowptr.extent(0));
    // copy the ordinals into the local view (type conversion)
    copy_view(Arowptr,Arowptr_int);


    m   = static_cast<int>(A.getNodeNumRows());
    n   = static_cast<int>(A.getNodeNumCols());
    nnz = static_cast<int>(Acolind.extent(0));
    vals   = reinterpret_cast<Scalar*>(Avals.data());
    cols   = const_cast<int*>(reinterpret_cast<const int*>(Acolind.data()));
    rowptr = reinterpret_cast<int*>(Arowptr_int.data());

    auto X_lcl = X.template getLocalView<device_type> ();
    auto Y_lcl = Y.template getLocalView<device_type> ();
    x = reinterpret_cast<Scalar*>(X_lcl.data());
    y = reinterpret_cast<Scalar*>(Y_lcl.data());

    magma_init ();
    int device;
    magma_getdevice( &device );
    magma_queue_create( device, &queue );

    // create the magma data container
    magma_dev_Acrs = {Magma_CSR};
    magma_Acrs     = {Magma_CSR};
    magma_dev_x    = {Magma_DENSE};
    magma_dev_y    = {Magma_DENSE};

    //
    magma_dvset_dev (m, 1,         // assume mx1 size
                     x,            // ptr to data on the device
                     &magma_dev_x, // magma vector to populate
                     queue);

    magma_dvset_dev (m, 1,         // assume mx1 size
                     y,            // ptr to data on the device
                     &magma_dev_y, // magma vector to populate
                     queue);
    magma_dcsrset( m, n, rowptr, cols, vals, &magma_Acrs, queue );
    magma_dmtransfer( magma_Acrs, &magma_dev_Acrs, Magma_DEV, Magma_DEV, queue );

  }

  ~MagmaSparse_SpmV_Pack()
  {
    // we dont' own the data... not sure what magma does here...
    magma_dmfree(&magma_dev_x, queue);
    magma_dmfree(&magma_dev_y, queue);
    magma_dmfree(&magma_dev_Acrs, queue);

    magma_finalize ();
  };

  bool spmv(const Scalar alpha, const Scalar beta)
  {
    magma_d_spmv( 1.0, magma_dev_Acrs, magma_dev_x, 0.0, magma_dev_y, queue );
  }

private:
  magma_d_matrix magma_dev_Acrs;
  magma_d_matrix magma_Acrs;
  magma_d_matrix magma_dev_x;
  magma_d_matrix magma_dev_y;
  magma_queue_t queue;

  int m = -1;
  int n = -1;
  int nnz = -1;
  Scalar * vals  = nullptr; // aliased
  int * cols     = nullptr; // aliased
  int * rowptr   = nullptr; // copied
  Scalar * x     = nullptr; // aliased
  Scalar * y     = nullptr; // aliased

  // handles to the copied data
  int_array_type Arowptr_int;
};
#endif


// =========================================================================
// CuSparse Testing
// =========================================================================
#if defined(HAVE_MUELU_CUSPARSE) && defined(HAVE_MUELU_TPETRA)
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class CuSparse_SpmV_Pack {
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> vector_type;

public:
  CuSparse_SpmV_Pack (const crs_matrix_type& A,
                      const vector_type& X,
                            vector_type& Y)
  {}

 ~CuSparse_SpmV_Pack() {};

  cusparseStatus_t spmv(const Scalar alpha, const Scalar beta) { return CUSPARSE_STATUS_SUCCESS; }
};

template<typename LocalOrdinal, typename GlobalOrdinal>
class CuSparse_SpmV_Pack<double,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosCudaWrapperNode> {
  // typedefs shared among other TPLs
  typedef double Scalar;
  typedef typename Kokkos::Compat::KokkosCudaWrapperNode Node;
  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> vector_type;
  typedef typename crs_matrix_type::local_matrix_type    KCRS;
  typedef typename KCRS::StaticCrsGraphType              graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type     c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::const_type     c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type     scalar_view_t;
  typedef typename Node::device_type device_type;

  typedef typename Kokkos::View<int*,
                                typename lno_nnz_view_t::array_layout,
                                typename lno_nnz_view_t::device_type> cusparse_int_type;
public:
  CuSparse_SpmV_Pack (const crs_matrix_type& A,
                      const vector_type& X,
                            vector_type& Y)
  {
    // data access common to other TPLs
    const KCRS & Amat = A.getLocalMatrix();
    c_lno_view_t Arowptr = Amat.graph.row_map;
    c_lno_nnz_view_t Acolind = Amat.graph.entries;
    const scalar_view_t Avals = Amat.values;


    Arowptr_cusparse = cusparse_int_type("Arowptr", Arowptr.extent(0));
    Acolind_cusparse = cusparse_int_type("Acolind", Acolind.extent(0));
    // copy the ordinals into the local view (type conversion)
    copy_view(Arowptr,Arowptr_cusparse);
    copy_view(Acolind,Acolind_cusparse);


    m   = static_cast<int>(A.getNodeNumRows());
    n   = static_cast<int>(A.getNodeNumCols());
    nnz = static_cast<int>(Acolind_cusparse.extent(0));
    vals   = reinterpret_cast<Scalar*>(Avals.data());
    cols   = reinterpret_cast<int*>(Acolind_cusparse.data());
    rowptr = reinterpret_cast<int*>(Arowptr_cusparse.data());

    auto X_lcl = X.template getLocalView<device_type> ();
    auto Y_lcl = Y.template getLocalView<device_type> ();
    x = reinterpret_cast<Scalar*>(X_lcl.data());
    y = reinterpret_cast<Scalar*>(Y_lcl.data());

    /* Get handle to the CUBLAS context */
    cublasStatus_t cublasStatus;
    cublasStatus = cublasCreate(&cublasHandle);

    //checkCudaErrors(cublasStatus);

    /* Get handle to the CUSPARSE context */
    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&cusparseHandle);

    //checkCudaErrors(cusparseStatus);
    cusparseStatus = cusparseCreateMatDescr(&descrA);

    //checkCudaErrors(cusparseStatus);

    cusparseSetMatType(descrA,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ZERO);
  }

  ~CuSparse_SpmV_Pack () {
    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);
  }

  cusparseStatus_t spmv(const Scalar alpha, const Scalar beta) {
    // compute: y = alpha*Ax + beta*y
    //cusparseDcsrmv(cusparseHandle_t handle,
    //               cusparseOperation_t transA,
    //               int m,
    //               int n,
    //               int nnz,
    //               const double *alpha,
    //               // CUSPARSE_MATRIX_TYPE_GENERAL
    //               const cusparseMatDescr_t descrA,
    //               const double *csrValA,
    //               const int    *csrRowPtrA,
    //               const int    *csrColIndA,
    //               const double *x,
    //               const double *beta,
    //                     double *y)
    cusparseStatus_t rc;
    rc = cusparseDcsrmv(
        cusparseHandle,
        transA,
        m, n, nnz,
        &alpha,
        descrA,
        vals,
        rowptr,
        cols,
        x,
        &beta,
        y);

    return (rc);
  }
private:
  cublasHandle_t     cublasHandle   = 0;
  cusparseHandle_t   cusparseHandle = 0;
  cusparseMatDescr_t descrA         = 0;

  // CUSPARSE_OPERATION_NON_TRANSPOSE
  // CUSPARSE_OPERATION_TRANSPOSE
  // CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE
  cusparseOperation_t transA = CUSPARSE_OPERATION_NON_TRANSPOSE;
  int m = -1;
  int n = -1;
  int nnz = -1;
  Scalar * vals  = nullptr; // aliased
  int * cols     = nullptr; // copied
  int * rowptr   = nullptr; // copied
  Scalar * x     = nullptr; // aliased
  Scalar * y     = nullptr; // aliased

  // handles to the copied data
  cusparse_int_type Arowptr_cusparse;
  cusparse_int_type Acolind_cusparse;
};
#endif

// =========================================================================
// MKL Testing
// =========================================================================
#if defined(HAVE_MUELU_MKL) && defined(HAVE_MUELU_TPETRA)
#include "mkl.h"


std::string mkl_error(sparse_status_t code) {
  switch(code) {
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

struct matrix_descr mkl_descr;

void MV_MKL(sparse_matrix_t & AMKL, double * x, double * y) {
  //sparse_status_t mkl_sparse_d_mv (sparse_operation_t operation, double alpha, const sparse_matrix_t A, struct matrix_descr descr, const double *x, double beta, double *y);

  mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,1.0,AMKL,mkl_descr,x,0.0,y);
}

#endif



// =========================================================================
// Tpetra Kernel Testing
// =========================================================================
#if defined(HAVE_MUELU_TPETRA)
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MV_Tpetra(const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,  const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &x,   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &y) {
  A.apply(x,y);
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MV_KK(const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,  const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &x,   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &y) {
  typedef typename Node::device_type device_type;
  const auto& AK = A.getLocalMatrix();
  auto X_lcl = x.template getLocalView<device_type> ();
  auto Y_lcl = y.template getLocalView<device_type> ();
  KokkosSparse::spmv(KokkosSparse::NoTranspose,Teuchos::ScalarTraits<Scalar>::one(),AK,X_lcl,Teuchos::ScalarTraits<Scalar>::zero(),Y_lcl);
}
#endif

// =========================================================================
// =========================================================================
// =========================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib& lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using std::endl;


  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    //    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> pOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *pOut;
    out.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    GO nx = 100, ny = 100, nz = 50;
    Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

    bool binaryFormat = false;  clp.setOption("binary", "ascii",  &binaryFormat,   "print timings to screen");
    std::string matrixFile;     clp.setOption("matrixfile",       &matrixFile,     "matrix market file containing matrix");
    std::string rowMapFile;     clp.setOption("rowmap",           &rowMapFile,     "map data file");
    std::string colMapFile;     clp.setOption("colmap",           &colMapFile,     "colmap data file");
    std::string domainMapFile;  clp.setOption("domainmap",        &domainMapFile,  "domainmap data file");
    std::string rangeMapFile;   clp.setOption("rangemap",         &rangeMapFile,   "rangemap data file");

    bool printTimings = true;   clp.setOption("timings", "notimings",  &printTimings, "print timings to screen");
    int  nrepeat      = 100;    clp.setOption("nrepeat",               &nrepeat,      "repeat the experiment N times");

    bool describeMatrix = true; clp.setOption("showmatrix", "noshowmatrix",  &describeMatrix, "describe matrix");
    bool useStackedTimer = false; clp.setOption("stackedtimer", "nostackedtimer",  &useStackedTimer, "use stacked timer");
    // the kernels
    bool do_mkl      = true;
    bool do_tpetra   = true;
    bool do_kk       = true;
    bool do_cusparse = true;
    bool do_magmasparse = true;

    bool report_error_norms = false;
    // handle the TPLs
    #if ! defined(HAVE_MUELU_MKL)
      do_mkl = false;
    #endif
    #if ! defined(HAVE_MUELU_CUSPARSE)
      do_cusparse = false;
    #endif
    #if ! defined(HAVE_MUELU_MAGMASPARSE)
      do_magmasparse = false;
    #endif

    clp.setOption("mkl",      "nomkl",      &do_mkl,        "Evaluate MKL");
    clp.setOption("tpetra",   "notpetra",   &do_tpetra,     "Evaluate Tpetra");
    clp.setOption("kk",       "nokk",       &do_kk,         "Evaluate KokkosKernels");
    clp.setOption("cusparse", "nocusparse", &do_cusparse,   "Evaluate CuSparse");
    clp.setOption("magamasparse", "nomagmasparse", &do_magmasparse,   "Evaluate MagmaSparse");

    clp.setOption("report_error_norms", "noreport_error_norms", &report_error_norms,   "Report L2 norms for the solution");

    std::ostringstream galeriStream;
    std::string rhsFile,coordFile,nullFile, materialFile; //unused
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
    typedef Xpetra::MultiVector<real_type,LO,GO,NO> RealValuedMultiVector;
    RCP<RealValuedMultiVector> coordinates;
    RCP<MultiVector> nullspace, material, x, b;
    RCP<Matrix> A;
    RCP<const Map> map;

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    if (do_kk && comm->getSize() > 1)
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"The Kokkos-Kernels matvec kernel cannot be run with more than one rank.");

    // Load the matrix off disk (or generate it via Galeri), assuming only one right hand side is loaded.
    MatrixLoad<SC,LO,GO,NO>(comm, lib, binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile, domainMapFile, rangeMapFile, coordFile, nullFile, materialFile, map, A, coordinates, nullspace, material, x, b, 1, galeriParameters, xpetraParameters, galeriStream);

    #ifndef HAVE_MUELU_MKL
    if (do_mkl) {
      out << "MKL was requested, but this kernel is not available. Disabling..." << endl;
      do_mkl = false;
    }
    #endif

    #if ! defined(HAVE_MUELU_CUSPARSE)
    if (do_cusparse) {
      out << "CuSparse was requested, but this kernel is not available. Disabling..." << endl;
      do_cusparse = false;
    }
    #else
    if(! std::is_same<NO, Kokkos::Compat::KokkosCudaWrapperNode>::value) do_cusparse = false;
    #endif

    #if ! defined(HAVE_MUELU_MAGMASPARSE)
    if (do_magmasparse) {
      out << "MagmaSparse was requested, but this kernel is not available. Disabling..." << endl;
      do_magmasparse = false;
    }
    #endif

    // simple hack to randomize order of experiments
    enum class Experiments { MKL=0, TPETRA, KK, CUSPARSE, MAGMASPARSE };
    const char * const experiment_id_to_string[] = {
        "MKL        ",
        "Tpetra     ",
        "KK         ",
        "CuSparse   ",
        "MagmaSparse"};
    std::vector<Experiments> my_experiments;
    // add the experiments we will run

    #ifdef HAVE_MUELU_MKL
    if (do_mkl) my_experiments.push_back(Experiments::MKL);   // MKL
    #endif

    #ifdef HAVE_MUELU_CUSPARSE
    if (do_cusparse) my_experiments.push_back(Experiments::CUSPARSE);   // CuSparse
    #endif

    #ifdef HAVE_MUELU_MAGMASPARSE
    if (do_magmasparse) my_experiments.push_back(Experiments::MAGMASPARSE);   // MagmaSparse
    #endif

    // assume these are available
    if (do_tpetra)  my_experiments.push_back(Experiments::TPETRA);     // Tpetra
    if (do_kk)     my_experiments.push_back(Experiments::KK);     // KK


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
        << "Vector:        " << Teuchos::demangleName(typeid(MultiVector).name()) << endl
        << "Hierarchy:     " << Teuchos::demangleName(typeid(Hierarchy).name()) << endl
        << "========================================================" << endl;

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_TPETRA_INST_OPENMP)
    out<< "Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency() = "<<Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency()<<endl
       << "========================================================" << endl;
#endif

    // =========================================================================
    // Problem construction
    // =========================================================================
    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
    RCP<TimeMonitor> globalTimeMonitor;
    if (useStackedTimer)
      stacked_timer = rcp(new Teuchos::StackedTimer("MueLu_MatvecKernelDriver"));
    else
      globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time")));

    comm->barrier();


    RCP<MultiVector> y = Xpetra::VectorFactory<SC,LO,GO,Node>::Build(A->getRowMap());
    RCP<MultiVector> y_baseline = Xpetra::VectorFactory<SC,LO,GO,Node>::Build(A->getRowMap());
    x->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    y->putScalar(Teuchos::ScalarTraits<Scalar>::nan());
    y_baseline->putScalar(Teuchos::ScalarTraits<Scalar>::nan());


#ifdef HAVE_MUELU_TPETRA
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> vector_type;
    RCP<const crs_matrix_type> At;
    vector_type xt;
    vector_type yt;

    At = Utilities::Op2TpetraCrs(A);
    const crs_matrix_type &Att = *At;
    xt = Xpetra::toTpetra(*x);
    yt = Xpetra::toTpetra(*y);

    size_t l_permutes = 0, g_permutes = 0;
    if(!At->getGraph()->getImporter().is_null()) {
      l_permutes = At->getGraph()->getImporter()->getNumPermuteIDs();
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM,1,&l_permutes,&g_permutes);
    }
    if(!comm->getRank()) printf("DEBUG: A's importer has %d total permutes globally\n",(int)g_permutes);

  #if defined(HAVE_MUELU_CUSPARSE)
    typedef CuSparse_SpmV_Pack<SC,LO,GO,Node> CuSparse_thing_t;
    CuSparse_thing_t cusparse_spmv(*At, xt, yt);
  #endif

  #if defined(HAVE_MUELU_MAGMASPARSE)
    typedef MagmaSparse_SpmV_Pack<SC,LO,GO,Node> MagmaSparse_thing_t;
    MagmaSparse_thing_t magmasparse_spmv(*At, xt, yt);
  #endif
  #if defined(HAVE_MUELU_MKL)
    // typedefs shared among other TPLs
    typedef typename crs_matrix_type::local_matrix_type    KCRS;
    typedef typename KCRS::StaticCrsGraphType              graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::row_map_type::const_type     c_lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename graph_t::entries_type::const_type     c_lno_nnz_view_t;
    typedef typename KCRS::values_type::non_const_type     scalar_view_t;
    typedef typename Node::device_type device_type;

    // data access common to other TPLs
    const KCRS & Amat = At->getLocalMatrix();
    c_lno_view_t Arowptr = Amat.graph.row_map;
    c_lno_nnz_view_t Acolind = Amat.graph.entries;
    const scalar_view_t Avals = Amat.values;

    // MKL specific things
    sparse_matrix_t mkl_A;
    typedef typename Kokkos::View<MKL_INT*,typename lno_nnz_view_t::array_layout,typename lno_nnz_view_t::device_type> mkl_int_type;

    mkl_int_type ArowptrMKL("Arowptr", Arowptr.extent(0));
    mkl_int_type AcolindMKL("Acolind", Acolind.extent(0));
    copy_view(Arowptr,ArowptrMKL);
    copy_view(Acolind,AcolindMKL);
    double * mkl_xdouble = nullptr;
    double * mkl_ydouble = nullptr;
    mkl_descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    if(Kokkos::Impl::is_same<Scalar,double>::value) {
      mkl_sparse_d_create_csr(&mkl_A,
                              SPARSE_INDEX_BASE_ZERO,
                              At->getNodeNumRows(),
                              At->getNodeNumCols(),
                              ArowptrMKL.data(),
                              ArowptrMKL.data()+1,
                              AcolindMKL.data(),
                              (double*)Avals.data());
      auto X_lcl = xt.template getLocalView<device_type> ();
      auto Y_lcl = yt.template getLocalView<device_type> ();
      mkl_xdouble = (double*)X_lcl.data();
      mkl_ydouble = (double*)Y_lcl.data();
    }
    else
      throw std::runtime_error("MKL Type Mismatch");

  #endif // end MKL
#endif


    globalTimeMonitor = Teuchos::null;
    comm->barrier();

    out << "Matrix Read complete." << endl;
    if (describeMatrix) {
        out << "Matrix A:" << endl
            << *A
            << "========================================================" << endl;
    }
    // save std::cout formating
    std::ios cout_default_fmt_flags(NULL);
    cout_default_fmt_flags.copyfmt(std::cout);

    // random source
    std::random_device rd;
    std::mt19937 random_source (rd());

    #ifdef HAVE_MUELU_TPETRA
    // compute the baseline
    vector_type yt_baseline = Xpetra::toTpetra(*y_baseline);
    if (report_error_norms) MV_Tpetra(*At,xt,yt_baseline);
    const bool error_check_y = true;
    #else
    const bool error_check_y = false;
    #endif
    std::vector<typename Teuchos::ScalarTraits< Scalar >::magnitudeType> dummy;
    dummy.resize(1);
    Teuchos::ArrayView< typename Teuchos::ScalarTraits< Scalar >::magnitudeType > y_norms (dummy.data(), 1);

    // no need for a barrier, because the randomization process uses a collective.
    if ( !my_experiments.empty() ) {


    for (int i=0; i<nrepeat; i++) {

      // randomize the experiments
      if (comm->getRank() == 0 ) {
        std::shuffle(my_experiments.begin(), my_experiments.end(), random_source);
      }

      // Broadcast this ordering to the other processes
      comm->broadcast(0,
                      static_cast<int>(sizeof(Experiments::TPETRA)*my_experiments.size()),
                      reinterpret_cast<char *>(my_experiments.data()));

      // loop over the randomized experiments
      for (const auto& experiment_id : my_experiments) {
        switch (experiment_id) {

        #ifdef HAVE_MUELU_MKL
        // MKL
        case Experiments::MKL:
        {
            TimeMonitor t(*TimeMonitor::getNewTimer("MV MKL: Total"));
            MV_MKL(mkl_A,mkl_xdouble,mkl_ydouble);
        }
          break;
        #endif

        #ifdef HAVE_MUELU_TPETRA
        // KK Algorithms
        case Experiments::KK:
        {
           TimeMonitor t(*TimeMonitor::getNewTimer("MV KK: Total"));
           MV_KK(Att,xt,yt);
        }
          break;
        // Tpetra
        case Experiments::TPETRA:
        {
           TimeMonitor t(*TimeMonitor::getNewTimer("MV Tpetra: Total"));
           MV_Tpetra(*At,xt,yt);
        }
          break;
        #endif

        #ifdef HAVE_MUELU_CUSPARSE
        // MKL
        case Experiments::CUSPARSE:
        {
           const Scalar alpha = 1.0;
           const Scalar beta = 0.0;
           TimeMonitor t(*TimeMonitor::getNewTimer("MV CuSparse: Total"));
           cusparse_spmv.spmv(alpha,beta);
        }
          break;
        #endif

        #ifdef HAVE_MUELU_MAGMASPARSE
        // Magma CSR
        case Experiments::MAGMASPARSE:
        {
           const Scalar alpha = 1.0;
           const Scalar beta = 0.0;
           TimeMonitor t(*TimeMonitor::getNewTimer("MV MagmaSparse: Total"));
           magmasparse_spmv.spmv(alpha,beta);
        }
          break;
        #endif
        default:
          std::cerr << "Unknown experiment ID encountered: " << (int) experiment_id << std::endl;
        }
        //TODO: add a correctness check
        // For now, all things alias x/y, so we can test yt (flakey and scary, but you only live once)
        if (error_check_y && report_error_norms) {
          // compute ||y||_2
          y_norms[0] = -1;
          y->norm2(y_norms);
          const auto y_norm2 = y_norms[0];

          y_norms[0] = -1;
          yt.norm2(y_norms);
          const auto y_mv_norm2 = y_norms[0];

          y->update(-Teuchos::ScalarTraits<Scalar>::one(), *y_baseline, Teuchos::ScalarTraits<Scalar>::one());

          y_norms[0] = -1;
          y->norm2(y_norms);
          const auto y_err =  y_norms[0];

          y->putScalar(Teuchos::ScalarTraits<Scalar>::nan());;

          y_norms[0] = -1;
          y_baseline->norm2(y_norms);
          const auto y_baseline_norm2 = y_norms[0];

          y_norms[0] = -1;
          yt.norm2(y_norms);
          const auto y_mv_norm2_next_itr = y_norms[0];

          std::cout << "ExperimentID: " << experiment_id_to_string[(int) experiment_id] << ", ||y-y_hat||_2 = "
                    << std::setprecision(std::numeric_limits<Scalar>::digits10 + 1)
                    << std::scientific << y_err
                    << ", ||y||_2 = " << y_norm2
                    << ", ||y_baseline||_2 = " << y_baseline_norm2
                    << ", ||y_ptr|| == ||y_mv||:  " << std::boolalpha << (y_mv_norm2 == y_norm2)
                    << ", setting y to nan ... ||y||_2 for next iter: " << y_mv_norm2_next_itr
                    << "\n";
        }
        comm->barrier();
      }// end random exp loop
    } // end repeat
    } // end ! my_experiments.empty()
    // restore the IO stream
    std::cout.copyfmt(cout_default_fmt_flags);

    if (useStackedTimer) {
      stacked_timer->stop("MueLu_MatvecKernelDriver");
      Teuchos::StackedTimer::OutputOptions options;
      options.output_fraction = options.output_histogram = options.output_minmax = true;
      stacked_timer->report(out, comm, options);
    } else {
      TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
    }

#if defined(HAVE_MUELU_MKL)
    mkl_sparse_destroy(mkl_A);
#endif

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} //main




//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}

