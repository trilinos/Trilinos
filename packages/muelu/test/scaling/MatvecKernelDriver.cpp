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
  auto AK = A.getLocalMatrix();
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

    #ifndef HAVE_MUELU_MKL
      do_mkl = false;
    #endif
    clp.setOption("mkl",      "nomkl",    &do_mkl,        "Evaluate MKL");
    clp.setOption("tpetra",   "notpetra", &do_tpetra,     "Evaluate Tpetra");
    clp.setOption("kk",       "nokk",     &do_kk,         "Evaluate KokkosKernels");

    std::ostringstream galeriStream;
    std::string rhsFile,coordFile,nullFile; //unused
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
    typedef Xpetra::MultiVector<real_type,LO,GO,NO> RealValuedMultiVector;
    RCP<RealValuedMultiVector> coordinates;
    RCP<MultiVector> nullspace, x, b;
    RCP<Matrix> A;
    RCP<const Map> map;
   
    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    if (do_kk && comm->getSize() > 1)
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"The Kokkos-Kernels matvec kernel cannot be run with more than one rank.");

    // Load the matrix off disk (or generate it via Galeri)
    MatrixLoad<SC,LO,GO,NO>(comm, lib, binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile, domainMapFile, rangeMapFile, coordFile, nullFile, map, A, coordinates, nullspace, x, b, galeriParameters, xpetraParameters, galeriStream);

#ifndef HAVE_MUELU_MKL
    if (do_mkl) {
      out << "MKL was requested, but this kernel is not available. Disabling..." << endl;
      do_mkl = false;
    }
#endif

    // simple hack to randomize order of experiments
    enum class Experiments { MKL=0, TPETRA, KK };
    std::vector<Experiments> my_experiments;
    // add the experiments we will run
  
    #ifdef HAVE_MUELU_MKL
    if (do_mkl) my_experiments.push_back(Experiments::MKL);   // MKL
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
    x->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    y->putScalar(Teuchos::ScalarTraits<Scalar>::zero());


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

#ifdef HAVE_MUELU_MKL
    sparse_matrix_t AMKL;
    typedef typename crs_matrix_type::local_matrix_type    KCRS;
    typedef typename KCRS::StaticCrsGraphType              graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::row_map_type::const_type     c_lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename graph_t::entries_type::const_type     c_lno_nnz_view_t;
    typedef typename KCRS::values_type::non_const_type     scalar_view_t;
    typedef typename Kokkos::View<MKL_INT*,typename lno_nnz_view_t::array_layout,typename lno_nnz_view_t::device_type> mkl_int_type;    
    typedef typename Node::device_type device_type;
    const KCRS & Amat = At->getLocalMatrix();  
    c_lno_view_t Arowptr = Amat.graph.row_map;
    c_lno_nnz_view_t Acolind = Amat.graph.entries;
    const scalar_view_t Avals = Amat.values;
    mkl_int_type ArowptrMKL("Arowptr",Arowptr.extent(0));
    mkl_int_type AcolindMKL("Acolind",Acolind.extent(0));
    copy_view(Arowptr,ArowptrMKL);
    copy_view(Acolind,AcolindMKL);
    double * xdouble=0, * ydouble=0;
    mkl_descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    if(Kokkos::Impl::is_same<Scalar,double>::value) {
      mkl_sparse_d_create_csr(&AMKL, SPARSE_INDEX_BASE_ZERO, At->getNodeNumRows(), At->getNodeNumCols(), ArowptrMKL.data(),ArowptrMKL.data()+1,AcolindMKL.data(),(double*)Avals.data());
      auto X_lcl = xt.template getLocalView<device_type> ();
      auto Y_lcl = yt.template getLocalView<device_type> ();
      xdouble = (double*)X_lcl.data();
      ydouble = (double*)Y_lcl.data();
    }
  else
    throw std::runtime_error("MKL Type Mismatch");

#endif
#endif


    globalTimeMonitor = Teuchos::null;
    comm->barrier();

    out << "Matrix Read complete." << endl;
    if (describeMatrix) {
        out << "Matrix A:" << endl
            << *A
            << "========================================================" << endl;
    }

    // random source
    std::random_device rd;
    std::mt19937 random_source (rd());

    

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
 
        // MKL 
        case Experiments::MKL:
          {
#ifdef HAVE_MUELU_MKL
            TimeMonitor t(*TimeMonitor::getNewTimer("MV MKL: Total"));            
            MV_MKL(AMKL,xdouble,ydouble);
#endif
          }
          break;
        // KK Algorithms
        case Experiments::KK:
          {
#ifdef HAVE_MUELU_TPETRA
            TimeMonitor t(*TimeMonitor::getNewTimer("MV KK: Total"));
            MV_KK(Att,xt,yt);
#endif
          }
          break;
        // Tpetra
        case Experiments::TPETRA:
          {
#ifdef HAVE_MUELU_TPETRA
            TimeMonitor t(*TimeMonitor::getNewTimer("MV Tpetra: Total"));
            MV_Tpetra(*At,xt,yt);
#endif
          }
          break;

        default:
          std::cerr << "Unknown experiment ID encountered: " << (int) experiment_id << std::endl;
        }
        comm->barrier();
      }// end random exp loop
    } // end repeat
    } // end ! my_experiments.empty()

    if (useStackedTimer) {
      stacked_timer->stop("MueLu_MatvecKernelDriver");
      Teuchos::StackedTimer::OutputOptions options;
      options.output_fraction = options.output_histogram = options.output_minmax = true;
      stacked_timer->report(out, comm, options);
    } else {
      TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
    }

#if defined(HAVE_MUELU_MKL) 
    mkl_sparse_destroy(AMKL);
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

