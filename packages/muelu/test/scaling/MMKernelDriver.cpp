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



// =========================================================================
// ViennaCL Testing
// =========================================================================
#ifdef HAVE_MUELU_VIENNACL
#define VIENNACL_WITH_OPENMP
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/host_based/common.hpp"

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Multiply_ViennaCL(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,  const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,  Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &C) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
  RCP<TimeMonitor> tm;

  // NOTE: ViennaCL matrices use "unsigned int" for rowptr and colind and are templated on scalar type (yay); which means (for us) POD only.

  if (lib == Xpetra::UseTpetra) {
#if defined(HAVE_MUELU_TPETRA)    
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
    typedef typename crs_matrix_type::local_matrix_type    KCRS;
    typedef typename KCRS::StaticCrsGraphType              graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::row_map_type::const_type     c_lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename graph_t::entries_type::const_type     c_lno_nnz_view_t;
    typedef typename KCRS::values_type::non_const_type     scalar_view_t;
    typedef typename Kokkos::View<unsigned int*,typename lno_nnz_view_t::array_layout,typename lno_nnz_view_t::device_type> vcl_size_t_type;

    RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A,false));
    RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B,false));
    RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C,false));

    const KCRS & Amat = Au->getLocalMatrix();
    const KCRS & Bmat = Bu->getLocalMatrix();
    KCRS Cmat = Cu->getLocalMatrix();
    
    c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
    lno_view_t Crowptr("Crowptr",C.getNodeNumRows()+1);
    c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
    lno_nnz_view_t Ccolind = Cmat.graph.entries;
    const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
    scalar_view_t Cvals = Cmat.values;

    // **********************************
    // Copy in the data for ViennaCL    
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM ViennaCL: CopyIn")));

    vcl_size_t_type ArowptrVCL("Arowptr",Arowptr.extent(0));
    vcl_size_t_type BrowptrVCL("Browptr",Browptr.extent(0));

    vcl_size_t_type AcolindVCL("Acolind",Acolind.extent(0));
    vcl_size_t_type BcolindVCL("Bcolind",Bcolind.extent(0));

    copy_view(Arowptr,ArowptrVCL);
    copy_view(Browptr,BrowptrVCL);
    copy_view(Acolind,AcolindVCL);
    copy_view(Bcolind,BcolindVCL);

    viennacl::compressed_matrix<Scalar> AVCL, BVCL; 
    AVCL.set(ArowptrVCL.ptr_on_device(), AcolindVCL.ptr_on_device(), Avals.ptr_on_device(),Au->getNodeNumRows(),Au->getNodeNumCols(),Au->getNodeNumEntries());
    BVCL.set(BrowptrVCL.ptr_on_device(), BcolindVCL.ptr_on_device(), Bvals.ptr_on_device(),Bu->getNodeNumRows(),Bu->getNodeNumCols(),Bu->getNodeNumEntries());

    tm = Teuchos::null;
    Au->getComm()->barrier();

    // **********************************
    // Multiply
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM ViennaCL: Multiply")));
    viennacl::compressed_matrix<Scalar> CVCL = viennacl::linalg::prod(AVCL, BVCL);
    KCRS::execution_space::fence();
      

    // **********************************
    // Copy out the data for ViennaCL    
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM ViennaCL: Copy Out")));
    size_t cnnz = (size_t)CVCL.nnz();
    Kokkos::resize(Ccolind,cnnz);
    Kokkos::resize(Cvals,cnnz);
#ifdef VIENNACL_WITH_CUDA
    const unsigned int * CrowptrVCL = viennacl::cuda_arg<unsigned int>(CVCL.handle1());
    const unsigned int * CcolindVCL = viennacl::cuda_arg<unsigned int>(CVCL.handle2());
    const Scalar * CvalsVCL       = viennacl::cuda_arg<Scalar>(CVCL.handle());
#else 
    const unsigned int * CrowptrVCL = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(CVCL.handle1());
    const unsigned int * CcolindVCL = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(CVCL.handle2());
    const Scalar * CvalsVCL   = viennacl::linalg::host_based::detail::extract_raw_pointer<Scalar>(CVCL.handle());
#endif

    copy_view_n(Crowptr.extent(0),CrowptrVCL,Crowptr);
    copy_view_n(cnnz,CcolindVCL,Ccolind);
    copy_view_n(cnnz,CvalsVCL,Cvals);
     
    Cmat.graph.row_map = Crowptr;
    Cmat.graph.entries = Ccolind;
    Cmat.values = Cvals;

#endif

  tm = Teuchos::null;
  Au->getComm()->barrier();
  }

}

#endif


// =========================================================================
// MKL Testing
// =========================================================================
#ifdef HAVE_MUELU_MKL
#include "mkl.h"

  // mkl_sparse_spmm
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Multiply_MKL_SPMM(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,  const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &C) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  // NOTE: MKL uses int as its columnn index type and the either double or float for its Scalar type

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
  RCP<TimeMonitor> tm;
#if defined(HAVE_MUELU_TPETRA)    
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
    typedef typename crs_matrix_type::local_matrix_type    KCRS;
    typedef typename KCRS::StaticCrsGraphType              graph_t;
    typedef typename graph_t::row_map_type::non_const_type lno_view_t;
    typedef typename graph_t::row_map_type::const_type     c_lno_view_t;
    typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
    typedef typename graph_t::entries_type::const_type     c_lno_nnz_view_t;
    typedef typename KCRS::values_type::non_const_type     scalar_view_t;

    typedef typename Kokkos::View<MKL_INT*,typename lno_nnz_view_t::array_layout,typename lno_nnz_view_t::device_type> mkl_int_type;

    RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A,false));
    RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B,false));
    RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C,false));

    const KCRS & Amat = Au->getLocalMatrix();
    const KCRS & Bmat = Bu->getLocalMatrix();
    KCRS Cmat = Cu->getLocalMatrix();
    if(A.getNodeNumRows()!=C.getNodeNumRows())  throw std::runtime_error("C is not sized correctly");

    c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
    lno_view_t Crowptr("Crowptr",C.getNodeNumRows()+1);
    c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
    lno_nnz_view_t Ccolind = Cmat.graph.entries;
    const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
    scalar_view_t Cvals = Cmat.values;

    sparse_matrix_t AMKL;
    sparse_matrix_t BMKL;
    sparse_matrix_t CMKL;

    // **********************************
    // Copy in the data for MKL
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM MKL: CopyIn")));

    mkl_int_type ArowptrMKL("Arowptr",Arowptr.extent(0));
    mkl_int_type BrowptrMKL("Browptr",Browptr.extent(0));

    mkl_int_type AcolindMKL("Acolind",Acolind.extent(0));
    mkl_int_type BcolindMKL("Bcolind",Bcolind.extent(0));

    copy_view(Arowptr,ArowptrMKL);
    copy_view(Browptr,BrowptrMKL);
    copy_view(Acolind,AcolindMKL);
    copy_view(Bcolind,BcolindMKL);

    if(Kokkos::Impl::is_same<Scalar,double>::value) {
      mkl_sparse_d_create_csr(&AMKL, SPARSE_INDEX_BASE_ZERO, Au->getNodeNumRows(), Au->getNodeNumCols(), ArowptrMKL.data(),ArowptrMKL.data()+1,AcolindMKL.data(),(double*)Avals.data());
      mkl_sparse_d_create_csr(&BMKL, SPARSE_INDEX_BASE_ZERO, Bu->getNodeNumRows(), Bu->getNodeNumCols(), BrowptrMKL.data(),BrowptrMKL.data()+1,BcolindMKL.data(),(double*)Bvals.data());
    }
    else
      throw std::runtime_error("MKL Type Mismatch");
    

    tm = Teuchos::null;
    Au->getComm()->barrier();

    // **********************************
    // Multiply
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM MKL: Multiply")));
    mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, AMKL, BMKL, &CMKL);
    KCRS::execution_space::fence();
      
    // **********************************
    // Copy out the data for MKL
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM MKL: Copy Out")));
    sparse_index_base_t c_indexing;
    MKL_INT c_rows, c_cols, *rows_start, *rows_end, *columns;
    double * values;
    mkl_sparse_d_export_csr(CMKL,&c_indexing, &c_rows, &c_cols, &rows_start, &rows_end, &columns, &values);                          
    size_t cnnz = rows_end[c_rows-1];
    Kokkos::resize(Ccolind,cnnz);
    Kokkos::resize(Cvals,cnnz);
    if(c_rows != A.getNodeNumRows() || c_rows+1 != Crowptr.extent(0)) throw std::runtime_error("C row size mismatch");
    copy_view_n(c_rows,rows_start,Crowptr); Crowptr(c_rows) = rows_end[c_rows-1];
    copy_view_n(cnnz,columns,Ccolind);
    copy_view_n(cnnz,values,Cvals);
     
    Cmat.graph.row_map = Crowptr;
    Cmat.graph.entries = Ccolind;
    Cmat.values = Cvals;

    mkl_sparse_destroy(AMKL);
    mkl_sparse_destroy(BMKL);
    mkl_sparse_destroy(CMKL);
#endif



  tm = Teuchos::null;
  Au->getComm()->barrier();

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

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> pOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *pOut;
    out.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    GO nx = 50, ny = 50, nz = 50;
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace3D"); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

    bool   printTimings = true;  clp.setOption("timings", "notimings",  &printTimings, "print timings to screen");
    int    nrepeat      = 100;   clp.setOption("nrepeat",               &nrepeat,      "repeat the experiment N times");

    std::string matrixFileNameA = "A.mm"; clp.setOption("matrixfileA", &matrixFileNameA, "matrix market file containing matrix");
    std::string matrixFileNameB = "B.mm"; clp.setOption("matrixfileB", &matrixFileNameB, "matrix market file containing matrix");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    out << "========================================================\n" << xpetraParameters << matrixParameters;

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time"))), tm;

    comm->barrier();


    RCP<Matrix> A = Xpetra::IO<SC,LO,GO,Node>::Read(std::string(matrixFileNameA), lib, comm);
    RCP<Matrix> B = Xpetra::IO<SC,LO,GO,Node>::Read(std::string(matrixFileNameB), lib, comm);
    RCP<Matrix> C;


    comm->barrier();
    tm = Teuchos::null;

    out << "Matrix Read complete.\n========================================================" << std::endl;
    comm->barrier();

    for (int i=0; i<nrepeat; i++) {

      // ViennaCL
#ifdef HAVE_MUELU_VIENNACL
        C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
        {
          TimeMonitor t(*TimeMonitor::getNewTimer("MM ViennaCL: Total"));
          Multiply_ViennaCL(*A,*B,*C);
        }
        comm->barrier();
#endif
      
        // MKL_SPMM
#ifdef HAVE_MUELU_MKL
        C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
        {
          TimeMonitor t(*TimeMonitor::getNewTimer("MM MKL: Total"));
          Multiply_MKL_SPMM(*A,*B,*C);
        }
        comm->barrier();
#endif

        // LTG

        
        // KK Algorithms

    }// endnrepeat


    tm = Teuchos::null;
    globalTimeMonitor = Teuchos::null;

    if (printTimings) {
      TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
    }

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

