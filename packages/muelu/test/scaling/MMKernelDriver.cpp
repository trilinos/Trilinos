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
    lno_view_t Crowptr("Crowptr",Cmat.graph.row_map.extent(0));// Because const
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

    Kokkos::parallel_for(Arowptr.extent(0),KOKKOS_LAMBDA(const size_t i) {
        ArowptrVCL(i) = Arowptr(i);
      });

    Kokkos::parallel_for(Browptr.extent(0),KOKKOS_LAMBDA(const size_t i) {
        BrowptrVCL(i) = Browptr(i);
      });

    Kokkos::parallel_for(Acolind.extent(0),KOKKOS_LAMBDA(const size_t i) {
        AcolindVCL(i) = Arowptr(i);
      });

    Kokkos::parallel_for(Bcolind.extent(0),KOKKOS_LAMBDA(const size_t i) {
        BcolindVCL(i) = Bcolind(i);
      });

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

    Kokkos::parallel_for(Crowptr.extent(0),KOKKOS_LAMBDA(const size_t i) {
        Crowptr(i) = CrowptrVCL[i];
      });


     Kokkos::parallel_for(Ccolind.extent(0),KOKKOS_LAMBDA(const size_t i) {
         Ccolind(i) = CcolindVCL[i];
         Cvals(i)  = CvalsVCL[i];
      });
     
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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Multiply_MKL(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,  const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &C) {

}

#endif



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
      
        // MKL
#ifdef HAVE_MUELU_MKL
        C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
        {
          TimeMonitor t(*TimeMonitor::getNewTimer("MM MKL: Total"));
          Multiply_MKL(*A,*B,*C);
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

