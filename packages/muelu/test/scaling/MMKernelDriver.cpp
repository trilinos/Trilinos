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
    AVCL.set(ArowptrVCL.data(), AcolindVCL.data(), Avals.data(),Au->getNodeNumRows(),Au->getNodeNumCols(),Au->getNodeNumEntries());
    BVCL.set(BrowptrVCL.data(), BcolindVCL.data(), Bvals.data(),Bu->getNodeNumRows(),Bu->getNodeNumCols(),Bu->getNodeNumEntries());

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
// Kokkos Kernels Testing
// =========================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Multiply_KokkosKernels(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,  const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,  Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &C, std::string algorithm_name, int team_work_size) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
  RCP<TimeMonitor> tm;


  std::string prefix = std::string("MM KokkosKernels ")+algorithm_name + std::string(": ");

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
    typedef typename KCRS::device_type device_t;

    RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A,false));
    RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B,false));
    RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C,false));

    const KCRS & Amat = Au->getLocalMatrix();
    const KCRS & Bmat = Bu->getLocalMatrix();
    KCRS Cmat = Cu->getLocalMatrix();

    c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
    lno_view_t Crowptr("Crowptr",A.getNodeNumRows()+1);
    c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
    lno_nnz_view_t Ccolind = Cmat.graph.entries;
    const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
    scalar_view_t Cvals = Cmat.values;

    // KokkosKernelsHandle
    typedef KokkosKernels::Experimental::KokkosKernelsHandle<
    typename lno_view_t::const_value_type,typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type,
      typename device_t::execution_space, typename device_t::memory_space,typename device_t::memory_space > KernelHandle;
    KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(algorithm_name);
    typename KernelHandle::nnz_lno_t AnumRows = Au->getNodeNumRows();
    typename KernelHandle::nnz_lno_t BnumRows = Bu->getNodeNumRows();
    typename KernelHandle::nnz_lno_t BnumCols = Bu->getNodeNumCols();

    // **********************************
    // Multiply
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("Multiply"))));
    KernelHandle kh;
    kh.create_spgemm_handle(alg_enum);
    kh.set_team_work_size(team_work_size);

    KokkosSparse::Experimental::spgemm_symbolic(&kh,AnumRows,BnumRows,BnumCols,Arowptr,Acolind,false,Browptr,Bcolind,false,Crowptr);

    size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size){
      Ccolind = lno_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      Cvals = scalar_view_t (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }
    KokkosSparse::Experimental::spgemm_numeric(&kh,AnumRows,BnumRows,BnumCols,Arowptr,Acolind,Avals,false,Browptr,Bcolind,Bvals,false,Crowptr,Ccolind,Cvals);
    kh.destroy_spgemm_handle();
    KCRS::execution_space::fence();

    tm = Teuchos::null;
    Au->getComm()->barrier();

    // **********************************
    // Copy out the data for KokkosKernels
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("Copy Out"))));
    Cmat.graph.row_map = Crowptr;
    Cmat.graph.entries = Ccolind;
    Cmat.values = Cvals;

  }
#endif
}


// =========================================================================
// Kokkos Kernels Testing
// =========================================================================
#ifdef HAVE_MUELU_TPETRA
#include "Tpetra_Import_Util2.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_decl.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_def.hpp"
#endif

//The LTG kernel is only defined for the Kokkos OpenMP node, so
//its test must only be enabled for OpenMP
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct LTG_Tests
{
  static void Multiply_LTG(
      const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&,
      const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&,
      Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>&)
  {}
};

#ifdef HAVE_TPETRA_INST_OPENMP

template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
struct LTG_Tests<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>
{
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;
  static void Multiply_LTG(const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A,  const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,  Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &C)
  {
#include <MueLu_UseShortNames.hpp>
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::TimeMonitor;

    Xpetra::UnderlyingLib lib = A.getRowMap()->lib();
    RCP<TimeMonitor> tm;

    if (lib == Xpetra::UseTpetra) {
      typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> crs_matrix_type;
      typedef Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>          import_type;
      typedef typename crs_matrix_type::local_matrix_type    KCRS;
      typedef typename KCRS::device_type device_t;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type lno_view_t;
      typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
      typedef typename KCRS::values_type::non_const_type scalar_view_t;
      typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;
      typedef Tpetra::Map<LO,GO,NO>                                     map_type;
      typedef typename map_type::local_map_type                         local_map_type;
      typedef typename Node::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
      LocalOrdinal LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
      RCP<const import_type> Cimport;
      RCP<const crs_matrix_type> Au = Utilities::Op2TpetraCrs(rcp(&A,false));
      RCP<const crs_matrix_type> Bu = Utilities::Op2TpetraCrs(rcp(&B,false));
      RCP<const crs_matrix_type> Cu = Utilities::Op2TpetraCrs(rcp(&C,false));
      RCP<crs_matrix_type> Cnc = Teuchos::rcp_const_cast<crs_matrix_type>(Cu);

      //    if(!Au->getComm()->getRank())
      //      std::cout<< "Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency() = "<<Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency()<<std::endl;



      // **********************************
      // Copy in the data for LTG
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM LTG: CopyIn")));

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
      lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"),Bview.colMap->getNodeNumElements()), Icol2Ccol;
      const LO colMapSize = static_cast<LO>(Bview.colMap->getNodeNumElements());
      Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::Bcol2Ccol_fill",
                           Kokkos::RangePolicy<execution_space, LO>(0, colMapSize),
                           [=](const LO i) {
                             Bcol2Ccol(i) = i;
                           });

      // Acol2Brow
      local_map_type Acolmap_local = Aview.colMap->getLocalMap();
      local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
      lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getNodeNumElements());
      lo_view_t targetMapToImportRow;
      Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::construct_tables",range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),[&](const LO i) {
        GO aidx = Acolmap_local.getGlobalElement(i);
        LO B_LID = Browmap_local.getLocalElement(aidx);
        if (B_LID != LO_INVALID) {
          targetMapToOrigRow(i)   = B_LID;
          //        targetMapToImportRow(i) = LO_INVALID;
        } else {
          // This shouldn't happen here
        }
      });
      tm = Teuchos::null;
      Au->getComm()->barrier();

      // **********************************
      // Multiply
      Teuchos::RCP<Teuchos::ParameterList> params;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MM LTG: Multiply")));
      Tpetra::MatrixMatrix::ExtraKernels::mult_A_B_newmatrix_LowThreadGustavsonKernel(Aview,Bview,targetMapToOrigRow,targetMapToImportRow,Bcol2Ccol,Icol2Ccol,*Cnc,Cimport,std::string("LTG Test"),params);

      tm = Teuchos::null;
      Au->getComm()->barrier();
    }
  }
};
#endif  //Tpetra and OpenMP

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
    //    GO nx = 50, ny = 50, nz = 50;
    // Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace3D"); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

    bool   printTimings = true;  clp.setOption("timings", "notimings",  &printTimings, "print timings to screen");
    int    nrepeat      = 100;   clp.setOption("nrepeat",               &nrepeat,      "repeat the experiment N times");
    // the kernels
    bool do_viennaCL = true;
    bool do_mkl      = true;
    bool do_kk_mem   = true;
    bool do_kk_dense = true;
    bool do_kk_default = true;
    bool do_ltg        = true;

    #ifndef HAVE_MUELU_VIENNACL
      do_viennaCL = false;
    #endif

    #ifndef HAVE_MUELU_MKL
      do_mkl = false;
    #endif
    clp.setOption("viennaCL",   "noviennaCL",   &do_viennaCL,   "Evaluate ViennaCL");
    clp.setOption("mkl",        "nomkl",        &do_mkl,        "Evaluate MKL SpMM");
    clp.setOption("kk_mem",     "nokk_mem",     &do_kk_mem,     "Evaluate KK Mem");
    clp.setOption("kk_dense",   "nokk_dense",   &do_kk_dense,   "Evaluate KK Dense");
    clp.setOption("kk_default", "nokk_default", &do_kk_default, "Evaluate KK Default");
    clp.setOption("ltg",        "noltg",        &do_ltg,        "Evaluate LTG");

    int kk_team_work_size=16;

    std::string matrixFileNameA = "A.mm"; clp.setOption("matrixfileA", &matrixFileNameA, "matrix market file containing matrix");
    std::string matrixFileNameB = "B.mm"; clp.setOption("matrixfileB", &matrixFileNameB, "matrix market file containing matrix");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    #ifndef HAVE_MUELU_VIENNACL
    if (do_viennaCL == true ) {
      out << "ViennaCL was requested, but this kernel is not available. Disabling..." << endl;
      do_viennaCL = false;
    }
    #endif

    #ifndef HAVE_MUELU_MKL
    if (do_mkl == true ) {
      out << "MKL was requested, but this kernel is not available. Disabling..." << endl;
      do_mkl = false;
    }
    #endif

    // simple hack to randomize order of experiments
    enum class Experiments { ViennaCL=0, MKL_SPMM, KK_MEM, KK_DENSE, KK_DEFAULT, LTG };
    std::vector<Experiments> my_experiments;
    // add the experiments we will run
    #ifdef HAVE_MUELU_VIENNACL
    if (do_viennaCL) my_experiments.push_back(Experiments::ViennaCL);   // ViennaCL
    #endif

    #ifdef HAVE_MUELU_MKL
    if (do_mkl) my_experiments.push_back(Experiments::MKL_SPMM);   // MKL SPMM
    #endif

    // assume these are available
    if (do_kk_mem)     my_experiments.push_back(Experiments::KK_MEM);     // KK Mem
    if (do_kk_dense)   my_experiments.push_back(Experiments::KK_DENSE);   // KK Dense
    if (do_kk_default) my_experiments.push_back(Experiments::KK_DEFAULT); // KK Default
    if (do_ltg)        my_experiments.push_back(Experiments::LTG);        // LTG


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

    // At the moment, this test only runs on one MPI rank
    if(comm->getSize() != 1) exit(1);

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("MatrixRead: S - Global Time")));

    comm->barrier();


    RCP<Matrix> A = Xpetra::IO<SC,LO,GO,Node>::Read(std::string(matrixFileNameA), lib, comm);
    RCP<Matrix> B = Xpetra::IO<SC,LO,GO,Node>::Read(std::string(matrixFileNameB), lib, comm);
    RCP<Matrix> C;

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
    std::mt19937 random_source (rd());

    // no need for a barrier, because the randomization process uses a collective.
    if ( my_experiments.empty() ) {

    for (int i=0; i<nrepeat; i++) {

      // randomize the experiments
      if (comm->getRank() == 0 ) {
        std::shuffle(my_experiments.begin(), my_experiments.end(), random_source);
      }
      // Broadcast this ordering to the other processes
      comm->broadcast(0,
                      static_cast<int>(sizeof(Experiments::LTG)*my_experiments.size()),
                      reinterpret_cast<char *>(my_experiments.data()));

      // loop over the randomized experiments
      for (const auto& experiment_id : my_experiments) {
        switch (experiment_id) {
        // ViennaCL
        case Experiments::ViennaCL:
          #ifdef HAVE_MUELU_VIENNACL
          C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
          {
            TimeMonitor t(*TimeMonitor::getNewTimer("MM ViennaCL: Total"));
            Multiply_ViennaCL(*A,*B,*C);
          }
          #endif
          break;
        // MKL_SPMM
        case Experiments::MKL_SPMM:
          #ifdef HAVE_MUELU_MKL
          C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
          {
            TimeMonitor t(*TimeMonitor::getNewTimer("MM MKL: Total"));
            Multiply_MKL_SPMM(*A,*B,*C);
          }
          #endif
          break;
        // KK Algorithms (KK Memory)
        case Experiments::KK_MEM:
          C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
          {
            TimeMonitor t(*TimeMonitor::getNewTimer("MM SPGEMM_KK_MEMORY: Total"));
            Multiply_KokkosKernels(*A,*B,*C,std::string("SPGEMM_KK_MEMORY"),kk_team_work_size);
          }
          break;
        // KK Algorithms (KK Dense)
        case Experiments::KK_DENSE:
          C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
          {
            TimeMonitor t(*TimeMonitor::getNewTimer("MM SPGEMM_KK_DENSE: Total"));
            Multiply_KokkosKernels(*A,*B,*C,std::string("SPGEMM_KK_DENSE"),kk_team_work_size);
          }
          break;
        // KK Algorithms (KK Default)
        case Experiments::KK_DEFAULT:
          C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
          {
            TimeMonitor t(*TimeMonitor::getNewTimer("MM SPGEMM_KK: Total"));
            Multiply_KokkosKernels(*A,*B,*C,std::string("SPGEMM_KK"),kk_team_work_size);
          }
          break;
        // LTG
        case Experiments::LTG:
          C = Xpetra::MatrixFactory<SC,LO,GO,Node>::Build(A->getRowMap(),0);
          {
            TimeMonitor t(*TimeMonitor::getNewTimer("MM LTG: Total"));
            LTG_Tests<SC,LO,GO,NO>::Multiply_LTG(*A,*B,*C);
          }
          break;
        default:
          std::cerr << "Unknown experiment ID encountered: " << (int) experiment_id << std::endl;
        }
        comm->barrier();
      }// end random exp loop
    } // end repeat
    } // end ! my_experiments.empty()

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

