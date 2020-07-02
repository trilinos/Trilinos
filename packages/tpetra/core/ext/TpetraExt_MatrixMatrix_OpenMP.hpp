// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_OPENMP_DEF_HPP
#define TPETRA_MATRIXMATRIX_OPENMP_DEF_HPP

#ifdef HAVE_TPETRA_INST_OPENMP
namespace Tpetra {
namespace MMdetails {

/*********************************************************************************************************/
// MMM KernelWrappers for Partial Specialization to OpenMP
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType> {
    static inline void mult_A_B_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                  const LocalOrdinalViewType & Acol2Brow,
                                                  const LocalOrdinalViewType & Acol2Irow,
                                                  const LocalOrdinalViewType & Bcol2Ccol,
                                                  const LocalOrdinalViewType & Icol2Ccol,
                                                  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                  const std::string& label = std::string(),
                                                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    static inline void mult_A_B_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                  const LocalOrdinalViewType & Acol2Brow,
                                                  const LocalOrdinalViewType & Acol2Irow,
                                                  const LocalOrdinalViewType & Bcol2Ccol,
                                                  const LocalOrdinalViewType & Icol2Ccol,
                                                  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                  const std::string& label = std::string(),
                                                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);


};


// Jacobi KernelWrappers for Partial Specialization to OpenMP
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers2<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType> {
    static inline void jacobi_A_B_newmatrix_kernel_wrapper(Scalar omega,
                                                           const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                           const LocalOrdinalViewType & Acol2Brow,
                                                           const LocalOrdinalViewType & Acol2Irow,
                                                           const LocalOrdinalViewType & Bcol2Ccol,
                                                           const LocalOrdinalViewType & Icol2Ccol,
                                                           CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                           Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                           const std::string& label = std::string(),
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    static inline void jacobi_A_B_reuse_kernel_wrapper(Scalar omega,
                                                       const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                       const LocalOrdinalViewType & Acol2Brow,
                                                       const LocalOrdinalViewType & Acol2Irow,
                                                       const LocalOrdinalViewType & Bcol2Ccol,
                                                       const LocalOrdinalViewType & Icol2Ccol,
                                                       CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                       Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                       const std::string& label = std::string(),
                                                       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    static inline void jacobi_A_B_newmatrix_KokkosKernels(Scalar omega,
                                                          const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                          CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                          CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                          const LocalOrdinalViewType & Acol2Brow,
                                                          const LocalOrdinalViewType & Acol2Irow,
                                                          const LocalOrdinalViewType & Bcol2Ccol,
                                                          const LocalOrdinalViewType & Icol2Ccol,
                                                          CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                          Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                          const std::string& label = std::string(),
                                                          const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

};


// Triple-Product KernelWrappers for Partial Specialization to OpenMP
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal, class LocalOrdinalViewType>
struct KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType> {
  static inline void mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Rview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                                         const LocalOrdinalViewType & Acol2Prow,
                                                         const LocalOrdinalViewType & Acol2PIrow,
                                                         const LocalOrdinalViewType & Pcol2Ccol,
                                                         const LocalOrdinalViewType & PIcol2Ccol,
                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                                         Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
                                                         const std::string& label = std::string(),
                                                         const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

static inline void mult_R_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Rview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                                         const LocalOrdinalViewType & Acol2Prow,
                                                         const LocalOrdinalViewType & Acol2PIrow,
                                                         const LocalOrdinalViewType & Pcol2Ccol,
                                                         const LocalOrdinalViewType & PIcol2Ccol,
                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                                         Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
                                                         const std::string& label = std::string(),
                                                         const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static inline void mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                                         const LocalOrdinalViewType & Acol2Prow,
                                                         const LocalOrdinalViewType & Acol2PIrow,
                                                         const LocalOrdinalViewType & Pcol2Ccol,
                                                         const LocalOrdinalViewType & PIcol2Ccol,
                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                                         Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
                                                         const std::string& label = std::string(),
                                                         const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  static inline void mult_PT_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                                         const LocalOrdinalViewType & Acol2Prow,
                                                         const LocalOrdinalViewType & Acol2PIrow,
                                                         const LocalOrdinalViewType & Pcol2Ccol,
                                                         const LocalOrdinalViewType & PIcol2Ccol,
                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                                         Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
                                                         const std::string& label = std::string(),
                                                         const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
};


/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::mult_A_B_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                                                               const LocalOrdinalViewType & Acol2Brow,
                                                                                               const LocalOrdinalViewType & Acol2Irow,
                                                                                               const LocalOrdinalViewType & Bcol2Ccol,
                                                                                               const LocalOrdinalViewType & Icol2Ccol,
                                                                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                                                               const std::string& label,
                                                                                               const Teuchos::RCP<Teuchos::ParameterList>& params) {

#ifdef HAVE_TPETRA_MMM_TIMINGS
    std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
    using Teuchos::TimeMonitor;
    Teuchos::RCP<TimeMonitor> MM;
#endif

  // Node-specific code
  std::string nodename("OpenMP");

  // Lots and lots of typedefs
  using Teuchos::RCP;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode>::local_matrix_type KCRS;
  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  // Options
  int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer
  std::string myalg("SPGEMM_KK_MEMORY");


  if(!params.is_null()) {
    if(params->isParameter("openmp: algorithm"))
      myalg = params->get("openmp: algorithm",myalg);
    if(params->isParameter("openmp: team work size"))
      team_work_size = params->get("openmp: team work size",team_work_size);
  }

  if(myalg == "LTG") {
    // Use the LTG kernel if requested
    ::Tpetra::MatrixMatrix::ExtraKernels::mult_A_B_newmatrix_LowThreadGustavsonKernel(Aview,Bview,Acol2Brow,Acol2Irow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);
  }
  else {
    // Use the Kokkos-Kernels OpenMP Kernel
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix OpenMPWrapper"))));
#endif
    // KokkosKernelsHandle
    typedef KokkosKernels::Experimental::KokkosKernelsHandle<
       typename lno_view_t::const_value_type,typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type,
       typename device_t::execution_space, typename device_t::memory_space,typename device_t::memory_space > KernelHandle;

    // Grab the  Kokkos::SparseCrsMatrices
    const KCRS & Ak = Aview.origMatrix->getLocalMatrix();
    // const KCRS & Bk = Bview.origMatrix->getLocalMatrix();

    // Get the algorithm mode
    std::string alg = nodename+std::string(" algorithm");
    //  printf("DEBUG: Using kernel: %s\n",myalg.c_str());
    if(!params.is_null() && params->isParameter(alg)) myalg = params->get(alg,myalg);
    KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(myalg);

    // Merge the B and Bimport matrices
    const KCRS Bmerged = Tpetra::MMdetails::merge_matrices(Aview,Bview,Acol2Brow,Acol2Irow,Bcol2Ccol,Icol2Ccol,C.getColMap()->getNodeNumElements());

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix OpenMPCore"))));
#endif

    // Do the multiply on whatever we've got
    typename KernelHandle::nnz_lno_t AnumRows = Ak.numRows();
                                                          //    typename KernelHandle::nnz_lno_t BnumRows = Bmerged->numRows();
                                                          //    typename KernelHandle::nnz_lno_t BnumCols = Bmerged->numCols();
    typename KernelHandle::nnz_lno_t BnumRows = Bmerged.numRows();
    typename KernelHandle::nnz_lno_t BnumCols = Bmerged.numCols();


    lno_view_t      row_mapC (Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), AnumRows + 1);
    lno_nnz_view_t  entriesC;
    scalar_view_t   valuesC;
    KernelHandle kh;
    kh.create_spgemm_handle(alg_enum);
    kh.set_team_work_size(team_work_size);
                                                          //    KokkosSparse::Experimental::spgemm_symbolic(&kh,AnumRows,BnumRows,BnumCols,Ak.graph.row_map,Ak.graph.entries,false,Bmerged->graph.row_map,Bmerged->graph.entries,false,row_mapC);
    KokkosSparse::Experimental::spgemm_symbolic(&kh,AnumRows,BnumRows,BnumCols,Ak.graph.row_map,Ak.graph.entries,false,Bmerged.graph.row_map,Bmerged.graph.entries,false,row_mapC);

    size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size){
      entriesC = lno_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC = scalar_view_t (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }
                                                          //    KokkosSparse::Experimental::spgemm_numeric(&kh,AnumRows,BnumRows,BnumCols,Ak.graph.row_map,Ak.graph.entries,Ak.values,false,Bmerged->graph.row_map,Bmerged->graph.entries,Bmerged->values,false,row_mapC,entriesC,valuesC);
    KokkosSparse::Experimental::spgemm_numeric(&kh,AnumRows,BnumRows,BnumCols,Ak.graph.row_map,Ak.graph.entries,Ak.values,false,Bmerged.graph.row_map,Bmerged.graph.entries,Bmerged.values,false,row_mapC,entriesC,valuesC);
    kh.destroy_spgemm_handle();

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix OpenMPSort"))));
#endif
    // Sort & set values
    if (params.is_null() || params->get("sort entries",true))
      Import_Util::sortCrsEntries(row_mapC, entriesC, valuesC);
    C.setAllValues(row_mapC,entriesC,valuesC);

  }// end OMP KokkosKernels loop

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix OpenMPESFC"))));
#endif

  // Final Fillcomplete
  RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
  labelList->set("Timer Label",label);
  if(!params.is_null()) labelList->set("compute global constants",params->get("compute global constants",true));
  RCP<const Export<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > dummyExport;
  C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport,dummyExport,labelList);

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


/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::mult_A_B_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                                                               const LocalOrdinalViewType & Acol2Brow,
                                                                                               const LocalOrdinalViewType & Acol2Irow,
                                                                                               const LocalOrdinalViewType & Bcol2Ccol,
                                                                                               const LocalOrdinalViewType & Icol2Ccol,
                                                                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
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
  int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer
  std::string myalg("LTG");
  if(!params.is_null()) {
    if(params->isParameter("openmp: algorithm"))
      myalg = params->get("openmp: algorithm",myalg);
    if(params->isParameter("openmp: team work size"))
      team_work_size = params->get("openmp: team work size",team_work_size);
  }

  if(myalg == "LTG") {
    // Use the LTG kernel if requested
    ::Tpetra::MatrixMatrix::ExtraKernels::mult_A_B_reuse_LowThreadGustavsonKernel(Aview,Bview,Acol2Brow,Acol2Irow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);
  }
  else {
    throw std::runtime_error("Tpetra::MatrixMatrix::MMM reuse unknown kernel");
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Reuse OpenMPESFC"))));
#endif
  C.fillComplete(C.getDomainMap(), C.getRangeMap());
}


/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers2<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::jacobi_A_B_newmatrix_kernel_wrapper(Scalar omega,
                                                                                               const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                                                               const LocalOrdinalViewType & Acol2Brow,
                                                                                               const LocalOrdinalViewType & Acol2Irow,
                                                                                               const LocalOrdinalViewType & Bcol2Ccol,
                                                                                               const LocalOrdinalViewType & Icol2Ccol,
                                                                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                                                               const std::string& label,
                                                                                               const Teuchos::RCP<Teuchos::ParameterList>& params) {

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
  if(!params.is_null()) {
    if(params->isParameter("openmp: jacobi algorithm"))
      myalg = params->get("openmp: jacobi algorithm",myalg);
    if(params->isParameter("openmp: team work size"))
      team_work_size = params->get("openmp: team work size",team_work_size);
  }

  if(myalg == "LTG") {
    // Use the LTG kernel if requested
    ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_newmatrix_LowThreadGustavsonKernel(omega,Dinv,Aview,Bview,Acol2Brow,Acol2Irow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);
  }
  else if(myalg == "MSAK") {
    ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_newmatrix_MultiplyScaleAddKernel(omega,Dinv,Aview,Bview,Acol2Brow,Acol2Irow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);
  }
  else if(myalg == "KK") {
    jacobi_A_B_newmatrix_KokkosKernels(omega,Dinv,Aview,Bview,Acol2Brow,Acol2Irow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);
  }
  else {
    throw std::runtime_error("Tpetra::MatrixMatrix::Jacobi newmatrix unknown kernel");
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix OpenMPESFC"))));
#endif

  // Final Fillcomplete
  RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
  labelList->set("Timer Label",label);
  if(!params.is_null()) labelList->set("compute global constants",params->get("compute global constants",true));

  // NOTE: MSAK already fillCompletes, so we have to check here
  if(!C.isFillComplete()) {
    RCP<const Export<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > dummyExport;
    C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport,dummyExport,labelList);
  }

}



/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers2<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::jacobi_A_B_reuse_kernel_wrapper(Scalar omega,
                                                                                               const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                                                               const LocalOrdinalViewType & Acol2Brow,
                                                                                               const LocalOrdinalViewType & Acol2Irow,
                                                                                               const LocalOrdinalViewType & Bcol2Ccol,
                                                                                               const LocalOrdinalViewType & Icol2Ccol,
                                                                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
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
  int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer
  std::string myalg("LTG");
  if(!params.is_null()) {
    if(params->isParameter("openmp: jacobi algorithm"))
      myalg = params->get("openmp: jacobi algorithm",myalg);
    if(params->isParameter("openmp: team work size"))
      team_work_size = params->get("openmp: team work size",team_work_size);
  }

  if(myalg == "LTG") {
    // Use the LTG kernel if requested
    ::Tpetra::MatrixMatrix::ExtraKernels::jacobi_A_B_reuse_LowThreadGustavsonKernel(omega,Dinv,Aview,Bview,Acol2Brow,Acol2Irow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);
  }
  else {
    throw std::runtime_error("Tpetra::MatrixMatrix::Jacobi reuse unknown kernel");
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Reuse OpenMPESFC"))));
#endif
  C.fillComplete(C.getDomainMap(), C.getRangeMap());

}

/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers2<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::jacobi_A_B_newmatrix_KokkosKernels(Scalar omega,
												const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
											        CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
												CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
												const LocalOrdinalViewType & Acol2Brow,
												const LocalOrdinalViewType & Acol2Irow,
												const LocalOrdinalViewType & Bcol2Ccol,
												const LocalOrdinalViewType & Icol2Ccol,
												CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
												Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
												const std::string& label,
												const Teuchos::RCP<Teuchos::ParameterList>& params) {

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM;
#endif

  // Check if the diagonal entries exist in debug mode
  const bool debug = Tpetra::Details::Behavior::debug();
  if(debug) {

    auto rowMap = Aview.origMatrix->getRowMap();
    Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode> diags(rowMap);
    Aview.origMatrix->getLocalDiagCopy(diags);
    size_t diagLength = rowMap->getNodeNumElements();
    Teuchos::Array<Scalar> diagonal(diagLength);
    diags.get1dCopy(diagonal());

    for(size_t i = 0; i < diagLength; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(diagonal[i] == Teuchos::ScalarTraits<Scalar>::zero(), 
				 std::runtime_error, 
				 "Matrix A has a zero/missing diagonal: " << diagonal[i] << std::endl <<
				 "KokkosKernels Jacobi-fused SpGEMM requires nonzero diagonal entries in A" << std::endl);
    }
  }

  // Usings
  using device_t = typename Kokkos::Compat::KokkosOpenMPWrapperNode::device_type;
  using matrix_t = typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode>::local_matrix_type;
  using graph_t = typename matrix_t::StaticCrsGraphType;
  using lno_view_t = typename graph_t::row_map_type::non_const_type;
  using c_lno_view_t = typename graph_t::row_map_type::const_type;
  using lno_nnz_view_t = typename graph_t::entries_type::non_const_type;
  using scalar_view_t = typename matrix_t::values_type::non_const_type;

  // KokkosKernels handle 
  using handle_t = typename KokkosKernels::Experimental::KokkosKernelsHandle<
    typename lno_view_t::const_value_type,typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type,
    typename device_t::execution_space, typename device_t::memory_space,typename device_t::memory_space >;

  // Get the rowPtr, colInd and vals of importMatrix
  c_lno_view_t Irowptr;
  lno_nnz_view_t Icolind;
  scalar_view_t Ivals;
  if(!Bview.importMatrix.is_null()) {
    Irowptr = Bview.importMatrix->getLocalMatrix().graph.row_map;
    Icolind = Bview.importMatrix->getLocalMatrix().graph.entries;
    Ivals   = Bview.importMatrix->getLocalMatrix().values;
  }

  // Merge the B and Bimport matrices
  const matrix_t Bmerged = Tpetra::MMdetails::merge_matrices(Aview,Bview,Acol2Brow,Acol2Irow,Bcol2Ccol,Icol2Ccol,C.getColMap()->getNodeNumElements());

  // Get the properties and arrays of input matrices
  const matrix_t & Amat = Aview.origMatrix->getLocalMatrix();
  const matrix_t & Bmat = Bview.origMatrix->getLocalMatrix();

  typename handle_t::nnz_lno_t AnumRows = Amat.numRows();
  typename handle_t::nnz_lno_t BnumRows = Bmerged.numRows();
  typename handle_t::nnz_lno_t BnumCols = Bmerged.numCols();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmerged.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmerged.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmerged.values;

  // Arrays of the output matrix
  lno_view_t row_mapC (Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), AnumRows + 1);
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  // Options
  int team_work_size = 16; 
  std::string myalg("SPGEMM_KK_MEMORY");
  if(!params.is_null()) {
    if(params->isParameter("cuda: algorithm"))
      myalg = params->get("cuda: algorithm",myalg);
    if(params->isParameter("cuda: team work size"))
      team_work_size = params->get("cuda: team work size",team_work_size);
  }

  // Get the algorithm mode
  std::string nodename("OpenMP");
  std::string alg = nodename + std::string(" algorithm");
  if(!params.is_null() && params->isParameter(alg)) myalg = params->get(alg,myalg);
  KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(myalg);


  // KokkosKernels call
  handle_t kh;
  kh.create_spgemm_handle(alg_enum);
  kh.set_team_work_size(team_work_size);

  KokkosSparse::Experimental::spgemm_symbolic(&kh, AnumRows, BnumRows, BnumCols,
					      Arowptr, Acolind, false,
					      Browptr, Bcolind, false,
					      row_mapC);

  size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
  if (c_nnz_size){
    entriesC = lno_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
    valuesC = scalar_view_t (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
  }

  KokkosSparse::Experimental::spgemm_jacobi(&kh, AnumRows, BnumRows, BnumCols,
					    Arowptr, Acolind, Avals, false,
					    Browptr, Bcolind, Bvals, false,
					    row_mapC, entriesC, valuesC,
					    omega, Dinv.getLocalViewDevice());
  kh.destroy_spgemm_handle();

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix OpenMPSort"))));
#endif

  // Sort & set values
  if (params.is_null() || params->get("sort entries",true))
    Import_Util::sortCrsEntries(row_mapC, entriesC, valuesC);
  C.setAllValues(row_mapC,entriesC,valuesC);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix OpenMPESFC"))));
#endif

  // Final Fillcomplete
  Teuchos::RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
  labelList->set("Timer Label",label);
  if(!params.is_null()) labelList->set("compute global constants",params->get("compute global constants",true));
  Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > dummyExport;
  C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport,dummyExport,labelList);
}


/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Rview,
                                                                                                                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                                                                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                                                                                                                                         const LocalOrdinalViewType & Acol2Prow,
                                                                                                                                                         const LocalOrdinalViewType & Acol2PIrow,
                                                                                                                                                         const LocalOrdinalViewType & Pcol2Accol,
                                                                                                                                                         const LocalOrdinalViewType & PIcol2Accol,
                                                                                                                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                                                                                                                                         Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
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

  if(!params.is_null()) {
    if(params->isParameter("openmp: rap algorithm"))
      myalg = params->get("openmp: rap algorithm",myalg);
  }

  if(myalg == "LTG") {
    // Use the LTG kernel if requested
    ::Tpetra::MatrixMatrix::ExtraKernels::mult_R_A_P_newmatrix_LowThreadGustavsonKernel(Rview,Aview,Pview,Acol2Prow,Acol2PIrow,Pcol2Accol,PIcol2Accol,Ac,Acimport,label,params);
  }
  else {
    throw std::runtime_error("Tpetra::MatrixMatrix::R_A_P newmatrix unknown kernel");
  }
}

/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::mult_R_A_P_reuse_kernel_wrapper(
                                                                                                                                                      CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Rview,
                                                                                                                                                      CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                                                                                                                      CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,

                                                                                               const LocalOrdinalViewType & Acol2Prow,
                                                                                               const LocalOrdinalViewType & Acol2Irow,
                                                                                               const LocalOrdinalViewType & Pcol2Ccol,
                                                                                               const LocalOrdinalViewType & Icol2Ccol,
                                                                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
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
  if(!params.is_null()) {
    if(params->isParameter("openmp: rap algorithm"))
      myalg = params->get("openmp: rap algorithm",myalg);
  }

   if(myalg == "LTG") {
    // Use the LTG kernel if requested
    ::Tpetra::MatrixMatrix::ExtraKernels::mult_R_A_P_reuse_LowThreadGustavsonKernel(Rview,Aview,Pview,Acol2Prow,Acol2Irow,Pcol2Ccol,Icol2Ccol,C,Cimport,label,params);
  }
  else {
    throw std::runtime_error("Tpetra::MatrixMatrix::R_A_P newmatrix unknown kernel");
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Reuse OpenMPESFC"))));
#endif
  C.fillComplete(C.getDomainMap(), C.getRangeMap());

}



/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,

                                                                                                                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                                                                                                                                         const LocalOrdinalViewType & Acol2Prow,
                                                                                                                                                         const LocalOrdinalViewType & Acol2PIrow,
                                                                                                                                                         const LocalOrdinalViewType & Pcol2Accol,
                                                                                                                                                         const LocalOrdinalViewType & PIcol2Accol,
                                                                                                                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                                                                                                                                         Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
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

  if(!params.is_null()) {
    if(params->isParameter("openmp: ptap algorithm"))
      myalg = params->get("openmp: ptap algorithm",myalg);
  }

  if(myalg == "LTG") {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose"))));
#endif

    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using SC = Scalar;

    // We don't need a kernel-level PTAP, we just transpose here
    using transposer_type =
      RowMatrixTransposer<SC, LO, GO, Kokkos::Compat::KokkosOpenMPWrapperNode>;
    transposer_type transposer (Pview.origMatrix, label + "XP: ");
    RCP<ParameterList> transposeParams (new ParameterList);
    if (! params.is_null ()) {
      transposeParams->set ("compute global constants",
                            params->get ("compute global constants: temporaries",
                                         false));
    }
    transposeParams->set ("sort", false);
    RCP<CrsMatrix<SC, LO, GO, Kokkos::Compat::KokkosOpenMPWrapperNode> > Ptrans =
      transposer.createTransposeLocal (transposeParams);
    CrsMatrixStruct<SC, LO, GO, Kokkos::Compat::KokkosOpenMPWrapperNode> Rview;
    Rview.origMatrix = Ptrans;

    using ::Tpetra::MatrixMatrix::ExtraKernels::mult_R_A_P_newmatrix_LowThreadGustavsonKernel;
    mult_R_A_P_newmatrix_LowThreadGustavsonKernel
      (Rview, Aview, Pview, Acol2Prow, Acol2PIrow, Pcol2Accol,
       PIcol2Accol, Ac, Acimport, label, params);
  }
  else {
    throw std::runtime_error("Tpetra::MatrixMatrix::PT_A_P newmatrix unknown kernel");
  }
}

/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode,LocalOrdinalViewType>::mult_PT_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,

                                                                                                                                                         CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                                                                                                                                         const LocalOrdinalViewType & Acol2Prow,
                                                                                                                                                         const LocalOrdinalViewType & Acol2PIrow,
                                                                                                                                                         const LocalOrdinalViewType & Pcol2Accol,
                                                                                                                                                         const LocalOrdinalViewType & PIcol2Accol,
                                                                                                                                                         CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                                                                                                                                         Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
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

  if(!params.is_null()) {
    if(params->isParameter("openmp: ptap algorithm"))
      myalg = params->get("openmp: ptap algorithm",myalg);
  }

  if(myalg == "LTG") {
#ifdef HAVE_TPETRA_MMM_TIMINGS
   MM = Teuchos::null;  MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose"))));
#endif

    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;
    using SC = Scalar;

    // We don't need a kernel-level PTAP, we just transpose here
    using transposer_type =
      RowMatrixTransposer<SC, LO, GO, Kokkos::Compat::KokkosOpenMPWrapperNode>;
    transposer_type transposer (Pview.origMatrix, label + "XP: ");
    RCP<ParameterList> transposeParams (new ParameterList);
    if (! params.is_null ()) {
      transposeParams->set ("compute global constants",
                            params->get ("compute global constants: temporaries",
                                         false));
    }
    transposeParams->set ("sort", false);
    RCP<CrsMatrix<SC, LO, GO, Kokkos::Compat::KokkosOpenMPWrapperNode> > Ptrans =
      transposer.createTransposeLocal (transposeParams);
    CrsMatrixStruct<SC, LO, GO, Kokkos::Compat::KokkosOpenMPWrapperNode> Rview;
    Rview.origMatrix = Ptrans;

    using ::Tpetra::MatrixMatrix::ExtraKernels::mult_R_A_P_reuse_LowThreadGustavsonKernel;
    mult_R_A_P_reuse_LowThreadGustavsonKernel
      (Rview, Aview, Pview, Acol2Prow, Acol2PIrow, Pcol2Accol,
       PIcol2Accol, Ac, Acimport, label, params);
  }
  else {
    throw std::runtime_error("Tpetra::MatrixMatrix::PT_A_P reuse unknown kernel");
  }
  Ac.fillComplete(Ac.getDomainMap(), Ac.getRangeMap());
}


}//MMdetails
}//Tpetra

#endif//OpenMP

#endif
