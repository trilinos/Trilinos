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

#ifndef TPETRA_MATRIXMATRIX_EXTRAKERNELS_DEF_HPP
#define TPETRA_MATRIXMATRIX_EXTRAKERNELS_DEF_HPP
#include "TpetraExt_MatrixMatrix_ExtraKernels_decl.hpp"

namespace Tpetra {

namespace MatrixMatrix{

namespace ExtraKernels{


template<class CrsMatrixType>
size_t C_estimate_nnz_per_row(CrsMatrixType & A, CrsMatrixType &B){
  // Follows the NZ estimate in ML's ml_matmatmult.c
  size_t Aest = 100, Best=100;
  if (A.getNodeNumEntries() > 0)
    Aest = (A.getNodeNumRows() > 0)?  A.getNodeNumEntries()/A.getNodeNumRows() : 100;
  if (B.getNodeNumEntries() > 0)
    Best = (B.getNodeNumRows() > 0) ? B.getNodeNumEntries()/B.getNodeNumRows() : 100;

  size_t nnzperrow = (size_t)(sqrt((double)Aest) + sqrt((double)Best) - 1);
  nnzperrow *= nnzperrow;

  return nnzperrow;
}


#if defined (HAVE_TPETRA_INST_OPENMP)
/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void mult_A_B_newmatrix_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                 const LocalOrdinalViewType & targetMapToOrigRow,
                                                 const LocalOrdinalViewType & targetMapToImportRow,
                                                 const LocalOrdinalViewType & Bcol2Ccol,
                                                 const LocalOrdinalViewType & Icol2Ccol,
                                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                 Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                 const std::string& label,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix LTGCore"))));
#endif

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;


  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode>::local_matrix_type KCRS;
  //  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  // Unmanaged versions of the above
  typedef UnmanagedView<lno_view_t> u_lno_view_t;
  typedef UnmanagedView<lno_nnz_view_t> u_lno_nnz_view_t;
  typedef UnmanagedView<scalar_view_t> u_scalar_view_t;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode NO;
  typedef Map<LO,GO,NO>                     map_type;

  // NOTE (mfh 15 Sep 2017) This is specifically only for
  // execution_space = Kokkos::OpenMP, so we neither need nor want
  // KOKKOS_LAMBDA (with its mandatory __device__ marking).
  typedef NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;

  // All of the invalid guys
  const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS & Amat = Aview.origMatrix->getLocalMatrix();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrix();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  size_t b_max_nnz_per_row = Bview.origMatrix->getNodeMaxNumRowEntries();

  c_lno_view_t  Irowptr;
  lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    Irowptr = Bview.importMatrix->getLocalMatrix().graph.row_map;
    Icolind = Bview.importMatrix->getLocalMatrix().graph.entries;
    Ivals   = Bview.importMatrix->getLocalMatrix().values;
    b_max_nnz_per_row = std::max(b_max_nnz_per_row,Bview.importMatrix->getNodeMaxNumRowEntries());
  }

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getNodeNumRows();
  size_t n = Ccolmap->getNodeNumElements();
  size_t Cest_nnz_per_row = 2*C_estimate_nnz_per_row(*Aview.origMatrix,*Bview.origMatrix);

  // Get my node / thread info (right from openmp or parameter list)
  size_t thread_max =  Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency();
  if(!params.is_null()) {
    if(params->isParameter("openmp: ltg thread max"))
      thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));
  }

  // Thread-local memory
  Kokkos::View<u_lno_view_t*> tl_rowptr("top_rowptr",thread_max);
  Kokkos::View<u_lno_nnz_view_t*> tl_colind("top_colind",thread_max);
  Kokkos::View<u_scalar_view_t*> tl_values("top_values",thread_max);

  double thread_chunk = (double)(m) / thread_max;

  // Run chunks of the matrix independently
  Kokkos::parallel_for("MMM::LTG::NewMatrix::ThreadLocal",range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid)
    {
      // Thread coordination stuff
      size_t my_thread_start =  tid * thread_chunk;
      size_t my_thread_stop  = tid == thread_max-1 ? m : (tid+1)*thread_chunk;
      size_t my_thread_m     = my_thread_stop - my_thread_start;

      // Size estimate
      size_t CSR_alloc = (size_t) (my_thread_m*Cest_nnz_per_row*0.75 + 100);

      // Allocations
      std::vector<size_t> c_status(n,INVALID);

      u_lno_view_t Crowptr((typename u_lno_view_t::data_type)malloc(u_lno_view_t::shmem_size(my_thread_m+1)),my_thread_m+1);
      u_lno_nnz_view_t Ccolind((typename u_lno_nnz_view_t::data_type)malloc(u_lno_nnz_view_t::shmem_size(CSR_alloc)),CSR_alloc);
      u_scalar_view_t Cvals((typename u_scalar_view_t::data_type)malloc(u_scalar_view_t::shmem_size(CSR_alloc)),CSR_alloc);

      // For each row of A/C
      size_t CSR_ip = 0, OLD_ip = 0;
      for (size_t i = my_thread_start; i < my_thread_stop; i++) {
        // mfh 27 Sep 2016: m is the number of rows in the input matrix A
        // on the calling process.
        Crowptr(i-my_thread_start) = CSR_ip;

        // mfh 27 Sep 2016: For each entry of A in the current row of A
        for (size_t k = Arowptr(i); k < Arowptr(i+1); k++) {
          LO Aik  = Acolind(k); // local column index of current entry of A
          const SC Aval = Avals(k);   // value of current entry of A
          if (Aval == SC_ZERO)
            continue; // skip explicitly stored zero values in A

          if (targetMapToOrigRow(Aik) != LO_INVALID) {
            // mfh 27 Sep 2016: If the entry of targetMapToOrigRow
            // corresponding to the current entry of A is populated, then
            // the corresponding row of B is in B_local (i.e., it lives on
            // the calling process).

            // Local matrix
            size_t Bk = Teuchos::as<size_t>(targetMapToOrigRow(Aik));

            // mfh 27 Sep 2016: Go through all entries in that row of B_local.
            for (size_t j = Browptr(Bk); j < Browptr(Bk+1); ++j) {
              LO Bkj = Bcolind(j);
              LO Cij = Bcol2Ccol(Bkj);

              if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip) {
                // New entry
                c_status[Cij]   = CSR_ip;
                Ccolind(CSR_ip) = Cij;
                Cvals(CSR_ip)   = Aval*Bvals(j);
                CSR_ip++;

              } else {
                Cvals(c_status[Cij]) += Aval*Bvals(j);
              }
            }

          } else {
            // mfh 27 Sep 2016: If the entry of targetMapToOrigRow
            // corresponding to the current entry of A NOT populated (has
            // a flag "invalid" value), then the corresponding row of B is
            // in B_local (i.e., it lives on the calling process).

            // Remote matrix
            size_t Ik = Teuchos::as<size_t>(targetMapToImportRow(Aik));
            for (size_t j = Irowptr(Ik); j < Irowptr(Ik+1); ++j) {
              LO Ikj = Icolind(j);
              LO Cij = Icol2Ccol(Ikj);

              if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip){
                // New entry
                c_status[Cij]   = CSR_ip;
                Ccolind(CSR_ip) = Cij;
                Cvals(CSR_ip)   = Aval*Ivals(j);
                CSR_ip++;

              } else {
                Cvals(c_status[Cij]) += Aval*Ivals(j);
              }
            }
          }
        }

        // Resize for next pass if needed
        if (i+1 < my_thread_stop && CSR_ip + std::min(n,(Arowptr(i+2)-Arowptr(i+1))*b_max_nnz_per_row) > CSR_alloc) {
          CSR_alloc *= 2;
          Ccolind = u_lno_nnz_view_t((typename u_lno_nnz_view_t::data_type)realloc(Ccolind.data(),u_lno_nnz_view_t::shmem_size(CSR_alloc)),CSR_alloc);
          Cvals = u_scalar_view_t((typename u_scalar_view_t::data_type)realloc(Cvals.data(),u_scalar_view_t::shmem_size(CSR_alloc)),CSR_alloc);
        }
        OLD_ip = CSR_ip;
      }

      tl_rowptr(tid) = Crowptr;
      tl_colind(tid) = Ccolind;
      tl_values(tid) = Cvals;
      Crowptr(my_thread_m) = CSR_ip;
  });

  // Do the copy out
  lno_view_t row_mapC("non_const_lnow_row", m + 1);
  lno_nnz_view_t  entriesC;
  scalar_view_t   valuesC;
  copy_out_from_thread_memory(tl_rowptr,tl_colind,tl_values,m,thread_chunk,row_mapC,entriesC,valuesC);

  //Free the unamanged views
  for(size_t i=0; i<thread_max; i++) {
    if(tl_rowptr(i).data()) free(tl_rowptr(i).data());
    if(tl_colind(i).data()) free(tl_colind(i).data());
    if(tl_values(i).data()) free(tl_values(i).data());
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix OpenMPSort"))));
#endif
    // Sort & set values
    if (params.is_null() || params->get("sort entries",true))
      Import_Util::sortCrsEntries(row_mapC, entriesC, valuesC);
    C.setAllValues(row_mapC,entriesC,valuesC);

}

/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void mult_A_B_reuse_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                 const LocalOrdinalViewType & targetMapToOrigRow,
                                                 const LocalOrdinalViewType & targetMapToImportRow,
                                                 const LocalOrdinalViewType & Bcol2Ccol,
                                                 const LocalOrdinalViewType & Icol2Ccol,
                                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                 Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                 const std::string& label,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Reuse LTGCore"))));
#endif

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode>::local_matrix_type KCRS;
  //  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode NO;
  typedef Map<LO,GO,NO>                     map_type;

  // NOTE (mfh 15 Sep 2017) This is specifically only for
  // execution_space = Kokkos::OpenMP, so we neither need nor want
  // KOKKOS_LAMBDA (with its mandatory __device__ marking).
  typedef NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;

  // All of the invalid guys
  const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS & Amat = Aview.origMatrix->getLocalMatrix();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrix();
  const KCRS & Cmat = C.getLocalMatrix();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map, Crowptr = Cmat.graph.row_map;
  const c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries, Ccolind = Cmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  scalar_view_t Cvals = Cmat.values;

  c_lno_view_t  Irowptr;
  c_lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    Irowptr = Bview.importMatrix->getLocalMatrix().graph.row_map;
    Icolind = Bview.importMatrix->getLocalMatrix().graph.entries;
    Ivals   = Bview.importMatrix->getLocalMatrix().values;
  }

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getNodeNumRows();
  size_t n = Ccolmap->getNodeNumElements();

  // Get my node / thread info (right from openmp or parameter list)
  size_t thread_max =  Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency();
  if(!params.is_null()) {
    if(params->isParameter("openmp: ltg thread max"))
      thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));
  }

  double thread_chunk = (double)(m) / thread_max;

  // Run chunks of the matrix independently
  Kokkos::parallel_for("MMM::LTG::Reuse::ThreadLocal",range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid)
    {
      // Thread coordination stuff
      size_t my_thread_start =  tid * thread_chunk;
      size_t my_thread_stop  = tid == thread_max-1 ? m : (tid+1)*thread_chunk;

      // Allocations
      std::vector<size_t> c_status(n,INVALID);

      // For each row of A/C
      size_t CSR_ip = 0, OLD_ip = 0;
      for (size_t i = my_thread_start; i < my_thread_stop; i++) {
        // First fill the c_status array w/ locations where we're allowed to
        // generate nonzeros for this row
        OLD_ip = Crowptr(i);
        CSR_ip = Crowptr(i+1);
        for (size_t k = OLD_ip; k < CSR_ip; k++) {
          c_status[Ccolind(k)] = k;
          // Reset values in the row of C
          Cvals(k) = SC_ZERO;
        }

        for (size_t k = Arowptr(i); k < Arowptr(i+1); k++) {
          LO Aik  = Acolind(k);
          const SC Aval = Avals(k);
          if (Aval == SC_ZERO)
            continue;

          if (targetMapToOrigRow(Aik) != LO_INVALID) {
            // Local matrix
            size_t Bk = Teuchos::as<size_t>(targetMapToOrigRow(Aik));

            for (size_t j = Browptr(Bk); j < Browptr(Bk+1); ++j) {
              LO Bkj = Bcolind(j);
              LO Cij = Bcol2Ccol(Bkj);

              TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
                                         std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
                                         "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

              Cvals(c_status[Cij]) += Aval * Bvals(j);
            }
          } else {
            // Remote matrix
            size_t Ik = Teuchos::as<size_t>(targetMapToImportRow(Aik));
            for (size_t j = Irowptr(Ik); j < Irowptr(Ik+1); ++j) {
              LO Ikj = Icolind(j);
              LO Cij = Icol2Ccol(Ikj);

              TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
                                         std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
                                         "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

              Cvals(c_status[Cij]) += Aval * Ivals(j);
            }
          }
        }
      }
    });

  // NOTE: No copy out or "set" of data is needed here, since we're working directly with Kokkos::Views
}

/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void jacobi_A_B_newmatrix_LowThreadGustavsonKernel(Scalar omega,
                                                   const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                   const LocalOrdinalViewType & targetMapToOrigRow,
                                                   const LocalOrdinalViewType & targetMapToImportRow,
                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                   const std::string& label,
                                                   const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix LTGCore"))));
#endif

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Kokkos::Compat::KokkosOpenMPWrapperNode Node;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type KCRS;
  //  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  // Unmanaged versions of the above
  typedef UnmanagedView<lno_view_t> u_lno_view_t;
  typedef UnmanagedView<lno_nnz_view_t> u_lno_nnz_view_t;
  typedef UnmanagedView<scalar_view_t> u_scalar_view_t;

  // Jacobi-specific
  typedef typename scalar_view_t::memory_space scalar_memory_space;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;
  typedef Map<LO,GO,NO>                     map_type;

  // NOTE (mfh 15 Sep 2017) This is specifically only for
  // execution_space = Kokkos::OpenMP, so we neither need nor want
  // KOKKOS_LAMBDA (with its mandatory __device__ marking).
  typedef NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;

  // All of the invalid guys
  const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS & Amat = Aview.origMatrix->getLocalMatrix();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrix();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  size_t b_max_nnz_per_row = Bview.origMatrix->getNodeMaxNumRowEntries();

  c_lno_view_t  Irowptr;
  lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    Irowptr = Bview.importMatrix->getLocalMatrix().graph.row_map;
    Icolind = Bview.importMatrix->getLocalMatrix().graph.entries;
    Ivals   = Bview.importMatrix->getLocalMatrix().values;
    b_max_nnz_per_row = std::max(b_max_nnz_per_row,Bview.importMatrix->getNodeMaxNumRowEntries());
  }

  // Jacobi-specific inner stuff
  auto Dvals = Dinv.template getLocalView<scalar_memory_space>();

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getNodeNumRows();
  size_t n = Ccolmap->getNodeNumElements();
  size_t Cest_nnz_per_row = 2*C_estimate_nnz_per_row(*Aview.origMatrix,*Bview.origMatrix);

  // Get my node / thread info (right from openmp)
  size_t thread_max =  Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency();
  if(!params.is_null()) {
    if(params->isParameter("openmp: ltg thread max"))
      thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));
  }

  // Thread-local memory
  Kokkos::View<u_lno_view_t*> tl_rowptr("top_rowptr",thread_max);
  Kokkos::View<u_lno_nnz_view_t*> tl_colind("top_colind",thread_max);
  Kokkos::View<u_scalar_view_t*> tl_values("top_values",thread_max);

  double thread_chunk = (double)(m) / thread_max;

  // Run chunks of the matrix independently
  Kokkos::parallel_for("Jacobi::LTG::NewMatrix::ThreadLocal",range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid)
    {
      // Thread coordination stuff
      size_t my_thread_start =  tid * thread_chunk;
      size_t my_thread_stop  = tid == thread_max-1 ? m : (tid+1)*thread_chunk;
      size_t my_thread_m     = my_thread_stop - my_thread_start;

      // Size estimate
      size_t CSR_alloc = (size_t) (my_thread_m*Cest_nnz_per_row*0.75 + 100);

      // Allocations
      std::vector<size_t> c_status(n,INVALID);

      u_lno_view_t Crowptr((typename u_lno_view_t::data_type)malloc(u_lno_view_t::shmem_size(my_thread_m+1)),my_thread_m+1);
      u_lno_nnz_view_t Ccolind((typename u_lno_nnz_view_t::data_type)malloc(u_lno_nnz_view_t::shmem_size(CSR_alloc)),CSR_alloc);
      u_scalar_view_t Cvals((typename u_scalar_view_t::data_type)malloc(u_scalar_view_t::shmem_size(CSR_alloc)),CSR_alloc);

      // For each row of A/C
      size_t CSR_ip = 0, OLD_ip = 0;
      for (size_t i = my_thread_start; i < my_thread_stop; i++) {
        //        printf("CMS: row %d CSR_alloc = %d\n",(int)i,(int)CSR_alloc);fflush(stdout);
        // mfh 27 Sep 2016: m is the number of rows in the input matrix A
        // on the calling process.
        Crowptr(i-my_thread_start) = CSR_ip;
        // NOTE: Vector::getLocalView returns a rank 2 view here
        SC minusOmegaDval = -omega*Dvals(i,0);

        // Entries of B
        for (size_t j = Browptr(i); j < Browptr(i+1); j++) {
          const SC Bval = Bvals(j);
          if (Bval == SC_ZERO)
            continue;
          LO Bij = Bcolind(j);
          LO Cij = Bcol2Ccol(Bij);

          // Assume no repeated entries in B
          c_status[Cij]   = CSR_ip;
          Ccolind(CSR_ip) = Cij;
          Cvals(CSR_ip)   = Bvals[j];
          CSR_ip++;
        }

        // Entries of -omega * Dinv * A * B
        // mfh 27 Sep 2016: For each entry of A in the current row of A
        for (size_t k = Arowptr(i); k < Arowptr(i+1); k++) {
          LO Aik  = Acolind(k); // local column index of current entry of A
          const SC Aval = Avals(k);   // value of current entry of A
          if (Aval == SC_ZERO)
            continue; // skip explicitly stored zero values in A

          if (targetMapToOrigRow(Aik) != LO_INVALID) {
            // mfh 27 Sep 2016: If the entry of targetMapToOrigRow
            // corresponding to the current entry of A is populated, then
            // the corresponding row of B is in B_local (i.e., it lives on
            // the calling process).

            // Local matrix
            size_t Bk = Teuchos::as<size_t>(targetMapToOrigRow(Aik));

            // mfh 27 Sep 2016: Go through all entries in that row of B_local.
            for (size_t j = Browptr(Bk); j < Browptr(Bk+1); ++j) {
              LO Bkj = Bcolind(j);
              LO Cij = Bcol2Ccol(Bkj);

              if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip) {
                // New entry
                c_status[Cij]   = CSR_ip;
                Ccolind(CSR_ip) = Cij;
                Cvals(CSR_ip)   = minusOmegaDval*Aval*Bvals(j);
                CSR_ip++;

              } else {
                Cvals(c_status[Cij]) += minusOmegaDval*Aval*Bvals(j);
              }
            }

          } else {
            // mfh 27 Sep 2016: If the entry of targetMapToOrigRow
            // corresponding to the current entry of A NOT populated (has
            // a flag "invalid" value), then the corresponding row of B is
            // in B_local (i.e., it lives on the calling process).

            // Remote matrix
            size_t Ik = Teuchos::as<size_t>(targetMapToImportRow(Aik));
            for (size_t j = Irowptr(Ik); j < Irowptr(Ik+1); ++j) {
              LO Ikj = Icolind(j);
              LO Cij = Icol2Ccol(Ikj);

              if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip){
                // New entry
                c_status[Cij]   = CSR_ip;
                Ccolind(CSR_ip) = Cij;
                Cvals(CSR_ip)   = minusOmegaDval*Aval*Ivals(j);
                CSR_ip++;

              } else {
                Cvals(c_status[Cij]) += minusOmegaDval*Aval*Ivals(j);
              }
            }
          }
        }

        // Resize for next pass if needed
        if (i+1 < my_thread_stop && CSR_ip + std::min(n,(Arowptr(i+2)-Arowptr(i+1)+1)*b_max_nnz_per_row) > CSR_alloc) {
          CSR_alloc *= 2;
          Ccolind = u_lno_nnz_view_t((typename u_lno_nnz_view_t::data_type)realloc(Ccolind.data(),u_lno_nnz_view_t::shmem_size(CSR_alloc)),CSR_alloc);
          Cvals = u_scalar_view_t((typename u_scalar_view_t::data_type)realloc(Cvals.data(),u_scalar_view_t::shmem_size(CSR_alloc)),CSR_alloc);
        }
        OLD_ip = CSR_ip;
      }

      tl_rowptr(tid) = Crowptr;
      tl_colind(tid) = Ccolind;
      tl_values(tid) = Cvals;
      Crowptr(my_thread_m) = CSR_ip;
  });



  // Do the copy out
  lno_view_t row_mapC("non_const_lnow_row", m + 1);
  lno_nnz_view_t  entriesC;
  scalar_view_t   valuesC;
  copy_out_from_thread_memory(tl_rowptr,tl_colind,tl_values,m,thread_chunk,row_mapC,entriesC,valuesC);

  //Free the unamanged views
  for(size_t i=0; i<thread_max; i++) {
    if(tl_rowptr(i).data()) free(tl_rowptr(i).data());
    if(tl_colind(i).data()) free(tl_colind(i).data());
    if(tl_values(i).data()) free(tl_values(i).data());
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix OpenMPSort"))));
#endif
    // Sort & set values
    if (params.is_null() || params->get("sort entries",true))
      Import_Util::sortCrsEntries(row_mapC, entriesC, valuesC);
    C.setAllValues(row_mapC,entriesC,valuesC);

}



/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void jacobi_A_B_reuse_LowThreadGustavsonKernel(Scalar omega,
                                                   const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode> & Dinv,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Bview,
                                                   const LocalOrdinalViewType & targetMapToOrigRow,
                                                   const LocalOrdinalViewType & targetMapToImportRow,
                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& C,
                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Cimport,
                                                   const std::string& label,
                                                   const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Reuse LTGCore"))));
#endif
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Kokkos::Compat::KokkosOpenMPWrapperNode Node;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type KCRS;
  //  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  // Jacobi-specific
  typedef typename scalar_view_t::memory_space scalar_memory_space;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;
  typedef Map<LO,GO,NO>                     map_type;

  // NOTE (mfh 15 Sep 2017) This is specifically only for
  // execution_space = Kokkos::OpenMP, so we neither need nor want
  // KOKKOS_LAMBDA (with its mandatory __device__ marking).
  typedef NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;

  // All of the invalid guys
  const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS & Amat = Aview.origMatrix->getLocalMatrix();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrix();
  const KCRS & Cmat = C.getLocalMatrix();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map, Crowptr = Cmat.graph.row_map;
  const c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries, Ccolind = Cmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  scalar_view_t Cvals = Cmat.values;

  c_lno_view_t  Irowptr;
  c_lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    Irowptr = Bview.importMatrix->getLocalMatrix().graph.row_map;
    Icolind = Bview.importMatrix->getLocalMatrix().graph.entries;
    Ivals   = Bview.importMatrix->getLocalMatrix().values;
  }

  // Jacobi-specific inner stuff
  auto Dvals = Dinv.template getLocalView<scalar_memory_space>();

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getNodeNumRows();
  size_t n = Ccolmap->getNodeNumElements();

  // Get my node / thread info (right from openmp or parameter list)
  size_t thread_max =  Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency();
  if(!params.is_null()) {
    if(params->isParameter("openmp: ltg thread max"))
      thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));
  }

  double thread_chunk = (double)(m) / thread_max;

  // Run chunks of the matrix independently
  Kokkos::parallel_for("Jacobi::LTG::Reuse::ThreadLocal",range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid)
    {
      // Thread coordination stuff
      size_t my_thread_start =  tid * thread_chunk;
      size_t my_thread_stop  = tid == thread_max-1 ? m : (tid+1)*thread_chunk;

      // Allocations
      std::vector<size_t> c_status(n,INVALID);

      // For each row of A/C
      size_t CSR_ip = 0, OLD_ip = 0;
      for (size_t i = my_thread_start; i < my_thread_stop; i++) {
        // First fill the c_status array w/ locations where we're allowed to
        // generate nonzeros for this row
        OLD_ip = Crowptr(i);
        CSR_ip = Crowptr(i+1);
        // NOTE: Vector::getLocalView returns a rank 2 view here
        SC minusOmegaDval = -omega*Dvals(i,0);

        for (size_t k = OLD_ip; k < CSR_ip; k++) {
          c_status[Ccolind(k)] = k;
          // Reset values in the row of C
          Cvals(k) = SC_ZERO;
        }

        // Entries of B
        for (size_t j = Browptr(i); j < Browptr(i+1); j++) {
          const SC Bval = Bvals(j);
          if (Bval == SC_ZERO)
            continue;
          LO Bij = Bcolind(j);
          LO Cij = Bcol2Ccol(Bij);

          // Assume no repeated entries in B
          Cvals(c_status[Cij]) += Bvals(j);
          CSR_ip++;
        }


        for (size_t k = Arowptr(i); k < Arowptr(i+1); k++) {
          LO Aik  = Acolind(k);
          const SC Aval = Avals(k);
          if (Aval == SC_ZERO)
            continue;

          if (targetMapToOrigRow(Aik) != LO_INVALID) {
            // Local matrix
            size_t Bk = Teuchos::as<size_t>(targetMapToOrigRow(Aik));

            for (size_t j = Browptr(Bk); j < Browptr(Bk+1); ++j) {
              LO Bkj = Bcolind(j);
              LO Cij = Bcol2Ccol(Bkj);

              TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
                                         std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
                                         "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

              Cvals(c_status[Cij]) += minusOmegaDval * Aval * Bvals(j);
            }
          } else {
            // Remote matrix
            size_t Ik = Teuchos::as<size_t>(targetMapToImportRow(Aik));
            for (size_t j = Irowptr(Ik); j < Irowptr(Ik+1); ++j) {
              LO Ikj = Icolind(j);
              LO Cij = Icol2Ccol(Ikj);

              TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
                                         std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
                                         "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

              Cvals(c_status[Cij]) += minusOmegaDval * Aval * Ivals(j);
            }
          }
        }
      }
    });

  // NOTE: No copy out or "set" of data is needed here, since we're working directly with Kokkos::Views
}


/*********************************************************************************************************/
template<class InRowptrArrayType, class InColindArrayType, class InValsArrayType,
         class OutRowptrType, class OutColindType, class OutValsType>
void copy_out_from_thread_memory(const InRowptrArrayType & Inrowptr, const InColindArrayType &Incolind, const InValsArrayType & Invalues,
                                   size_t m, double thread_chunk,
                                   OutRowptrType & row_mapC, OutColindType &entriesC, OutValsType & valuesC ) {
  typedef OutRowptrType lno_view_t;
  typedef OutColindType lno_nnz_view_t;
  typedef OutValsType scalar_view_t;
  typedef typename lno_view_t::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;

  // Generate the starting nnz number per thread
  size_t thread_max =  Inrowptr.size();
  size_t c_nnz_size=0;
  lno_view_t thread_start_nnz("thread_nnz",thread_max+1);
  Kokkos::parallel_scan("LTG::Scan",range_type(0,thread_max).set_chunk_size(1), [=] (const size_t i, size_t& update, const bool final) {
      size_t mynnz = Inrowptr(i)(Inrowptr(i).extent(0)-1);
      if(final) thread_start_nnz(i) = update;
      update+=mynnz;
      if(final && i+1==thread_max) thread_start_nnz(i+1)=update;
    });
  c_nnz_size = thread_start_nnz(thread_max);

  // Allocate output
  lno_nnz_view_t  entriesC_(Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size); entriesC = entriesC_;
  scalar_view_t   valuesC_(Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);  valuesC = valuesC_;

  // Copy out
  Kokkos::parallel_for("LTG::CopyOut", range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid) {
      size_t my_thread_start =  tid * thread_chunk;
      size_t my_thread_stop  = tid == thread_max-1 ? m : (tid+1)*thread_chunk;
      size_t nnz_thread_start = thread_start_nnz(tid);

      for (size_t i = my_thread_start; i < my_thread_stop; i++) {
        size_t ii = i - my_thread_start;
        // Rowptr
        row_mapC(i) = nnz_thread_start + Inrowptr(tid)(ii);
        if (i==m-1) {
          row_mapC(m) = nnz_thread_start + Inrowptr(tid)(ii+1);
        }

        // Colind / Values
        for(size_t j = Inrowptr(tid)(ii); j<Inrowptr(tid)(ii+1); j++) {
          entriesC(nnz_thread_start + j) = Incolind(tid)(j);
          valuesC(nnz_thread_start + j)  = Invalues(tid)(j);
        }
      }
    });
}//end copy_out

#endif // OpenMP


/*********************************************************************************************************/
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalOrdinalViewType>
void jacobi_A_B_newmatrix_MultiplyScaleAddKernel(Scalar omega,
                                                  const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Node> & Dinv,
                                                  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                  const LocalOrdinalViewType & Acol2Brow,
                                                  const LocalOrdinalViewType & Acol2Irow,
                                                  const LocalOrdinalViewType & Bcol2Ccol,
                                                  const LocalOrdinalViewType & Icol2Ccol,
                                                  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                  Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                  const std::string& label,
                                                  const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix MSAK"))));
  Teuchos::RCP<Teuchos::TimeMonitor> MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix MSAK Multiply"))));
#endif
  typedef  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix_t;

  // This kernel computes (I-omega Dinv A) B the slow way (for generality's sake, you must understand)

  // 1) Multiply A*B
  Teuchos::RCP<Matrix_t> AB = Teuchos::rcp(new Matrix_t(C.getRowMap(),0));
  Tpetra::MMdetails::mult_A_B_newmatrix(Aview,Bview,*AB,label+std::string(" MSAK"),params);

#ifdef HAVE_TPETRA_MMM_TIMINGS
MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix MSAK Scale"))));
#endif

  // 2) Scale A by Dinv
  AB->leftScale(Dinv);

#ifdef HAVE_TPETRA_MMM_TIMINGS
MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix MSAK Add"))));
#endif

  // 3) Add [-omega Dinv A] + B
 Teuchos::ParameterList jparams;
  if(!params.is_null()) {
    jparams = *params;
    jparams.set("label",label+std::string(" MSAK Add"));
  }
  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  Tpetra::MatrixMatrix::add(one,false,*Bview.origMatrix,Scalar(-omega),false,*AB,C,AB->getDomainMap(),AB->getRangeMap(),Teuchos::rcp(&jparams,false));

 }// jacobi_A_B_newmatrix_MultiplyScaleAddKernel



}//ExtraKernels
}//MatrixMatrix
}//Tpetra


#endif
