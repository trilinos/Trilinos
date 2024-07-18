// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_EXTRAKERNELS_DEF_HPP
#define TPETRA_MATRIXMATRIX_EXTRAKERNELS_DEF_HPP
#include "TpetraExt_MatrixMatrix_ExtraKernels_decl.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

namespace Tpetra {

namespace MatrixMatrix{

namespace ExtraKernels{


// Double product version
template<class CrsMatrixType>
size_t C_estimate_nnz_per_row(CrsMatrixType & A, CrsMatrixType &B){
  // Follows the NZ estimate in ML's ml_matmatmult.c
  size_t Aest = 100, Best=100;
  if (A.getLocalNumEntries() > 0)
    Aest = (A.getLocalNumRows() > 0)?  A.getLocalNumEntries()/A.getLocalNumRows() : 100;
  if (B.getLocalNumEntries() > 0)
    Best = (B.getLocalNumRows() > 0) ? B.getLocalNumEntries()/B.getLocalNumRows() : 100;

  size_t nnzperrow = (size_t)(sqrt((double)Aest) + sqrt((double)Best) - 1);
  nnzperrow *= nnzperrow;

  return nnzperrow;
}


// Triple product version
template<class CrsMatrixType>
size_t Ac_estimate_nnz(CrsMatrixType & A, CrsMatrixType &P){
  size_t nnzPerRowA = 100, Pcols = 100;
  if (A.getLocalNumEntries() > 0)
    nnzPerRowA = (A.getLocalNumRows() > 0)?  A.getLocalNumEntries()/A.getLocalNumRows() : 9;
  if (P.getLocalNumEntries() > 0)
    Pcols = (P.getLocalNumCols() > 0) ? P.getLocalNumCols() : 100;
  return (size_t)(Pcols*nnzPerRowA + 5*nnzPerRowA + 300);
}

#if defined (HAVE_TPETRA_INST_OPENMP)
/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class LocalOrdinalViewType>
void mult_A_B_newmatrix_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                 const LocalOrdinalViewType & targetMapToOrigRow,
                                                 const LocalOrdinalViewType & targetMapToImportRow,
                                                 const LocalOrdinalViewType & Bcol2Ccol,
                                                 const LocalOrdinalViewType & Icol2Ccol,
                                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                 Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                 const std::string& label,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix LTGCore"))));
#endif
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;


  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>::local_matrix_device_type KCRS;
  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  // Unmanaged versions of the above
  //typedef UnmanagedView<lno_view_t> u_lno_view_t; // unused
  typedef UnmanagedView<lno_nnz_view_t> u_lno_nnz_view_t;
  typedef UnmanagedView<scalar_view_t> u_scalar_view_t;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode NO;
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
  const KCRS & Amat = Aview.origMatrix->getLocalMatrixDevice();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrixDevice();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  size_t b_max_nnz_per_row = Bview.origMatrix->getLocalMaxNumRowEntries();

  c_lno_view_t  Irowptr;
  lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixDevice();
    Irowptr = lclB.graph.row_map;
    Icolind = lclB.graph.entries;
    Ivals   = lclB.values;
    b_max_nnz_per_row = std::max(b_max_nnz_per_row,Bview.importMatrix->getLocalMaxNumRowEntries());
  }

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getLocalNumRows();
  size_t n = Ccolmap->getLocalNumElements();
  size_t Cest_nnz_per_row = 2*C_estimate_nnz_per_row(*Aview.origMatrix,*Bview.origMatrix);

  // Get my node / thread info (right from openmp or parameter list)
  size_t thread_max =  Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency();
  if(!params.is_null()) {
    if(params->isParameter("openmp: ltg thread max"))
      thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));
  }

  // 2019 Apr 10 jje: We can do rowptr in place, and no need to inialize since we can fault as we go
  lno_view_t row_mapC(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
  // we will not touch these until the final copyout step
  lno_nnz_view_t entriesC;
  scalar_view_t  valuesC;

  // add this, since we do the rowptr in place, we could figure this out
  // using the rowptr, but it requires an unusual loop (that jumps all over memory)
  lno_nnz_view_t thread_total_nnz("thread_total_nnz",thread_max+1);

  // Thread-local memory
  Kokkos::View<u_lno_nnz_view_t*, device_t> tl_colind("top_colind",thread_max);
  Kokkos::View<u_scalar_view_t*, device_t> tl_values("top_values",thread_max);

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

      u_lno_nnz_view_t Ccolind((typename u_lno_nnz_view_t::data_type)malloc(u_lno_nnz_view_t::shmem_size(CSR_alloc)),CSR_alloc);
      u_scalar_view_t Cvals((typename u_scalar_view_t::data_type)malloc(u_scalar_view_t::shmem_size(CSR_alloc)),CSR_alloc);

      // For each row of A/C
      size_t CSR_ip = 0, OLD_ip = 0;
      for (size_t i = my_thread_start; i < my_thread_stop; i++) {
        // mfh 27 Sep 2016: m is the number of rows in the input matrix A
        // on the calling process.
        // JJE 10 Apr 2019 index directly into the rowptr
        row_mapC(i) = CSR_ip;

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
          Cvals = u_scalar_view_t((typename u_scalar_view_t::data_type)realloc((void*) Cvals.data(),u_scalar_view_t::shmem_size(CSR_alloc)),CSR_alloc);
        }
        OLD_ip = CSR_ip;
      }
      thread_total_nnz(tid) = CSR_ip;
      tl_colind(tid) = Ccolind;
      tl_values(tid) = Cvals;
  });

  // Do the copy out
  copy_out_from_thread_memory(thread_total_nnz,tl_colind,tl_values,m,thread_chunk,row_mapC,entriesC,valuesC);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix OpenMPSort"))));
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
void mult_A_B_reuse_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                 const LocalOrdinalViewType & targetMapToOrigRow,
                                                 const LocalOrdinalViewType & targetMapToImportRow,
                                                 const LocalOrdinalViewType & Bcol2Ccol,
                                                 const LocalOrdinalViewType & Icol2Ccol,
                                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                 Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                 const std::string& label,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Reuse LTGCore"))));
#endif

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>::local_matrix_device_type KCRS;
  //  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode NO;
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
  const KCRS & Amat = Aview.origMatrix->getLocalMatrixDevice();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrixDevice();
  const KCRS & Cmat = C.getLocalMatrixDevice();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map, Crowptr = Cmat.graph.row_map;
  const c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries, Ccolind = Cmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  scalar_view_t Cvals = Cmat.values;

  c_lno_view_t  Irowptr;
  c_lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixDevice();
    Irowptr = lclB.graph.row_map;
    Icolind = lclB.graph.entries;
    Ivals   = lclB.values;
  }

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getLocalNumRows();
  size_t n = Ccolmap->getLocalNumElements();

  // Get my node / thread info (right from openmp or parameter list)
  size_t thread_max =  Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency();
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
                                                   const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> & Dinv,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                   const LocalOrdinalViewType & targetMapToOrigRow,
                                                   const LocalOrdinalViewType & targetMapToImportRow,
                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                   const std::string& label,
                                                   const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": "); 
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix LTGCore"))));
#endif

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type KCRS;
  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  // Unmanaged versions of the above
  //typedef UnmanagedView<lno_view_t> u_lno_view_t; // unused
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
  const KCRS & Amat = Aview.origMatrix->getLocalMatrixDevice();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrixDevice();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  size_t b_max_nnz_per_row = Bview.origMatrix->getLocalMaxNumRowEntries();

  c_lno_view_t  Irowptr;
  lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixDevice();
    Irowptr = lclB.graph.row_map;
    Icolind = lclB.graph.entries;
    Ivals   = lclB.values;
    b_max_nnz_per_row = std::max(b_max_nnz_per_row,Bview.importMatrix->getLocalMaxNumRowEntries());
  }

  // Jacobi-specific inner stuff
  auto Dvals = 
       Dinv.template getLocalView<scalar_memory_space>(Access::ReadOnly);

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getLocalNumRows();
  size_t n = Ccolmap->getLocalNumElements();
  size_t Cest_nnz_per_row = 2*C_estimate_nnz_per_row(*Aview.origMatrix,*Bview.origMatrix);

  // Get my node / thread info (right from openmp)
  size_t thread_max =  Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency();
  if(!params.is_null()) {
    if(params->isParameter("openmp: ltg thread max"))
      thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));
  }

  // 2019 Apr 10 jje: We can do rowptr in place, and no need to inialize since we can fault as we go
  lno_view_t row_mapC(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
  // we will not touch these until the final copyout step
  lno_nnz_view_t entriesC;
  scalar_view_t  valuesC;

  // add this, since we do the rowptr in place, we could figure this out
  // using the rowptr, but it requires an unusual loop (that jumps all over memory)
  lno_nnz_view_t thread_total_nnz("thread_total_nnz",thread_max+1);

  // Thread-local memory
  Kokkos::View<u_lno_nnz_view_t*, device_t> tl_colind("top_colind",thread_max);
  Kokkos::View<u_scalar_view_t*, device_t> tl_values("top_values",thread_max);

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

      u_lno_nnz_view_t Ccolind((typename u_lno_nnz_view_t::data_type)malloc(u_lno_nnz_view_t::shmem_size(CSR_alloc)),CSR_alloc);
      u_scalar_view_t Cvals((typename u_scalar_view_t::data_type)malloc(u_scalar_view_t::shmem_size(CSR_alloc)),CSR_alloc);

      // For each row of A/C
      size_t CSR_ip = 0, OLD_ip = 0;
      for (size_t i = my_thread_start; i < my_thread_stop; i++) {
        //        printf("CMS: row %d CSR_alloc = %d\n",(int)i,(int)CSR_alloc);fflush(stdout);
        // mfh 27 Sep 2016: m is the number of rows in the input matrix A
        // on the calling process.
        // JJE: directly access the rowptr here (indexed using our loop var)
        row_mapC(i) = CSR_ip;
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
          Cvals = u_scalar_view_t((typename u_scalar_view_t::data_type)realloc((void*) Cvals.data(),u_scalar_view_t::shmem_size(CSR_alloc)),CSR_alloc);
        }
        OLD_ip = CSR_ip;
      }
      thread_total_nnz(tid) = CSR_ip;
      tl_colind(tid) = Ccolind;
      tl_values(tid) = Cvals;
  });


  // Do the copy out (removed the tl_rowptr!)
  copy_out_from_thread_memory(thread_total_nnz,tl_colind,tl_values,m,thread_chunk,row_mapC,entriesC,valuesC);


#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix OpenMPSort"))));
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
                                                   const Vector<Scalar,LocalOrdinal,GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> & Dinv,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                   CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Bview,
                                                   const LocalOrdinalViewType & targetMapToOrigRow,
                                                   const LocalOrdinalViewType & targetMapToImportRow,
                                                   const LocalOrdinalViewType & Bcol2Ccol,
                                                   const LocalOrdinalViewType & Icol2Ccol,
                                                   CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& C,
                                                   Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Cimport,
                                                   const std::string& label,
                                                   const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Reuse LTGCore"))));
#endif
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type KCRS;
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
  const KCRS & Amat = Aview.origMatrix->getLocalMatrixDevice();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrixDevice();
  const KCRS & Cmat = C.getLocalMatrixDevice();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map, Crowptr = Cmat.graph.row_map;
  const c_lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries, Ccolind = Cmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  scalar_view_t Cvals = Cmat.values;

  c_lno_view_t  Irowptr;
  c_lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixDevice();
    Irowptr = lclB.graph.row_map;
    Icolind = lclB.graph.entries;
    Ivals   = lclB.values;
  }

  // Jacobi-specific inner stuff
  auto Dvals = 
       Dinv.template getLocalView<scalar_memory_space>(Access::ReadOnly);

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getLocalNumRows();
  size_t n = Ccolmap->getLocalNumElements();

  // Get my node / thread info (right from openmp or parameter list)
  size_t thread_max =  Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency();
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
template<class InColindArrayType, class InValsArrayType,
         class OutRowptrType, class OutColindType, class OutValsType>
void copy_out_from_thread_memory(const OutColindType& thread_total_nnz,
                                 const InColindArrayType& Incolind,
                                 const InValsArrayType& Invalues,
                                 const size_t m,
                                 const double thread_chunk,
                                 OutRowptrType& row_mapC,
                                 OutColindType& entriesC,
                                 OutValsType& valuesC ) {
  typedef OutRowptrType lno_view_t;
  typedef OutColindType lno_nnz_view_t;
  typedef OutValsType scalar_view_t;
  typedef typename lno_view_t::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;

  // Generate the starting nnz number per thread
  size_t thread_max =  Incolind.size();
  size_t c_nnz_size=0;
  // since this kernel is 'low thread count' is very likely, that we could
  // sum over the thread_total_nnz and not go parallel, but since it is a view
  // we don't know where the memory actually lives... so we kinda need to go parallel

  lno_view_t thread_start_nnz("thread_nnz",thread_max+1);
  Kokkos::parallel_scan("LTG::Scan",range_type(0,thread_max).set_chunk_size(1), [=] (const size_t i, size_t& update, const bool final) {
      size_t mynnz = thread_total_nnz(i);
      if(final) thread_start_nnz(i) = update;
      update+=mynnz;
      if(final && i+1==thread_max) thread_start_nnz(i+1)=update;
    });
  c_nnz_size = thread_start_nnz(thread_max);

  // 2019 Apr 10 JJE: update the rowptr's final entry here
  row_mapC(m) = thread_start_nnz(thread_max);

  // Allocate output
  lno_nnz_view_t  entriesC_(Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size); entriesC = entriesC_;
  scalar_view_t   valuesC_(Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);  valuesC = valuesC_;

  // Copy out
  Kokkos::parallel_for("LTG::CopyOut", range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid) {
    const size_t my_thread_start =  tid * thread_chunk;
    const size_t my_thread_stop  = tid == thread_max-1 ? m : (tid+1)*thread_chunk;
    const size_t nnz_thread_start = thread_start_nnz(tid);
    // we need this info, since we did the rowptr in place
    const size_t nnz_thread_stop  = thread_start_nnz(tid+1);

    // There are two fundamental operations:
    // * Updateing the rowptr with the correct offset
    // * copying entries and values

    // First, update the rowptr, it since we did it inplace, it is a += operation
    // in the paper, I experimented with a coloring scheme that had threads do
    // do their copies in different orders. It wasn't obvious it was beneficial
    // but it can be replicated, by choosing the array to copy first based on your
    // thread id % 3
    if (my_thread_start != 0 ) {
      for (size_t i = my_thread_start; i < my_thread_stop; i++) {
        row_mapC(i) += nnz_thread_start;
      }
    }

    // try to Kokkos::single() the alloc here. It should implicitly barrier
    // thread 0 doesn't copy the rowptr, so it could hit the single first
    // in the paper, I played a game that effectively got LTG down to a single
    // OpenMP barrier. But doing that requires the ability to spawn a single
    // parallel region.  The Scan above was implemented using atomic adds
    // and the barrier was only needed so you could allocate
    //
    // Since we can't spawn a single region, we could move the allocations
    // here, and using 'single'. Most likely, thread 0 will hit it first allowing
    // the other threads to update the rowptr while it allocates.


    // next, bulk copy the vals/colind arrays
    const size_t my_num_nnz = nnz_thread_stop - nnz_thread_start;
    for (size_t i = 0; i < my_num_nnz; ++i) {
      entriesC(nnz_thread_start + i) = Incolind(tid)(i);
      valuesC(nnz_thread_start + i)  = Invalues(tid)(i);
    }

    //Free the unamanged views, let each thread deallocate its memory
    // May need to cast away const here..
    if(Incolind(tid).data()) free(Incolind(tid).data());
    if(Invalues(tid).data()) free(Invalues(tid).data());
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
  Teuchos::RCP<TimeMonitor> MM = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix MSAK"))));
  Teuchos::RCP<TimeMonitor> MM2 = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix MSAK Multiply"))));
  using Teuchos::rcp;
#endif
  typedef  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix_t;

  // This kernel computes (I-omega Dinv A) B the slow way (for generality's sake, you must understand)

  // 1) Multiply A*B
  Teuchos::RCP<Matrix_t> AB = Teuchos::rcp(new Matrix_t(C.getRowMap(),0));
  Tpetra::MMdetails::mult_A_B_newmatrix(Aview,Bview,*AB,label+std::string(" MSAK"),params);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM2=Teuchos::null; MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix MSAK Scale"))));
#endif

  // 2) Scale A by Dinv
  AB->leftScale(Dinv);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM2=Teuchos::null; MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix MSAK Add"))));
#endif

  // 3) Add [-omega Dinv A] + B
 Teuchos::ParameterList jparams;
  if(!params.is_null()) {
    jparams = *params;
    jparams.set("label",label+std::string(" MSAK Add"));
  }
  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  Tpetra::MatrixMatrix::add(one,false,*Bview.origMatrix,Scalar(-omega),false,*AB,C,AB->getDomainMap(),AB->getRangeMap(),Teuchos::rcp(&jparams,false));
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM2=Teuchos::null;
#endif
 }// jacobi_A_B_newmatrix_MultiplyScaleAddKernel



#if defined (HAVE_TPETRA_INST_OPENMP)
/*********************************************************************************************************/
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
static inline void mult_R_A_P_newmatrix_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Rview,
                                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                                 const LocalOrdinalViewType & Acol2PRow,
                                                                 const LocalOrdinalViewType & Acol2PRowImport,
                                                                 const LocalOrdinalViewType & Pcol2Accol,
                                                                 const LocalOrdinalViewType & PIcol2Accol,
                                                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                                 Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                                 const std::string& label,
                                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {

        using Teuchos::RCP;
        using Tpetra::MatrixMatrix::UnmanagedView;
  #ifdef HAVE_TPETRA_MMM_TIMINGS
        std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
        using Teuchos::rcp;
        using Teuchos::TimeMonitor;
        RCP<TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix LTGCore"))));
  #endif

        typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;
        typedef Scalar        SC;
        typedef LocalOrdinal  LO;
        typedef GlobalOrdinal GO;
        typedef Node          NO;
        typedef Map<LO,GO,NO> map_type;
        typedef typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_matrix_device_type KCRS;
        typedef typename KCRS::StaticCrsGraphType graph_t;
        typedef typename graph_t::row_map_type::non_const_type lno_view_t;
        typedef typename graph_t::row_map_type::const_type c_lno_view_t;
        typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
        typedef typename KCRS::values_type::non_const_type scalar_view_t;
        typedef typename KCRS::device_type device_t;
        typedef typename device_t::execution_space execution_space;
        typedef Kokkos::RangePolicy<execution_space, size_t> range_type;

        // Unmanaged versions of the above
        typedef UnmanagedView<lno_view_t> u_lno_view_t;
        typedef UnmanagedView<lno_nnz_view_t> u_lno_nnz_view_t;
        typedef UnmanagedView<scalar_view_t> u_scalar_view_t;

        const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
        const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

        // Sizes
        RCP<const map_type> Accolmap = Ac.getColMap();
        size_t m = Rview.origMatrix->getLocalNumRows();
        size_t n = Accolmap->getLocalNumElements();

        // Get raw Kokkos matrices, and the raw CSR views
        const KCRS & Rmat = Rview.origMatrix->getLocalMatrixDevice();
        const KCRS & Amat = Aview.origMatrix->getLocalMatrixDevice();
        const KCRS & Pmat = Pview.origMatrix->getLocalMatrixDevice();

        c_lno_view_t Rrowptr = Rmat.graph.row_map, 
                     Arowptr = Amat.graph.row_map, 
                     Prowptr = Pmat.graph.row_map, Irowptr;
        const lno_nnz_view_t Rcolind = Rmat.graph.entries, 
                             Acolind = Amat.graph.entries, 
                             Pcolind = Pmat.graph.entries;
        lno_nnz_view_t Icolind;
        const scalar_view_t Rvals = Rmat.values, 
                            Avals = Amat.values, 
                            Pvals = Pmat.values;
        scalar_view_t Ivals;

        if (!Pview.importMatrix.is_null())
        {
          const KCRS& Imat = Pview.importMatrix->getLocalMatrixDevice();
          Irowptr = Imat.graph.row_map;
          Icolind = Imat.graph.entries;
          Ivals = Imat.values;
        }

        // Classic csr assembly (low memory edition)
        //
        // mfh 27 Sep 2016: Ac_estimate_nnz does not promise an upper bound.
        // The method loops over rows of R, and may resize after processing
        // each row.  Chris Siefert says that this reflects experience in
        // ML; for the non-threaded case, ML found it faster to spend less
        // effort on estimation and risk an occasional reallocation.

        size_t Acest_nnz_per_row = std::ceil(Ac_estimate_nnz(*Aview.origMatrix, *Pview.origMatrix) / m);
        size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

        // Get my node / thread info (right from openmp or parameter list)
        size_t thread_max =  Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency();
        if(!params.is_null()) {
          if(params->isParameter("openmp: ltg thread max"))
            thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));
        }

        double thread_chunk = (double)(m) / thread_max;

        // we can construct the rowptr inplace, allocate here and fault in parallel
        lno_view_t rowmapAc(Kokkos::ViewAllocateWithoutInitializing("non_const_lnow_row"), m + 1);
        // we do not touch these until copy out
        lno_nnz_view_t entriesAc;
        scalar_view_t valuesAc;

        // add this, since we do the rowptr in place, we could figure this out
        // using the rowptr, but it requires an unusual loop (that jumps all over memory)
        lno_nnz_view_t thread_total_nnz("thread_total_nnz",thread_max+1);

        // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
        // routine.  The routine computes Ac := R * A * (P_local + P_remote).
        //
        // For column index Aik in row i of A, Acol2PRow[Aik] tells
        // you whether the corresponding row of P belongs to P_local
        // ("orig") or P_remote ("Import").

        // Thread-local memory
        Kokkos::View<u_lno_nnz_view_t*, device_t> tl_colind("top_colind", thread_max);
        Kokkos::View<u_scalar_view_t*, device_t> tl_values("top_values", thread_max);

        // For each row of R
        Kokkos::parallel_for("MMM::RAP::NewMatrix::LTG::ThreadLocal",range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid)
        {
          // Thread coordination stuff
          size_t my_thread_start = tid * thread_chunk;
          size_t my_thread_stop  = tid == thread_max-1 ? m : (tid+1)*thread_chunk;
          size_t my_thread_m     = my_thread_stop - my_thread_start;

          size_t nnzAllocated = (size_t) (my_thread_m * Acest_nnz_per_row + 100);

          std::vector<size_t> ac_status(n, INVALID);

          //manually allocate the thread-local storage for Ac
          u_lno_view_t Acrowptr((typename u_lno_view_t::data_type) malloc(u_lno_view_t::shmem_size(my_thread_m+1)), my_thread_m + 1);
          u_lno_nnz_view_t Accolind((typename u_lno_nnz_view_t::data_type) malloc(u_lno_nnz_view_t::shmem_size(nnzAllocated)), nnzAllocated);
          u_scalar_view_t Acvals((typename u_scalar_view_t::data_type) malloc(u_scalar_view_t::shmem_size(nnzAllocated)), nnzAllocated);

          //count entries as they are added to Ac
          size_t nnz = 0, nnz_old = 0;
          // bmk: loop over the rows of R which are assigned to thread tid
          for (size_t i = my_thread_start; i < my_thread_stop; i++) {
            // directly index into the rowptr
            rowmapAc(i) = nnz;
            // mfh 27 Sep 2016: For each entry of R in the current row of R
            for (size_t kk = Rrowptr(i); kk < Rrowptr(i+1); kk++) {
              LO k  = Rcolind(kk); // local column index of current entry of R
              const SC Rik = Rvals(kk);   // value of current entry of R
              if (Rik == SC_ZERO)
                continue; // skip explicitly stored zero values in R
              // For each entry of A in the current row of A
              for (size_t ll = Arowptr(k); ll < Arowptr(k+1); ll++) {
                LO l = Acolind(ll); // local column index of current entry of A
                const SC Akl = Avals(ll);   // value of current entry of A
                if (Akl == SC_ZERO)
                  continue; // skip explicitly stored zero values in A

                if (Acol2PRow[l] != LO_INVALID) {
                  // mfh 27 Sep 2016: If the entry of Acol2PRow
                  // corresponding to the current entry of A is populated, then
                  // the corresponding row of P is in P_local (i.e., it lives on
                  // the calling process).

                  // Local matrix
                  size_t Pl = Teuchos::as<size_t>(Acol2PRow(l));

                  // mfh 27 Sep 2016: Go through all entries in that row of P_local.
                  for (size_t jj = Prowptr(Pl); jj < Prowptr(Pl+1); jj++) {
                    LO j = Pcolind(jj);
                    LO Acj = Pcol2Accol(j);
                    SC Plj = Pvals(jj);

                    if (ac_status[Acj] == INVALID || ac_status[Acj] < nnz_old) {
    #ifdef HAVE_TPETRA_DEBUG
                      // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                      TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.size()),
                                                 std::runtime_error,
                                                 label << " ERROR, not enough memory allocated for matrix product. Allocated: " << Accolind.size() << std::endl);
    #endif
                      // New entry
                      ac_status[Acj] = nnz;
                      Accolind(nnz) = Acj;
                      Acvals(nnz) = Rik*Akl*Plj;
                      nnz++;
                    } else {
                      Acvals(ac_status[Acj]) += Rik*Akl*Plj;
                    }
                  }
                } else {
                  // mfh 27 Sep 2016: If the entry of Acol2PRow
                  // corresponding to the current entry of A is NOT populated (has
                  // a flag "invalid" value), then the corresponding row of P is
                  // in P_remote (i.e., it does not live on the calling process).

                  // Remote matrix
                  size_t Il = Teuchos::as<size_t>(Acol2PRowImport(l));
                  for (size_t jj = Irowptr(Il); jj < Irowptr(Il+1); jj++) {
                    LO j = Icolind(jj);
                    LO Acj = PIcol2Accol(j);
                    SC Plj = Ivals(jj);

                    if (ac_status[Acj] == INVALID || ac_status[Acj] < nnz_old) {
    #ifdef HAVE_TPETRA_DEBUG
                      // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                      TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.extent(0)),
                                                 std::runtime_error,
                                                 label << " ERROR, not enough memory allocated for matrix product. Allocated: "  << Accolind.size() << std::endl);
    #endif
                      // New entry
                      ac_status[Acj] = nnz;
                      Accolind(nnz) = Acj;
                      Acvals(nnz) = Rik*Akl*Plj;
                      nnz++;
                    } else {
                      Acvals(ac_status[Acj]) += Rik*Akl*Plj;
                    }
                  }
                }
              }
            }
            // Resize for next pass if needed
            if (nnz + n > nnzAllocated) {
              nnzAllocated *= 2;
              Accolind = u_lno_nnz_view_t((typename u_lno_nnz_view_t::data_type) realloc(Accolind.data(), u_lno_nnz_view_t::shmem_size(nnzAllocated)), nnzAllocated);
              Acvals = u_scalar_view_t((typename u_scalar_view_t::data_type) realloc((void*) Acvals.data(), u_scalar_view_t::shmem_size(nnzAllocated)), nnzAllocated);
            }
            nnz_old = nnz;
          }
          thread_total_nnz(tid) = nnz;
          tl_colind(tid) = Accolind;
          tl_values(tid) = Acvals;
        });
  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix copy from thread local"))));
  #endif

        copy_out_from_thread_memory(thread_total_nnz,tl_colind, tl_values, m, thread_chunk, rowmapAc, entriesAc, valuesAc);

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix Final Sort"))));
  #endif

        // Final sort & set of CRS arrays
        Import_Util::sortCrsEntries(rowmapAc, entriesAc, valuesAc);
        // mfh 27 Sep 2016: This just sets pointers.
        Ac.setAllValues(rowmapAc, entriesAc, valuesAc);

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix ESFC"))));
  #endif

        // Final FillComplete
        //
        // mfh 27 Sep 2016: So-called "expert static fill complete" bypasses
        // Import (from domain Map to column Map) construction (which costs
        // lots of communication) by taking the previously constructed
        // Import object.  We should be able to do this without interfering
        // with the implementation of the local part of sparse matrix-matrix
        // multiply above.
        RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
        labelList->set("Timer Label",label);
        if(!params.is_null()) labelList->set("compute global constants",params->get("compute global constants",true));
        RCP<const Export<LO,GO,NO> > dummyExport;
        Ac.expertStaticFillComplete(Pview.origMatrix->getDomainMap(),
                                    Rview.origMatrix->getRangeMap(),
                                    Acimport,
                                    dummyExport,
                                    labelList);
}



/*********************************************************************************************************/
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class LocalOrdinalViewType>
static inline void mult_R_A_P_reuse_LowThreadGustavsonKernel(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Rview,
                                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Aview,
                                                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Pview,
                                                                 const LocalOrdinalViewType & Acol2PRow,
                                                                 const LocalOrdinalViewType & Acol2PRowImport,
                                                                 const LocalOrdinalViewType & Pcol2Accol,
                                                                 const LocalOrdinalViewType & PIcol2Accol,
                                                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>& Ac,
                                                                 Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode> > Acimport,
                                                                 const std::string& label,
                                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {

        using Teuchos::RCP;
        using Tpetra::MatrixMatrix::UnmanagedView;
  #ifdef HAVE_TPETRA_MMM_TIMINGS
        std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
        using Teuchos::TimeMonitor;
        using Teuchos::rcp;
        RCP<TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Reuse LTGCore"))));
  #endif

        typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;
        typedef Scalar        SC;
        typedef LocalOrdinal  LO;
        typedef GlobalOrdinal GO;
        typedef Node          NO;
        typedef Map<LO,GO,NO> map_type;
        typedef typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_matrix_device_type KCRS;
        typedef typename KCRS::StaticCrsGraphType graph_t;
        typedef typename graph_t::row_map_type::const_type c_lno_view_t;
        typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
        typedef typename KCRS::values_type::non_const_type scalar_view_t;
        typedef typename KCRS::device_type device_t;
        typedef typename device_t::execution_space execution_space;
        typedef Kokkos::RangePolicy<execution_space, size_t> range_type;

        const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
        const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();
        const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

        // Sizes
        RCP<const map_type> Accolmap = Ac.getColMap();
        size_t m = Rview.origMatrix->getLocalNumRows();
        size_t n = Accolmap->getLocalNumElements();

        // Get raw Kokkos matrices, and the raw CSR views
        const KCRS & Rmat = Rview.origMatrix->getLocalMatrixDevice();
        const KCRS & Amat = Aview.origMatrix->getLocalMatrixDevice();
        const KCRS & Pmat = Pview.origMatrix->getLocalMatrixDevice();
        const KCRS & Cmat = Ac.getLocalMatrixDevice();

        c_lno_view_t Rrowptr = Rmat.graph.row_map, Arowptr = Amat.graph.row_map, Prowptr = Pmat.graph.row_map, Crowptr = Cmat.graph.row_map, Irowptr;
        const lno_nnz_view_t Rcolind = Rmat.graph.entries, Acolind = Amat.graph.entries, Pcolind = Pmat.graph.entries, Ccolind = Cmat.graph.entries;
        lno_nnz_view_t Icolind;
        const scalar_view_t Rvals = Rmat.values, Avals = Amat.values, Pvals = Pmat.values;
        scalar_view_t Cvals = Cmat.values;
        scalar_view_t Ivals;

        if (!Pview.importMatrix.is_null())
        {
          const KCRS& Imat = Pview.importMatrix->getLocalMatrixDevice();
          Irowptr = Imat.graph.row_map;
          Icolind = Imat.graph.entries;
          Ivals = Imat.values;
        }

        // Get my node / thread info (right from openmp or parameter list)
        size_t thread_max =  Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency();
        if(!params.is_null()) {
          if(params->isParameter("openmp: ltg thread max"))
            thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));
        }

        double thread_chunk = (double)(m) / thread_max;

        // For each row of R
        Kokkos::parallel_for("MMM::RAP::Reuse::LTG::ThreadLocal",range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid)
        {
          // Thread coordination stuff
          size_t my_thread_start = tid * thread_chunk;
          size_t my_thread_stop  = tid == thread_max-1 ? m : (tid+1)*thread_chunk;

          std::vector<size_t> c_status(n, INVALID);

          //count entries as they are added to Ac
          size_t OLD_ip = 0, CSR_ip = 0;
          // bmk: loop over the rows of R which are assigned to thread tid
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

            // mfh 27 Sep 2016: For each entry of R in the current row of R
            for (size_t kk = Rrowptr(i); kk < Rrowptr(i+1); kk++) {
              LO k  = Rcolind(kk); // local column index of current entry of R
              const SC Rik = Rvals(kk);   // value of current entry of R
              if (Rik == SC_ZERO)
                continue; // skip explicitly stored zero values in R
              // For each entry of A in the current row of A
              for (size_t ll = Arowptr(k); ll < Arowptr(k+1); ll++) {
                LO l = Acolind(ll); // local column index of current entry of A
                const SC Akl = Avals(ll);   // value of current entry of A
                if (Akl == SC_ZERO)
                  continue; // skip explicitly stored zero values in A

                if (Acol2PRow[l] != LO_INVALID) {
                  // mfh 27 Sep 2016: If the entry of Acol2PRow
                  // corresponding to the current entry of A is populated, then
                  // the corresponding row of P is in P_local (i.e., it lives on
                  // the calling process).

                  // Local matrix
                  size_t Pl = Teuchos::as<size_t>(Acol2PRow(l));

                  // mfh 27 Sep 2016: Go through all entries in that row of P_local.
                  for (size_t jj = Prowptr(Pl); jj < Prowptr(Pl+1); jj++) {
                    LO j = Pcolind(jj);
                    LO Cij = Pcol2Accol(j);
                    SC Plj = Pvals(jj);
                    TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
                                         std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
                                         "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");


                    Cvals(c_status[Cij]) += Rik*Akl*Plj;
                  }
                } else {
                  // mfh 27 Sep 2016: If the entry of Acol2PRow
                  // corresponding to the current entry of A is NOT populated (has
                  // a flag "invalid" value), then the corresponding row of P is
                  // in P_remote (i.e., it does not live on the calling process).

                  // Remote matrix
                  size_t Il = Teuchos::as<size_t>(Acol2PRowImport(l));
                  for (size_t jj = Irowptr(Il); jj < Irowptr(Il+1); jj++) {
                    LO j = Icolind(jj);
                    LO Cij = PIcol2Accol(j);
                    SC Plj = Ivals(jj);
                    TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
                                               std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
                                               "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

                    Cvals(c_status[Cij]) += Rik*Akl*Plj;
                  }
                }
              }
            }
          }
        });
        // NOTE: No copy out or "set" of data is needed here, since we're working directly with Kokkos::Views
}



#endif


}//ExtraKernels
}//MatrixMatrix
}//Tpetra


#endif
