// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
 * \file   Amesos2_EpetraCrsMatrix_MatrixAdapter_decl.hpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Tue Jun 14 17:17:00 MDT 2011
 * 
 * \brief Specialization of the ConcreteMatrixAdapter for
 * Epetra_CrsMatrix.  Inherits all its functionality from the
 * Epetra_RowMatrix specialization of \c AbstractConcreteMatrixAdapter.
 */

#ifndef AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP
#define AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP

#include "Amesos2_config.h"

#include <Epetra_CrsMatrix.h>
#ifdef HAVE_AMESOS2_EPETRAEXT
#include <EpetraExt_Reindex_CrsMatrix.h>
#endif
#include <Epetra_Comm.h>
#include <Epetra_SerialComm.h>
#ifdef HAVE_MPI
#  include <Epetra_MpiComm.h>
#endif

#include "Amesos2_EpetraRowMatrix_AbstractMatrixAdapter_decl.hpp"
#include "Amesos2_MatrixAdapter_decl.hpp"

namespace Amesos2 {

  // template <class M, class D> class AbstractConcreteMatrixAdapter;

  /**
   * \brief MatrixAdapter definitions for Epetra_CrsMatrix objects.
   *
   * Defines only the get_impl() method, which returns an instance of
   * a Amesos2::MatrixAdapter whose underlying matrix has the given
   * distribution based on the Tpetra::Map.
   *
   * All other significant functionality is inherited from this
   * class's superclass.
   *
   * \ingroup amesos2_matrix_adapters
   */
  template <>
  class ConcreteMatrixAdapter< Epetra_CrsMatrix >
    : public AbstractConcreteMatrixAdapter< Epetra_RowMatrix, Epetra_CrsMatrix >
  {
    // Give our matrix adapter class access to our private
    // implementation functions
    friend class MatrixAdapter< Epetra_RowMatrix >;
  public:
    typedef Epetra_CrsMatrix                               matrix_t;
  private:
    typedef AbstractConcreteMatrixAdapter<Epetra_RowMatrix,
                                          Epetra_CrsMatrix> super_t;
  public:
    // 'import' superclass types
    typedef super_t::scalar_t                                 scalar_t;
    typedef super_t::local_ordinal_t                   local_ordinal_t;
    typedef super_t::global_ordinal_t                 global_ordinal_t;
    typedef super_t::node_t                                     node_t;
    typedef super_t::global_size_t                       global_size_t;

    typedef Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> map_t;
    typedef ConcreteMatrixAdapter<matrix_t>                       type;
    
    ConcreteMatrixAdapter(RCP<matrix_t> m);
    
    // get & reindex
    RCP<const MatrixAdapter<matrix_t> > get_impl(const Teuchos::Ptr<const map_t> map, EDistribution distribution = ROOTED) const;
    RCP<const MatrixAdapter<matrix_t> > reindex_impl(Teuchos::RCP<const map_t> &contigRowMap,
                                                     Teuchos::RCP<const map_t> &contigColMap) const;

    // gather
    template<typename KV_S, typename KV_GO, typename KV_GS>
    local_ordinal_t gather_impl(KV_S& nzvals, KV_GO& indices, KV_GS& pointers, bool column_major, EPhase current_phase) const
    {
      local_ordinal_t ret = -1;
      //TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "ConcreteMatrixAdapter<Epetra_CrsMatrix> has not been implemented gather_impl.");
      {
#ifdef HAVE_MPI
        auto Comm = this->mat_->Comm();
        printf( " %d/%d\n",Comm.MyPID(),Comm.NumProc());
        #if 0
        Comm.PrintInfo(std::cout);
        Teuchos::RCP<const Epetra_MpiComm>
          mpiEpetraComm = Teuchos::rcp_dynamic_cast<const Epetra_MpiComm>(Teuchos::rcpFromRef(Comm));
        if( !mpiEpetraComm.get() ) {
          printf( " ERROR\n" ); fflush(stdout);
          return ret;
        }
        RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
          rawMpiComm = Teuchos::opaqueWrapper(mpiEpetraComm->Comm());
        RCP<Teuchos::Comm<int>> comm = Teuchos::createMpiComm<int>(rawMpiComm);
        #else
        auto comm = Util::to_teuchos_comm(Teuchos::rcpFromRef(Comm));
        #endif
        auto nRanks = comm->getSize();
        auto myRank = comm->getRank();
#else
        RCP<Teuchos::Comm<int>> comm = Teuchos::createSerialComm<int>();
        int nRanks = 1;
        int myRank = 0;
#endif // HAVE_MPI

        int *lclRowptr = nullptr;
        int *lclColind = nullptr;
        double *lclNzvals = nullptr;
        this->mat_->ExtractCrsDataPointers(lclRowptr, lclColind, lclNzvals);

        // workspace for column major
        KV_GS pointers_t;
        KV_S nzvals_t;
        KV_GO indices_t;

        // gether rowptr
        int nRows = this->mat_->NumGlobalRows();
        int myNRows = this->mat_->NumMyRows();
        int myNnz = this->mat_->NumMyNonzeros();
        KV_GS recvCounts("recvCounts",nRanks);
        KV_GS recvDispls("recvDispls",nRanks+1);
        {
#ifdef HAVE_AMESOS2_TIMERS
          Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(rowptr)");
          Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
          Teuchos::gather<int, int> (&myNRows, 1, recvCounts.data(), 1, 0, *comm);
          if (myRank == 0) {
            Kokkos::resize(recvDispls, nRanks+1);
            recvDispls(0) = 0;
            for (int p = 1; p <= nRanks; p++) {
              recvDispls(p) = recvDispls(p-1) + recvCounts(p-1);
            }
            fflush(stdout);
            if (recvDispls(nRanks) != nRows) {
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Amesos2_TpetraCrsMatrix_MatrixAdapter::gather_impl : mismatch between gathered(local nrows) and global nrows.");
            }
          } else {
            for (int p = 0; p < nRanks; p++) {
              recvCounts(p) = 0;
              recvDispls(p) = 0;
            }
            recvDispls(nRanks) = 0;
          }
          if (myRank == 0 && column_major) {
            Kokkos::resize(pointers_t, nRows+1);
          } else if (myRank != 0) {
            Kokkos::resize(pointers_t, 2);
          }
          int *pointers_ = (myRank != 0 || column_major ? pointers_t.data() : pointers.data());
          Teuchos::gatherv<int, local_ordinal_t> (&lclRowptr[1], myNRows, &pointers_[1], 
                                                  recvCounts.data(), recvDispls.data(),
                                                  0, *comm);
          if (myRank == 0) {
            // shift to global pointers
            pointers_[0] = 0;
            recvCounts(0) = pointers_[recvDispls(1)];
            local_ordinal_t displs = recvCounts(0);
            for (int p = 1; p < nRanks; p++) {
              // save recvCounts from pth MPI
              recvCounts(p) = pointers_[recvDispls(p+1)];
              // shift pointers for pth MPI to global
              for (int i = 1+recvDispls(p); i <= recvDispls(p+1); i++) {
                pointers_[i] += displs;
              }
              displs += recvCounts(p);
            }
            ret = pointers_[nRows];
          }
        }
        {
#ifdef HAVE_AMESOS2_TIMERS
          Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(colind)");
          Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
          // gather colinds & nzvals
          if (myRank == 0) {
            recvDispls(0) = 0;
            for (int p = 0; p < nRanks; p++) {
              recvDispls(p+1) = recvDispls(p) + recvCounts(p);
            }
          }
          if (column_major) {
            Kokkos::resize(indices_t, indices.extent(0));
            Teuchos::gatherv<int, local_ordinal_t> (lclColind, myNnz, indices_t.data(), 
                                                    recvCounts.data(), recvDispls.data(),
                                                    0, *comm);
          } else {
            Teuchos::gatherv<int, local_ordinal_t> (lclColind, myNnz, indices.data(), 
                                                    recvCounts.data(), recvDispls.data(),
                                                    0, *comm);
          }
        }
        {
#ifdef HAVE_AMESOS2_TIMERS
          Teuchos::RCP< Teuchos::Time > gatherTime = Teuchos::TimeMonitor::getNewCounter ("Amesos2::gather(nzvals)");
          Teuchos::TimeMonitor GatherTimer(*gatherTime);
#endif
          if (column_major) {
            Kokkos::resize(nzvals_t, nzvals.extent(0));
            Teuchos::gatherv<int, scalar_t> (lclNzvals, myNnz, nzvals_t.data(), 
                                             recvCounts.data(), recvDispls.data(),
                                             0, *comm);
          } else {
            Teuchos::gatherv<int, scalar_t> (lclNzvals, myNnz, nzvals.data(), 
                                             recvCounts.data(), recvDispls.data(),
                                             0, *comm);
          }
        }
      }
      return ret;
    }

    //! Print a description of this adapter to the given output stream
    void
    describe (Teuchos::FancyOStream& os,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;
#ifdef HAVE_AMESOS2_EPETRAEXT
  private:
    mutable RCP<EpetraExt::CrsMatrix_Reindex> StdIndex_;
    mutable RCP<Epetra_CrsMatrix> ContigMat_;
#endif
  };

} // end namespace Amesos2

#endif // AMESOS2_EPETRACRSMATRIX_MATRIXADAPTER_DECL_HPP
