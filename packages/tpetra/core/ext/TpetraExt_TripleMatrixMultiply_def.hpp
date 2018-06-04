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
#ifndef TPETRA_TRIPLEMATRIXMULTIPLY_DEF_HPP
#define TPETRA_TRIPLEMATRIXMULTIPLY_DEF_HPP

#include "TpetraExt_MatrixMatrix_def.hpp"
#include "TpetraExt_MatrixMatrix_ExtraKernels_decl.hpp" //for UnmanagedView
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Array.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "TpetraExt_MMHelpers_def.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import_Util.hpp"
#include "Tpetra_Import_Util2.hpp"
#include <algorithm>
#include <cmath>
#include "Teuchos_FancyOStream.hpp"
// #include "KokkosSparse_spgemm.hpp"


/*! \file TpetraExt_TripleMatrixMultiply_def.hpp

  The implementations for the members of class Tpetra::TripleMatrixMultiply and related non-member constructors.
*/

namespace Tpetra {

  namespace TripleMatrixMultiply{

    //
    // This method forms the matrix-matrix product Ac = op(R) * op(A) * op(P), where
    // op(A) == A   if transposeA is false,
    // op(A) == A^T if transposeA is true,
    // and similarly for op(R) and op(P).
    //
    template <class Scalar,
              class LocalOrdinal,
              class GlobalOrdinal,
              class Node>
    void MultiplyRAP(
                     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& R,
                     bool transposeR,
                     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                     bool transposeA,
                     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& P,
                     bool transposeP,
                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                     bool call_FillComplete_on_result,
                     const std::string& label,
                     const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      using Teuchos::null;
      using Teuchos::RCP;
      typedef Scalar                            SC;
      typedef LocalOrdinal                      LO;
      typedef GlobalOrdinal                     GO;
      typedef Node                              NO;
      typedef CrsMatrix<SC,LO,GO,NO>            crs_matrix_type;
      typedef Import<LO,GO,NO>                  import_type;
      typedef Export<LO,GO,NO>                  export_type;
      typedef CrsMatrixStruct<SC,LO,GO,NO>      crs_matrix_struct_type;
      typedef Map<LO,GO,NO>                     map_type;
      typedef RowMatrixTransposer<SC,LO,GO,NO>  transposer_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP All Setup"))));
#endif

      const std::string prefix = "TpetraExt::TripleMatrixMultiply::MultiplyRAP(): ";

      // TEUCHOS_FUNC_TIME_MONITOR_DIFF("My Matrix Mult", mmm_multiply);

      // The input matrices R, A and P must both be fillComplete.
      TEUCHOS_TEST_FOR_EXCEPTION(!R.isFillComplete(), std::runtime_error, prefix << "Matrix R is not fill complete.");
      TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error, prefix << "Matrix A is not fill complete.");
      TEUCHOS_TEST_FOR_EXCEPTION(!P.isFillComplete(), std::runtime_error, prefix << "Matrix P is not fill complete.");

      // If transposeA is true, then Rprime will be the transpose of R
      // (computed explicitly via RowMatrixTransposer).  Otherwise, Rprime
      // will just be a pointer to R.
      RCP<const crs_matrix_type> Rprime = null;
      // If transposeA is true, then Aprime will be the transpose of A
      // (computed explicitly via RowMatrixTransposer).  Otherwise, Aprime
      // will just be a pointer to A.
      RCP<const crs_matrix_type> Aprime = null;
      // If transposeB is true, then Pprime will be the transpose of P
      // (computed explicitly via RowMatrixTransposer).  Otherwise, Pprime
      // will just be a pointer to P.
      RCP<const crs_matrix_type> Pprime = null;

      // Is this a "clean" matrix?
      //
      // mfh 27 Sep 2016: Historically, if Epetra_CrsMatrix was neither
      // locally nor globally indexed, then it was empty.  I don't like
      // this, because the most straightforward implementation presumes
      // lazy allocation of indices.  However, historical precedent
      // demands that we keep around this predicate as a way to test
      // whether the matrix is empty.
      const bool newFlag = !Ac.getGraph()->isLocallyIndexed() && !Ac.getGraph()->isGloballyIndexed();

      if (transposeR && &R != &P) {
        transposer_type transposer(rcpFromRef (R));
        Rprime = transposer.createTranspose();
      } else {
        Rprime = rcpFromRef(R);
      }

      if (transposeA) {
        transposer_type transposer(rcpFromRef (A));
        Aprime = transposer.createTranspose();
      } else {
        Aprime = rcpFromRef(A);
      }

      if (transposeP) {
        transposer_type transposer(rcpFromRef (P));
        Pprime = transposer.createTranspose();
      } else {
        Pprime = rcpFromRef(P);
      }

      // Check size compatibility
      global_size_t numRCols = R.getDomainMap()->getGlobalNumElements();
      global_size_t numACols = A.getDomainMap()->getGlobalNumElements();
      global_size_t numPCols = P.getDomainMap()->getGlobalNumElements();
      global_size_t Rleft    = transposeR ? numRCols             : R.getGlobalNumRows();
      global_size_t Rright   = transposeR ? R.getGlobalNumRows() : numRCols;
      global_size_t Aleft    = transposeA ? numACols             : A.getGlobalNumRows();
      global_size_t Aright   = transposeA ? A.getGlobalNumRows() : numACols;
      global_size_t Pleft    = transposeP ? numPCols             : P.getGlobalNumRows();
      global_size_t Pright   = transposeP ? P.getGlobalNumRows() : numPCols;
      TEUCHOS_TEST_FOR_EXCEPTION(Rright != Aleft, std::runtime_error,
                                 prefix << "ERROR, inner dimensions of op(R) and op(A) "
                                 "must match for matrix-matrix product. op(R) is "
                                 << Rleft << "x" << Rright << ", op(A) is "<< Aleft << "x" << Aright);

      TEUCHOS_TEST_FOR_EXCEPTION(Aright != Pleft, std::runtime_error,
                                 prefix << "ERROR, inner dimensions of op(A) and op(P) "
                                 "must match for matrix-matrix product. op(A) is "
                                 << Aleft << "x" << Aright << ", op(P) is "<< Pleft << "x" << Pright);

      // The result matrix Ac must at least have a row-map that reflects the correct
      // row-size. Don't check the number of columns because rectangular matrices
      // which were constructed with only one map can still end up having the
      // correct capacity and dimensions when filled.
      TEUCHOS_TEST_FOR_EXCEPTION(Rleft > Ac.getGlobalNumRows(), std::runtime_error,
                                 prefix << "ERROR, dimensions of result Ac must "
                                 "match dimensions of op(R) * op(A) * op(P). Ac has " << Ac.getGlobalNumRows()
                                 << " rows, should have at least " << Rleft << std::endl);

      // It doesn't matter whether Ac is already Filled or not. If it is already
      // Filled, it must have space allocated for the positions that will be
      // referenced in forming Ac = op(R)*op(A)*op(P). If it doesn't have enough space,
      // we'll error out later when trying to store result values.

      // CGB: However, matrix must be in active-fill
      TEUCHOS_TEST_FOR_EXCEPT( Ac.isFillActive() == false );

      // We're going to need to import remotely-owned sections of P if
      // more than one processor is performing this run, depending on the scenario.
      int numProcs = P.getComm()->getSize();

      // Declare a couple of structs that will be used to hold views of the data
      // of R, A and P, to be used for fast access during the matrix-multiplication.
      crs_matrix_struct_type Rview;
      crs_matrix_struct_type Aview;
      crs_matrix_struct_type Pview;

      RCP<const map_type> targetMap_R = Rprime->getRowMap();
      RCP<const map_type> targetMap_A = Aprime->getRowMap();
      RCP<const map_type> targetMap_P = Pprime->getRowMap();

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP All I&X"))));
#endif

      // Now import any needed remote rows and populate the Aview struct
      // NOTE: We assert that an import isn't needed --- since we do the transpose
      // above to handle that.
      RCP<const import_type> dummyImporter;

      if (!(transposeR && &R == &P))
        MMdetails::import_and_extract_views(*Rprime, targetMap_R, Rview, dummyImporter, true, label, params);

      MMdetails::import_and_extract_views(*Aprime, targetMap_A, Aview, dummyImporter, true, label, params);

      // We will also need local access to all rows of P that correspond to the
      // column-map of op(A).
      if (numProcs > 1)
        targetMap_P = Aprime->getColMap();

      // Import any needed remote rows and populate the Pview struct.
      MMdetails::import_and_extract_views(*Pprime, targetMap_P, Pview, Aprime->getGraph()->getImporter(), false, label, params);


      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Actemp;

      bool needs_final_export = !Pprime->getGraph()->getImporter().is_null();
      if (needs_final_export)
        Actemp = rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Pprime->getColMap(),0));
      else
        Actemp = rcp(&Ac,false);// don't allow deallocation

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP All Multiply"))));
#endif

      // Call the appropriate method to perform the actual multiplication.
      if (call_FillComplete_on_result && newFlag) {
        if (transposeR && &R == &P)
          MMdetails::mult_PT_A_P_newmatrix(Aview, Pview, *Actemp, label, params);
        else
          MMdetails::mult_R_A_P_newmatrix(Rview, Aview, Pview, *Actemp, label, params);
      } else if (call_FillComplete_on_result) {
        // MMdetails::mult_A__reuse(Aview, Bview, C, label,params);
        // Not implemented
        if (transposeR && &R == &P)
          MMdetails::mult_PT_A_P_newmatrix(Aview, Pview, *Actemp, label, params);
        else
          MMdetails::mult_R_A_P_newmatrix(Rview, Aview, Pview, *Actemp, label, params);
      } else {
        // mfh 27 Sep 2016: Is this the "slow" case?  This
        // "CrsWrapper_CrsMatrix" thing could perhaps be made to support
        // thread-parallel inserts, but that may take some effort.
        // CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsmat(Ac);

        //     MMdetails::mult_A_B(Aview, Bview, crsmat, label,params);

        // #ifdef HAVE_TPETRA_MMM_TIMINGS
        //     MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP All FillComplete"))));
        // #endif
        //     if (call_FillComplete_on_result) {
        //       // We'll call FillComplete on the C matrix before we exit, and give it a
        //       // domain-map and a range-map.
        //       // The domain-map will be the domain-map of B, unless
        //       // op(B)==transpose(B), in which case the range-map of B will be used.
        //       // The range-map will be the range-map of A, unless op(A)==transpose(A),
        //       // in which case the domain-map of A will be used.
        //       if (!C.isFillComplete())
        //         C.fillComplete(Bprime->getDomainMap(), Aprime->getRangeMap());
        //     }
        // Not implemented
        if (transposeR && &R == &P)
          MMdetails::mult_PT_A_P_newmatrix(Aview, Pview, *Actemp, label, params);
        else
          MMdetails::mult_R_A_P_newmatrix(Rview, Aview, Pview, *Actemp, label, params);
      }

      if (needs_final_export) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP exportAndFillComplete"))));
#endif
        Teuchos::ParameterList labelList;
        labelList.set("Timer Label", label);
        RCP<crs_matrix_type> Acprime = rcpFromRef(Ac);
        if(!params.is_null()) labelList.set("compute global constants",params->get("compute global constants",true));
        export_type exporter = export_type(*Pprime->getGraph()->getImporter());
        Actemp->exportAndFillComplete(Acprime,
                                      exporter,
                                      Acprime->getDomainMap(),
                                      Acprime->getRangeMap(),
                                      rcp(&labelList,false));

      }
#ifdef HAVE_TPETRA_MMM_STATISTICS
      printMultiplicationStatistics(Actemp->getGraph()->getExporter(), label+std::string(" RAP MMM"));
#endif

    }


  } //End namespace TripleMatrixMultiply

  namespace MMdetails{


    template<class CrsMatrixType>
    size_t Ac_estimate_nnz(CrsMatrixType & A, CrsMatrixType &P){
      size_t nnzPerRowA = 100, Pcols = 100;
      if (A.getNodeNumEntries() > 0)
        nnzPerRowA = (A.getNodeNumRows() > 0)?  A.getNodeNumEntries()/A.getNodeNumRows() : 9;
      if (P.getNodeNumEntries() > 0)
        Pcols = (P.getNodeNumCols() > 0) ? P.getNodeNumCols() : 100;
      return (size_t)(Pcols*nnzPerRowA + 5*nnzPerRowA + 300);
    }



    // Kernel method for computing the local portion of Ac = R*A*P
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void mult_R_A_P_newmatrix(
                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                              const std::string& label,
                              const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;

      //typedef Scalar            SC; // unused
      typedef LocalOrdinal      LO;
      typedef GlobalOrdinal     GO;
      typedef Node              NO;

      typedef Import<LO,GO,NO>  import_type;
      typedef Map<LO,GO,NO>     map_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP M5 Cmap")))));
#endif
      LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

      // Build the final importer / column map, hash table lookups for Ac
      RCP<const import_type> Acimport;
      RCP<const map_type>    Accolmap;
      RCP<const import_type> Pimport = Pview.origMatrix->getGraph()->getImporter();
      RCP<const import_type> Iimport = Pview.importMatrix.is_null() ?
        Teuchos::null : Pview.importMatrix->getGraph()->getImporter();

      // mfh 27 Sep 2016: Pcol2Accol is a table that maps from local column
      // indices of P, to local column indices of Ac.  (P and Ac have the
      // same number of columns.)  The kernel uses this, instead of
      // copying the entire input matrix P and converting its column
      // indices to those of Ac.
      Array<LO> Pcol2Accol(Pview.colMap->getNodeNumElements()), PIcol2Accol;

      if (Pview.importMatrix.is_null()) {
        // mfh 27 Sep 2016: P has no "remotes," so P and Ac have the same column Map.
        Acimport = Pimport;
        Accolmap = Pview.colMap;
        // Pcol2Accol is trivial
        for (size_t i = 0; i < Pview.colMap->getNodeNumElements(); i++)
          Pcol2Accol[i] = Teuchos::as<LO>(i);

      } else {
        // mfh 27 Sep 2016: P has "remotes," so we need to build the
        // column Map of Ac, as well as Ac's Import object (from its domain
        // Map to its column Map).  Ac's column Map is the union of the
        // column Maps of (the local part of) P, and the "remote" part of
        // P.  Ditto for the Import.  We have optimized this "setUnion"
        // operation on Import objects and Maps.

        // Choose the right variant of setUnion
        if (!Pimport.is_null() && !Iimport.is_null())
          Acimport = Pimport->setUnion(*Iimport);

        else if (!Pimport.is_null() && Iimport.is_null())
          Acimport = Pimport->setUnion();

        else if (Pimport.is_null() && !Iimport.is_null())
          Acimport = Iimport->setUnion();

        else
          throw std::runtime_error("TpetraExt::MMM status of matrix importers is nonsensical");

        Accolmap = Acimport->getTargetMap();

        // FIXME (mfh 27 Sep 2016) This error check requires an all-reduce
        // in general.  We should get rid of it in order to reduce
        // communication costs of sparse matrix-matrix multiply.
        TEUCHOS_TEST_FOR_EXCEPTION(!Acimport->getSourceMap()->isSameAs(*Pview.origMatrix->getDomainMap()),
                                   std::runtime_error, "Tpetra::MMM: Import setUnion messed with the DomainMap in an unfortunate way");

        // NOTE: This is not efficient and should be folded into setUnion
        //
        // mfh 27 Sep 2016: What the above comment means, is that the
        // setUnion operation on Import objects could also compute these
        // local index - to - local index look-up tables.
        PIcol2Accol.resize(Pview.importMatrix->getColMap()->getNodeNumElements());
        ArrayView<const GO> Pgid = Pview.origMatrix->getColMap()->getNodeElementList();
        ArrayView<const GO> Igid = Pview.importMatrix->getColMap()->getNodeElementList();

        for (size_t i = 0; i < Pview.origMatrix->getColMap()->getNodeNumElements(); i++)
          Pcol2Accol[i] = Accolmap->getLocalElement(Pgid[i]);
        for (size_t i = 0; i < Pview.importMatrix->getColMap()->getNodeNumElements(); i++)
          PIcol2Accol[i] = Accolmap->getLocalElement(Igid[i]);
      }

      // Replace the column map
      //
      // mfh 27 Sep 2016: We do this because Ac was originally created
      // without a column Map.  Now we have its column Map.
      Ac.replaceColMap(Accolmap);

      // mfh 27 Sep 2016: Construct tables that map from local column
      // indices of A, to local row indices of either P_local (the locally
      // owned part of P), or P_remote (the "imported" remote part of P).
      //
      // For column index Aik in row i of A, if the corresponding row of P
      // exists in the local part of P ("orig") (which I'll call P_local),
      // then Acol2PRow[Aik] is the local index of that row of P.
      // Otherwise, Acol2PRow[Aik] is "invalid" (a flag value).
      //
      // For column index Aik in row i of A, if the corresponding row of P
      // exists in the remote part of P ("Import") (which I'll call
      // P_remote), then Acol2PRowImport[Aik] is the local index of
      // that row of B.  Otherwise, Acol2PRowImport[Aik] is "invalid"
      // (a flag value).

      // Run through all the hash table lookups once and for all
      Array<LO> Acol2PRow  (Aview.colMap->getNodeNumElements(), LO_INVALID);
      Array<LO> Acol2PRowImport(Aview.colMap->getNodeNumElements(), LO_INVALID);

      for (LO i = Aview.colMap->getMinLocalIndex(); i <= Aview.colMap->getMaxLocalIndex(); i++) {
        LO P_LID = Pview.origMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
        if (P_LID != LO_INVALID) {
          Acol2PRow[i] = P_LID;

        } else {
          LO I_LID = Pview.importMatrix->getRowMap()->getLocalElement(Aview.colMap->getGlobalElement(i));
          Acol2PRowImport[i] = I_LID;
        }
      }

      // Call the actual kernel.  We'll rely on partial template specialization to call the correct one ---
      // Either the straight-up Tpetra code (SerialNode) or the KokkosKernels one (other NGP node types)
      KernelWrappers3MMM_Specialized<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
        mult_R_A_P_newmatrix_kernel_wrapper(Rview, Aview, Pview,
                                            Acol2PRow, Acol2PRowImport, Pcol2Accol, PIcol2Accol,
                                            Ac, Acimport, label, params);
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
      {
        lno_view_t thread_nnz_count("thread_nnz_counts", thread_max);
        for(size_t i = 0; i < thread_max; i++)
          thread_nnz_count(i) = Inrowptr(i)(Inrowptr(i).dimension(0) - 1);
        Tpetra::Details::computeOffsetsFromCounts(thread_start_nnz, thread_nnz_count);
      }
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

    /*********************************************************************************************************/
    // RAP NewMatrix Kernel wrappers (Default non-threaded version)
    // Computes R * A * P -> Ac using classic Gustavson approach
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void KernelWrappers3MMM_Specialized<Scalar,LocalOrdinal,GlobalOrdinal,Node>::mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                                                                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                                                                                     const Teuchos::Array<LocalOrdinal> & Acol2PRow,
                                                                                                                     const Teuchos::Array<LocalOrdinal> & Acol2PRowImport,
                                                                                                                     const Teuchos::Array<LocalOrdinal> & Pcol2Accol,
                                                                                                                     const Teuchos::Array<LocalOrdinal> & PIcol2Accol,
                                                                                                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                                                                                     Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                                                                                     const std::string& label,
                                                                                                                     const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix SerialCore"))));
#endif

      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;

      typedef Scalar            SC;
      typedef LocalOrdinal      LO;
      typedef GlobalOrdinal     GO;
      typedef Node              NO;
      typedef Map<LO,GO,NO>     map_type;
      const size_t ST_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
      const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
      const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

      // Sizes
      RCP<const map_type> Accolmap = Ac.getColMap();
      size_t m = Rview.origMatrix->getNodeNumRows();
      size_t n = Accolmap->getNodeNumElements();

      // Get Data Pointers
      ArrayRCP<const size_t> Rrowptr_RCP, Arowptr_RCP, Prowptr_RCP, Irowptr_RCP;
      ArrayRCP<size_t> Acrowptr_RCP;
      ArrayRCP<const LO> Rcolind_RCP, Acolind_RCP, Pcolind_RCP, Icolind_RCP;
      ArrayRCP<LO> Accolind_RCP;
      ArrayRCP<const Scalar> Rvals_RCP, Avals_RCP, Pvals_RCP, Ivals_RCP;
      ArrayRCP<SC> Acvals_RCP;

      // mfh 27 Sep 2016: "getAllValues" just gets the three CSR arrays
      // out of the CrsMatrix.  This code computes R * A * (P_local +
      // P_remote), where P_local contains the locally owned rows of P,
      // and P_remote the (previously Import'ed) remote rows of P.

      Rview.origMatrix->getAllValues(Rrowptr_RCP, Rcolind_RCP, Rvals_RCP);
      Aview.origMatrix->getAllValues(Arowptr_RCP, Acolind_RCP, Avals_RCP);
      Pview.origMatrix->getAllValues(Prowptr_RCP, Pcolind_RCP, Pvals_RCP);

      if (!Pview.importMatrix.is_null())
        Pview.importMatrix->getAllValues(Irowptr_RCP, Icolind_RCP, Ivals_RCP);

      // mfh 27 Sep 2016: Remark below "For efficiency" refers to an issue
      // where Teuchos::ArrayRCP::operator[] may be slower than
      // Teuchos::ArrayView::operator[].

      // For efficiency
      ArrayView<const size_t>   Rrowptr, Arowptr, Prowptr, Irowptr;
      ArrayView<const LO>       Rcolind, Acolind, Pcolind, Icolind;
      ArrayView<const SC>       Rvals, Avals, Pvals, Ivals;
      ArrayView<size_t>         Acrowptr;
      ArrayView<LO> Accolind;
      ArrayView<SC> Acvals;
      Rrowptr = Rrowptr_RCP();  Rcolind = Rcolind_RCP();  Rvals = Rvals_RCP();
      Arowptr = Arowptr_RCP();  Acolind = Acolind_RCP();  Avals = Avals_RCP();
      Prowptr = Prowptr_RCP();  Pcolind = Pcolind_RCP();  Pvals = Pvals_RCP();
      if (!Pview.importMatrix.is_null()) {
        Irowptr = Irowptr_RCP(); Icolind = Icolind_RCP(); Ivals = Ivals_RCP();
      }

      // Classic csr assembly (low memory edition)
      //
      // mfh 27 Sep 2016: Ac_estimate_nnz does not promise an upper bound.
      // The method loops over rows of R, and may resize after processing
      // each row.  Chris Siefert says that this reflects experience in
      // ML; for the non-threaded case, ML found it faster to spend less
      // effort on estimation and risk an occasional reallocation.
      size_t nnzAllocated = std::max(Ac_estimate_nnz(*Aview.origMatrix, *Pview.origMatrix), n);
      Acrowptr_RCP.resize(m+1);          Acrowptr = Acrowptr_RCP();
      Accolind_RCP.resize(nnzAllocated); Accolind = Accolind_RCP();
      Acvals_RCP.resize(nnzAllocated);   Acvals   = Acvals_RCP();

      // mfh 27 Sep 2016: The ac_status array is an implementation detail
      // of the local sparse matrix-matrix multiply routine.

      // The status array will contain the index into colind where this entry was last deposited.
      //   ac_status[i] <  nnz - not in the row yet
      //   ac_status[i] >= nnz - this is the entry where you can find the data
      // We start with this filled with INVALID's indicating that there are no entries yet.
      // Sadly, this complicates the code due to the fact that size_t's are unsigned.
      const size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
      Array<size_t> ac_status(n, ST_INVALID);

      // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
      // routine.  The routine computes Ac := R * A * (P_local + P_remote).
      //
      // For column index Aik in row i of A, Acol2PRow[Aik] tells
      // you whether the corresponding row of P belongs to P_local
      // ("orig") or P_remote ("Import").

      // For each row of R
      size_t nnz = 0, nnz_old = 0;
      for (size_t i = 0; i < m; i++) {
        // mfh 27 Sep 2016: m is the number of rows in the input matrix R
        // on the calling process.
        Acrowptr[i] = nnz;

        // mfh 27 Sep 2016: For each entry of R in the current row of R
        for (size_t kk = Rrowptr[i]; kk < Rrowptr[i+1]; kk++) {
          LO k  = Rcolind[kk]; // local column index of current entry of R
          const SC Rik = Rvals[kk];   // value of current entry of R
          if (Rik == SC_ZERO)
            continue; // skip explicitly stored zero values in R
          // For each entry of A in the current row of A
          for (size_t ll = Arowptr[k]; ll < Arowptr[k+1]; ll++) {
            LO l = Acolind[ll]; // local column index of current entry of A
            const SC Akl = Avals[ll];   // value of current entry of A
            if (Akl == SC_ZERO)
              continue; // skip explicitly stored zero values in A


            if (Acol2PRow[l] != LO_INVALID) {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is populated, then
              // the corresponding row of P is in P_local (i.e., it lives on
              // the calling process).

              // Local matrix
              size_t Pl = Teuchos::as<size_t>(Acol2PRow[l]);

              // mfh 27 Sep 2016: Go through all entries in that row of P_local.
              for (size_t jj = Prowptr[Pl]; jj < Prowptr[Pl+1]; jj++) {
                LO j = Pcolind[jj];
                LO Acj = Pcol2Accol[j];
                SC Plj = Pvals[jj];

                if (ac_status[Acj] == INVALID || ac_status[Acj] < nnz_old) {
#ifdef HAVE_TPETRA_DEBUG
                  // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                  TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.size()),
                                             std::runtime_error,
                                             label << " ERROR, not enough memory allocated for matrix product. Allocated: " << Accolind.size() << std::endl);
#endif
                  // New entry
                  ac_status[Acj] = nnz;
                  Accolind[nnz] = Acj;
                  Acvals[nnz] = Rik*Akl*Plj;
                  nnz++;
                } else {
                  Acvals[ac_status[Acj]] += Rik*Akl*Plj;
                }
              }
            } else {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is NOT populated (has
              // a flag "invalid" value), then the corresponding row of P is
              // in P_remote (i.e., it does not live on the calling process).

              // Remote matrix
              size_t Il = Teuchos::as<size_t>(Acol2PRowImport[l]);
              for (size_t jj = Irowptr[Il]; jj < Irowptr[Il+1]; jj++) {
                LO j = Icolind[jj];
                LO Acj = PIcol2Accol[j];
                SC Plj = Ivals[jj];

                if (ac_status[Acj] == INVALID || ac_status[Acj] < nnz_old) {
#ifdef HAVE_TPETRA_DEBUG
                  // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                  TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.size()),
                                             std::runtime_error,
                                             label << " ERROR, not enough memory allocated for matrix product. Allocated: "  << Accolind.size() << std::endl);
#endif
                  // New entry
                  ac_status[Acj] = nnz;
                  Accolind[nnz] = Acj;
                  Acvals[nnz] = Rik*Akl*Plj;
                  nnz++;
                } else {
                  Acvals[ac_status[Acj]] += Rik*Akl*Plj;
                }
              }
            }
          }
        }
        // Resize for next pass if needed
        if (nnz + n > nnzAllocated) {
          nnzAllocated *= 2;
          Accolind_RCP.resize(nnzAllocated); Accolind = Accolind_RCP();
          Acvals_RCP.resize(nnzAllocated);   Acvals   = Acvals_RCP();
        }
        nnz_old = nnz;
      }

      Acrowptr[m] = nnz;

      // Downward resize
      Acvals_RCP.resize(nnz);
      Accolind_RCP.resize(nnz);

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix Final Sort"))));
#endif

      // Final sort & set of CRS arrays
      //
      // TODO (mfh 27 Sep 2016) Will the thread-parallel "local" sparse
      // matrix-matrix multiply routine sort the entries for us?
      Import_Util::sortCrsEntries(Acrowptr_RCP(), Accolind_RCP(), Acvals_RCP());
      // mfh 27 Sep 2016: This just sets pointers.
      Ac.setAllValues(Acrowptr_RCP, Accolind_RCP, Acvals_RCP);

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix ESFC"))));
#endif

      // Final FillComplete
      //
      // mfh 27 Sep 2016: So-called "expert static fill complete" bypasses
      // Import (from domain Map to column Map) construction (which costs
      // lots of communication) by taking the previously constructed
      // Import object.  We should be able to do this without interfering
      // with the implementation of the local part of sparse matrix-matrix
      // multply above.
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

#ifdef HAVE_TPETRA_INST_OPENMP
    /*********************************************************************************************************/
    // RAP NewMatrix Kernel wrappers (Threaded LTG version, the OpenMP specialization)
    // Computes Ac = R * A * P
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal>
    struct KernelWrappers3MMM_Specialized<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode>
    {
      static inline void mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Rview,
                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                               const Teuchos::Array<LocalOrdinal> & Acol2PRow,
                                               const Teuchos::Array<LocalOrdinal> & Acol2PRowImport,
                                               const Teuchos::Array<LocalOrdinal> & Pcol2Accol,
                                               const Teuchos::Array<LocalOrdinal> & PIcol2Accol,
                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
                                               const std::string& label,
                                               const Teuchos::RCP<Teuchos::ParameterList>& params) {
        using Teuchos::RCP;
        using Tpetra::MatrixMatrix::UnmanagedView;
  #ifdef HAVE_TPETRA_MMM_TIMINGS
        std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
        using Teuchos::TimeMonitor;
        Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix SerialCore"))));
  #endif

        typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;
        typedef Scalar        SC;
        typedef LocalOrdinal  LO;
        typedef GlobalOrdinal GO;
        typedef Node          NO;
        typedef Map<LO,GO,NO> map_type;
        typedef typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_matrix_type KCRS;
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
        size_t m = Rview.origMatrix->getNodeNumRows();
        size_t n = Accolmap->getNodeNumElements();

        // Get raw Kokkos matrices, and the raw CSR views
        const KCRS & Rmat = Rview.origMatrix->getLocalMatrix();
        const KCRS & Amat = Aview.origMatrix->getLocalMatrix();
        const KCRS & Pmat = Pview.origMatrix->getLocalMatrix();

        c_lno_view_t Rrowptr = Rmat.graph.row_map, Arowptr = Amat.graph.row_map, Prowptr = Pmat.graph.row_map, Irowptr;
        const lno_nnz_view_t Rcolind = Rmat.graph.entries, Acolind = Amat.graph.entries, Pcolind = Pmat.graph.entries;
        lno_nnz_view_t Icolind;
        const scalar_view_t Rvals = Rmat.values, Avals = Amat.values, Pvals = Pmat.values;
        scalar_view_t Ivals;

        if (!Pview.importMatrix.is_null())
        {
          const KCRS& Imat = Pview.importMatrix->getLocalMatrix();
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
        size_t thread_max =  Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency();
        if(!params.is_null()) {
          if(params->isParameter("openmp: ltg thread max"))
            thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));    
        }

        double thread_chunk = (double)(m) / thread_max;

        // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
        // routine.  The routine computes Ac := R * A * (P_local + P_remote).
        //
        // For column index Aik in row i of A, Acol2PRow[Aik] tells
        // you whether the corresponding row of P belongs to P_local
        // ("orig") or P_remote ("Import").

        // Thread-local memory
        Kokkos::View<u_lno_view_t*> tl_rowptr("top_rowptr", thread_max);
        Kokkos::View<u_lno_nnz_view_t*> tl_colind("top_colind", thread_max);
        Kokkos::View<u_scalar_view_t*> tl_values("top_values", thread_max);

        // For each row of R
        Kokkos::parallel_for("MMM::RAP::NewMatrix::ThreadLocal",range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid)
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
            Acrowptr(i - my_thread_start) = nnz;
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
                  size_t Pl = Teuchos::as<size_t>(Acol2PRow[l]);

                  // mfh 27 Sep 2016: Go through all entries in that row of P_local.
                  for (size_t jj = Prowptr(Pl); jj < Prowptr(Pl+1); jj++) {
                    LO j = Pcolind(jj);
                    LO Acj = Pcol2Accol[j];
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
                  size_t Il = Teuchos::as<size_t>(Acol2PRowImport[l]);
                  for (size_t jj = Irowptr(Il); jj < Irowptr(Il+1); jj++) {
                    LO j = Icolind(jj);
                    LO Acj = PIcol2Accol[j];
                    SC Plj = Ivals(jj);

                    if (ac_status[Acj] == INVALID || ac_status[Acj] < nnz_old) {
    #ifdef HAVE_TPETRA_DEBUG
                      // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                      TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.dimension_0()),
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
              Acvals = u_scalar_view_t((typename u_scalar_view_t::data_type) realloc(Acvals.data(), u_scalar_view_t::shmem_size(nnzAllocated)), nnzAllocated);
            }
            nnz_old = nnz;
          }
          Acrowptr(my_thread_m) = nnz;
          tl_rowptr(tid) = Acrowptr;
          tl_colind(tid) = Accolind;
          tl_values(tid) = Acvals;
        });
  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix copy from thread local"))));
  #endif

        lno_view_t rowmapAc("non_const_lnow_row", m + 1);
        lno_nnz_view_t entriesAc;
        scalar_view_t valuesAc;
        copy_out_from_thread_memory(tl_rowptr, tl_colind, tl_values, m, thread_chunk, rowmapAc, entriesAc, valuesAc);

        for(size_t i=0; i<thread_max; i++) {
          if(tl_rowptr(i).data()) free(tl_rowptr(i).data());
          if(tl_colind(i).data()) free(tl_colind(i).data());
          if(tl_values(i).data()) free(tl_values(i).data());
        }   

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix Final Sort"))));
  #endif

        // Final sort & set of CRS arrays
        Import_Util::sortCrsEntries(rowmapAc, entriesAc, valuesAc);
        // mfh 27 Sep 2016: This just sets pointers.
        Ac.setAllValues(rowmapAc, entriesAc, valuesAc);

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix ESFC"))));
  #endif

        // Final FillComplete
        //
        // mfh 27 Sep 2016: So-called "expert static fill complete" bypasses
        // Import (from domain Map to column Map) construction (which costs
        // lots of communication) by taking the previously constructed
        // Import object.  We should be able to do this without interfering
        // with the implementation of the local part of sparse matrix-matrix
        // multply above.
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
      // PT_A_P NewMatrix Kernel wrappers (LTG OpenMP specialization)
      // Computes P.T * A * P -> Ac
      static inline void mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Aview,
                                                              CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Pview,
                                                              const Teuchos::Array<LocalOrdinal> & Acol2PRow,
                                                              const Teuchos::Array<LocalOrdinal> & Acol2PRowImport,
                                                              const Teuchos::Array<LocalOrdinal> & Pcol2Accol,
                                                              const Teuchos::Array<LocalOrdinal> & PIcol2Accol,
                                                              CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::Compat::KokkosOpenMPWrapperNode>& Ac,
                                                              Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosOpenMPWrapperNode> > Acimport,
                                                              const std::string& label,
                                                              const Teuchos::RCP<Teuchos::ParameterList>& params) {
  #ifdef HAVE_TPETRA_MMM_TIMINGS
        std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
        using Teuchos::TimeMonitor;
        Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix OpenMP"))));
  #endif

        using Teuchos::RCP;
        using Tpetra::MatrixMatrix::UnmanagedView;

        typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;
        typedef Scalar        SC;
        typedef LocalOrdinal  LO;
        typedef GlobalOrdinal GO;
        typedef Node          NO;
        typedef RowMatrixTransposer<SC,LO,GO,NO>  transposer_type;
        typedef typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_matrix_type KCRS;
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

        // number of rows on the process of the fine matrix
        // size_t m = Pview.origMatrix->getNodeNumRows();
        // number of rows on the process of the coarse matrix
        size_t n = Ac.getRowMap()->getNodeNumElements();
        LO maxAccol = Ac.getColMap()->getMaxLocalIndex();

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose"))));
  #endif

        //////////////////////////////////////////////////////////////////////
        transposer_type transposer (Pview.origMatrix,label+std::string("XP: "));
        RCP<Teuchos::ParameterList> transposeParams = Teuchos::rcp(new Teuchos::ParameterList);
        if (!params.is_null())
          transposeParams->set("compute global constants",
                               params->get("compute global constants: temporaries",
                                           false));
        RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ptrans = transposer.createTransposeLocal(transposeParams);

        // Get raw Kokkos matrices, and the raw CSR views
        const KCRS & Rmat = Ptrans->getLocalMatrix();
        const KCRS & Amat = Aview.origMatrix->getLocalMatrix();
        const KCRS & Pmat = Pview.origMatrix->getLocalMatrix();

        c_lno_view_t Rrowptr = Rmat.graph.row_map, Arowptr = Amat.graph.row_map, Prowptr = Pmat.graph.row_map, Irowptr;
        const lno_nnz_view_t Rcolind = Rmat.graph.entries, Acolind = Amat.graph.entries, Pcolind = Pmat.graph.entries;
        lno_nnz_view_t Icolind;
        const scalar_view_t Rvals = Rmat.values, Avals = Amat.values, Pvals = Pmat.values;
        scalar_view_t Ivals;

        if (!Pview.importMatrix.is_null())
        {
          const KCRS & Imat = Pview.importMatrix->getLocalMatrix();
          Irowptr = Imat.graph.row_map;
          Icolind = Imat.graph.entries;
          Ivals = Imat.values;
        }

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix Alloc"))));
  #endif

        size_t Acest_total_nnz = std::max(Ac_estimate_nnz(*Aview.origMatrix, *Pview.origMatrix),
                                    Ac.getColMap()->getNodeNumElements());
        size_t Acest_nnz_per_row = std::ceil(Acest_total_nnz / n);
        size_t nnzPerRowA = 100;
        if (Aview.origMatrix->getNodeNumEntries() > 0)
          nnzPerRowA = Aview.origMatrix->getNodeNumEntries()/Aview.origMatrix->getNodeNumRows();

        // Classic csr assembly (low memory edition)
        //
        // mfh 27 Sep 2016: Ac_estimate_nnz does not promise an upper bound.
        // The method loops over rows of R, and may resize after processing
        // each row.  Chris Siefert says that this reflects experience in
        // ML; for the non-threaded case, ML found it faster to spend less
        // effort on estimation and risk an occasional reallocation.

        const size_t ST_INVALID = Teuchos::OrdinalTraits<size_t>::invalid();

        // Get my node / thread info (right from openmp or parameter list)
        size_t thread_max =  Kokkos::Compat::KokkosOpenMPWrapperNode::execution_space::concurrency();
        if(!params.is_null()) {
          if(params->isParameter("openmp: ltg thread max"))
            thread_max = std::max((size_t)1,std::min(thread_max,params->get("openmp: ltg thread max",thread_max)));    
        }

        double thread_chunk = (double)(n) / thread_max;


  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix Fill Matrix"))));
  #endif

        // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
        // routine.  The routine computes Ac := R * A * (P_local + P_remote).
        //
        // For column index Aik in row i of A, Acol2PRow[Aik] tells
        // you whether the corresponding row of P belongs to P_local
        // ("orig") or P_remote ("Import").

        // Thread-local memory
        Kokkos::View<u_lno_view_t*> tl_rowptr("top_rowptr", thread_max);
        Kokkos::View<u_lno_nnz_view_t*> tl_colind("top_colind", thread_max);
        Kokkos::View<u_scalar_view_t*> tl_values("top_values", thread_max);

        Kokkos::parallel_for("MMM::PTAP::NewMatrix::ThreadLocal",range_type(0, thread_max).set_chunk_size(1),[=](const size_t tid)
        {
          // Thread coordination stuff
          size_t my_thread_start =  tid * thread_chunk;
          size_t my_thread_stop  = tid == thread_max-1 ? n : (tid+1)*thread_chunk;
          size_t my_thread_n     = my_thread_stop - my_thread_start;

          size_t nnzAllocated = (size_t) (my_thread_n * Acest_nnz_per_row + 100);

          std::vector<size_t> ac_status(maxAccol + 1, ST_INVALID);

          //manually allocate the thread-local storage for Ac
          u_lno_view_t Acrowptr((typename u_lno_view_t::data_type) malloc(u_lno_view_t::shmem_size(my_thread_n+1)), my_thread_n + 1);
          u_lno_nnz_view_t Accolind((typename u_lno_nnz_view_t::data_type) malloc(u_lno_nnz_view_t::shmem_size(nnzAllocated)), nnzAllocated);
          u_scalar_view_t Acvals((typename u_scalar_view_t::data_type) malloc(u_scalar_view_t::shmem_size(nnzAllocated)), nnzAllocated);
          // mfh 27 Sep 2016: m is the number of rows in the input matrix R
          // on the calling process.
          size_t nnz = 0, nnz_old = 0;
          // bmk: loop over the rows of R which are assigned to thread tid
          for (size_t i = my_thread_start; i < my_thread_stop; i++) {
            Acrowptr(i - my_thread_start) = nnz;
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
                  size_t Pl = Teuchos::as<size_t>(Acol2PRow[l]);

                  // mfh 27 Sep 2016: Go through all entries in that row of P_local.
                  for (size_t jj = Prowptr(Pl); jj < Prowptr(Pl+1); jj++) {
                    LO j = Pcolind(jj);
                    SC Plj = Pvals(jj);
                    if (Plj == SC_ZERO)
                      continue;
                    LO Acj = Pcol2Accol[j];

                    if (ac_status[Acj] == ST_INVALID || ac_status[Acj] < nnz_old) {
    #ifdef HAVE_TPETRA_DEBUG
                      // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                      TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.dimension_0()),
                                                 std::runtime_error,
                                                 label << " ERROR, not enough memory allocated for matrix product. Allocated: "  << Accolind.dimension_0() << std::endl);
    #endif
                      // New entry
                      ac_status[Acj]   = nnz;
                      Accolind(nnz) = Acj;
                      Acvals(nnz)   = Rik*Akl*Plj;
                      nnz++;
                    } else {
                      Acvals(ac_status[Acj]) += Rik*Akl*Plj;
                    }
                  }
                } else {
                  // mfh 27 Sep 2016: If the entry of Acol2PRow
                  // corresponding to the current entry of A is NOT populated (has
                  // a flag "invalid" value), then the corresponding row of P is
                  // in P_reote (i.e., it does not live on the calling process).

                  // Remote matrix
                  size_t Il = Teuchos::as<size_t>(Acol2PRowImport[l]);
                  for (size_t jj = Irowptr[Il]; jj < Irowptr[Il+1]; jj++) {
                    LO j = Icolind(jj);
                    SC Plj = Ivals(jj);
                    if (Plj == SC_ZERO)
                      continue;
                    LO Acj = PIcol2Accol[j];

                    if (ac_status[Acj] == ST_INVALID || ac_status[Acj] < nnz_old){
    #ifdef HAVE_TPETRA_DEBUG
                      // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                      TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.size()),
                                                 std::runtime_error,
                                                 label << " ERROR, not enough memory allocated for matrix product. Allocated: "  << Accolind.size() << std::endl);
    #endif
                      // New entry
                      ac_status[Acj]   = nnz;
                      Accolind(nnz) = Acj;
                      Acvals(nnz)   = Rik*Akl*Plj;
                      nnz++;
                    } else {
                      Acvals(ac_status[Acj]) += Rik*Akl*Plj;
                    }
                  }
                }
              }
            }
            // Resize for next pass if needed
            if (nnz + std::max(5*nnzPerRowA, n) > nnzAllocated) {
              nnzAllocated *= 2;
              nnzAllocated = std::max(nnzAllocated, nnz + std::max(5*nnzPerRowA, n));
              Accolind = u_lno_nnz_view_t((typename u_lno_nnz_view_t::data_type) realloc(Accolind.data(), u_lno_nnz_view_t::shmem_size(nnzAllocated)), nnzAllocated);
              Acvals = u_scalar_view_t((typename u_scalar_view_t::data_type) realloc(Acvals.data(), u_scalar_view_t::shmem_size(nnzAllocated)), nnzAllocated);
            }
            nnz_old = nnz;
          }
          tl_rowptr(tid) = Acrowptr;
          tl_colind(tid) = Accolind;
          tl_values(tid) = Acvals;      
          Acrowptr(my_thread_n) = nnz;
        });

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP copy from thread local"))));
  #endif
        lno_view_t rowmapAc("non_const_lnow_row", n + 1);
        lno_nnz_view_t entriesAc;
        scalar_view_t valuesAc;
        copy_out_from_thread_memory(tl_rowptr, tl_colind, tl_values, n, thread_chunk, rowmapAc, entriesAc, valuesAc);

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP sort"))));
  #endif

        // Final sort & set of CRS arrays
        //
        // TODO (mfh 27 Sep 2016) Will the thread-parallel "local" sparse
        // matrix-matrix multiply routine sort the entries for us?
        Import_Util::sortCrsEntries(rowmapAc, entriesAc, valuesAc);

        // mfh 27 Sep 2016: This just sets pointers.
        Ac.setAllValues(rowmapAc, entriesAc, valuesAc);

  #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix ESFC"))));
  #endif

        // Final FillComplete
        //
        // mfh 27 Sep 2016: So-called "expert static fill complete" bypasses
        // Import (from domain Map to column Map) construction (which costs
        // lots of communication) by taking the previously constructed
        // Import object.  We should be able to do this without interfering
        // with the implementation of the local part of sparse matrix-matrix
        // multply above.
        RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
        labelList->set("Timer Label",label);
        // labelList->set("Sort column Map ghost GIDs")
        if(!params.is_null()) labelList->set("compute global constants",params->get("compute global constants",true));
        RCP<const Export<LO,GO,NO> > dummyExport;
        Ac.expertStaticFillComplete(Pview.origMatrix->getDomainMap(),
                                    Pview.origMatrix->getDomainMap(),
                                    Acimport,
                                    dummyExport, labelList);
      }
    };
#endif  //OpenMP


    // Kernel method for computing the local portion of Ac = PT*A*P
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void mult_PT_A_P_newmatrix(
                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                               const std::string& label,
                               const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;

      //typedef Scalar            SC; // unused
      typedef LocalOrdinal      LO;
      typedef GlobalOrdinal     GO;
      typedef Node              NO;

      typedef Import<LO,GO,NO>  import_type;
      typedef Map<LO,GO,NO>     map_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP")))));
#endif

      LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

      // Build the final importer / column map, hash table lookups for Ac
      RCP<const import_type> Acimport;
      RCP<const map_type>    Accolmap;
      RCP<const import_type> Pimport = Pview.origMatrix->getGraph()->getImporter();
      RCP<const import_type> Iimport = Pview.importMatrix.is_null() ?
        Teuchos::null : Pview.importMatrix->getGraph()->getImporter();

      // mfh 27 Sep 2016: Pcol2Accol is a table that maps from local column
      // indices of P, to local column indices of Ac.  (P and Ac have the
      // same number of columns.)  The kernel uses this, instead of
      // copying the entire input matrix P and converting its column
      // indices to those of Ac.
      Array<LO> Pcol2Accol(Pview.colMap->getNodeNumElements()), PIcol2Accol;

      if (Pview.importMatrix.is_null()) {
        // mfh 27 Sep 2016: P has no "remotes," so P and Ac have the same column Map.
        Acimport = Pimport;
        Accolmap = Pview.colMap;
        // Pcol2Accol is trivial
        for (size_t i = 0; i < Pview.colMap->getNodeNumElements(); i++)
          Pcol2Accol[i] = Teuchos::as<LO>(i);

      } else {
        // mfh 27 Sep 2016: P has "remotes," so we need to build the
        // column Map of Ac, as well as Ac's Import object (from its domain
        // Map to its column Map).  Ac's column Map is the union of the
        // column Maps of (the local part of) P, and the "remote" part of
        // P.  Ditto for the Import.  We have optimized this "setUnion"
        // operation on Import objects and Maps.

        // Choose the right variant of setUnion
        if (!Pimport.is_null() && !Iimport.is_null())
          Acimport = Pimport->setUnion(*Iimport);

        else if (!Pimport.is_null() && Iimport.is_null())
          Acimport = Pimport->setUnion();

        else if (Pimport.is_null() && !Iimport.is_null())
          Acimport = Iimport->setUnion();

        else
          throw std::runtime_error("TpetraExt::PTAP status of matrix importers is nonsensical");

        Accolmap = Acimport->getTargetMap();

        // FIXME (mfh 27 Sep 2016) This error check requires an all-reduce
        // in general.  We should get rid of it in order to reduce
        // communication costs of sparse matrix-matrix multiply.
        TEUCHOS_TEST_FOR_EXCEPTION(!Acimport->getSourceMap()->isSameAs(*Pview.origMatrix->getDomainMap()),
                                   std::runtime_error, "Tpetra::PTAP: Import setUnion messed with the DomainMap in an unfortunate way");

        // NOTE: This is not efficient and should be folded into setUnion
        //
        // mfh 27 Sep 2016: What the above comment means, is that the
        // setUnion operation on Import objects could also compute these
        // local index - to - local index look-up tables.
        PIcol2Accol.resize(Pview.importMatrix->getColMap()->getNodeNumElements());
        ArrayView<const GO> Pgid = Pview.origMatrix->getColMap()->getNodeElementList();
        ArrayView<const GO> Igid = Pview.importMatrix->getColMap()->getNodeElementList();

        for (size_t i = 0; i < Pview.origMatrix->getColMap()->getNodeNumElements(); i++)
          Pcol2Accol[i] = Accolmap->getLocalElement(Pgid[i]);

        for (size_t i = 0; i < Pview.importMatrix->getColMap()->getNodeNumElements(); i++)
          PIcol2Accol[i] = Accolmap->getLocalElement(Igid[i]);

      }

      // Replace the column map
      //
      // mfh 27 Sep 2016: We do this because Ac was originally created
      // without a column Map.  Now we have its column Map.
      Ac.replaceColMap(Accolmap);

      // mfh 27 Sep 2016: Construct tables that map from local column
      // indices of A, to local row indices of either P_local (the locally
      // owned part of P), or P_remote (the "imported" remote part of P).
      //
      // For column index Aik in row i of A, if the corresponding row of P
      // exists in the local part of P ("orig") (which I'll call P_local),
      // then Acol2PRow[Aik] is the local index of that row of P.
      // Otherwise, Acol2PRow[Aik] is "invalid" (a flag value).
      //
      // For column index Aik in row i of A, if the corresponding row of P
      // exists in the remote part of P ("Import") (which I'll call
      // P_remote), then Acol2PRowImport[Aik] is the local index of
      // that row of P.  Otherwise, Acol2PRowImport[Aik] is "invalid"
      // (a flag value).

      // Run through all the hash table lookups once and for all

      // These two are needed for the product A*P
      Array<LO> Acol2PRow (Aview.colMap->getNodeNumElements(), LO_INVALID);
      Array<LO> Acol2PRowImport(Aview.colMap->getNodeNumElements(), LO_INVALID);
      for (LO A_LID = Aview.colMap->getMinLocalIndex(); A_LID <= Aview.colMap->getMaxLocalIndex(); A_LID++) {
        GO A_GID = Aview.colMap->getGlobalElement(A_LID);
        LO P_LID = Pview.origMatrix->getRowMap()->getLocalElement(A_GID);
        if (P_LID != LO_INVALID) {
          Acol2PRow[A_LID] = P_LID;
        } else {
          LO I_LID = Pview.importMatrix->getRowMap()->getLocalElement(A_GID);
          Acol2PRowImport[A_LID] = I_LID;
        }
      }

      // Call the actual kernel.  We'll rely on partial template specialization to call the correct one ---
      // Either the straight-up Tpetra code (SerialNode) or the KokkosKernels one (other NGP node types)
      KernelWrappers3MMM_Specialized<Scalar,LocalOrdinal,GlobalOrdinal,Node>::mult_PT_A_P_newmatrix_kernel_wrapper(Aview,
                                                                                                                   Pview,
                                                                                                                   Acol2PRow,
                                                                                                                   Acol2PRowImport,
                                                                                                                   Pcol2Accol,
                                                                                                                   PIcol2Accol,
                                                                                                                   Ac,
                                                                                                                   Acimport,
                                                                                                                   label,
                                                                                                                   params);

    }


    /*********************************************************************************************************/
    // PT_A_P NewMatrix Kernel wrappers (Default, general, non-threaded version)
    // Computes P.T * A * P -> Ac
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void KernelWrappers3MMM_Specialized<Scalar,LocalOrdinal,GlobalOrdinal,Node>::mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                                          CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                                                                          const Teuchos::Array<LocalOrdinal> & Acol2PRow,
                                                                                                          const Teuchos::Array<LocalOrdinal> & Acol2PRowImport,
                                                                                                          const Teuchos::Array<LocalOrdinal> & Pcol2Accol,
                                                                                                          const Teuchos::Array<LocalOrdinal> & PIcol2Accol,
                                                                                                          CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                                                                          Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                                                                          const std::string& label,
                                                                                                          const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix SerialCore"))));
#endif

      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;

      typedef Scalar            SC;
      typedef LocalOrdinal      LO;
      typedef GlobalOrdinal     GO;
      typedef Node              NO;
      typedef RowMatrixTransposer<SC,LO,GO,NO>  transposer_type;
      const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
      const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

      // number of rows on the process of the fine matrix
      // size_t m = Pview.origMatrix->getNodeNumRows();
      // number of rows on the process of the coarse matrix
      size_t n = Ac.getRowMap()->getNodeNumElements();
      LO maxAccol = Ac.getColMap()->getMaxLocalIndex();

      // Get Data Pointers
      ArrayRCP<const size_t> Arowptr_RCP, Prowptr_RCP, Irowptr_RCP;
      ArrayRCP<size_t> Acrowptr_RCP;
      ArrayRCP<const LO> Acolind_RCP, Pcolind_RCP, Icolind_RCP;
      ArrayRCP<LO> Accolind_RCP;
      ArrayRCP<const Scalar> Avals_RCP, Pvals_RCP, Ivals_RCP;
      ArrayRCP<SC> Acvals_RCP;

      // mfh 27 Sep 2016: "getAllValues" just gets the three CSR arrays
      // out of the CrsMatrix.  This code computes R * A * (P_local +
      // P_remote), where P_local contains the locally owned rows of P,
      // and P_remote the (previously Import'ed) remote rows of P.

      Aview.origMatrix->getAllValues(Arowptr_RCP, Acolind_RCP, Avals_RCP);
      Pview.origMatrix->getAllValues(Prowptr_RCP, Pcolind_RCP, Pvals_RCP);

      if (!Pview.importMatrix.is_null())
        Pview.importMatrix->getAllValues(Irowptr_RCP, Icolind_RCP, Ivals_RCP);

      // mfh 27 Sep 2016: Remark below "For efficiency" refers to an issue
      // where Teuchos::ArrayRCP::operator[] may be slower than
      // Teuchos::ArrayView::operator[].

      // For efficiency
      ArrayView<const size_t>   Arowptr, Prowptr, Irowptr;
      ArrayView<const LO>       Acolind, Pcolind, Icolind;
      ArrayView<const SC>       Avals, Pvals, Ivals;
      ArrayView<size_t>         Acrowptr;
      ArrayView<LO> Accolind;
      ArrayView<SC> Acvals;
      Arowptr = Arowptr_RCP();  Acolind = Acolind_RCP();  Avals = Avals_RCP();
      Prowptr = Prowptr_RCP();  Pcolind = Pcolind_RCP();  Pvals = Pvals_RCP();
      if (!Pview.importMatrix.is_null()) {
        Irowptr = Irowptr_RCP(); Icolind = Icolind_RCP(); Ivals = Ivals_RCP();
      }

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose"))));
#endif

      //////////////////////////////////////////////////////////////////////
      // Get the local transpose of the graph of P by locally transposing
      // all of P

      ArrayRCP<const size_t>  Rrowptr_RCP;
      ArrayRCP<const LO>      Rcolind_RCP;
      ArrayRCP<const Scalar>  Rvals_RCP;
      ArrayView<const size_t> Rrowptr;
      ArrayView<const LO>     Rcolind;
      ArrayView<const SC>     Rvals;

      transposer_type transposer (Pview.origMatrix,label+std::string("XP: "));
      RCP<Teuchos::ParameterList> transposeParams = Teuchos::rcp(new Teuchos::ParameterList);
      if (!params.is_null())
        transposeParams->set("compute global constants",
                             params->get("compute global constants: temporaries",
                                         false));
      RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ptrans = transposer.createTransposeLocal(transposeParams);

      Ptrans->getAllValues(Rrowptr_RCP, Rcolind_RCP, Rvals_RCP);
      Rrowptr = Rrowptr_RCP();
      Rcolind = Rcolind_RCP();
      Rvals = Rvals_RCP();

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix Alloc"))));
#endif

      size_t nnz_alloc = std::max(Ac_estimate_nnz(*Aview.origMatrix, *Pview.origMatrix),
                                  Ac.getColMap()->getNodeNumElements());
      size_t nnzPerRowA = 100;
      if (Aview.origMatrix->getNodeNumEntries() > 0)
        nnzPerRowA = Aview.origMatrix->getNodeNumEntries()/Aview.origMatrix->getNodeNumRows();
      Acrowptr_RCP.resize(n+1);       Acrowptr = Acrowptr_RCP();
      Accolind_RCP.resize(nnz_alloc); Accolind = Accolind_RCP();
      Acvals_RCP.resize(nnz_alloc);   Acvals   = Acvals_RCP();

      // Classic csr assembly (low memory edition)
      //
      // mfh 27 Sep 2016: Ac_estimate_nnz does not promise an upper bound.
      // The method loops over rows of R, and may resize after processing
      // each row.  Chris Siefert says that this reflects experience in
      // ML; for the non-threaded case, ML found it faster to spend less
      // effort on estimation and risk an occasional reallocation.

      // mfh 27 Sep 2016: The ac_status array is an implementation detail
      // of the local sparse matrix-matrix multiply routine.

      // The status array will contain the index into colind where this entry was last deposited.
      //   ac_status[i] <  nnz - not in the row yet
      //   ac_status[i] >= nnz - this is the entry where you can find the data
      // We start with this filled with INVALID's indicating that there are no entries yet.
      // Sadly, this complicates the code due to the fact that size_t's are unsigned.
      const size_t ST_INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
      Array<size_t> ac_status(maxAccol+1, ST_INVALID);


#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix Fill Matrix"))));
#endif

      // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
      // routine.  The routine computes Ac := R * A * (P_local + P_remote).
      //
      // For column index Aik in row i of A, Acol2PRow[Aik] tells
      // you whether the corresponding row of P belongs to P_local
      // ("orig") or P_remote ("Import").

      size_t nnz = 0, nnz_old = 0;

      // For each row of R
      for (size_t i = 0; i < n; i++) {
        std::vector<size_t> ac_status(maxAccol + 1, ST_INVALID);
        // mfh 27 Sep 2016: m is the number of rows in the input matrix R
        // on the calling process.
        Acrowptr[i] = nnz;

        // mfh 27 Sep 2016: For each entry of R in the current row of R
        for (size_t kk = Rrowptr[i]; kk < Rrowptr[i+1]; kk++) {
          LO k  = Rcolind[kk]; // local column index of current entry of R
          const SC Rik = Rvals[kk];   // value of current entry of R
          if (Rik == SC_ZERO)
            continue; // skip explicitly stored zero values in R
          // For each entry of A in the current row of A
          for (size_t ll = Arowptr[k]; ll < Arowptr[k+1]; ll++) {
            LO l = Acolind[ll]; // local column index of current entry of A
            const SC Akl = Avals[ll];   // value of current entry of A
            if (Akl == SC_ZERO)
              continue; // skip explicitly stored zero values in A

            if (Acol2PRow[l] != LO_INVALID) {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is populated, then
              // the corresponding row of P is in P_local (i.e., it lives on
              // the calling process).

              // Local matrix
              size_t Pl = Teuchos::as<size_t>(Acol2PRow[l]);

              // mfh 27 Sep 2016: Go through all entries in that row of P_local.
              for (size_t jj = Prowptr[Pl]; jj < Prowptr[Pl+1]; jj++) {
                LO j = Pcolind[jj];
                SC Plj = Pvals[jj];
                if (Plj == SC_ZERO)
                  continue;
                LO Acj = Pcol2Accol[j];

                if (ac_status[Acj] == ST_INVALID || ac_status[Acj] < nnz_old) {
#ifdef HAVE_TPETRA_DEBUG
                  // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                  TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.size()),
                                             std::runtime_error,
                                             label << " ERROR, not enough memory allocated for matrix product. Allocated: "  << Accolind.size() << std::endl);
#endif
                  // New entry
                  ac_status[Acj]   = nnz;
                  Accolind[nnz] = Acj;
                  Acvals[nnz]   = Rik*Akl*Plj;
                  nnz++;
                } else {
                  Acvals[ac_status[Acj]] += Rik*Akl*Plj;
                }
              }
            } else {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is NOT populated (has
              // a flag "invalid" value), then the corresponding row of P is
              // in P_reote (i.e., it does not live on the calling process).

              // Remote matrix
              size_t Il = Teuchos::as<size_t>(Acol2PRowImport[l]);
              for (size_t jj = Irowptr[Il]; jj < Irowptr[Il+1]; jj++) {
                LO j = Icolind[jj];
                SC Plj = Ivals[jj];
                if (Plj == SC_ZERO)
                  continue;
                LO Acj = PIcol2Accol[j];

                if (ac_status[Acj] == ST_INVALID || ac_status[Acj] < nnz_old){
#ifdef HAVE_TPETRA_DEBUG
                  // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                  TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Accolind.size()),
                                             std::runtime_error,
                                             label << " ERROR, not enough memory allocated for matrix product. Allocated: "  << Accolind.size() << std::endl);
#endif
                  // New entry
                  ac_status[Acj]   = nnz;
                  Accolind[nnz] = Acj;
                  Acvals[nnz]   = Rik*Akl*Plj;
                  nnz++;
                } else {
                  Acvals[ac_status[Acj]] += Rik*Akl*Plj;
                }
              }
            }
          }
        }
        // Resize for next pass if needed
        if (nnz + std::max(5*nnzPerRowA, n) > nnz_alloc) {
          nnz_alloc *= 2;
          nnz_alloc = std::max(nnz_alloc, nnz + std::max(5*nnzPerRowA, n));
          Accolind_RCP.resize(nnz_alloc); Accolind = Accolind_RCP();
          Acvals_RCP.resize(nnz_alloc);   Acvals   = Acvals_RCP();
        }
        nnz_old = nnz;
      }

      Acrowptr[n] = nnz;

      // Downward resize
      Acvals_RCP.resize(nnz);
      Accolind_RCP.resize(nnz);

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP sort"))));
#endif

      // Final sort & set of CRS arrays
      //
      // TODO (mfh 27 Sep 2016) Will the thread-parallel "local" sparse
      // matrix-matrix multiply routine sort the entries for us?
      Import_Util::sortCrsEntries(Acrowptr_RCP(), Accolind_RCP(), Acvals_RCP());

      // mfh 27 Sep 2016: This just sets pointers.
      Ac.setAllValues(Acrowptr_RCP, Accolind_RCP, Acvals_RCP);

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix ESFC"))));
#endif

      // Final FillComplete
      //
      // mfh 27 Sep 2016: So-called "expert static fill complete" bypasses
      // Import (from domain Map to column Map) construction (which costs
      // lots of communication) by taking the previously constructed
      // Import object.  We should be able to do this without interfering
      // with the implementation of the local part of sparse matrix-matrix
      // multply above.
      RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
      labelList->set("Timer Label",label);
      // labelList->set("Sort column Map ghost GIDs")
      if(!params.is_null()) labelList->set("compute global constants",params->get("compute global constants",true));
      RCP<const Export<LO,GO,NO> > dummyExport;
      Ac.expertStaticFillComplete(Pview.origMatrix->getDomainMap(),
                                  Pview.origMatrix->getDomainMap(),
                                  Acimport,
                                  dummyExport, labelList);
    }


    /*********************************************************************************************************/
    // PT_A_P NewMatrix Kernel wrappers (Default non-threaded version)
    // Computes P.T * A * P -> Ac using a 2-pass algorithm.
    // This turned out to be slower on SerialNode, but it might still be helpful when going to Kokkos, so I left it in.
    // Currently, this implementation never gets called.
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void KernelWrappers3MMM<Scalar,LocalOrdinal,GlobalOrdinal,Node>::mult_PT_A_P_newmatrix_kernel_wrapper_2pass(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                                                CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                                                                                const Teuchos::Array<LocalOrdinal> & Acol2PRow,
                                                                                                                const Teuchos::Array<LocalOrdinal> & Acol2PRowImport,
                                                                                                                const Teuchos::Array<LocalOrdinal> & Pcol2Accol,
                                                                                                                const Teuchos::Array<LocalOrdinal> & PIcol2Accol,
                                                                                                                CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                                                                                Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                                                                                const std::string& label,
                                                                                                                const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix SerialCore"))));
#endif

      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;

      typedef Scalar            SC;
      typedef LocalOrdinal      LO;
      typedef GlobalOrdinal     GO;
      typedef Node              NO;
      typedef RowMatrixTransposer<SC,LO,GO,NO>  transposer_type;
      const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
      const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

      // number of rows on the process of the fine matrix
      // size_t m = Pview.origMatrix->getNodeNumRows();
      // number of rows on the process of the coarse matrix
      size_t n = Ac.getRowMap()->getNodeNumElements();
      LO maxAccol = Ac.getColMap()->getMaxLocalIndex();

      // Get Data Pointers
      ArrayRCP<const size_t> Arowptr_RCP, Prowptr_RCP, Irowptr_RCP;
      ArrayRCP<size_t> Acrowptr_RCP;
      ArrayRCP<const LO> Acolind_RCP, Pcolind_RCP, Icolind_RCP;
      ArrayRCP<LO> Accolind_RCP;
      ArrayRCP<const Scalar> Avals_RCP, Pvals_RCP, Ivals_RCP;
      ArrayRCP<SC> Acvals_RCP;

      // mfh 27 Sep 2016: "getAllValues" just gets the three CSR arrays
      // out of the CrsMatrix.  This code computes R * A * (P_local +
      // P_remote), where P_local contains the locally owned rows of P,
      // and P_remote the (previously Import'ed) remote rows of P.

      Aview.origMatrix->getAllValues(Arowptr_RCP, Acolind_RCP, Avals_RCP);
      Pview.origMatrix->getAllValues(Prowptr_RCP, Pcolind_RCP, Pvals_RCP);

      if (!Pview.importMatrix.is_null())
        Pview.importMatrix->getAllValues(Irowptr_RCP, Icolind_RCP, Ivals_RCP);

      // mfh 27 Sep 2016: Remark below "For efficiency" refers to an issue
      // where Teuchos::ArrayRCP::operator[] may be slower than
      // Teuchos::ArrayView::operator[].

      // For efficiency
      ArrayView<const size_t>   Arowptr, Prowptr, Irowptr;
      ArrayView<const LO>       Acolind, Pcolind, Icolind;
      ArrayView<const SC>       Avals, Pvals, Ivals;
      ArrayView<size_t>         Acrowptr;
      ArrayView<LO> Accolind;
      ArrayView<SC> Acvals;
      Arowptr = Arowptr_RCP();  Acolind = Acolind_RCP();  Avals = Avals_RCP();
      Prowptr = Prowptr_RCP();  Pcolind = Pcolind_RCP();  Pvals = Pvals_RCP();
      if (!Pview.importMatrix.is_null()) {
        Irowptr = Irowptr_RCP(); Icolind = Icolind_RCP(); Ivals = Ivals_RCP();
      }

      //////////////////////////////////////////////////////////////////////
      // In a first pass, determine the graph of Ac.
      //////////////////////////////////////////////////////////////////////

      //////////////////////////////////////////////////////////////////////
      // Get the graph of Ac. This gets the local transpose of P,
      // then loops over R, A, P to get the graph of Ac.
      //////////////////////////////////////////////////////////////////////

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose"))));
#endif

      //////////////////////////////////////////////////////////////////////
      // Get the local transpose of the graph of P by locally transposing
      // all of P

      ArrayRCP<const size_t>  Rrowptr_RCP;
      ArrayRCP<const LO>      Rcolind_RCP;
      ArrayRCP<const Scalar>  Rvals_RCP;
      ArrayView<const size_t> Rrowptr;
      ArrayView<const LO>     Rcolind;
      ArrayView<const SC>     Rvals;

      transposer_type transposer (Pview.origMatrix,label+std::string("XP: "));
      RCP<Teuchos::ParameterList> transposeParams = Teuchos::rcp(new Teuchos::ParameterList);
      if (!params.is_null())
        transposeParams->set("compute global constants",
                             params->get("compute global constants: temporaries",
                                         false));
      RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ptrans = transposer.createTransposeLocal(transposeParams);

      Ptrans->getAllValues(Rrowptr_RCP, Rcolind_RCP, Rvals_RCP);
      Rrowptr = Rrowptr_RCP();
      Rcolind = Rcolind_RCP();
      Rvals = Rvals_RCP();

      //////////////////////////////////////////////////////////////////////
      // Construct graph

      #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP graph"))));
      #endif

      const size_t ST_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
      Array<size_t> ac_status(maxAccol + 1, ST_INVALID);

      size_t nnz_alloc = std::max(Ac_estimate_nnz(*Aview.origMatrix, *Pview.origMatrix), n);
      size_t nnzPerRowA = 100;
      if (Aview.origMatrix->getNodeNumEntries() > 0)
        nnzPerRowA = Aview.origMatrix->getNodeNumEntries()/Aview.origMatrix->getNodeNumRows();
      Acrowptr_RCP.resize(n+1);
      Acrowptr = Acrowptr_RCP();
      Accolind_RCP.resize(nnz_alloc);
      Accolind = Accolind_RCP();

      size_t nnz = 0, nnz_old = 0;
      for (size_t i = 0; i < n; i++) {
        // mfh 27 Sep 2016: m is the number of rows in the input matrix R
        // on the calling process.
        Acrowptr[i] = nnz;

        // mfh 27 Sep 2016: For each entry of R in the current row of R
        for (size_t kk = Rrowptr[i]; kk < Rrowptr[i+1]; kk++) {
          LO k  = Rcolind[kk]; // local column index of current entry of R
          // For each entry of A in the current row of A
          for (size_t ll = Arowptr[k]; ll < Arowptr[k+1]; ll++) {
            LO l = Acolind[ll]; // local column index of current entry of A

            if (Acol2PRow[l] != LO_INVALID) {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is populated, then
              // the corresponding row of P is in P_local (i.e., it lives on
              // the calling process).

              // Local matrix
              size_t Pl = Teuchos::as<size_t>(Acol2PRow[l]);

              // mfh 27 Sep 2016: Go through all entries in that row of P_local.
              for (size_t jj = Prowptr[Pl]; jj < Prowptr[Pl+1]; jj++) {
                LO j = Pcolind[jj];
                LO Acj = Pcol2Accol[j];

                if (ac_status[Acj] == ST_INVALID || ac_status[Acj] < nnz_old) {
                  // New entry
                  ac_status[Acj]   = nnz;
                  Accolind[nnz] = Acj;
                  nnz++;
                }
              }
            } else {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is NOT populated (has
              // a flag "invalid" value), then the corresponding row of P is
              // in P_remote (i.e., it does not live on the calling process).

              // Remote matrix
              size_t Il = Teuchos::as<size_t>(Acol2PRowImport[l]);
              for (size_t jj = Irowptr[Il]; jj < Irowptr[Il+1]; jj++) {
                LO j = Icolind[jj];
                LO Acj = PIcol2Accol[j];

                if (ac_status[Acj] == ST_INVALID || ac_status[Acj] < nnz_old){
                  // New entry
                  ac_status[Acj]   = nnz;
                  Accolind[nnz] = Acj;
                  nnz++;
                }
              }
            }
          }
        }
        // Resize for next pass if needed
        // cag: Maybe we can do something more subtle here, and not double
        //      the size right away.
        if (nnz + std::max(5*nnzPerRowA, n) > nnz_alloc) {
          nnz_alloc *= 2;
          nnz_alloc = std::max(nnz_alloc, nnz + std::max(5*nnzPerRowA, n));
          Accolind_RCP.resize(nnz_alloc); Accolind = Accolind_RCP();
          Acvals_RCP.resize(nnz_alloc);   Acvals   = Acvals_RCP();
        }
        nnz_old = nnz;
      }
      Acrowptr[n] = nnz;

      // Downward resize
      Accolind_RCP.resize(nnz);
      Accolind = Accolind_RCP();

      // Allocate Acvals
      Acvals_RCP.resize(nnz, SC_ZERO);
      Acvals = Acvals_RCP();


      //////////////////////////////////////////////////////////////////////
      // In a second pass, enter the values into Acvals.
      //////////////////////////////////////////////////////////////////////

      #ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix Fill Matrix"))));
      #endif


      for (size_t k = 0; k < n; k++) {
        for (size_t ii = Prowptr[k]; ii < Prowptr[k+1]; ii++) {
          LO i = Pcolind[ii];
          const SC Pki = Pvals[ii];
          for (size_t ll = Arowptr[k]; ll < Arowptr[k+1]; ll++) {
            LO l = Acolind[ll];
            const SC Akl = Avals[ll];
            if (Akl == 0.)
              continue;
            if (Acol2PRow[l] != LO_INVALID) {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is populated, then
              // the corresponding row of P is in P_local (i.e., it lives on
              // the calling process).

              // Local matrix
              size_t Pl = Teuchos::as<size_t>(Acol2PRow[l]);
              for (size_t jj = Prowptr[Pl]; jj < Prowptr[Pl+1]; jj++) {
                LO j = Pcolind[jj];
                LO Acj = Pcol2Accol[j];
                size_t pp;
                for (pp = Acrowptr[i]; pp < Acrowptr[i+1]; pp++)
                  if (Accolind[pp] == Acj)
                    break;
                // TEUCHOS_TEST_FOR_EXCEPTION(Accolind[pp] != Acj,
                //                            std::runtime_error, "problem with Ac column indices");
                Acvals[pp] += Pki*Akl*Pvals[jj];
              }
            } else {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A NOT populated (has
              // a flag "invalid" value), then the corresponding row of P is
              // in P_remote (i.e., it does not live on the calling process).

              // Remote matrix
              size_t Il = Teuchos::as<size_t>(Acol2PRowImport[l]);
              for (size_t jj = Irowptr[Il]; jj < Irowptr[Il+1]; jj++) {
                LO j = Icolind[jj];
                LO Acj = PIcol2Accol[j];
                size_t pp;
                for (pp = Acrowptr[i]; pp < Acrowptr[i+1]; pp++)
                  if (Accolind[pp] == Acj)
                    break;
                // TEUCHOS_TEST_FOR_EXCEPTION(Accolind[pp] != Acj,
                //                            std::runtime_error, "problem with Ac column indices");
                Acvals[pp] += Pki*Akl*Ivals[jj];
              }
            }
          }
        }
      }


#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP sort"))));
#endif

      // Final sort & set of CRS arrays
      //
      // TODO (mfh 27 Sep 2016) Will the thread-parallel "local" sparse
      // matrix-matrix multiply routine sort the entries for us?
      Import_Util::sortCrsEntries(Acrowptr_RCP(), Accolind_RCP(), Acvals_RCP());

      // mfh 27 Sep 2016: This just sets pointers.
      Ac.setAllValues(Acrowptr_RCP, Accolind_RCP, Acvals_RCP);

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP Newmatrix ESFC"))));
#endif

      // Final FillComplete
      //
      // mfh 27 Sep 2016: So-called "expert static fill complete" bypasses
      // Import (from domain Map to column Map) construction (which costs
      // lots of communication) by taking the previously constructed
      // Import object.  We should be able to do this without interfering
      // with the implementation of the local part of sparse matrix-matrix
      // multply above.
      RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
      labelList->set("Timer Label",label);
      // labelList->set("Sort column Map ghost GIDs")
      if(!params.is_null()) labelList->set("compute global constants",params->get("compute global constants",true));
      RCP<const Export<LO,GO,NO> > dummyExport;
      Ac.expertStaticFillComplete(Pview.origMatrix->getDomainMap(),
                                  Pview.origMatrix->getDomainMap(),
                                  Acimport,
                                  dummyExport, labelList);
    }



  } //End namepsace MMdetails

} //End namespace Tpetra
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_TRIPLEMATRIXMULTIPLY_INSTANT(SCALAR,LO,GO,NODE)                  \
                                                                        \
  template                                                              \
  void TripleMatrixMultiply::MultiplyRAP(                               \
                                         const CrsMatrix< SCALAR , LO , GO , NODE >& R, \
                                         bool transposeR,               \
                                         const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
                                         bool transposeA,               \
                                         const CrsMatrix< SCALAR , LO , GO , NODE >& P, \
                                         bool transposeP,               \
                                         CrsMatrix< SCALAR , LO , GO , NODE >& Ac, \
                                         bool call_FillComplete_on_result, \
                                         const std::string & label,     \
                                         const Teuchos::RCP<Teuchos::ParameterList>& params); \
                                                                        \


#endif // TPETRA_TRIPLEMATRIXMULTIPLY_DEF_HPP
