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



/*********************************************************************************************************/
// Include the architecture-specific kernel partial specializations here
// NOTE: This needs to be outside all namespaces
#include "TpetraExt_MatrixMatrix_OpenMP.hpp"
#include "TpetraExt_MatrixMatrix_Cuda.hpp"

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

      using Teuchos::ParameterList;
      RCP<ParameterList> transposeParams (new ParameterList);
      transposeParams->set ("sort", false);

      if (transposeR && &R != &P) {
        transposer_type transposer(rcpFromRef (R));
        Rprime = transposer.createTranspose (transposeParams);
      } else {
        Rprime = rcpFromRef(R);
      }

      if (transposeA) {
        transposer_type transposer(rcpFromRef (A));
        Aprime = transposer.createTranspose (transposeParams);
      } else {
        Aprime = rcpFromRef(A);
      }

      if (transposeP) {
        transposer_type transposer(rcpFromRef (P));
        Pprime = transposer.createTranspose (transposeParams);
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
      MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP All I&X"))));
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
      MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP All Multiply"))));
#endif

      // Call the appropriate method to perform the actual multiplication.
      if (call_FillComplete_on_result && newFlag) {
        if (transposeR && &R == &P)
          MMdetails::mult_PT_A_P_newmatrix(Aview, Pview, *Actemp, label, params);
        else
          MMdetails::mult_R_A_P_newmatrix(Rview, Aview, Pview, *Actemp, label, params);
      } else if (call_FillComplete_on_result) {
        if (transposeR && &R == &P)
          MMdetails::mult_PT_A_P_reuse(Aview, Pview, *Actemp, label, params);
        else
          MMdetails::mult_R_A_P_reuse(Rview, Aview, Pview, *Actemp, label, params);
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
        MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP exportAndFillComplete"))));
#endif
        Teuchos::ParameterList labelList;
        labelList.set("Timer Label", label);
        Teuchos::ParameterList& labelList_subList = labelList.sublist("matrixmatrix: kernel params",false);

        RCP<crs_matrix_type> Acprime = rcpFromRef(Ac);
        bool isMM = true;
        bool overrideAllreduce = false;
        int mm_optimization_core_count=::Tpetra::Details::Behavior::TAFC_OptimizationCoreCount();
        if(!params.is_null()) {
          Teuchos::ParameterList& params_sublist = params->sublist("matrixmatrix: kernel params",false);
          mm_optimization_core_count = ::Tpetra::Details::Behavior::TAFC_OptimizationCoreCount();
          mm_optimization_core_count = params_sublist.get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
          int mm_optimization_core_count2 = params->get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
          if(mm_optimization_core_count2<mm_optimization_core_count) mm_optimization_core_count=mm_optimization_core_count2;
          isMM = params_sublist.get("isMatrixMatrix_TransferAndFillComplete",false);
          overrideAllreduce = params_sublist.get("MM_TAFC_OverrideAllreduceCheck",false);

          labelList.set("compute global constants",params->get("compute global constants",true));
        }
        labelList_subList.set("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count,"Core Count above which the optimized neighbor discovery is used");

        labelList_subList.set("isMatrixMatrix_TransferAndFillComplete",isMM,
                                "This parameter should be set to true only for MatrixMatrix operations: the optimization in Epetra that was ported to Tpetra does _not_ take into account the possibility that for any given source PID, a particular GID may not exist on the target PID: i.e. a transfer operation. A fix for this general case is in development.");
        labelList_subList.set("MM_TAFC_OverrideAllreduceCheck",overrideAllreduce);

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

      // Kokkos typedefs
      typedef typename map_type::local_map_type local_map_type;
      typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type KCRS;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type lno_view_t;
      typedef typename NO::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
      typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;

#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP M5 Cmap")))));
#endif
      LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

      // Build the final importer / column map, hash table lookups for Ac
      RCP<const import_type> Cimport;
      RCP<const map_type>    Ccolmap;
      RCP<const import_type> Pimport = Pview.origMatrix->getGraph()->getImporter();
      RCP<const import_type> Iimport = Pview.importMatrix.is_null() ?  Teuchos::null : Pview.importMatrix->getGraph()->getImporter();
      local_map_type Acolmap_local = Aview.colMap->getLocalMap();
      local_map_type Prowmap_local = Pview.origMatrix->getRowMap()->getLocalMap();
      local_map_type Irowmap_local;  if(!Pview.importMatrix.is_null()) Irowmap_local = Pview.importMatrix->getRowMap()->getLocalMap();
      local_map_type Pcolmap_local = Pview.origMatrix->getColMap()->getLocalMap();
      local_map_type Icolmap_local;  if(!Pview.importMatrix.is_null()) Icolmap_local = Pview.importMatrix->getColMap()->getLocalMap();


      // mfh 27 Sep 2016: Pcol2Ccol is a table that maps from local column
      // indices of B, to local column indices of Ac.  (B and Ac have the
      // same number of columns.)  The kernel uses this, instead of
      // copying the entire input matrix B and converting its column
      // indices to those of C.
      lo_view_t Pcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Pcol2Ccol"),Pview.colMap->getNodeNumElements()), Icol2Ccol;

      if (Pview.importMatrix.is_null()) {
        // mfh 27 Sep 2016: B has no "remotes," so P and C have the same column Map.
        Cimport = Pimport;
        Ccolmap = Pview.colMap;
        const LO colMapSize = static_cast<LO>(Pview.colMap->getNodeNumElements());
        // Pcol2Ccol is trivial
        Kokkos::parallel_for("Tpetra::mult_R_A_P_newmatrix::Pcol2Ccol_fill",
                             Kokkos::RangePolicy<execution_space, LO>(0, colMapSize),
                             KOKKOS_LAMBDA(const LO i) {
                               Pcol2Ccol(i) = i;
                             });
      }
      else {
        // mfh 27 Sep 2016: P has "remotes," so we need to build the
        // column Map of C, as well as C's Import object (from its domain
        // Map to its column Map).  C's column Map is the union of the
        // column Maps of (the local part of) P, and the "remote" part of
        // P.  Ditto for the Import.  We have optimized this "setUnion"
        // operation on Import objects and Maps.

        // Choose the right variant of setUnion
        if (!Pimport.is_null() && !Iimport.is_null()) {
          Cimport = Pimport->setUnion(*Iimport);
        }
        else if (!Pimport.is_null() && Iimport.is_null()) {
          Cimport = Pimport->setUnion();
        }
        else if (Pimport.is_null() && !Iimport.is_null()) {
          Cimport = Iimport->setUnion();
        }
        else {
          throw std::runtime_error("TpetraExt::RAP status of matrix importers is nonsensical");
        }
        Ccolmap = Cimport->getTargetMap();

        // FIXME (mfh 27 Sep 2016) This error check requires an all-reduce
        // in general.  We should get rid of it in order to reduce
        // communication costs of sparse matrix-matrix multiply.
        TEUCHOS_TEST_FOR_EXCEPTION(!Cimport->getSourceMap()->isSameAs(*Pview.origMatrix->getDomainMap()),
                                   std::runtime_error, "Tpetra::RAP: Import setUnion messed with the DomainMap in an unfortunate way");

        // NOTE: This is not efficient and should be folded into setUnion
        //
        // mfh 27 Sep 2016: What the above comment means, is that the
        // setUnion operation on Import objects could also compute these
        // local index - to - local index look-up tables.
        Kokkos::resize(Icol2Ccol,Pview.importMatrix->getColMap()->getNodeNumElements());
        local_map_type Ccolmap_local = Ccolmap->getLocalMap();
        Kokkos::parallel_for("Tpetra::mult_R_A_P_newmatrix::Pcol2Ccol_getGlobalElement",range_type(0,Pview.origMatrix->getColMap()->getNodeNumElements()),KOKKOS_LAMBDA(const LO i) {
            Pcol2Ccol(i) = Ccolmap_local.getLocalElement(Pcolmap_local.getGlobalElement(i));
          });
        Kokkos::parallel_for("Tpetra::mult_R_A_P_newmatrix::Icol2Ccol_getGlobalElement",range_type(0,Pview.importMatrix->getColMap()->getNodeNumElements()),KOKKOS_LAMBDA(const LO i) {
            Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
          });

      }

      // Replace the column map
      //
      // mfh 27 Sep 2016: We do this because C was originally created
      // without a column Map.  Now we have its column Map.
      Ac.replaceColMap(Ccolmap);

      // mfh 27 Sep 2016: Construct tables that map from local column
      // indices of A, to local row indices of either B_local (the locally
      // owned part of B), or B_remote (the "imported" remote part of B).
      //
      // For column index Aik in row i of A, if the corresponding row of B
      // exists in the local part of B ("orig") (which I'll call B_local),
      // then targetMapToOrigRow[Aik] is the local index of that row of B.
      // Otherwise, targetMapToOrigRow[Aik] is "invalid" (a flag value).
      //
      // For column index Aik in row i of A, if the corresponding row of B
      // exists in the remote part of B ("Import") (which I'll call
      // B_remote), then targetMapToImportRow[Aik] is the local index of
      // that row of B.  Otherwise, targetMapToOrigRow[Aik] is "invalid"
      // (a flag value).

      // Run through all the hash table lookups once and for all
      lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getNodeNumElements());
      lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),Aview.colMap->getNodeNumElements());
      Kokkos::parallel_for("Tpetra::mult_R_A_P_newmatrix::construct_tables",range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),KOKKOS_LAMBDA(const LO i) {
          GO aidx = Acolmap_local.getGlobalElement(i);
          LO P_LID = Prowmap_local.getLocalElement(aidx);
          if (P_LID != LO_INVALID) {
            targetMapToOrigRow(i)   = P_LID;
            targetMapToImportRow(i) = LO_INVALID;
          } else {
            LO I_LID = Irowmap_local.getLocalElement(aidx);
            targetMapToOrigRow(i)   = LO_INVALID;
            targetMapToImportRow(i) = I_LID;
          }
        });

      // Call the actual kernel.  We'll rely on partial template specialization to call the correct one ---
      // Either the straight-up Tpetra code (SerialNode) or the KokkosKernels one (other NGP node types)
      KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Node,lo_view_t>::
        mult_R_A_P_newmatrix_kernel_wrapper(Rview, Aview, Pview,
                                            targetMapToOrigRow,targetMapToImportRow, Pcol2Ccol, Icol2Ccol,
                                            Ac, Cimport, label, params);
    }


    // Kernel method for computing the local portion of Ac = R*A*P (reuse mode)
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void mult_R_A_P_reuse(
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

      // Kokkos typedefs
      typedef typename map_type::local_map_type local_map_type;
      typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type KCRS;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type lno_view_t;
      typedef typename NO::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
      typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;

#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP M5 Cmap")))));
#endif
      LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

      // Build the final importer / column map, hash table lookups for Ac
      RCP<const import_type> Cimport = Ac.getGraph()->getImporter();
      RCP<const map_type>    Ccolmap = Ac.getColMap();
      RCP<const import_type> Pimport = Pview.origMatrix->getGraph()->getImporter();
      RCP<const import_type> Iimport = Pview.importMatrix.is_null() ?  Teuchos::null : Pview.importMatrix->getGraph()->getImporter();
      local_map_type Acolmap_local = Aview.colMap->getLocalMap();
      local_map_type Prowmap_local = Pview.origMatrix->getRowMap()->getLocalMap();
      local_map_type Irowmap_local;  if(!Pview.importMatrix.is_null()) Irowmap_local = Pview.importMatrix->getRowMap()->getLocalMap();
      local_map_type Pcolmap_local = Pview.origMatrix->getColMap()->getLocalMap();
      local_map_type Icolmap_local;  if(!Pview.importMatrix.is_null()) Icolmap_local = Pview.importMatrix->getColMap()->getLocalMap();
      local_map_type Ccolmap_local = Ccolmap->getLocalMap();

      // Build the final importer / column map, hash table lookups for C
      lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"),Pview.colMap->getNodeNumElements()), Icol2Ccol;
      {
        // Bcol2Col may not be trivial, as Ccolmap is compressed during fillComplete in newmatrix
        // So, column map of C may be a strict subset of the column map of B
        Kokkos::parallel_for(range_type(0,Pview.origMatrix->getColMap()->getNodeNumElements()),KOKKOS_LAMBDA(const LO i) {
            Bcol2Ccol(i) = Ccolmap_local.getLocalElement(Pcolmap_local.getGlobalElement(i));
          });

        if (!Pview.importMatrix.is_null()) {
          TEUCHOS_TEST_FOR_EXCEPTION(!Cimport->getSourceMap()->isSameAs(*Pview.origMatrix->getDomainMap()),
                                     std::runtime_error, "Tpetra::MMM: Import setUnion messed with the DomainMap in an unfortunate way");

          Kokkos::resize(Icol2Ccol,Pview.importMatrix->getColMap()->getNodeNumElements());
          Kokkos::parallel_for(range_type(0,Pview.importMatrix->getColMap()->getNodeNumElements()),KOKKOS_LAMBDA(const LO i) {
              Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
            });
        }
      }

      // Run through all the hash table lookups once and for all
      lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getNodeNumElements());
      lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),Aview.colMap->getNodeNumElements());
      Kokkos::parallel_for(range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),KOKKOS_LAMBDA(const LO i) {
          GO aidx = Acolmap_local.getGlobalElement(i);
          LO B_LID = Prowmap_local.getLocalElement(aidx);
          if (B_LID != LO_INVALID) {
            targetMapToOrigRow(i)   = B_LID;
            targetMapToImportRow(i) = LO_INVALID;
          } else {
            LO I_LID = Irowmap_local.getLocalElement(aidx);
            targetMapToOrigRow(i)   = LO_INVALID;
            targetMapToImportRow(i) = I_LID;

          }
        });

      // Call the actual kernel.  We'll rely on partial template specialization to call the correct one ---
      // Either the straight-up Tpetra code (SerialNode) or the KokkosKernels one (other NGP node types)
      KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Node,lo_view_t>::
        mult_R_A_P_reuse_kernel_wrapper(Rview, Aview, Pview,
                                            targetMapToOrigRow,targetMapToImportRow, Bcol2Ccol, Icol2Ccol,
                                            Ac, Cimport, label, params);
    }


   // Kernel method for computing the local portion of Ac = R*A*P
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

      // Kokkos typedefs
      typedef typename map_type::local_map_type local_map_type;
      typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type KCRS;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type lno_view_t;
      typedef typename NO::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
      typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;

#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP M5 Cmap")))));
#endif
      LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

      // Build the final importer / column map, hash table lookups for Ac
      RCP<const import_type> Cimport;
      RCP<const map_type>    Ccolmap;
      RCP<const import_type> Pimport = Pview.origMatrix->getGraph()->getImporter();
      RCP<const import_type> Iimport = Pview.importMatrix.is_null() ?  Teuchos::null : Pview.importMatrix->getGraph()->getImporter();
      local_map_type Acolmap_local = Aview.colMap->getLocalMap();
      local_map_type Prowmap_local = Pview.origMatrix->getRowMap()->getLocalMap();
      local_map_type Irowmap_local;  if(!Pview.importMatrix.is_null()) Irowmap_local = Pview.importMatrix->getRowMap()->getLocalMap();
      local_map_type Pcolmap_local = Pview.origMatrix->getColMap()->getLocalMap();
      local_map_type Icolmap_local;  if(!Pview.importMatrix.is_null()) Icolmap_local = Pview.importMatrix->getColMap()->getLocalMap();


      // mfh 27 Sep 2016: Pcol2Ccol is a table that maps from local column
      // indices of B, to local column indices of Ac.  (B and Ac have the
      // same number of columns.)  The kernel uses this, instead of
      // copying the entire input matrix B and converting its column
      // indices to those of C.
      lo_view_t Pcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Pcol2Ccol"),Pview.colMap->getNodeNumElements()), Icol2Ccol;

      if (Pview.importMatrix.is_null()) {
        // mfh 27 Sep 2016: B has no "remotes," so P and C have the same column Map.
        Cimport = Pimport;
        Ccolmap = Pview.colMap;
        const LO colMapSize = static_cast<LO>(Pview.colMap->getNodeNumElements());
        // Pcol2Ccol is trivial
        Kokkos::parallel_for("Tpetra::mult_R_A_P_newmatrix::Pcol2Ccol_fill",
                             Kokkos::RangePolicy<execution_space, LO>(0, colMapSize),
                             KOKKOS_LAMBDA(const LO i) {
                               Pcol2Ccol(i) = i;
                             });
      }
      else {
        // mfh 27 Sep 2016: P has "remotes," so we need to build the
        // column Map of C, as well as C's Import object (from its domain
        // Map to its column Map).  C's column Map is the union of the
        // column Maps of (the local part of) P, and the "remote" part of
        // P.  Ditto for the Import.  We have optimized this "setUnion"
        // operation on Import objects and Maps.

        // Choose the right variant of setUnion
        if (!Pimport.is_null() && !Iimport.is_null()) {
          Cimport = Pimport->setUnion(*Iimport);
        }
        else if (!Pimport.is_null() && Iimport.is_null()) {
          Cimport = Pimport->setUnion();
        }
        else if (Pimport.is_null() && !Iimport.is_null()) {
          Cimport = Iimport->setUnion();
        }
        else {
          throw std::runtime_error("TpetraExt::RAP status of matrix importers is nonsensical");
        }
        Ccolmap = Cimport->getTargetMap();

        // FIXME (mfh 27 Sep 2016) This error check requires an all-reduce
        // in general.  We should get rid of it in order to reduce
        // communication costs of sparse matrix-matrix multiply.
        TEUCHOS_TEST_FOR_EXCEPTION(!Cimport->getSourceMap()->isSameAs(*Pview.origMatrix->getDomainMap()),
                                   std::runtime_error, "Tpetra::RAP: Import setUnion messed with the DomainMap in an unfortunate way");

        // NOTE: This is not efficient and should be folded into setUnion
        //
        // mfh 27 Sep 2016: What the above comment means, is that the
        // setUnion operation on Import objects could also compute these
        // local index - to - local index look-up tables.
        Kokkos::resize(Icol2Ccol,Pview.importMatrix->getColMap()->getNodeNumElements());
        local_map_type Ccolmap_local = Ccolmap->getLocalMap();
        Kokkos::parallel_for("Tpetra::mult_R_A_P_newmatrix::Pcol2Ccol_getGlobalElement",range_type(0,Pview.origMatrix->getColMap()->getNodeNumElements()),KOKKOS_LAMBDA(const LO i) {
            Pcol2Ccol(i) = Ccolmap_local.getLocalElement(Pcolmap_local.getGlobalElement(i));
          });
        Kokkos::parallel_for("Tpetra::mult_R_A_P_newmatrix::Icol2Ccol_getGlobalElement",range_type(0,Pview.importMatrix->getColMap()->getNodeNumElements()),KOKKOS_LAMBDA(const LO i) {
            Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
          });

      }

      // Replace the column map
      //
      // mfh 27 Sep 2016: We do this because C was originally created
      // without a column Map.  Now we have its column Map.
      Ac.replaceColMap(Ccolmap);

      // mfh 27 Sep 2016: Construct tables that map from local column
      // indices of A, to local row indices of either B_local (the locally
      // owned part of B), or B_remote (the "imported" remote part of B).
      //
      // For column index Aik in row i of A, if the corresponding row of B
      // exists in the local part of B ("orig") (which I'll call B_local),
      // then targetMapToOrigRow[Aik] is the local index of that row of B.
      // Otherwise, targetMapToOrigRow[Aik] is "invalid" (a flag value).
      //
      // For column index Aik in row i of A, if the corresponding row of B
      // exists in the remote part of B ("Import") (which I'll call
      // B_remote), then targetMapToImportRow[Aik] is the local index of
      // that row of B.  Otherwise, targetMapToOrigRow[Aik] is "invalid"
      // (a flag value).

      // Run through all the hash table lookups once and for all
      lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getNodeNumElements());
      lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),Aview.colMap->getNodeNumElements());

      Kokkos::parallel_for("Tpetra::mult_R_A_P_newmatrix::construct_tables",range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),KOKKOS_LAMBDA(const LO i) {
          GO aidx = Acolmap_local.getGlobalElement(i);
          LO P_LID = Prowmap_local.getLocalElement(aidx);
          if (P_LID != LO_INVALID) {
            targetMapToOrigRow(i)   = P_LID;
            targetMapToImportRow(i) = LO_INVALID;
          } else {
            LO I_LID = Irowmap_local.getLocalElement(aidx);
            targetMapToOrigRow(i)   = LO_INVALID;
            targetMapToImportRow(i) = I_LID;
          }
        });

      // Call the actual kernel.  We'll rely on partial template specialization to call the correct one ---
      // Either the straight-up Tpetra code (SerialNode) or the KokkosKernels one (other NGP node types)
      KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Node,lo_view_t>::
        mult_PT_A_P_newmatrix_kernel_wrapper(Aview, Pview,
                                            targetMapToOrigRow,targetMapToImportRow, Pcol2Ccol, Icol2Ccol,
                                            Ac, Cimport, label, params);
    }

    // Kernel method for computing the local portion of Ac = R*A*P (reuse mode)
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node>
    void mult_PT_A_P_reuse(
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

      // Kokkos typedefs
      typedef typename map_type::local_map_type local_map_type;
      typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type KCRS;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::non_const_type lno_view_t;
      typedef typename NO::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
      typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;

#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP M5 Cmap")))));
#endif
      LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

      // Build the final importer / column map, hash table lookups for Ac
      RCP<const import_type> Cimport = Ac.getGraph()->getImporter();
      RCP<const map_type>    Ccolmap = Ac.getColMap();
      RCP<const import_type> Pimport = Pview.origMatrix->getGraph()->getImporter();
      RCP<const import_type> Iimport = Pview.importMatrix.is_null() ?  Teuchos::null : Pview.importMatrix->getGraph()->getImporter();
      local_map_type Acolmap_local = Aview.colMap->getLocalMap();
      local_map_type Prowmap_local = Pview.origMatrix->getRowMap()->getLocalMap();
      local_map_type Irowmap_local;  if(!Pview.importMatrix.is_null()) Irowmap_local = Pview.importMatrix->getRowMap()->getLocalMap();
      local_map_type Pcolmap_local = Pview.origMatrix->getColMap()->getLocalMap();
      local_map_type Icolmap_local;  if(!Pview.importMatrix.is_null()) Icolmap_local = Pview.importMatrix->getColMap()->getLocalMap();
      local_map_type Ccolmap_local = Ccolmap->getLocalMap();

      // Build the final importer / column map, hash table lookups for C
      lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"),Pview.colMap->getNodeNumElements()), Icol2Ccol;
      {
        // Bcol2Col may not be trivial, as Ccolmap is compressed during fillComplete in newmatrix
        // So, column map of C may be a strict subset of the column map of B
        Kokkos::parallel_for(range_type(0,Pview.origMatrix->getColMap()->getNodeNumElements()),KOKKOS_LAMBDA(const LO i) {
            Bcol2Ccol(i) = Ccolmap_local.getLocalElement(Pcolmap_local.getGlobalElement(i));
          });

        if (!Pview.importMatrix.is_null()) {
          TEUCHOS_TEST_FOR_EXCEPTION(!Cimport->getSourceMap()->isSameAs(*Pview.origMatrix->getDomainMap()),
                                     std::runtime_error, "Tpetra::MMM: Import setUnion messed with the DomainMap in an unfortunate way");

          Kokkos::resize(Icol2Ccol,Pview.importMatrix->getColMap()->getNodeNumElements());
          Kokkos::parallel_for(range_type(0,Pview.importMatrix->getColMap()->getNodeNumElements()),KOKKOS_LAMBDA(const LO i) {
              Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
            });
        }
      }

      // Run through all the hash table lookups once and for all
      lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getNodeNumElements());
      lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),Aview.colMap->getNodeNumElements());
      Kokkos::parallel_for(range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),KOKKOS_LAMBDA(const LO i) {
          GO aidx = Acolmap_local.getGlobalElement(i);
          LO B_LID = Prowmap_local.getLocalElement(aidx);
          if (B_LID != LO_INVALID) {
            targetMapToOrigRow(i)   = B_LID;
            targetMapToImportRow(i) = LO_INVALID;
          } else {
            LO I_LID = Irowmap_local.getLocalElement(aidx);
            targetMapToOrigRow(i)   = LO_INVALID;
            targetMapToImportRow(i) = I_LID;

          }
        });

      // Call the actual kernel.  We'll rely on partial template specialization to call the correct one ---
      // Either the straight-up Tpetra code (SerialNode) or the KokkosKernels one (other NGP node types)
      KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Node,lo_view_t>::
        mult_PT_A_P_reuse_kernel_wrapper(Aview, Pview,
                                            targetMapToOrigRow,targetMapToImportRow, Bcol2Ccol, Icol2Ccol,
                                            Ac, Cimport, label, params);
    }


    /*********************************************************************************************************/
    // RAP NewMatrix Kernel wrappers (Default non-threaded version)
    // Computes R * A * P -> Ac using classic Gustavson approach
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node,
             class LocalOrdinalViewType>
    void KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalOrdinalViewType>::mult_R_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                                                                                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                                                                                           const LocalOrdinalViewType & Acol2Prow,
                                                                                                                           const LocalOrdinalViewType & Acol2PIrow,
                                                                                                                           const LocalOrdinalViewType & Pcol2Accol,
                                                                                                                           const LocalOrdinalViewType & PIcol2Accol,
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

      // Lots and lots of typedefs
      typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type KCRS;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::const_type c_lno_view_t;
      typedef typename graph_t::row_map_type::non_const_type lno_view_t;
      typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
      typedef typename KCRS::values_type::non_const_type scalar_view_t;

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
      size_t p_max_nnz_per_row = Pview.origMatrix->getNodeMaxNumRowEntries();

      // Grab the  Kokkos::SparseCrsMatrices & inner stuff
      const KCRS & Amat = Aview.origMatrix->getLocalMatrix();
      const KCRS & Pmat = Pview.origMatrix->getLocalMatrix();
      const KCRS & Rmat = Rview.origMatrix->getLocalMatrix();

      c_lno_view_t Arowptr = Amat.graph.row_map, Prowptr = Pmat.graph.row_map,  Rrowptr = Rmat.graph.row_map;
      const lno_nnz_view_t Acolind = Amat.graph.entries, Pcolind = Pmat.graph.entries , Rcolind = Rmat.graph.entries;
      const scalar_view_t Avals = Amat.values, Pvals = Pmat.values, Rvals = Rmat.values;

      c_lno_view_t  Irowptr;
      lno_nnz_view_t  Icolind;
      scalar_view_t  Ivals;
      if(!Pview.importMatrix.is_null()) {
        Irowptr = Pview.importMatrix->getLocalMatrix().graph.row_map;
        Icolind = Pview.importMatrix->getLocalMatrix().graph.entries;
        Ivals   = Pview.importMatrix->getLocalMatrix().values;
        p_max_nnz_per_row = std::max(p_max_nnz_per_row,Pview.importMatrix->getNodeMaxNumRowEntries());
      }

#ifdef HAVE_TPETRA_MMM_TIMINGS
      RCP<TimeMonitor> MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix SerialCore - Compare"))));
#endif

      // Classic csr assembly (low memory edition)
      //
      // mfh 27 Sep 2016: Ac_estimate_nnz does not promise an upper bound.
      // The method loops over rows of R, and may resize after processing
      // each row.  Chris Siefert says that this reflects experience in
      // ML; for the non-threaded case, ML found it faster to spend less
      // effort on estimation and risk an occasional reallocation.
      size_t CSR_alloc = std::max(C_estimate_nnz(*Aview.origMatrix, *Pview.origMatrix), n);
      lno_view_t Crowptr(Kokkos::ViewAllocateWithoutInitializing("Crowptr"),m+1);
      lno_nnz_view_t Ccolind(Kokkos::ViewAllocateWithoutInitializing("Ccolind"),CSR_alloc);
      scalar_view_t Cvals(Kokkos::ViewAllocateWithoutInitializing("Cvals"),CSR_alloc);

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
      // For column index Aik in row i of A, Acol2Prow[Aik] tells
      // you whether the corresponding row of P belongs to P_local
      // ("orig") or P_remote ("Import").

      // For each row of R
      size_t nnz = 0, nnz_old = 0;
      for (size_t i = 0; i < m; i++) {
        // mfh 27 Sep 2016: m is the number of rows in the input matrix R
        // on the calling process.
        Crowptr[i] = nnz;

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


            if (Acol2Prow[l] != LO_INVALID) {
              // mfh 27 Sep 2016: If the entry of Acol2Prow
              // corresponding to the current entry of A is populated, then
              // the corresponding row of P is in P_local (i.e., it lives on
              // the calling process).

              // Local matrix
              size_t Pl = Teuchos::as<size_t>(Acol2Prow[l]);

              // mfh 27 Sep 2016: Go through all entries in that row of P_local.
              for (size_t jj = Prowptr[Pl]; jj < Prowptr[Pl+1]; jj++) {
                LO j = Pcolind[jj];
                LO Acj = Pcol2Accol[j];
                SC Plj = Pvals[jj];

                if (ac_status[Acj] == INVALID || ac_status[Acj] < nnz_old) {
#ifdef HAVE_TPETRA_DEBUG
                  // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                  TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Ccolind.size()),
                                             std::runtime_error,
                                             label << " ERROR, not enough memory allocated for matrix product. Allocated: " << Ccolind.extent(0) << std::endl);
#endif
                  // New entry
                  ac_status[Acj] = nnz;
                  Ccolind[nnz] = Acj;
                  Cvals[nnz] = Rik*Akl*Plj;
                  nnz++;
                } else {
                  Cvals[ac_status[Acj]] += Rik*Akl*Plj;
                }
              }
            } else {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is NOT populated (has
              // a flag "invalid" value), then the corresponding row of P is
              // in P_remote (i.e., it does not live on the calling process).

              // Remote matrix
              size_t Il = Teuchos::as<size_t>(Acol2PIrow[l]);
              for (size_t jj = Irowptr[Il]; jj < Irowptr[Il+1]; jj++) {
                LO j = Icolind[jj];
                LO Acj = PIcol2Accol[j];
                SC Plj = Ivals[jj];

                if (ac_status[Acj] == INVALID || ac_status[Acj] < nnz_old) {
#ifdef HAVE_TPETRA_DEBUG
                  // Ac_estimate_nnz() is probably not perfect yet. If this happens, we need to allocate more memory..
                  TEUCHOS_TEST_FOR_EXCEPTION(nnz >= Teuchos::as<size_t>(Ccolind.size()),
                                             std::runtime_error,
                                             label << " ERROR, not enough memory allocated for matrix product. Allocated: "  << Ccolind.extent(0) << std::endl);
#endif
                  // New entry
                  ac_status[Acj] = nnz;
                  Ccolind[nnz] = Acj;
                  Cvals[nnz] = Rik*Akl*Plj;
                  nnz++;
                } else {
                  Cvals[ac_status[Acj]] += Rik*Akl*Plj;
                }
              }
            }
          }
        }
        // Resize for next pass if needed
        if (nnz + n > CSR_alloc) {
          CSR_alloc *= 2;
          Kokkos::resize(Ccolind,CSR_alloc);
          Kokkos::resize(Cvals,CSR_alloc);
        }
        nnz_old = nnz;
      }

      Crowptr[m] = nnz;

      // Downward resize
      Kokkos::resize(Ccolind,nnz);
      Kokkos::resize(Cvals,nnz);

#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = Teuchos::null; MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix Final Sort"))));
#endif

      // Final sort & set of CRS arrays
      if (params.is_null() || params->get("sort entries",true))
        Import_Util::sortCrsEntries(Crowptr,Ccolind, Cvals);
      Ac.setAllValues(Crowptr, Ccolind, Cvals);

#ifdef HAVE_TPETRA_MMM_TIMINGS
     MM = Teuchos::null;  MM = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Newmatrix ESFC"))));
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
    // RAP Reuse Kernel wrappers (Default non-threaded version)
    // Computes R * A * P -> Ac using reuse Gustavson
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node,
             class LocalOrdinalViewType>
    void KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalOrdinalViewType>::mult_R_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Rview,
                                                                                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                                                                                           const LocalOrdinalViewType & Acol2Prow,
                                                                                                                           const LocalOrdinalViewType & Acol2PIrow,
                                                                                                                           const LocalOrdinalViewType & Pcol2Accol,
                                                                                                                           const LocalOrdinalViewType & PIcol2Accol,
                                                                                                                           CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                                                                                           Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                                                                                           const std::string& label,
                                                                                                                           const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Reuse SerialCore"))));
#endif

      using Teuchos::Array;
      using Teuchos::ArrayRCP;
      using Teuchos::ArrayView;
      using Teuchos::RCP;
      using Teuchos::rcp;

      // Lots and lots of typedefs
      typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type KCRS;
      typedef typename KCRS::StaticCrsGraphType graph_t;
      typedef typename graph_t::row_map_type::const_type c_lno_view_t;
      typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
      typedef typename KCRS::values_type::non_const_type scalar_view_t;

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
      size_t p_max_nnz_per_row = Pview.origMatrix->getNodeMaxNumRowEntries();

      // Grab the  Kokkos::SparseCrsMatrices & inner stuff
      const KCRS & Amat = Aview.origMatrix->getLocalMatrix();
      const KCRS & Pmat = Pview.origMatrix->getLocalMatrix();
      const KCRS & Rmat = Rview.origMatrix->getLocalMatrix();
      const KCRS & Cmat = Ac.getLocalMatrix();

      c_lno_view_t Arowptr = Amat.graph.row_map, Prowptr = Pmat.graph.row_map,  Rrowptr = Rmat.graph.row_map, Crowptr =  Cmat.graph.row_map;
      const lno_nnz_view_t Acolind = Amat.graph.entries, Pcolind = Pmat.graph.entries , Rcolind = Rmat.graph.entries, Ccolind = Cmat.graph.entries;
      const scalar_view_t Avals = Amat.values, Pvals = Pmat.values, Rvals = Rmat.values;
      scalar_view_t Cvals = Cmat.values;

      c_lno_view_t  Irowptr;
      lno_nnz_view_t  Icolind;
      scalar_view_t  Ivals;
      if(!Pview.importMatrix.is_null()) {
        Irowptr = Pview.importMatrix->getLocalMatrix().graph.row_map;
        Icolind = Pview.importMatrix->getLocalMatrix().graph.entries;
        Ivals   = Pview.importMatrix->getLocalMatrix().values;
        p_max_nnz_per_row = std::max(p_max_nnz_per_row,Pview.importMatrix->getNodeMaxNumRowEntries());
      }

#ifdef HAVE_TPETRA_MMM_TIMINGS
      RCP<TimeMonitor> MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Reuse SerialCore - Compare"))));
#endif

      // mfh 27 Sep 2016: The ac_status array is an implementation detail
      // of the local sparse matrix-matrix multiply routine.

      // The status array will contain the index into colind where this entry was last deposited.
      //   ac_status[i] <  nnz - not in the row yet
      //   ac_status[i] >= nnz - this is the entry where you can find the data
      // We start with this filled with INVALID's indicating that there are no entries yet.
      // Sadly, this complicates the code due to the fact that size_t's are unsigned.
      Array<size_t> ac_status(n, ST_INVALID);

      // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
      // routine.  The routine computes Ac := R * A * (P_local + P_remote).
      //
      // For column index Aik in row i of A, Acol2Prow[Aik] tells
      // you whether the corresponding row of P belongs to P_local
      // ("orig") or P_remote ("Import").

      // Necessary until following UVM host accesses are changed - for example Crowptr
      // Also probably needed in mult_R_A_P_newmatrix_kernel_wrapper - did not demonstrate this in test failure yet
      Kokkos::fence();

      // For each row of R
      size_t OLD_ip = 0, CSR_ip = 0;
      for (size_t i = 0; i < m; i++) {
        // First fill the c_status array w/ locations where we're allowed to
        // generate nonzeros for this row
        OLD_ip = Crowptr[i];
        CSR_ip = Crowptr[i+1];
        for (size_t k = OLD_ip; k < CSR_ip; k++) {
          ac_status[Ccolind[k]] = k;

          // Reset values in the row of C
          Cvals[k] = SC_ZERO;
        }

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


            if (Acol2Prow[l] != LO_INVALID) {
              // mfh 27 Sep 2016: If the entry of Acol2Prow
              // corresponding to the current entry of A is populated, then
              // the corresponding row of P is in P_local (i.e., it lives on
              // the calling process).

              // Local matrix
              size_t Pl = Teuchos::as<size_t>(Acol2Prow[l]);

              // mfh 27 Sep 2016: Go through all entries in that row of P_local.
              for (size_t jj = Prowptr[Pl]; jj < Prowptr[Pl+1]; jj++) {
                LO j = Pcolind[jj];
                LO Cij = Pcol2Accol[j];
                SC Plj = Pvals[jj];

                TEUCHOS_TEST_FOR_EXCEPTION(ac_status[Cij] < OLD_ip || ac_status[Cij] >= CSR_ip,
                                           std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
                                           "(c_status = " << ac_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

                Cvals[ac_status[Cij]] += Rik*Akl*Plj;
              }
            } else {
              // mfh 27 Sep 2016: If the entry of Acol2PRow
              // corresponding to the current entry of A is NOT populated (has
              // a flag "invalid" value), then the corresponding row of P is
              // in P_remote (i.e., it does not live on the calling process).

              // Remote matrix
              size_t Il = Teuchos::as<size_t>(Acol2PIrow[l]);
              for (size_t jj = Irowptr[Il]; jj < Irowptr[Il+1]; jj++) {
                LO j = Icolind[jj];
                LO Cij = PIcol2Accol[j];
                SC Plj = Ivals[jj];

                TEUCHOS_TEST_FOR_EXCEPTION(ac_status[Cij] < OLD_ip || ac_status[Cij] >= CSR_ip,
                                           std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
                                           "(c_status = " << ac_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

                Cvals[ac_status[Cij]] += Rik*Akl*Plj;
              }
            }
          }
        }
      }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  auto MM3 = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("RAP Reuse ESFC"))));
#endif

  Ac.fillComplete(Ac.getDomainMap(), Ac.getRangeMap());

}


    /*********************************************************************************************************/
    // PT_A_P NewMatrix Kernel wrappers (Default, general, non-threaded version)
    // Computes P.T * A * P -> Ac
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node,
             class LocalOrdinalViewType>
    void KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalOrdinalViewType>::mult_PT_A_P_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                                                                       const LocalOrdinalViewType & Acol2Prow,
                                                                                                       const LocalOrdinalViewType & Acol2PIrow,
                                                                                                       const LocalOrdinalViewType & Pcol2Accol,
                                                                                                       const LocalOrdinalViewType & PIcol2Accol,
                                                                                                       CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                                                                       Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                                                                       const std::string& label,
                                                                                                       const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      Teuchos::TimeMonitor MM(*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose")));
#endif

      // We don't need a kernel-level PTAP, we just transpose here
      typedef RowMatrixTransposer<Scalar,LocalOrdinal,GlobalOrdinal, Node>  transposer_type;
      transposer_type transposer (Pview.origMatrix,label+std::string("XP: "));

      using Teuchos::ParameterList;
      using Teuchos::RCP;
      RCP<ParameterList> transposeParams (new ParameterList);
      transposeParams->set ("sort", false);

      if (! params.is_null ()) {
        transposeParams->set ("compute global constants",
                              params->get ("compute global constants: temporaries",
                                           false));
      }
      RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ptrans =
        transposer.createTransposeLocal (transposeParams);
      CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node> Rview;
      Rview.origMatrix = Ptrans;

      mult_R_A_P_newmatrix_kernel_wrapper(Rview,Aview,Pview,Acol2Prow,Acol2PIrow,Pcol2Accol,PIcol2Accol,Ac,Acimport,label,params);
    }

   /*********************************************************************************************************/
    // PT_A_P Reuse Kernel wrappers (Default, general, non-threaded version)
    // Computes P.T * A * P -> Ac
    template<class Scalar,
             class LocalOrdinal,
             class GlobalOrdinal,
             class Node,
             class LocalOrdinalViewType>
    void KernelWrappers3<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalOrdinalViewType>::mult_PT_A_P_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                                       CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Pview,
                                                                                                       const LocalOrdinalViewType & Acol2Prow,
                                                                                                       const LocalOrdinalViewType & Acol2PIrow,
                                                                                                       const LocalOrdinalViewType & Pcol2Accol,
                                                                                                       const LocalOrdinalViewType & PIcol2Accol,
                                                                                                       CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Ac,
                                                                                                       Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Acimport,
                                                                                                       const std::string& label,
                                                                                                       const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
      using Teuchos::TimeMonitor;
      Teuchos::TimeMonitor MM(*TimeMonitor::getNewTimer(prefix_mmm + std::string("PTAP local transpose")));
#endif

      // We don't need a kernel-level PTAP, we just transpose here
      typedef RowMatrixTransposer<Scalar,LocalOrdinal,GlobalOrdinal, Node>  transposer_type;
      transposer_type transposer (Pview.origMatrix,label+std::string("XP: "));

      using Teuchos::ParameterList;
      using Teuchos::RCP;
      RCP<ParameterList> transposeParams (new ParameterList);
      transposeParams->set ("sort", false);

      if (! params.is_null ()) {
        transposeParams->set ("compute global constants",
                              params->get ("compute global constants: temporaries",
                                           false));
      }
      RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ptrans =
        transposer.createTransposeLocal (transposeParams);
      CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node> Rview;
      Rview.origMatrix = Ptrans;

      mult_R_A_P_reuse_kernel_wrapper(Rview,Aview,Pview,Acol2Prow,Acol2PIrow,Pcol2Accol,PIcol2Accol,Ac,Acimport,label,params);
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

      using Teuchos::ParameterList;
      RCP<ParameterList> transposeParams (new ParameterList);
      transposeParams->set ("sort", false);  
      if (! params.is_null ()) {
        transposeParams->set ("compute global constants",
                              params->get ("compute global constants: temporaries",
                                           false));
      }
      RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ptrans =
        transposer.createTransposeLocal (transposeParams);

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
