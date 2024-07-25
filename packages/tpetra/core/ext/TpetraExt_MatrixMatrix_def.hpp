// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MATRIXMATRIX_DEF_HPP
#define TPETRA_MATRIXMATRIX_DEF_HPP
#include "TpetraExt_MatrixMatrix_fwd.hpp"
#include "KokkosSparse_Utils.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "TpetraExt_MatrixMatrix_decl.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Array.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "TpetraExt_MMHelpers_def.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import_Util.hpp"
#include "Tpetra_Import_Util2.hpp"
#include <algorithm>
#include <type_traits>
#include "Teuchos_FancyOStream.hpp"

#include "TpetraExt_MatrixMatrix_ExtraKernels_def.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"

#include "KokkosSparse_spgemm.hpp"
#include "KokkosSparse_spadd.hpp"
#include "Kokkos_Bitset.hpp"

#include <MatrixMarket_Tpetra.hpp>

/*! \file TpetraExt_MatrixMatrix_def.hpp

    The implementations for the members of class Tpetra::MatrixMatrixMultiply and related non-member constructors.
 */



/*********************************************************************************************************/
// Include the architecture-specific kernel partial specializations here
// NOTE: This needs to be outside all namespaces
#include "TpetraExt_MatrixMatrix_OpenMP.hpp"
#include "TpetraExt_MatrixMatrix_Cuda.hpp"
#include "TpetraExt_MatrixMatrix_HIP.hpp"
#include "TpetraExt_MatrixMatrix_SYCL.hpp"

namespace Tpetra {

namespace MatrixMatrix{

//
// This method forms the matrix-matrix product C = op(A) * op(B), where
// op(A) == A   if transposeA is false,
// op(A) == A^T if transposeA is true,
// and similarly for op(B).
//
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Multiply(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  bool transposeB,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  bool call_FillComplete_on_result,
  const std::string& label,
  const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Scalar                            SC;
  typedef LocalOrdinal                      LO;
  typedef GlobalOrdinal                     GO;
  typedef Node                              NO;
  typedef CrsMatrix<SC,LO,GO,NO>            crs_matrix_type;
  typedef Import<LO,GO,NO>                  import_type;
  typedef CrsMatrixStruct<SC,LO,GO,NO>      crs_matrix_struct_type;
  typedef Map<LO,GO,NO>                     map_type;
  typedef RowMatrixTransposer<SC,LO,GO,NO>  transposer_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  //MM is used to time setup, and then multiply.
  
  RCP<TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM All Setup"))));
#endif

  const std::string prefix = "TpetraExt::MatrixMatrix::Multiply(): ";

  // TEUCHOS_FUNC_TIME_MONITOR_DIFF("My Matrix Mult", mmm_multiply);

  // The input matrices A and B must both be fillComplete.
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error, prefix << "Matrix A is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(), std::runtime_error, prefix << "Matrix B is not fill complete.");

  // If transposeA is true, then Aprime will be the transpose of A
  // (computed explicitly via RowMatrixTransposer).  Otherwise, Aprime
  // will just be a pointer to A.
  RCP<const crs_matrix_type> Aprime = null;
  // If transposeB is true, then Bprime will be the transpose of B
  // (computed explicitly via RowMatrixTransposer).  Otherwise, Bprime
  // will just be a pointer to B.
  RCP<const crs_matrix_type> Bprime = null;

  // Is this a "clean" matrix?
  //
  // mfh 27 Sep 2016: Historically, if Epetra_CrsMatrix was neither
  // locally nor globally indexed, then it was empty.  I don't like
  // this, because the most straightforward implementation presumes
  // lazy allocation of indices.  However, historical precedent
  // demands that we keep around this predicate as a way to test
  // whether the matrix is empty.
  const bool newFlag = !C.getGraph()->isLocallyIndexed() && !C.getGraph()->isGloballyIndexed();

  bool use_optimized_ATB = false;
  if (transposeA && !transposeB && call_FillComplete_on_result && newFlag)
    use_optimized_ATB = true;

#ifdef USE_OLD_TRANSPOSE // NOTE: For Grey Ballard's use.  Remove this later.
  use_optimized_ATB = false;
#endif

  using Teuchos::ParameterList;
  RCP<ParameterList> transposeParams (new ParameterList);
  transposeParams->set ("sort", true); // Kokkos Kernels spgemm requires inputs to be sorted

  if (!use_optimized_ATB && transposeA) {
    transposer_type transposer (rcpFromRef (A));
    Aprime = transposer.createTranspose (transposeParams);
  }
  else {
    Aprime = rcpFromRef(A);
  }

  if (transposeB) {
    transposer_type transposer (rcpFromRef (B));
    Bprime = transposer.createTranspose (transposeParams);
  }
  else {
    Bprime = rcpFromRef(B);
  }

  // Check size compatibility
  global_size_t numACols = A.getDomainMap()->getGlobalNumElements();
  global_size_t numBCols = B.getDomainMap()->getGlobalNumElements();
  global_size_t Aouter   = transposeA ? numACols             : A.getGlobalNumRows();
  global_size_t Bouter   = transposeB ? B.getGlobalNumRows() : numBCols;
  global_size_t Ainner   = transposeA ? A.getGlobalNumRows() : numACols;
  global_size_t Binner   = transposeB ? numBCols             : B.getGlobalNumRows();
  TEUCHOS_TEST_FOR_EXCEPTION(Ainner != Binner, std::runtime_error,
    prefix << "ERROR, inner dimensions of op(A) and op(B) "
    "must match for matrix-matrix product. op(A) is "
    << Aouter << "x" << Ainner << ", op(B) is "<< Binner << "x" << Bouter);

  // The result matrix C must at least have a row-map that reflects the correct
  // row-size. Don't check the number of columns because rectangular matrices
  // which were constructed with only one map can still end up having the
  // correct capacity and dimensions when filled.
  TEUCHOS_TEST_FOR_EXCEPTION(Aouter > C.getGlobalNumRows(), std::runtime_error,
    prefix << "ERROR, dimensions of result C must "
    "match dimensions of op(A) * op(B). C has " << C.getGlobalNumRows()
     << " rows, should have at least " << Aouter << std::endl);

  // It doesn't matter whether C is already Filled or not. If it is already
  // Filled, it must have space allocated for the positions that will be
  // referenced in forming C = op(A)*op(B). If it doesn't have enough space,
  // we'll error out later when trying to store result values.

  // CGB: However, matrix must be in active-fill
  if (!C.isFillActive()) C.resumeFill();

  // We're going to need to import remotely-owned sections of A and/or B if
  // more than one processor is performing this run, depending on the scenario.
  int numProcs = A.getComm()->getSize();

  // Declare a couple of structs that will be used to hold views of the data
  // of A and B, to be used for fast access during the matrix-multiplication.
  crs_matrix_struct_type Aview;
  crs_matrix_struct_type Bview;

  RCP<const map_type> targetMap_A = Aprime->getRowMap();
  RCP<const map_type> targetMap_B = Bprime->getRowMap();

#ifdef HAVE_TPETRA_MMM_TIMINGS
  {
  TimeMonitor MM_importExtract(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM All I&X")));
#endif

  // Now import any needed remote rows and populate the Aview struct
  // NOTE: We assert that an import isn't needed --- since we do the transpose
  // above to handle that.
  if (!use_optimized_ATB) {
    RCP<const import_type> dummyImporter;
    MMdetails::import_and_extract_views(*Aprime, targetMap_A, Aview, dummyImporter, true, label, params);
  }

  // We will also need local access to all rows of B that correspond to the
  // column-map of op(A).
  if (numProcs > 1)
    targetMap_B = Aprime->getColMap();

  // Import any needed remote rows and populate the Bview struct.
  if (!use_optimized_ATB)
    MMdetails::import_and_extract_views(*Bprime, targetMap_B, Bview, Aprime->getGraph()->getImporter(), Aprime->getGraph()->getImporter().is_null(), label, params);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  } //stop MM_importExtract here
  //stop the setup timer, and start the multiply timer
  MM = Teuchos::null; MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM All Multiply"))));
#endif

  // Call the appropriate method to perform the actual multiplication.
  if (use_optimized_ATB) {
    MMdetails::mult_AT_B_newmatrix(A, B, C, label,params);

  } else if (call_FillComplete_on_result && newFlag) {
    MMdetails::mult_A_B_newmatrix(Aview, Bview, C, label,params);

  } else if (call_FillComplete_on_result) {
    MMdetails::mult_A_B_reuse(Aview, Bview, C, label,params);

  } else {
    // mfh 27 Sep 2016: Is this the "slow" case?  This
    // "CrsWrapper_CrsMatrix" thing could perhaps be made to support
    // thread-parallel inserts, but that may take some effort.
    CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsmat(C);

    MMdetails::mult_A_B(Aview, Bview, crsmat, label,params);
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  TimeMonitor MM4(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM All FillComplete")));
#endif
  if (call_FillComplete_on_result && !C.isFillComplete()) {
    // We'll call FillComplete on the C matrix before we exit, and give it a
    // domain-map and a range-map.
    // The domain-map will be the domain-map of B, unless
    // op(B)==transpose(B), in which case the range-map of B will be used.
    // The range-map will be the range-map of A, unless op(A)==transpose(A),
    // in which case the domain-map of A will be used.
    C.fillComplete(Bprime->getDomainMap(), Aprime->getRangeMap());
  }
}

//
// This method forms the matrix-matrix product C = op(A) * op(B), where
// op(A) == A and similarly for op(B). op(A) = A^T is not yet implemented.
//
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Multiply(
  const Teuchos::RCP<const BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
  bool transposeA,
  const Teuchos::RCP<const BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& B,
  bool transposeB,
  Teuchos::RCP<BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& C,
  const std::string& label)
{
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Scalar                            SC;
  typedef LocalOrdinal                      LO;
  typedef GlobalOrdinal                     GO;
  typedef Node                              NO;
  typedef BlockCrsMatrixStruct<SC,LO,GO,NO> blockcrs_matrix_struct_type;
  typedef Map<LO,GO,NO>                     map_type;
  typedef Import<LO,GO,NO>                  import_type;

  std::string prefix = std::string("TpetraExt ") + label + std::string(": ");

  TEUCHOS_TEST_FOR_EXCEPTION(transposeA==true, std::runtime_error, prefix << "Matrix A cannot be transposed.");
  TEUCHOS_TEST_FOR_EXCEPTION(transposeB==true, std::runtime_error, prefix << "Matrix B cannot be transposed.");

  // Check size compatibility
  global_size_t numACols = A->getGlobalNumCols();
  global_size_t numBCols = B->getGlobalNumCols();
  global_size_t numARows = A->getGlobalNumRows();
  global_size_t numBRows = B->getGlobalNumRows();

  global_size_t Aouter   = numARows;
  global_size_t Bouter   = numBCols;
  global_size_t Ainner   = numACols;
  global_size_t Binner   = numBRows;
  TEUCHOS_TEST_FOR_EXCEPTION(Ainner != Binner, std::runtime_error,
    prefix << "ERROR, inner dimensions of op(A) and op(B) "
    "must match for matrix-matrix product. op(A) is "
    << Aouter << "x" << Ainner << ", op(B) is "<< Binner << "x" << Bouter);

  // We're going to need to import remotely-owned sections of A and/or B if
  // more than one processor is performing this run, depending on the scenario.
  int numProcs = A->getComm()->getSize();

  const LO blocksize = A->getBlockSize();
  TEUCHOS_TEST_FOR_EXCEPTION(blocksize != B->getBlockSize(), std::runtime_error,
    prefix << "ERROR, Blocksizes do not match. A.blocksize = " <<
              blocksize << ", B.blocksize = " << B->getBlockSize() );
  
  // Declare a couple of structs that will be used to hold views of the data
  // of A and B, to be used for fast access during the matrix-multiplication.
  blockcrs_matrix_struct_type Aview(blocksize);
  blockcrs_matrix_struct_type Bview(blocksize);

  RCP<const map_type> targetMap_A = A->getRowMap();
  RCP<const map_type> targetMap_B = B->getRowMap();

  // Populate the Aview struct. No remotes are needed.
  RCP<const import_type> dummyImporter;
  MMdetails::import_and_extract_views(*A, targetMap_A, Aview, dummyImporter, true);

  // We will also need local access to all rows of B that correspond to the
  // column-map of op(A).
  if (numProcs > 1)
    targetMap_B = A->getColMap();

  // Import any needed remote rows and populate the Bview struct.
  MMdetails::import_and_extract_views(*B, targetMap_B, Bview, A->getGraph()->getImporter(),
                                      A->getGraph()->getImporter().is_null());

  // Call the appropriate method to perform the actual multiplication.
  MMdetails::mult_A_B_newmatrix(Aview, Bview, C);
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Jacobi(Scalar omega,
            const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
            const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
            const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
            CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
            bool call_FillComplete_on_result,
                 const std::string& label,
            const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::RCP;
  typedef Scalar                            SC;
  typedef LocalOrdinal                      LO;
  typedef GlobalOrdinal                     GO;
  typedef Node                              NO;
  typedef Import<LO,GO,NO>                  import_type;
  typedef CrsMatrixStruct<SC,LO,GO,NO>      crs_matrix_struct_type;
  typedef Map<LO,GO,NO>                     map_type;
  typedef CrsMatrix<SC,LO,GO,NO>            crs_matrix_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ")+ label + std::string(": ");
  using Teuchos::TimeMonitor;
  TimeMonitor MM(*TimeMonitor::getNewTimer(prefix_mmm+std::string("Jacobi All Setup")));
#endif

  const std::string prefix = "TpetraExt::MatrixMatrix::Jacobi(): ";

  // A and B should already be Filled.
  // Should we go ahead and call FillComplete() on them if necessary or error
  // out? For now, we choose to error out.
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(),  std::runtime_error, prefix << "Matrix A is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION(!B.isFillComplete(),  std::runtime_error, prefix << "Matrix B is not fill complete.");

  RCP<const crs_matrix_type> Aprime = rcpFromRef(A);
  RCP<const crs_matrix_type> Bprime = rcpFromRef(B);

  // Now check size compatibility
  global_size_t numACols = A.getDomainMap()->getGlobalNumElements();
  global_size_t numBCols = B.getDomainMap()->getGlobalNumElements();
  global_size_t Aouter   = A.getGlobalNumRows();
  global_size_t Bouter   = numBCols;
  global_size_t Ainner   = numACols;
  global_size_t Binner   = B.getGlobalNumRows();
  TEUCHOS_TEST_FOR_EXCEPTION(Ainner != Binner, std::runtime_error,
    prefix << "ERROR, inner dimensions of op(A) and op(B) "
    "must match for matrix-matrix product. op(A) is "
    << Aouter << "x" << Ainner << ", op(B) is "<< Binner << "x" << Bouter);

  // The result matrix C must at least have a row-map that reflects the correct
  // row-size. Don't check the number of columns because rectangular matrices
  // which were constructed with only one map can still end up having the
  // correct capacity and dimensions when filled.
  TEUCHOS_TEST_FOR_EXCEPTION(Aouter > C.getGlobalNumRows(), std::runtime_error,
    prefix << "ERROR, dimensions of result C must "
    "match dimensions of op(A) * op(B). C has "<< C.getGlobalNumRows()
     << " rows, should have at least "<< Aouter << std::endl);

  // It doesn't matter whether C is already Filled or not. If it is already
  // Filled, it must have space allocated for the positions that will be
  // referenced in forming C = op(A)*op(B). If it doesn't have enough space,
  // we'll error out later when trying to store result values.

  // CGB: However, matrix must be in active-fill
  TEUCHOS_TEST_FOR_EXCEPT( C.isFillActive() == false );

  // We're going to need to import remotely-owned sections of A and/or B if
  // more than one processor is performing this run, depending on the scenario.
  int numProcs = A.getComm()->getSize();

  // Declare a couple of structs that will be used to hold views of the data of
  // A and B, to be used for fast access during the matrix-multiplication.
  crs_matrix_struct_type Aview;
  crs_matrix_struct_type Bview;

  RCP<const map_type> targetMap_A = Aprime->getRowMap();
  RCP<const map_type> targetMap_B = Bprime->getRowMap();

#ifdef HAVE_TPETRA_MMM_TIMINGS
  TimeMonitor MM2(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi All I&X")));
  {
#endif

  // Enable globalConstants by default
  // NOTE: the I&X routine sticks an importer on the paramlist as output, so we have to use a unique guy here
  RCP<Teuchos::ParameterList> importParams1 = Teuchos::rcp(new Teuchos::ParameterList);
  if(!params.is_null()) {
      importParams1->set("compute global constants",params->get("compute global constants: temporaries",false));
      int mm_optimization_core_count=0;
      auto slist = params->sublist("matrixmatrix: kernel params",false);
      mm_optimization_core_count = slist.get("MM_TAFC_OptimizationCoreCount",::Tpetra::Details::Behavior::TAFC_OptimizationCoreCount ());
      int mm_optimization_core_count2 = params->get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
      if(mm_optimization_core_count2<mm_optimization_core_count) mm_optimization_core_count=mm_optimization_core_count2;
      bool isMM = slist.get("isMatrixMatrix_TransferAndFillComplete",false);
      bool overrideAllreduce = slist.get("MM_TAFC_OverrideAllreduceCheck",false);
      auto & ip1slist = importParams1->sublist("matrixmatrix: kernel params",false);
      ip1slist.set("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
      ip1slist.set("isMatrixMatrix_TransferAndFillComplete",isMM);
      ip1slist.set("MM_TAFC_OverrideAllreduceCheck",overrideAllreduce);
  }

  //Now import any needed remote rows and populate the Aview struct.
  RCP<const import_type> dummyImporter;
  MMdetails::import_and_extract_views(*Aprime, targetMap_A, Aview, dummyImporter, true, label,importParams1);

  // We will also need local access to all rows of B that correspond to the
  // column-map of op(A).
  if (numProcs > 1)
    targetMap_B = Aprime->getColMap();

  // Now import any needed remote rows and populate the Bview struct.
  // Enable globalConstants by default
  // NOTE: the I&X routine sticks an importer on the paramlist as output, so we have to use a unique guy here
  RCP<Teuchos::ParameterList> importParams2 = Teuchos::rcp(new Teuchos::ParameterList);
  if(!params.is_null()) {
      importParams2->set("compute global constants",params->get("compute global constants: temporaries",false));

      auto slist = params->sublist("matrixmatrix: kernel params",false);
      int mm_optimization_core_count = slist.get("MM_TAFC_OptimizationCoreCount",::Tpetra::Details::Behavior::TAFC_OptimizationCoreCount () );
      bool isMM = slist.get("isMatrixMatrix_TransferAndFillComplete",false);
      bool overrideAllreduce = slist.get("MM_TAFC_OverrideAllreduceCheck",false);
      auto & ip2slist = importParams2->sublist("matrixmatrix: kernel params",false);
      int mm_optimization_core_count2 = params->get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
      if(mm_optimization_core_count2<mm_optimization_core_count) mm_optimization_core_count=mm_optimization_core_count2;
      ip2slist.set("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
      ip2slist.set("isMatrixMatrix_TransferAndFillComplete",isMM);
      ip2slist.set("MM_TAFC_OverrideAllreduceCheck",overrideAllreduce);
  }

  MMdetails::import_and_extract_views(*Bprime, targetMap_B, Bview, Aprime->getGraph()->getImporter(), Aprime->getGraph()->getImporter().is_null(), label,importParams2);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  }
  TimeMonitor MM3(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi All Multiply")));
#endif

  // Now call the appropriate method to perform the actual multiplication.
  CrsWrapper_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsmat(C);

  // Is this a "clean" matrix
  bool newFlag = !C.getGraph()->isLocallyIndexed() && !C.getGraph()->isGloballyIndexed();

  if (call_FillComplete_on_result && newFlag) {
    MMdetails::jacobi_A_B_newmatrix(omega, Dinv, Aview, Bview, C, label, params);

  } else if (call_FillComplete_on_result) {
    MMdetails::jacobi_A_B_reuse(omega, Dinv, Aview, Bview, C, label, params);

  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "jacobi_A_B_general not implemented");
  }

  if(!params.is_null()) {
    bool removeZeroEntries = params->get("remove zeros", false);
    if (removeZeroEntries) {
      typedef Teuchos::ScalarTraits<Scalar> STS;
      typename STS::magnitudeType threshold = params->get("remove zeros threshold", STS::magnitude(STS::zero()));
      removeCrsMatrixZeros(C, threshold);
    }
  }
}


template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  Scalar scalarA,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  Scalar scalarB )
{
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::null;
  typedef Scalar                            SC;
  typedef LocalOrdinal                      LO;
  typedef GlobalOrdinal                     GO;
  typedef Node                              NO;
  typedef CrsMatrix<SC,LO,GO,NO>            crs_matrix_type;
  typedef RowMatrixTransposer<SC,LO,GO,NO>  transposer_type;

  const std::string prefix = "TpetraExt::MatrixMatrix::Add(): ";

  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::runtime_error,
    prefix << "ERROR, input matrix A.isFillComplete() is false; it is required to be true. "
    "(Result matrix B is not required to be isFillComplete()).");
  TEUCHOS_TEST_FOR_EXCEPTION(B.isFillComplete() , std::runtime_error,
    prefix << "ERROR, input matrix B must not be fill complete!");
  TEUCHOS_TEST_FOR_EXCEPTION(B.isStaticGraph() , std::runtime_error,
    prefix << "ERROR, input matrix B must not have static graph!");
  TEUCHOS_TEST_FOR_EXCEPTION(B.isLocallyIndexed() , std::runtime_error,
    prefix << "ERROR, input matrix B must not be locally indexed!");

  using Teuchos::ParameterList;
  RCP<ParameterList> transposeParams (new ParameterList);
  transposeParams->set ("sort", false);

  RCP<const crs_matrix_type> Aprime = null;
  if (transposeA) {
    transposer_type transposer (rcpFromRef (A));
    Aprime = transposer.createTranspose (transposeParams);
  }
  else {
    Aprime = rcpFromRef(A);
  }

  size_t a_numEntries;
  typename crs_matrix_type::nonconst_global_inds_host_view_type a_inds("a_inds",A.getLocalMaxNumRowEntries());
  typename crs_matrix_type::nonconst_values_host_view_type a_vals("a_vals",A.getLocalMaxNumRowEntries());
  GO row;

  if (scalarB != Teuchos::ScalarTraits<SC>::one())
    B.scale(scalarB);

  size_t numMyRows = B.getLocalNumRows();
  if (scalarA != Teuchos::ScalarTraits<SC>::zero()) {
    for (LO i = 0; (size_t)i < numMyRows; ++i) {
      row = B.getRowMap()->getGlobalElement(i);
      Aprime->getGlobalRowCopy(row, a_inds, a_vals, a_numEntries);

      if (scalarA != Teuchos::ScalarTraits<SC>::one()) {
        for (size_t j = 0; j < a_numEntries; ++j)
          a_vals[j] *= scalarA;
      }
      B.insertGlobalValues(row,  a_numEntries, reinterpret_cast<Scalar *>(a_vals.data()), a_inds.data());
    }
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
add (const Scalar& alpha,
     const bool transposeA,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
     const Scalar& beta,
     const bool transposeB,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap,
     const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  typedef CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>  crs_matrix_type;
  if(!params.is_null())
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        params->isParameter("Call fillComplete") && !params->get<bool>("Call fillComplete"),
        std::invalid_argument,
        "Tpetra::MatrixMatrix::add(): this version of add() always calls fillComplete\n"
        "on the result, but you explicitly set 'Call fillComplete' = false in the parameter list. Don't set this explicitly.");
    params->set("Call fillComplete", true);
  }
  //If transposeB, must compute B's explicit transpose to
  //get the correct row map for C.
  RCP<const crs_matrix_type> Brcp = rcpFromRef(B);
  if(transposeB)
  {
    RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node> transposer(Brcp);
    Brcp = transposer.createTranspose();
  }
  //Check that A,B are fillComplete before getting B's column map
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete () || ! Brcp->isFillComplete (), std::invalid_argument,
     "TpetraExt::MatrixMatrix::add(): A and B must both be fill complete.");
  RCP<crs_matrix_type> C = rcp(new crs_matrix_type(Brcp->getRowMap(), 0));
  //this version of add() always fill completes the result, no matter what is in params on input
  add(alpha, transposeA, A, beta, false, *Brcp, *C, domainMap, rangeMap, params);
  return C;
}

//This functor does the same thing as CrsGraph::convertColumnIndicesFromGlobalToLocal,
//but since the spadd() output is always packed there is no need for a separate
//numRowEntries here.
//
template<class LO, class GO, class LOView, class GOView, class LocalMap>
struct ConvertGlobalToLocalFunctor
{
  ConvertGlobalToLocalFunctor(LOView& lids_, const GOView& gids_, const LocalMap localColMap_)
    : lids(lids_), gids(gids_), localColMap(localColMap_)
  {}

  KOKKOS_FUNCTION void operator() (const GO i) const
  {
    lids(i) = localColMap.getLocalElement(gids(i));
  }

  LOView lids;
  const GOView gids;
  const LocalMap localColMap;
};

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
add (const Scalar& alpha,
     const bool transposeA,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
     const Scalar& beta,
     const bool transposeB,
     const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap,
     const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap,
     const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::TimeMonitor;
  using SC = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NO = Node;
  using crs_matrix_type = CrsMatrix<SC,LO,GO,NO>;
  using crs_graph_type  = CrsGraph<LO,GO,NO>;
  using map_type        = Map<LO,GO,NO>;
  using transposer_type = RowMatrixTransposer<SC,LO,GO,NO>;
  using import_type     = Import<LO,GO,NO>;
  using export_type     = Export<LO,GO,NO>;
  using exec_space      = typename crs_graph_type::execution_space;
  using AddKern         = MMdetails::AddKernels<SC,LO,GO,NO>;
  const char* prefix_mmm = "TpetraExt::MatrixMatrix::add: ";
  constexpr bool debug = false;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("transpose"))));
#endif

  if (debug) {
    std::ostringstream os;
    os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
       << "TpetraExt::MatrixMatrix::add" << std::endl;
    std::cerr << os.str ();
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (C.isLocallyIndexed() || C.isGloballyIndexed(), std::invalid_argument,
     prefix_mmm << "C must be a 'new' matrix (neither locally nor globally indexed).");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete () || ! B.isFillComplete (), std::invalid_argument,
     prefix_mmm << "A and B must both be fill complete.");
#ifdef HAVE_TPETRA_DEBUG
  // The matrices don't have domain or range Maps unless they are fill complete.
  if (A.isFillComplete () && B.isFillComplete ()) {
    const bool domainMapsSame =
      (! transposeA && ! transposeB &&
       ! A.getDomainMap()->locallySameAs (*B.getDomainMap ())) ||
      (! transposeA &&   transposeB &&
       ! A.getDomainMap()->isSameAs (*B.getRangeMap  ())) ||
      (  transposeA && ! transposeB &&
       ! A.getRangeMap ()->isSameAs (*B.getDomainMap ()));
    TEUCHOS_TEST_FOR_EXCEPTION(domainMapsSame, std::invalid_argument,
      prefix_mmm << "The domain Maps of Op(A) and Op(B) are not the same.");

    const bool rangeMapsSame =
      (! transposeA && ! transposeB &&
       ! A.getRangeMap ()->isSameAs (*B.getRangeMap ())) ||
      (! transposeA &&   transposeB &&
       ! A.getRangeMap ()->isSameAs (*B.getDomainMap())) ||
      (  transposeA && ! transposeB &&
       ! A.getDomainMap()->isSameAs (*B.getRangeMap ()));
    TEUCHOS_TEST_FOR_EXCEPTION(rangeMapsSame, std::invalid_argument,
      prefix_mmm << "The range Maps of Op(A) and Op(B) are not the same.");
  }
#endif // HAVE_TPETRA_DEBUG

  using Teuchos::ParameterList;
  // Form the explicit transpose of A if necessary.
  RCP<const crs_matrix_type> Aprime = rcpFromRef(A);
  if (transposeA) {
    transposer_type transposer (Aprime);
    Aprime = transposer.createTranspose ();
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION
    (Aprime.is_null (), std::logic_error,
     prefix_mmm << "Failed to compute Op(A). "
     "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  // Form the explicit transpose of B if necessary.
  RCP<const crs_matrix_type> Bprime = rcpFromRef(B);
  if (transposeB) {
    if (debug) {
      std::ostringstream os;
      os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
         << "Form explicit xpose of B" << std::endl;
      std::cerr << os.str ();
    }
    transposer_type transposer (Bprime);
    Bprime = transposer.createTranspose ();
  }
#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Bprime.is_null (), std::logic_error,
    prefix_mmm << "Failed to compute Op(B). Please report this bug to the Tpetra developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    !Aprime->isFillComplete () || !Bprime->isFillComplete (), std::invalid_argument,
    prefix_mmm << "Aprime and Bprime must both be fill complete.  "
    "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  RCP<const map_type> CDomainMap = domainMap;
  RCP<const map_type> CRangeMap = rangeMap;
  if(CDomainMap.is_null())
  {
    CDomainMap = Bprime->getDomainMap();
  }
  if(CRangeMap.is_null())
  {
    CRangeMap = Bprime->getRangeMap();
  }
  assert(!(CDomainMap.is_null()));
  assert(!(CRangeMap.is_null()));
  typedef typename AddKern::values_array values_array;
  typedef typename AddKern::row_ptrs_array row_ptrs_array;
  typedef typename AddKern::col_inds_array col_inds_array;
  bool AGraphSorted = Aprime->getCrsGraph()->isSorted();
  bool BGraphSorted = Bprime->getCrsGraph()->isSorted();
  values_array vals;
  row_ptrs_array rowptrs;
  col_inds_array colinds;
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("rowmap check/import"))));
#endif
  if(!(Aprime->getRowMap()->isSameAs(*(Bprime->getRowMap()))))
  {
    //import Aprime into Bprime's row map so the local matrices have same # of rows
    auto import = rcp(new import_type(Aprime->getRowMap(), Bprime->getRowMap()));
    // cbl do not set
    // parameterlist "isMatrixMatrix_TransferAndFillComplete" true here as
    // this import _may_ take the form of a transfer. In practice it would be unlikely,
    // but the general case is not so forgiving.
    Aprime = importAndFillCompleteCrsMatrix<crs_matrix_type>(Aprime, *import, Bprime->getDomainMap(), Bprime->getRangeMap());
  }
  bool matchingColMaps = Aprime->getColMap()->isSameAs(*(Bprime->getColMap()));
  bool sorted = AGraphSorted && BGraphSorted;
  RCP<const import_type> Cimport = Teuchos::null;
  RCP<export_type> Cexport = Teuchos::null;
  bool doFillComplete = true;
  if(Teuchos::nonnull(params) && params->isParameter("Call fillComplete"))
  {
    doFillComplete = params->get<bool>("Call fillComplete");
  }
  auto Alocal = Aprime->getLocalMatrixDevice();
  auto Blocal = Bprime->getLocalMatrixDevice();
  LO numLocalRows = Alocal.numRows();
  if(numLocalRows == 0)
  {
    //KokkosKernels spadd assumes rowptrs.extent(0) + 1 == nrows,
    //but an empty Tpetra matrix is allowed to have rowptrs.extent(0) == 0.
    //Handle this case now
    //(without interfering with collective operations, since it's possible for
    //some ranks to have 0 local rows and others not).
    rowptrs = row_ptrs_array("C rowptrs", 0);
  }
  auto Acolmap = Aprime->getColMap();
  auto Bcolmap = Bprime->getColMap();
  if(!matchingColMaps)
  {
    using global_col_inds_array = typename AddKern::global_col_inds_array;
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("mismatched col map full kernel"))));
#endif
    //use kernel that converts col indices in both A and B to common domain map before adding
    auto AlocalColmap = Acolmap->getLocalMap();
    auto BlocalColmap = Bcolmap->getLocalMap();
    global_col_inds_array globalColinds;
    if (debug) {
      std::ostringstream os;
      os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
         << "Call AddKern::convertToGlobalAndAdd(...)" << std::endl;
      std::cerr << os.str ();
    }
    AddKern::convertToGlobalAndAdd(
      Alocal, alpha, Blocal, beta, AlocalColmap, BlocalColmap,
      vals, rowptrs, globalColinds);
    if (debug) {
      std::ostringstream os;
      os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
         << "Finished AddKern::convertToGlobalAndAdd(...)" << std::endl;
      std::cerr << os.str ();
    }
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Constructing graph"))));
#endif
    RCP<const map_type> CcolMap;
    Tpetra::Details::makeColMap<LocalOrdinal, GlobalOrdinal, Node>
      (CcolMap, CDomainMap, globalColinds);
    C.replaceColMap(CcolMap);
    col_inds_array localColinds("C colinds", globalColinds.extent(0));
    Kokkos::parallel_for(Kokkos::RangePolicy<exec_space>(0, globalColinds.extent(0)),
        ConvertGlobalToLocalFunctor<LocalOrdinal, GlobalOrdinal,
                              col_inds_array, global_col_inds_array,
                              typename map_type::local_map_type>
        (localColinds, globalColinds, CcolMap->getLocalMap()));
    KokkosSparse::sort_crs_matrix<exec_space, row_ptrs_array, col_inds_array, values_array>(rowptrs, localColinds, vals);
    C.setAllValues(rowptrs, localColinds, vals);
    C.fillComplete(CDomainMap, CRangeMap, params);
    if(!doFillComplete)
      C.resumeFill();
  }
  else
  {
    //Aprime, Bprime and C all have the same column maps
    auto Avals = Alocal.values;
    auto Bvals = Blocal.values;
    auto Arowptrs = Alocal.graph.row_map;
    auto Browptrs = Blocal.graph.row_map;
    auto Acolinds = Alocal.graph.entries;
    auto Bcolinds = Blocal.graph.entries;
    if(sorted)
    {
      //use sorted kernel
#ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = Teuchos::null;
        MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("sorted entries full kernel"))));
#endif
      if (debug) {
        std::ostringstream os;
        os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
           << "Call AddKern::addSorted(...)" << std::endl;
        std::cerr << os.str ();
      }
#if KOKKOSKERNELS_VERSION >= 40299
      AddKern::addSorted(Avals, Arowptrs, Acolinds, alpha, Bvals, Browptrs, Bcolinds, beta, Aprime->getGlobalNumCols(), vals, rowptrs, colinds);
#else
      AddKern::addSorted(Avals, Arowptrs, Acolinds, alpha, Bvals, Browptrs, Bcolinds, beta, vals, rowptrs, colinds);
#endif
    }
    else
    {
      //use unsorted kernel
#ifdef HAVE_TPETRA_MMM_TIMINGS
        MM = Teuchos::null;
        MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("mm add unsorted entries full kernel"))));
#endif
      if (debug) {
        std::ostringstream os;
        os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
           << "Call AddKern::addUnsorted(...)" << std::endl;
        std::cerr << os.str ();
      }
      AddKern::addUnsorted(Avals, Arowptrs, Acolinds, alpha, Bvals, Browptrs, Bcolinds, beta, Aprime->getGlobalNumCols(), vals, rowptrs, colinds);
    }
    //Bprime col map works as C's row map, since Aprime and Bprime have the same colmaps.
    RCP<const map_type> Ccolmap = Bcolmap;
    C.replaceColMap(Ccolmap);
    C.setAllValues(rowptrs, colinds, vals);
    if(doFillComplete)
    {
    #ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = Teuchos::null;
      MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Tpetra::Crs expertStaticFillComplete"))));
    #endif
      if(!CDomainMap->isSameAs(*Ccolmap))
      {
        if (debug) {
          std::ostringstream os;
          os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
             << "Create Cimport" << std::endl;
          std::cerr << os.str ();
        }
        Cimport = rcp(new import_type(CDomainMap, Ccolmap));
      }
      if(!C.getRowMap()->isSameAs(*CRangeMap))
      {
        if (debug) {
          std::ostringstream os;
          os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
             << "Create Cexport" << std::endl;
          std::cerr << os.str ();
        }
        Cexport = rcp(new export_type(C.getRowMap(), CRangeMap));
      }

      if (debug) {
        std::ostringstream os;
        os << "Proc " << A.getMap ()->getComm ()->getRank () << ": "
           << "Call C->expertStaticFillComplete(...)" << std::endl;
        std::cerr << os.str ();
      }
      C.expertStaticFillComplete(CDomainMap, CRangeMap, Cimport, Cexport, params);
    }
  }
}

// This version of Add takes C as RCP&, so C may be null on input (in this case,
// it is allocated and constructed in this function).
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  Scalar scalarA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  bool transposeB,
  Scalar scalarB,
  Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& C)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;
  using std::endl;
  //  typedef typename ArrayView<const Scalar>::size_type size_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef Map<LocalOrdinal, GlobalOrdinal, Node>                            map_type;
  //  typedef Import<LocalOrdinal, GlobalOrdinal, Node>                         import_type;
  //  typedef RowGraph<LocalOrdinal, GlobalOrdinal, Node>                       row_graph_type;
  //  typedef CrsGraph<LocalOrdinal, GlobalOrdinal, Node>                       crs_graph_type;
  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>              crs_matrix_type;
  typedef RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>    transposer_type;

  std::string prefix = "TpetraExt::MatrixMatrix::Add(): ";

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.isFillComplete () || ! B.isFillComplete (), std::invalid_argument,
    prefix << "A and B must both be fill complete before calling this function.");

  if(C.is_null()) {
    TEUCHOS_TEST_FOR_EXCEPTION(!A.haveGlobalConstants(), std::logic_error,
        prefix << "C is null (must be allocated), but A.haveGlobalConstants() is false. "
        "Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(!B.haveGlobalConstants(), std::logic_error,
        prefix << "C is null (must be allocated), but B.haveGlobalConstants() is false. "
        "Please report this bug to the Tpetra developers.");
  }

#ifdef HAVE_TPETRA_DEBUG
  {
    const bool domainMapsSame =
      (! transposeA && ! transposeB && ! A.getDomainMap ()->isSameAs (* (B.getDomainMap ()))) ||
      (! transposeA && transposeB && ! A.getDomainMap ()->isSameAs (* (B.getRangeMap ()))) ||
      (transposeA && ! transposeB && ! A.getRangeMap ()->isSameAs (* (B.getDomainMap ())));
    TEUCHOS_TEST_FOR_EXCEPTION(domainMapsSame, std::invalid_argument,
      prefix << "The domain Maps of Op(A) and Op(B) are not the same.");

    const bool rangeMapsSame =
      (! transposeA && ! transposeB && ! A.getRangeMap ()->isSameAs (* (B.getRangeMap ()))) ||
      (! transposeA && transposeB && ! A.getRangeMap ()->isSameAs (* (B.getDomainMap ()))) ||
      (transposeA && ! transposeB && ! A.getDomainMap ()->isSameAs (* (B.getRangeMap ())));
    TEUCHOS_TEST_FOR_EXCEPTION(rangeMapsSame, std::invalid_argument,
      prefix << "The range Maps of Op(A) and Op(B) are not the same.");
  }
#endif // HAVE_TPETRA_DEBUG

  using Teuchos::ParameterList;
  RCP<ParameterList> transposeParams (new ParameterList);
  transposeParams->set ("sort", false);

  // Form the explicit transpose of A if necessary.
  RCP<const crs_matrix_type> Aprime;
  if (transposeA) {
    transposer_type theTransposer (rcpFromRef (A));
    Aprime = theTransposer.createTranspose (transposeParams);
  }
  else {
    Aprime = rcpFromRef (A);
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Aprime.is_null (), std::logic_error,
    prefix << "Failed to compute Op(A). Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  // Form the explicit transpose of B if necessary.
  RCP<const crs_matrix_type> Bprime;
  if (transposeB) {
    transposer_type theTransposer (rcpFromRef (B));
    Bprime = theTransposer.createTranspose (transposeParams);
  }
  else {
    Bprime = rcpFromRef (B);
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Bprime.is_null (), std::logic_error,
    prefix << "Failed to compute Op(B). Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  bool CwasFillComplete = false;

  // Allocate or zero the entries of the result matrix.
  if (! C.is_null ()) {
    CwasFillComplete = C->isFillComplete();
    if(CwasFillComplete)
      C->resumeFill();
    C->setAllToScalar (STS::zero ());
  } else {
    // FIXME (mfh 08 May 2013) When I first looked at this method, I
    // noticed that C was being given the row Map of Aprime (the
    // possibly transposed version of A).  Is this what we want?

    // It is a precondition that Aprime and Bprime have the same domain and range maps.
    // However, they may have different row maps. In this case, it's difficult to
    // get a precise upper bound on the number of entries in each local row of C, so
    // just use the looser upper bound based on the max number of entries in any row of Aprime and Bprime.
    if(Aprime->getRowMap()->isSameAs(*Bprime->getRowMap())) {
      LocalOrdinal numLocalRows = Aprime->getLocalNumRows();
      Array<size_t> CmaxEntriesPerRow(numLocalRows);
      for(LocalOrdinal i = 0; i < numLocalRows; i++) {
        CmaxEntriesPerRow[i] = Aprime->getNumEntriesInLocalRow(i) + Bprime->getNumEntriesInLocalRow(i);
      }
      C = rcp (new crs_matrix_type (Aprime->getRowMap (), CmaxEntriesPerRow()));
    }
    else {
      // Note: above we checked that Aprime and Bprime have global constants, so it's safe to ask for max entries per row.
      C = rcp (new crs_matrix_type (Aprime->getRowMap (), Aprime->getGlobalMaxNumRowEntries() + Bprime->getGlobalMaxNumRowEntries()));
    }
  }

#ifdef HAVE_TPETRA_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(Aprime.is_null (), std::logic_error,
    prefix << "At this point, Aprime is null. Please report this bug to the Tpetra developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(Bprime.is_null (), std::logic_error,
    prefix << "At this point, Bprime is null. Please report this bug to the Tpetra developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
    prefix << "At this point, C is null. Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

  Array<RCP<const crs_matrix_type> > Mat =
    tuple<RCP<const crs_matrix_type> > (Aprime, Bprime);
  Array<Scalar> scalar = tuple<Scalar> (scalarA, scalarB);

  // do a loop over each matrix to add: A reordering might be more efficient
  for (int k = 0; k < 2; ++k) {
    typename crs_matrix_type::nonconst_global_inds_host_view_type Indices;
    typename crs_matrix_type::nonconst_values_host_view_type Values;

    // Loop over each locally owned row of the current matrix (either
    // Aprime or Bprime), and sum its entries into the corresponding
    // row of C.  This works regardless of whether Aprime or Bprime
    // has the same row Map as C, because both sumIntoGlobalValues and
    // insertGlobalValues allow summing resp. inserting into nonowned
    // rows of C.
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(Mat[k].is_null (), std::logic_error,
      prefix << "At this point, curRowMap is null. Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    RCP<const map_type> curRowMap = Mat[k]->getRowMap ();
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(curRowMap.is_null (), std::logic_error,
      prefix << "At this point, curRowMap is null. Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG

    const size_t localNumRows = Mat[k]->getLocalNumRows ();
    for (size_t i = 0; i < localNumRows; ++i) {
      const GlobalOrdinal globalRow = curRowMap->getGlobalElement (i);
      size_t numEntries = Mat[k]->getNumEntriesInGlobalRow (globalRow);
      if (numEntries > 0) {
        if(numEntries > Indices.extent(0)) {
          Kokkos::resize(Indices, numEntries);
          Kokkos::resize(Values, numEntries);
        }
        Mat[k]->getGlobalRowCopy (globalRow, Indices, Values, numEntries);

        if (scalar[k] != STS::one ()) {
          for (size_t j = 0; j < numEntries; ++j) {
            Values[j] *= scalar[k];
          }
        }

        if (CwasFillComplete) {
          size_t result = C->sumIntoGlobalValues (globalRow, numEntries, 
                                  reinterpret_cast<Scalar *>(Values.data()), Indices.data());
          TEUCHOS_TEST_FOR_EXCEPTION(result != numEntries, std::logic_error,
              prefix << "sumIntoGlobalValues failed to add entries from A or B into C.");
        } else {
          C->insertGlobalValues (globalRow,  numEntries, 
                                 reinterpret_cast<Scalar *>(Values.data()), Indices.data());
        }
      }
    }
  }
  if(CwasFillComplete) {
    C->fillComplete(C->getDomainMap (),
                    C->getRangeMap ());
  }
}

// This version of Add takes C as const RCP&, so C must not be null on input. Otherwise, its behavior is identical
// to the above version where C is RCP&.
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void Add(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  bool transposeA,
  Scalar scalarA,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  bool transposeB,
  Scalar scalarB,
  const Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& C)
{
  std::string prefix = "TpetraExt::MatrixMatrix::Add(): ";

  TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::invalid_argument,
    prefix << "C must not be null");

  Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > C_ = C;
  Add(A, transposeA, scalarA, B, transposeB, scalarB, C_);
}

} //End namespace MatrixMatrix

namespace MMdetails{

/*********************************************************************************************************/
//// Prints MMM-style statistics on communication done with an Import or Export object
//template <class TransferType>
//void printMultiplicationStatistics(Teuchos::RCP<TransferType > Transfer, const std::string &label) {
//  if (Transfer.is_null())
//    return;
//
//  const Distributor & Distor                   = Transfer->getDistributor();
//  Teuchos::RCP<const Teuchos::Comm<int> > Comm = Transfer->getSourceMap()->getComm();
//
//  size_t rows_send   = Transfer->getNumExportIDs();
//  size_t rows_recv   = Transfer->getNumRemoteIDs();
//
//  size_t round1_send = Transfer->getNumExportIDs() * sizeof(size_t);
//  size_t round1_recv = Transfer->getNumRemoteIDs() * sizeof(size_t);
//  size_t num_send_neighbors = Distor.getNumSends();
//  size_t num_recv_neighbors = Distor.getNumReceives();
//  size_t round2_send, round2_recv;
//  Distor.getLastDoStatistics(round2_send,round2_recv);
//
//  int myPID    = Comm->getRank();
//  int NumProcs = Comm->getSize();
//
//  // Processor by processor statistics
//  //    printf("[%d] %s Statistics: neigh[s/r]=%d/%d rows[s/r]=%d/%d r1bytes[s/r]=%d/%d r2bytes[s/r]=%d/%d\n",
//  //    myPID, label.c_str(),num_send_neighbors,num_recv_neighbors,rows_send,rows_recv,round1_send,round1_recv,round2_send,round2_recv);
//
//  // Global statistics
//  size_t lstats[8] = {num_send_neighbors,num_recv_neighbors,rows_send,rows_recv,round1_send,round1_recv,round2_send,round2_recv};
//  size_t gstats_min[8], gstats_max[8];
//
//  double lstats_avg[8], gstats_avg[8];
//  for(int i=0; i<8; i++)
//    lstats_avg[i] = ((double)lstats[i])/NumProcs;
//
//  Teuchos::reduceAll(*Comm(),Teuchos::REDUCE_MIN,8,lstats,gstats_min);
//  Teuchos::reduceAll(*Comm(),Teuchos::REDUCE_MAX,8,lstats,gstats_max);
//  Teuchos::reduceAll(*Comm(),Teuchos::REDUCE_SUM,8,lstats_avg,gstats_avg);
//
//  if(!myPID) {
//    printf("%s Send Statistics[min/avg/max]: neigh=%d/%4.1f/%d rows=%d/%4.1f/%d round1=%d/%4.1f/%d round2=%d/%4.1f/%d\n", label.c_str(),
//           (int)gstats_min[0],gstats_avg[0],(int)gstats_max[0], (int)gstats_min[2],gstats_avg[2],(int)gstats_max[2],
//           (int)gstats_min[4],gstats_avg[4],(int)gstats_max[4], (int)gstats_min[6],gstats_avg[6],(int)gstats_max[6]);
//    printf("%s Recv Statistics[min/avg/max]: neigh=%d/%4.1f/%d rows=%d/%4.1f/%d round1=%d/%4.1f/%d round2=%d/%4.1f/%d\n", label.c_str(),
//           (int)gstats_min[1],gstats_avg[1],(int)gstats_max[1], (int)gstats_min[3],gstats_avg[3],(int)gstats_max[3],
//           (int)gstats_min[5],gstats_avg[5],(int)gstats_max[5], (int)gstats_min[7],gstats_avg[7],(int)gstats_max[7]);
//  }
//}

// Kernel method for computing the local portion of C = A*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_AT_B_newmatrix(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string & label,
  const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Scalar                            SC;
  typedef LocalOrdinal                      LO;
  typedef GlobalOrdinal                     GO;
  typedef Node                              NO;
  typedef CrsMatrixStruct<SC,LO,GO,NO>      crs_matrix_struct_type;
  typedef RowMatrixTransposer<SC,LO,GO,NO>  transposer_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> MM = rcp (new TimeMonitor
    (*TimeMonitor::getNewTimer (prefix_mmm + "MMM-T Transpose")));
#endif

  /*************************************************************/
  /* 1) Local Transpose of A                                   */
  /*************************************************************/
  transposer_type transposer (rcpFromRef (A), label + std::string("XP: "));

  using Teuchos::ParameterList;
  RCP<ParameterList> transposeParams (new ParameterList);
  transposeParams->set ("sort", true); // Kokkos Kernels spgemm requires inputs to be sorted
  if(! params.is_null ()) {
    transposeParams->set ("compute global constants",
                          params->get ("compute global constants: temporaries",
                                       false));
  }
  RCP<Tpetra::CrsMatrix<SC, LO, GO, NO>> Atrans =
    transposer.createTransposeLocal (transposeParams);

  /*************************************************************/
  /* 2/3) Call mult_A_B_newmatrix w/ fillComplete              */
  /*************************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp (new TimeMonitor
    (*TimeMonitor::getNewTimer (prefix_mmm + std::string ("MMM-T I&X"))));
#endif

  // Get views, asserting that no import is required to speed up computation
  crs_matrix_struct_type Aview;
  crs_matrix_struct_type Bview;
  RCP<const Import<LO, GO, NO> > dummyImporter;

  // NOTE: the I&X routine sticks an importer on the paramlist as output, so we have to use a unique guy here
  RCP<Teuchos::ParameterList> importParams1 (new ParameterList);
  if (! params.is_null ()) {
    importParams1->set ("compute global constants",
                        params->get ("compute global constants: temporaries",
                                     false));
    auto slist = params->sublist ("matrixmatrix: kernel params", false);
    bool isMM = slist.get ("isMatrixMatrix_TransferAndFillComplete", false);
    bool overrideAllreduce = slist.get("MM_TAFC_OverrideAllreduceCheck", false);
    int mm_optimization_core_count =
      ::Tpetra::Details::Behavior::TAFC_OptimizationCoreCount();
    mm_optimization_core_count =
      slist.get ("MM_TAFC_OptimizationCoreCount", mm_optimization_core_count);
    int mm_optimization_core_count2 =
      params->get ("MM_TAFC_OptimizationCoreCount", mm_optimization_core_count);
    if (mm_optimization_core_count2 < mm_optimization_core_count) {
      mm_optimization_core_count = mm_optimization_core_count2;
    }
    auto & sip1 = importParams1->sublist ("matrixmatrix: kernel params", false);
    sip1.set ("MM_TAFC_OptimizationCoreCount", mm_optimization_core_count);
    sip1.set ("isMatrixMatrix_TransferAndFillComplete", isMM);
    sip1.set ("MM_TAFC_OverrideAllreduceCheck", overrideAllreduce);
  }

  MMdetails::import_and_extract_views (*Atrans, Atrans->getRowMap (),
                                       Aview, dummyImporter, true,
                                       label, importParams1);

  RCP<ParameterList> importParams2 (new ParameterList);
  if (! params.is_null ()) {
    importParams2->set ("compute global constants",
                        params->get ("compute global constants: temporaries",
                                     false));
    auto slist = params->sublist ("matrixmatrix: kernel params", false);
    bool isMM = slist.get ("isMatrixMatrix_TransferAndFillComplete", false);
    bool overrideAllreduce = slist.get ("MM_TAFC_OverrideAllreduceCheck", false);
    int mm_optimization_core_count =
      ::Tpetra::Details::Behavior::TAFC_OptimizationCoreCount();
    mm_optimization_core_count =
      slist.get ("MM_TAFC_OptimizationCoreCount",
                 mm_optimization_core_count);
    int mm_optimization_core_count2 =
      params->get ("MM_TAFC_OptimizationCoreCount",
                   mm_optimization_core_count);
    if (mm_optimization_core_count2 < mm_optimization_core_count) {
      mm_optimization_core_count = mm_optimization_core_count2;
    }
    auto & sip2 = importParams2->sublist ("matrixmatrix: kernel params", false);
    sip2.set ("MM_TAFC_OptimizationCoreCount", mm_optimization_core_count);
    sip2.set ("isMatrixMatrix_TransferAndFillComplete", isMM);
    sip2.set ("MM_TAFC_OverrideAllreduceCheck", overrideAllreduce);
  }

  if(B.getRowMap()->isSameAs(*Atrans->getColMap())){
    MMdetails::import_and_extract_views(B, B.getRowMap(), Bview, dummyImporter,true, label,importParams2);
  }
  else {
    MMdetails::import_and_extract_views(B, Atrans->getColMap(), Bview, dummyImporter,false, label,importParams2);
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM-T AB-core"))));
#endif

  RCP<Tpetra::CrsMatrix<SC, LO, GO, NO>> Ctemp;

  // If Atrans has no Exporter, we can use C instead of having to create a temp matrix
  bool needs_final_export = ! Atrans->getGraph ()->getExporter ().is_null();
  if (needs_final_export) {
    Ctemp = rcp (new Tpetra::CrsMatrix<SC, LO, GO, NO> (Atrans->getRowMap (), 0));
  }
  else {
    Ctemp = rcp (&C, false);
  }

  mult_A_B_newmatrix(Aview, Bview, *Ctemp, label,params);

  /*************************************************************/
  /* 4) exportAndFillComplete matrix                           */
  /*************************************************************/
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM-T exportAndFillComplete"))));
#endif

  RCP<Tpetra::CrsMatrix<SC, LO, GO, NO>> Crcp (&C, false);

  if (needs_final_export) {
    ParameterList labelList;
    labelList.set("Timer Label", label);
    if(!params.is_null()) {
      ParameterList& params_sublist = params->sublist("matrixmatrix: kernel params",false);
      ParameterList& labelList_subList = labelList.sublist("matrixmatrix: kernel params",false);
      int mm_optimization_core_count = ::Tpetra::Details::Behavior::TAFC_OptimizationCoreCount();
      mm_optimization_core_count = params_sublist.get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
      int mm_optimization_core_count2 = params->get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
      if(mm_optimization_core_count2<mm_optimization_core_count) mm_optimization_core_count=mm_optimization_core_count2;
      labelList_subList.set("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count,"Core Count above which the optimized neighbor discovery is used");
      bool isMM = params_sublist.get("isMatrixMatrix_TransferAndFillComplete",false);
      bool overrideAllreduce = params_sublist.get("MM_TAFC_OverrideAllreduceCheck",false);

      labelList_subList.set ("isMatrixMatrix_TransferAndFillComplete", isMM,
                             "This parameter should be set to true only for MatrixMatrix operations: the optimization in Epetra that was ported to Tpetra does _not_ take into account the possibility that for any given source PID, a particular GID may not exist on the target PID: i.e. a transfer operation. A fix for this general case is in development.");
      labelList.set("compute global constants",params->get("compute global constants",true));
      labelList.set("MM_TAFC_OverrideAllreduceCheck",overrideAllreduce);
    }
    Ctemp->exportAndFillComplete (Crcp,
                                  *Ctemp->getGraph ()->getExporter (),
                                  B.getDomainMap (),
                                  A.getDomainMap (),
                                  rcp (&labelList, false));
  }
#ifdef HAVE_TPETRA_MMM_STATISTICS
  printMultiplicationStatistics(Ctemp->getGraph()->getExporter(), label+std::string(" AT_B MMM"));
#endif
}

/*********************************************************************************************************/
// Kernel method for computing the local portion of C = A*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsWrapper<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& /* label */,
  const Teuchos::RCP<Teuchos::ParameterList>& /* params */)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::OrdinalTraits;
  using Teuchos::null;

  typedef Teuchos::ScalarTraits<Scalar> STS;
  // TEUCHOS_FUNC_TIME_MONITOR_DIFF("mult_A_B", mult_A_B);
  LocalOrdinal C_firstCol = Bview.colMap->getMinLocalIndex();
  LocalOrdinal C_lastCol  = Bview.colMap->getMaxLocalIndex();

  LocalOrdinal C_firstCol_import = OrdinalTraits<LocalOrdinal>::zero();
  LocalOrdinal C_lastCol_import  = OrdinalTraits<LocalOrdinal>::invalid();

  ArrayView<const GlobalOrdinal> bcols = Bview.colMap->getLocalElementList();
  ArrayView<const GlobalOrdinal> bcols_import = null;
  if (Bview.importColMap != null) {
    C_firstCol_import = Bview.importColMap->getMinLocalIndex();
    C_lastCol_import  = Bview.importColMap->getMaxLocalIndex();

    bcols_import = Bview.importColMap->getLocalElementList();
  }

  size_t C_numCols = C_lastCol - C_firstCol +
                        OrdinalTraits<LocalOrdinal>::one();
  size_t C_numCols_import = C_lastCol_import - C_firstCol_import +
                                OrdinalTraits<LocalOrdinal>::one();

  if (C_numCols_import > C_numCols)
    C_numCols = C_numCols_import;

  Array<Scalar> dwork = Array<Scalar>(C_numCols);
  Array<GlobalOrdinal> iwork = Array<GlobalOrdinal>(C_numCols);
  Array<size_t> iwork2 = Array<size_t>(C_numCols);

  Array<Scalar> C_row_i = dwork;
  Array<GlobalOrdinal> C_cols = iwork;
  Array<size_t> c_index = iwork2;
  Array<GlobalOrdinal> combined_index = Array<GlobalOrdinal>(2*C_numCols);
  Array<Scalar> combined_values = Array<Scalar>(2*C_numCols);

  size_t C_row_i_length, j, k, last_index;

  // Run through all the hash table lookups once and for all
  LocalOrdinal LO_INVALID = OrdinalTraits<LocalOrdinal>::invalid();
  Array<LocalOrdinal> Acol2Brow(Aview.colMap->getLocalNumElements(),LO_INVALID);
  Array<LocalOrdinal> Acol2Irow(Aview.colMap->getLocalNumElements(),LO_INVALID);
  if(Aview.colMap->isSameAs(*Bview.origMatrix->getRowMap())){
    // Maps are the same: Use local IDs as the hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <=
            Aview.colMap->getMaxLocalIndex(); i++)
      Acol2Brow[i]=i;
  }
  else {
    // Maps are not the same:  Use the map's hash
    for(LocalOrdinal i=Aview.colMap->getMinLocalIndex(); i <=
          Aview.colMap->getMaxLocalIndex(); i++) {
      GlobalOrdinal GID = Aview.colMap->getGlobalElement(i);
      LocalOrdinal BLID = Bview.origMatrix->getRowMap()->getLocalElement(GID);
      if(BLID != LO_INVALID) Acol2Brow[i] = BLID;
      else Acol2Irow[i] = Bview.importMatrix->getRowMap()->getLocalElement(GID);
    }
  }

  // To form C = A*B we're going to execute this expression:
  //
  //  C(i,j) = sum_k( A(i,k)*B(k,j) )
  //
  // Our goal, of course, is to navigate the data in A and B once, without
  // performing searches for column-indices, etc.
  auto Arowptr = Aview.origMatrix->getLocalRowPtrsHost();
  auto Acolind = Aview.origMatrix->getLocalIndicesHost();
  auto Avals   = Aview.origMatrix->getLocalValuesHost(Tpetra::Access::ReadOnly);
  auto Browptr = Bview.origMatrix->getLocalRowPtrsHost();
  auto Bcolind = Bview.origMatrix->getLocalIndicesHost();
  auto Bvals   = Bview.origMatrix->getLocalValuesHost(Tpetra::Access::ReadOnly);
  decltype(Browptr) Irowptr;
  decltype(Bcolind) Icolind;
  decltype(Bvals) Ivals;
  if(!Bview.importMatrix.is_null()) {
    Irowptr = Bview.importMatrix->getLocalRowPtrsHost();
    Icolind = Bview.importMatrix->getLocalIndicesHost();
    Ivals   = Bview.importMatrix->getLocalValuesHost(Tpetra::Access::ReadOnly);
  }

  bool C_filled = C.isFillComplete();

  for (size_t i = 0; i < C_numCols; i++)
      c_index[i] = OrdinalTraits<size_t>::invalid();

  // Loop over the rows of A.
  size_t Arows = Aview.rowMap->getLocalNumElements();
  for(size_t i=0; i<Arows; ++i) {

    // Only navigate the local portion of Aview... which is, thankfully, all of
    // A since this routine doesn't do transpose modes
    GlobalOrdinal global_row = Aview.rowMap->getGlobalElement(i);

    // Loop across the i-th row of A and for each corresponding row in B, loop
    // across columns and accumulate product A(i,k)*B(k,j) into our partial sum
    // quantities C_row_i. In other words, as we stride across B(k,:) we're
    // calculating updates for row i of the result matrix C.
    C_row_i_length = OrdinalTraits<size_t>::zero();

    for (k = Arowptr[i]; k < Arowptr[i+1]; ++k) {
      LocalOrdinal Ak = Acol2Brow[Acolind[k]];
      const Scalar Aval = Avals[k];
      if (Aval == STS::zero())
        continue;

      if (Ak == LO_INVALID)
        continue;

      for (j = Browptr[Ak]; j < Browptr[Ak+1]; ++j) {
          LocalOrdinal col = Bcolind[j];
          //assert(col >= 0 && col < C_numCols);

          if (c_index[col] == OrdinalTraits<size_t>::invalid()){
          //assert(C_row_i_length >= 0 && C_row_i_length < C_numCols);
            // This has to be a +=  so insertGlobalValue goes out
            C_row_i[C_row_i_length] = Aval*Bvals[j];
            C_cols[C_row_i_length] = col;
            c_index[col] = C_row_i_length;
            C_row_i_length++;

          } else {
            // static cast from impl_scalar_type to Scalar needed for complex
            C_row_i[c_index[col]] += Aval * static_cast<Scalar>(Bvals[j]);
          }
        }
    }

    for (size_t ii = 0; ii < C_row_i_length; ii++) {
      c_index[C_cols[ii]] = OrdinalTraits<size_t>::invalid();
      C_cols[ii] = bcols[C_cols[ii]];
      combined_index[ii] = C_cols[ii];
      combined_values[ii] = C_row_i[ii];
    }
    last_index = C_row_i_length;

    //
    //Now put the C_row_i values into C.
    //
    // We might have to revamp this later.
    C_row_i_length = OrdinalTraits<size_t>::zero();

    for (k = Arowptr[i]; k < Arowptr[i+1]; ++k) {
      LocalOrdinal Ak = Acol2Brow[Acolind[k]];
      const Scalar Aval = Avals[k];
      if (Aval == STS::zero())
        continue;

      if (Ak!=LO_INVALID) continue;

      Ak = Acol2Irow[Acolind[k]];
      for (j = Irowptr[Ak]; j < Irowptr[Ak+1]; ++j) {
          LocalOrdinal col = Icolind[j];
          //assert(col >= 0 && col < C_numCols);

          if (c_index[col] == OrdinalTraits<size_t>::invalid()) {
            //assert(C_row_i_length >= 0 && C_row_i_length < C_numCols);
            // This has to be a +=  so insertGlobalValue goes out
            C_row_i[C_row_i_length] = Aval*Ivals[j];
            C_cols[C_row_i_length] = col;
            c_index[col] = C_row_i_length;
            C_row_i_length++;

            } else {
              // This has to be a +=  so insertGlobalValue goes out
              // static cast from impl_scalar_type to Scalar needed for complex
              C_row_i[c_index[col]] += Aval * static_cast<Scalar>(Ivals[j]);
            }
        }
    }

    for (size_t ii = 0; ii < C_row_i_length; ii++) {
      c_index[C_cols[ii]] = OrdinalTraits<size_t>::invalid();
      C_cols[ii] = bcols_import[C_cols[ii]];
      combined_index[last_index] = C_cols[ii];
      combined_values[last_index] = C_row_i[ii];
      last_index++;
    }

    // Now put the C_row_i values into C.
    // We might have to revamp this later.
    C_filled ?
      C.sumIntoGlobalValues(
          global_row,
          combined_index.view(OrdinalTraits<size_t>::zero(), last_index),
          combined_values.view(OrdinalTraits<size_t>::zero(), last_index))
      :
      C.insertGlobalValues(
          global_row,
          combined_index.view(OrdinalTraits<size_t>::zero(), last_index),
          combined_values.view(OrdinalTraits<size_t>::zero(), last_index));

  }
}

/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void setMaxNumEntriesPerRow(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Mview) {
  typedef typename Teuchos::Array<Teuchos::ArrayView<const LocalOrdinal> >::size_type local_length_size;
  Mview.maxNumRowEntries = Teuchos::OrdinalTraits<local_length_size>::zero();

  if (Mview.indices.size() > Teuchos::OrdinalTraits<local_length_size>::zero()) {
    Mview.maxNumRowEntries = Mview.indices[0].size();

    for (local_length_size i = 1; i < Mview.indices.size(); ++i)
      if (Mview.indices[i].size() > Mview.maxNumRowEntries)
        Mview.maxNumRowEntries = Mview.indices[i].size();
  }
}

/*********************************************************************************************************/
template<class CrsMatrixType>
size_t C_estimate_nnz(CrsMatrixType & A, CrsMatrixType &B){
  // Follows the NZ estimate in ML's ml_matmatmult.c
  size_t Aest = 100, Best=100;
  if (A.getLocalNumEntries() >= A.getLocalNumRows())
    Aest = (A.getLocalNumRows() > 0) ? A.getLocalNumEntries()/A.getLocalNumRows() : 100;
  if (B.getLocalNumEntries() >= B.getLocalNumRows())
    Best = (B.getLocalNumRows() > 0) ? B.getLocalNumEntries()/B.getLocalNumRows() : 100;

  size_t nnzperrow = (size_t)(sqrt((double)Aest) + sqrt((double)Best) - 1);
  nnzperrow *= nnzperrow;

  return (size_t)(A.getLocalNumRows()*nnzperrow*0.75 + 100);
}


/*********************************************************************************************************/
// Kernel method for computing the local portion of C = A*B for CrsMatrix
//
// mfh 27 Sep 2016: Currently, mult_AT_B_newmatrix() also calls this
// function, so this is probably the function we want to
// thread-parallelize.
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B_newmatrix(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& label,
  const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Tpetra typedefs
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;
  typedef Import<LO,GO,NO>  import_type;
  typedef Map<LO,GO,NO>     map_type;

  // Kokkos typedefs
  typedef typename map_type::local_map_type local_map_type;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
  typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> MM = rcp(new TimeMonitor(*(TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM M5 Cmap")))));
#endif
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // Build the final importer / column map, hash table lookups for C
  RCP<const import_type> Cimport;
  RCP<const map_type>    Ccolmap;
  RCP<const import_type> Bimport = Bview.origMatrix->getGraph()->getImporter();
  RCP<const import_type> Iimport = Bview.importMatrix.is_null() ?
      Teuchos::null : Bview.importMatrix->getGraph()->getImporter();
  local_map_type Acolmap_local = Aview.colMap->getLocalMap();
  local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
  local_map_type Irowmap_local;  if(!Bview.importMatrix.is_null()) Irowmap_local = Bview.importMatrix->getRowMap()->getLocalMap();
  local_map_type Bcolmap_local = Bview.origMatrix->getColMap()->getLocalMap();
  local_map_type Icolmap_local;  if(!Bview.importMatrix.is_null()) Icolmap_local = Bview.importMatrix->getColMap()->getLocalMap();

  // mfh 27 Sep 2016: Bcol2Ccol is a table that maps from local column
  // indices of B, to local column indices of C.  (B and C have the
  // same number of columns.)  The kernel uses this, instead of
  // copying the entire input matrix B and converting its column
  // indices to those of C.
  lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"),Bview.colMap->getLocalNumElements()), Icol2Ccol;

  if (Bview.importMatrix.is_null()) {
    // mfh 27 Sep 2016: B has no "remotes," so B and C have the same column Map.
    Cimport = Bimport;
    Ccolmap = Bview.colMap;
    const LO colMapSize = static_cast<LO>(Bview.colMap->getLocalNumElements());
    // Bcol2Ccol is trivial
    Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::Bcol2Ccol_fill",
      Kokkos::RangePolicy<execution_space, LO>(0, colMapSize),
      KOKKOS_LAMBDA(const LO i) {
        Bcol2Ccol(i) = i;
      });
  }
  else {
    // mfh 27 Sep 2016: B has "remotes," so we need to build the
    // column Map of C, as well as C's Import object (from its domain
    // Map to its column Map).  C's column Map is the union of the
    // column Maps of (the local part of) B, and the "remote" part of
    // B.  Ditto for the Import.  We have optimized this "setUnion"
    // operation on Import objects and Maps.

    // Choose the right variant of setUnion
    if (!Bimport.is_null() && !Iimport.is_null()) {
      Cimport = Bimport->setUnion(*Iimport);
    }
    else if (!Bimport.is_null() && Iimport.is_null()) {
      Cimport = Bimport->setUnion();
    }
    else if (Bimport.is_null() && !Iimport.is_null()) {
      Cimport = Iimport->setUnion();
    }
    else {
      throw std::runtime_error("TpetraExt::MMM status of matrix importers is nonsensical");
    }
    Ccolmap = Cimport->getTargetMap();

    // FIXME (mfh 27 Sep 2016) This error check requires an all-reduce
    // in general.  We should get rid of it in order to reduce
    // communication costs of sparse matrix-matrix multiply.
    TEUCHOS_TEST_FOR_EXCEPTION(!Cimport->getSourceMap()->isSameAs(*Bview.origMatrix->getDomainMap()),
      std::runtime_error, "Tpetra::MMM: Import setUnion messed with the DomainMap in an unfortunate way");

    // NOTE: This is not efficient and should be folded into setUnion
    //
    // mfh 27 Sep 2016: What the above comment means, is that the
    // setUnion operation on Import objects could also compute these
    // local index - to - local index look-up tables.
    Kokkos::resize(Icol2Ccol,Bview.importMatrix->getColMap()->getLocalNumElements());
    local_map_type Ccolmap_local = Ccolmap->getLocalMap();
    Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::Bcol2Ccol_getGlobalElement",range_type(0,Bview.origMatrix->getColMap()->getLocalNumElements()),KOKKOS_LAMBDA(const LO i) {
        Bcol2Ccol(i) = Ccolmap_local.getLocalElement(Bcolmap_local.getGlobalElement(i));
      });
    Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::Icol2Ccol_getGlobalElement",range_type(0,Bview.importMatrix->getColMap()->getLocalNumElements()),KOKKOS_LAMBDA(const LO i) {
        Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
      });

  }

  // Replace the column map
  //
  // mfh 27 Sep 2016: We do this because C was originally created
  // without a column Map.  Now we have its column Map.
  C.replaceColMap(Ccolmap);

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
  lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getLocalNumElements());
  lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),Aview.colMap->getLocalNumElements());

  Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::construct_tables",range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),KOKKOS_LAMBDA(const LO i) {
      GO aidx = Acolmap_local.getGlobalElement(i);
      LO B_LID = Browmap_local.getLocalElement(aidx);
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
  KernelWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Node,lo_view_t>::mult_A_B_newmatrix_kernel_wrapper(Aview,Bview,targetMapToOrigRow,targetMapToImportRow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);

}

/*********************************************************************************************************/
// Kernel method for computing the local portion of C = A*B for BlockCrsMatrix
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B_newmatrix(BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                        BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                        Teuchos::RCP<BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &C)
{
  using Teuchos::null;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Tpetra typedefs
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;
  typedef Import<LO,GO,NO>  import_type;
  typedef Map<LO,GO,NO>     map_type;
  typedef BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> block_crs_matrix_type;
  typedef typename block_crs_matrix_type::crs_graph_type graph_t;

  // Kokkos typedefs
  typedef typename map_type::local_map_type local_map_type;
  typedef typename block_crs_matrix_type::local_matrix_device_type KBSR;
  typedef typename KBSR::device_type device_t;
  typedef typename KBSR::StaticCrsGraphType static_graph_t;
  typedef typename static_graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename static_graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KBSR::values_type::non_const_type scalar_view_t;
  typedef typename NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
  typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;

  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // Build the final importer / column map, hash table lookups for C
  RCP<const import_type> Cimport;
  RCP<const map_type>    Ccolmap;
  RCP<const import_type> Bimport = Bview.origMatrix->getGraph()->getImporter();
  RCP<const import_type> Iimport = Bview.importMatrix.is_null() ?
                                   Teuchos::null : Bview.importMatrix->getGraph()->getImporter();
  local_map_type Acolmap_local = Aview.colMap->getLocalMap();
  local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
  local_map_type Irowmap_local;  if(!Bview.importMatrix.is_null()) Irowmap_local = Bview.importMatrix->getRowMap()->getLocalMap();
  local_map_type Bcolmap_local = Bview.origMatrix->getColMap()->getLocalMap();
  local_map_type Icolmap_local;  if(!Bview.importMatrix.is_null()) Icolmap_local = Bview.importMatrix->getColMap()->getLocalMap();

  // Bcol2Ccol is a table that maps from local column
  // indices of B, to local column indices of C.  (B and C have the
  // same number of columns.)  The kernel uses this, instead of
  // copying the entire input matrix B and converting its column
  // indices to those of C.
  lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"),Bview.colMap->getLocalNumElements()), Icol2Ccol;

  if (Bview.importMatrix.is_null()) {
    // mfh 27 Sep 2016: B has no "remotes," so B and C have the same column Map.
    Cimport = Bimport;
    Ccolmap = Bview.colMap;
    const LO colMapSize = static_cast<LO>(Bview.colMap->getLocalNumElements());
    // Bcol2Ccol is trivial
    Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::Bcol2Ccol_fill",
      Kokkos::RangePolicy<execution_space, LO>(0, colMapSize),
      KOKKOS_LAMBDA(const LO i) {
        Bcol2Ccol(i) = i;
      });
  }
  else {
    // B has "remotes," so we need to build the
    // column Map of C, as well as C's Import object (from its domain
    // Map to its column Map).  C's column Map is the union of the
    // column Maps of (the local part of) B, and the "remote" part of
    // B.  Ditto for the Import.  We have optimized this "setUnion"
    // operation on Import objects and Maps.

    // Choose the right variant of setUnion
    if (!Bimport.is_null() && !Iimport.is_null()) {
      Cimport = Bimport->setUnion(*Iimport);
    }
    else if (!Bimport.is_null() && Iimport.is_null()) {
      Cimport = Bimport->setUnion();
    }
    else if (Bimport.is_null() && !Iimport.is_null()) {
      Cimport = Iimport->setUnion();
    }
    else {
      throw std::runtime_error("TpetraExt::MMM status of matrix importers is nonsensical");
    }
    Ccolmap = Cimport->getTargetMap();

    // NOTE: This is not efficient and should be folded into setUnion
    //
    // What the above comment means, is that the
    // setUnion operation on Import objects could also compute these
    // local index - to - local index look-up tables.
    Kokkos::resize(Icol2Ccol,Bview.importMatrix->getColMap()->getLocalNumElements());
    local_map_type Ccolmap_local = Ccolmap->getLocalMap();
    Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::Bcol2Ccol_getGlobalElement",
                         range_type(0,Bview.origMatrix->getColMap()->getLocalNumElements()),
                         KOKKOS_LAMBDA(const LO i) {
        Bcol2Ccol(i) = Ccolmap_local.getLocalElement(Bcolmap_local.getGlobalElement(i));
      });
    Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::Icol2Ccol_getGlobalElement",
                         range_type(0,Bview.importMatrix->getColMap()->getLocalNumElements()),
                         KOKKOS_LAMBDA(const LO i) {
        Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
      });
  }

  // Construct tables that map from local column
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
  lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),
                               Aview.colMap->getLocalNumElements());
  lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),
                                 Aview.colMap->getLocalNumElements());

  Kokkos::parallel_for("Tpetra::mult_A_B_newmatrix::construct_tables",
                       range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),
                       KOKKOS_LAMBDA(const LO i) {
      GO aidx = Acolmap_local.getGlobalElement(i);
      LO B_LID = Browmap_local.getLocalElement(aidx);
      if (B_LID != LO_INVALID) {
        targetMapToOrigRow(i)   = B_LID;
        targetMapToImportRow(i) = LO_INVALID;
      } else {
        LO I_LID = Irowmap_local.getLocalElement(aidx);
        targetMapToOrigRow(i)   = LO_INVALID;
        targetMapToImportRow(i) = I_LID;
      }
    });

  // Create the KernelHandle
  using KernelHandle =
    KokkosKernels::Experimental::KokkosKernelsHandle<typename lno_view_t::const_value_type,
                                                     typename lno_nnz_view_t::const_value_type,
                                                     typename scalar_view_t::const_value_type,
                                                     typename device_t::execution_space,
                                                     typename device_t::memory_space,
                                                     typename device_t::memory_space>;
  int team_work_size = 16;  // Defaults to 16 as per Deveci 12/7/16 - csiefer
  std::string myalg("SPGEMM_KK_MEMORY");
  KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(myalg);

  KernelHandle kh;
  kh.create_spgemm_handle(alg_enum);
  kh.set_team_work_size(team_work_size);

  // Get KokkosSparse::BsrMatrix for A and Bmerged (B and BImport)
  const KBSR Amat = Aview.origMatrix->getLocalMatrixDevice();
  const KBSR Bmerged = Tpetra::MMdetails::merge_matrices(Aview,Bview,
                                                         targetMapToOrigRow,targetMapToImportRow,
                                                         Bcol2Ccol,Icol2Ccol,
                                                         Ccolmap.getConst()->getLocalNumElements());

  RCP<graph_t> graphC;
  typename KBSR::values_type values;
  {
    // Call KokkosSparse routines to calculate Amat*Bmerged on device.
    // NOTE: Need to scope guard this since the BlockCrs constructor will need to copy the host graph
    KBSR Cmat;
    KokkosSparse::block_spgemm_symbolic(kh, Amat, false, Bmerged, false, Cmat);
    KokkosSparse::block_spgemm_numeric (kh, Amat, false, Bmerged, false, Cmat);
    kh.destroy_spgemm_handle();
    
    // Build Tpetra::BlockCrsMatrix from KokkosSparse::BsrMatrix
    graphC = rcp(new graph_t(Cmat.graph, Aview.origMatrix->getRowMap(), Ccolmap.getConst()));
    values = Cmat.values;
  }
  C = rcp (new block_crs_matrix_type (*graphC, values, Aview.blocksize));

}

/*********************************************************************************************************/
// AB NewMatrix Kernel wrappers (Default non-threaded version for CrsMatrix)
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class LocalOrdinalViewType>
void KernelWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalOrdinalViewType>::mult_A_B_newmatrix_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                                                               const LocalOrdinalViewType & targetMapToOrigRow,
                                                                                               const LocalOrdinalViewType & targetMapToImportRow,
                                                                                               const LocalOrdinalViewType & Bcol2Ccol,
                                                                                               const LocalOrdinalViewType & Icol2Ccol,
                                                                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                                                               const std::string& label,
                                                                                               const Teuchos::RCP<Teuchos::ParameterList>& params) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix SerialCore"))));
#endif

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_host_type KCRS;
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
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getLocalNumRows();
  size_t n = Ccolmap->getLocalNumElements();
  size_t b_max_nnz_per_row = Bview.origMatrix->getLocalMaxNumRowEntries();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS & Amat = Aview.origMatrix->getLocalMatrixHost();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrixHost();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;

  c_lno_view_t  Irowptr;
  lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixHost();
    Irowptr = lclB.graph.row_map;
    Icolind = lclB.graph.entries;
    Ivals   = lclB.values;
    b_max_nnz_per_row = std::max(b_max_nnz_per_row,Bview.importMatrix->getLocalMaxNumRowEntries());
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  RCP<TimeMonitor> MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix SerialCore - Compare"))));
#endif

  // Classic csr assembly (low memory edition)
  //
  // mfh 27 Sep 2016: C_estimate_nnz does not promise an upper bound.
  // The method loops over rows of A, and may resize after processing
  // each row.  Chris Siefert says that this reflects experience in
  // ML; for the non-threaded case, ML found it faster to spend less
  // effort on estimation and risk an occasional reallocation.
  size_t CSR_alloc = std::max(C_estimate_nnz(*Aview.origMatrix, *Bview.origMatrix), n);
  lno_view_t Crowptr(Kokkos::ViewAllocateWithoutInitializing("Crowptr"),m+1);
  lno_nnz_view_t Ccolind(Kokkos::ViewAllocateWithoutInitializing("Ccolind"),CSR_alloc);
  scalar_view_t Cvals(Kokkos::ViewAllocateWithoutInitializing("Cvals"),CSR_alloc);

  // mfh 27 Sep 2016: The c_status array is an implementation detail
  // of the local sparse matrix-matrix multiply routine.

  // The status array will contain the index into colind where this entry was last deposited.
  //   c_status[i] <  CSR_ip - not in the row yet
  //   c_status[i] >= CSR_ip - this is the entry where you can find the data
  // We start with this filled with INVALID's indicating that there are no entries yet.
  // Sadly, this complicates the code due to the fact that size_t's are unsigned.
  size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
  std::vector<size_t> c_status(n, ST_INVALID);

  // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
  // routine.  The routine computes C := A * (B_local + B_remote).
  //
  // For column index Aik in row i of A, targetMapToOrigRow[Aik] tells
  // you whether the corresponding row of B belongs to B_local
  // ("orig") or B_remote ("Import").

  // For each row of A/C
  size_t CSR_ip = 0, OLD_ip = 0;
  for (size_t i = 0; i < m; i++) {
    // mfh 27 Sep 2016: m is the number of rows in the input matrix A
    // on the calling process.
    Crowptr[i] = CSR_ip;

    // mfh 27 Sep 2016: For each entry of A in the current row of A
    for (size_t k = Arowptr[i]; k < Arowptr[i+1]; k++) {
      LO Aik  = Acolind[k]; // local column index of current entry of A
      const SC Aval = Avals[k];   // value of current entry of A
      if (Aval == SC_ZERO)
        continue; // skip explicitly stored zero values in A

      if (targetMapToOrigRow[Aik] != LO_INVALID) {
        // mfh 27 Sep 2016: If the entry of targetMapToOrigRow
        // corresponding to the current entry of A is populated, then
        // the corresponding row of B is in B_local (i.e., it lives on
        // the calling process).

        // Local matrix
        size_t Bk = static_cast<size_t> (targetMapToOrigRow[Aik]);

        // mfh 27 Sep 2016: Go through all entries in that row of B_local.
        for (size_t j = Browptr[Bk]; j < Browptr[Bk+1]; ++j) {
          LO Bkj = Bcolind[j];
          LO Cij = Bcol2Ccol[Bkj];

          if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip) {
            // New entry
            c_status[Cij]   = CSR_ip;
            Ccolind[CSR_ip] = Cij;
            Cvals[CSR_ip]   = Aval*Bvals[j];
            CSR_ip++;

          } else {
            Cvals[c_status[Cij]] += Aval*Bvals[j];
          }
        }

      } else {
        // mfh 27 Sep 2016: If the entry of targetMapToOrigRow
        // corresponding to the current entry of A NOT populated (has
        // a flag "invalid" value), then the corresponding row of B is
        // in B_local (i.e., it lives on the calling process).

        // Remote matrix
        size_t Ik = static_cast<size_t> (targetMapToImportRow[Aik]);
        for (size_t j = Irowptr[Ik]; j < Irowptr[Ik+1]; ++j) {
          LO Ikj = Icolind[j];
          LO Cij = Icol2Ccol[Ikj];

          if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip){
            // New entry
            c_status[Cij]   = CSR_ip;
            Ccolind[CSR_ip] = Cij;
            Cvals[CSR_ip]   = Aval*Ivals[j];
            CSR_ip++;
          } else {
            Cvals[c_status[Cij]] += Aval*Ivals[j];
          }
        }
      }
    }

    // Resize for next pass if needed
    if (i+1 < m && CSR_ip + std::min(n,(Arowptr[i+2]-Arowptr[i+1])*b_max_nnz_per_row) > CSR_alloc) {
      CSR_alloc *= 2;
      Kokkos::resize(Ccolind,CSR_alloc);
      Kokkos::resize(Cvals,CSR_alloc);
    }
    OLD_ip = CSR_ip;
  }

  Crowptr[m] = CSR_ip;

  // Downward resize
  Kokkos::resize(Ccolind,CSR_ip);
  Kokkos::resize(Cvals,CSR_ip);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  {
  auto MM3(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix Final Sort")));
#endif

  // Final sort & set of CRS arrays
  if (params.is_null() || params->get("sort entries",true))
    Import_Util::sortCrsEntries(Crowptr,Ccolind, Cvals);
  C.setAllValues(Crowptr,Ccolind, Cvals);


#ifdef HAVE_TPETRA_MMM_TIMINGS
  }
  auto MM4(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix ESFC")));
  {
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
  C.expertStaticFillComplete(Bview. origMatrix->getDomainMap(), Aview. origMatrix->getRangeMap(), Cimport,dummyExport,labelList);
#ifdef HAVE_TPETRA_MMM_TIMINGS
  }
  MM2 = Teuchos::null;
#endif

}
/*********************************************************************************************************/
// Kernel method for computing the local portion of C = A*B (reuse)
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void mult_A_B_reuse(
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& label,
  const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Tpetra typedefs
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;
  typedef Import<LO,GO,NO>  import_type;
  typedef Map<LO,GO,NO>     map_type;

  // Kokkos typedefs
  typedef typename map_type::local_map_type local_map_type;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
  typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Reuse Cmap"))));
#endif
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // Grab all the maps
  RCP<const import_type> Cimport = C.getGraph()->getImporter();
  RCP<const map_type>    Ccolmap = C.getColMap();
  local_map_type Acolmap_local = Aview.colMap->getLocalMap();
  local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
  local_map_type Irowmap_local;  if(!Bview.importMatrix.is_null()) Irowmap_local = Bview.importMatrix->getRowMap()->getLocalMap();
  local_map_type Bcolmap_local = Bview.origMatrix->getColMap()->getLocalMap();
  local_map_type Icolmap_local;  if(!Bview.importMatrix.is_null()) Icolmap_local = Bview.importMatrix->getColMap()->getLocalMap();
  local_map_type Ccolmap_local = Ccolmap->getLocalMap();

  // Build the final importer / column map, hash table lookups for C
  lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"),Bview.colMap->getLocalNumElements()), Icol2Ccol;
  {
    // Bcol2Col may not be trivial, as Ccolmap is compressed during fillComplete in newmatrix
    // So, column map of C may be a strict subset of the column map of B
    Kokkos::parallel_for(range_type(0,Bview.origMatrix->getColMap()->getLocalNumElements()),KOKKOS_LAMBDA(const LO i) {
        Bcol2Ccol(i) = Ccolmap_local.getLocalElement(Bcolmap_local.getGlobalElement(i));
      });

    if (!Bview.importMatrix.is_null()) {
      TEUCHOS_TEST_FOR_EXCEPTION(!Cimport->getSourceMap()->isSameAs(*Bview.origMatrix->getDomainMap()),
                                 std::runtime_error, "Tpetra::MMM: Import setUnion messed with the DomainMap in an unfortunate way");

      Kokkos::resize(Icol2Ccol,Bview.importMatrix->getColMap()->getLocalNumElements());
      Kokkos::parallel_for(range_type(0,Bview.importMatrix->getColMap()->getLocalNumElements()),KOKKOS_LAMBDA(const LO i) {
          Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
        });
    }
  }

  // Run through all the hash table lookups once and for all
  lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getLocalNumElements());
  lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),Aview.colMap->getLocalNumElements());
  Kokkos::parallel_for(range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),KOKKOS_LAMBDA(const LO i) {
      GO aidx = Acolmap_local.getGlobalElement(i);
      LO B_LID = Browmap_local.getLocalElement(aidx);
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
  KernelWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Node,lo_view_t>::mult_A_B_reuse_kernel_wrapper(Aview,Bview,targetMapToOrigRow,targetMapToImportRow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);
}

/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class LocalOrdinalViewType>
void KernelWrappers<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalOrdinalViewType>::mult_A_B_reuse_kernel_wrapper(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                                                               const LocalOrdinalViewType & targetMapToOrigRow,
                                                                                               const LocalOrdinalViewType & targetMapToImportRow,
                                                                                               const LocalOrdinalViewType & Bcol2Ccol,
                                                                                               const LocalOrdinalViewType & Icol2Ccol,
                                                                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > /* Cimport */,
                                                                                               const std::string& label,
                                                                                               const Teuchos::RCP<Teuchos::ParameterList>& /* params */) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Reuse SerialCore"))));
  Teuchos::RCP<Teuchos::TimeMonitor> MM2;
#else
  (void)label;
#endif
  using Teuchos::RCP;
  using Teuchos::rcp;


  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_host_type KCRS;
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
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getLocalNumRows();
  size_t n = Ccolmap->getLocalNumElements();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS & Amat = Aview.origMatrix->getLocalMatrixHost();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrixHost();
  const KCRS & Cmat = C.getLocalMatrixHost();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map, Crowptr = Cmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries, Ccolind = Cmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  scalar_view_t Cvals = Cmat.values;

  c_lno_view_t  Irowptr;
  lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixHost();
    Irowptr = lclB.graph.row_map;
    Icolind = lclB.graph.entries;
    Ivals   = lclB.values;
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Newmatrix SerialCore - Compare"))));
#endif

  // Classic csr assembly (low memory edition)
  // mfh 27 Sep 2016: The c_status array is an implementation detail
  // of the local sparse matrix-matrix multiply routine.

  // The status array will contain the index into colind where this entry was last deposited.
  //   c_status[i] <  CSR_ip - not in the row yet
  //   c_status[i] >= CSR_ip - this is the entry where you can find the data
  // We start with this filled with INVALID's indicating that there are no entries yet.
  // Sadly, this complicates the code due to the fact that size_t's are unsigned.
  std::vector<size_t> c_status(n, ST_INVALID);

  // For each row of A/C
  size_t CSR_ip = 0, OLD_ip = 0;
  for (size_t i = 0; i < m; i++) {
    // First fill the c_status array w/ locations where we're allowed to
    // generate nonzeros for this row
    OLD_ip = Crowptr[i];
    CSR_ip = Crowptr[i+1];
    for (size_t k = OLD_ip; k < CSR_ip; k++) {
      c_status[Ccolind[k]] = k;

      // Reset values in the row of C
      Cvals[k] = SC_ZERO;
    }

    for (size_t k = Arowptr[i]; k < Arowptr[i+1]; k++) {
      LO Aik  = Acolind[k];
      const SC Aval = Avals[k];
      if (Aval == SC_ZERO)
        continue;

      if (targetMapToOrigRow[Aik] != LO_INVALID) {
        // Local matrix
        size_t Bk = static_cast<size_t> (targetMapToOrigRow[Aik]);

        for (size_t j = Browptr[Bk]; j < Browptr[Bk+1]; ++j) {
          LO Bkj = Bcolind[j];
          LO Cij = Bcol2Ccol[Bkj];

          TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
            std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
            "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

          Cvals[c_status[Cij]] += Aval * Bvals[j];
        }

      } else {
        // Remote matrix
        size_t Ik = static_cast<size_t> (targetMapToImportRow[Aik]);
        for (size_t j = Irowptr[Ik]; j < Irowptr[Ik+1]; ++j) {
          LO Ikj = Icolind[j];
          LO Cij = Icol2Ccol[Ikj];

          TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
            std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph " <<
            "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

          Cvals[c_status[Cij]] += Aval * Ivals[j];
        }
      }
    }
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  auto MM3 = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM Reuse ESFC"))));
#endif

  C.fillComplete(C.getDomainMap(), C.getRangeMap());
}


/*********************************************************************************************************/
// Kernel method for computing the local portion of C = (I-omega D^{-1} A)*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void jacobi_A_B_newmatrix(
  Scalar omega,
  const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& label,
  const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  //  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;

  typedef Import<LO,GO,NO>  import_type;
  typedef Map<LO,GO,NO>     map_type;
  typedef typename map_type::local_map_type local_map_type;

  // All of the Kokkos typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
  typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;


#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi M5 Cmap"))));
#endif
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // Build the final importer / column map, hash table lookups for C
  RCP<const import_type> Cimport;
  RCP<const map_type>    Ccolmap;
  RCP<const import_type> Bimport = Bview.origMatrix->getGraph()->getImporter();
  RCP<const import_type> Iimport = Bview.importMatrix.is_null() ?
      Teuchos::null : Bview.importMatrix->getGraph()->getImporter();
  local_map_type Acolmap_local = Aview.colMap->getLocalMap();
  local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
  local_map_type Irowmap_local;  if(!Bview.importMatrix.is_null()) Irowmap_local = Bview.importMatrix->getRowMap()->getLocalMap();
  local_map_type Bcolmap_local = Bview.origMatrix->getColMap()->getLocalMap();
  local_map_type Icolmap_local;  if(!Bview.importMatrix.is_null()) Icolmap_local = Bview.importMatrix->getColMap()->getLocalMap();

  // mfh 27 Sep 2016: Bcol2Ccol is a table that maps from local column
  // indices of B, to local column indices of C.  (B and C have the
  // same number of columns.)  The kernel uses this, instead of
  // copying the entire input matrix B and converting its column
  // indices to those of C.
  lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"),Bview.colMap->getLocalNumElements()), Icol2Ccol;

  if (Bview.importMatrix.is_null()) {
    // mfh 27 Sep 2016: B has no "remotes," so B and C have the same column Map.
    Cimport = Bimport;
    Ccolmap = Bview.colMap;
    // Bcol2Ccol is trivial
    // Bcol2Ccol is trivial

    Kokkos::RangePolicy<execution_space, LO> range (0, static_cast<LO> (Bview.colMap->getLocalNumElements ()));
    Kokkos::parallel_for (range, KOKKOS_LAMBDA (const size_t i) {
        Bcol2Ccol(i) = static_cast<LO> (i);
      });
  } else {
    // mfh 27 Sep 2016: B has "remotes," so we need to build the
    // column Map of C, as well as C's Import object (from its domain
    // Map to its column Map).  C's column Map is the union of the
    // column Maps of (the local part of) B, and the "remote" part of
    // B.  Ditto for the Import.  We have optimized this "setUnion"
    // operation on Import objects and Maps.

    // Choose the right variant of setUnion
    if (!Bimport.is_null() && !Iimport.is_null()){
      Cimport = Bimport->setUnion(*Iimport);
      Ccolmap = Cimport->getTargetMap();

    } else if (!Bimport.is_null() && Iimport.is_null()) {
      Cimport = Bimport->setUnion();

    } else if(Bimport.is_null() && !Iimport.is_null()) {
      Cimport = Iimport->setUnion();

    } else
      throw std::runtime_error("TpetraExt::Jacobi status of matrix importers is nonsensical");

    Ccolmap = Cimport->getTargetMap();

    TEUCHOS_TEST_FOR_EXCEPTION(!Cimport->getSourceMap()->isSameAs(*Bview.origMatrix->getDomainMap()),
      std::runtime_error, "Tpetra:Jacobi Import setUnion messed with the DomainMap in an unfortunate way");

    // NOTE: This is not efficient and should be folded into setUnion
    //
    // mfh 27 Sep 2016: What the above comment means, is that the
    // setUnion operation on Import objects could also compute these
    // local index - to - local index look-up tables.
    Kokkos::resize(Icol2Ccol,Bview.importMatrix->getColMap()->getLocalNumElements());
    local_map_type Ccolmap_local = Ccolmap->getLocalMap();
    Kokkos::parallel_for(range_type(0,Bview.origMatrix->getColMap()->getLocalNumElements()),KOKKOS_LAMBDA(const LO i) {
        Bcol2Ccol(i) = Ccolmap_local.getLocalElement(Bcolmap_local.getGlobalElement(i));
      });
    Kokkos::parallel_for(range_type(0,Bview.importMatrix->getColMap()->getLocalNumElements()),KOKKOS_LAMBDA(const LO i) {
        Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
      });

  }

  // Replace the column map
  //
  // mfh 27 Sep 2016: We do this because C was originally created
  // without a column Map.  Now we have its column Map.
  C.replaceColMap(Ccolmap);

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
  lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getLocalNumElements());
  lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),Aview.colMap->getLocalNumElements());
  Kokkos::parallel_for(range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),KOKKOS_LAMBDA(const LO i) {
      GO aidx = Acolmap_local.getGlobalElement(i);
      LO B_LID = Browmap_local.getLocalElement(aidx);
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
  KernelWrappers2<Scalar,LocalOrdinal,GlobalOrdinal,Node,lo_view_t>::jacobi_A_B_newmatrix_kernel_wrapper(omega,Dinv,Aview,Bview,targetMapToOrigRow,targetMapToImportRow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);

}


/*********************************************************************************************************/
// Jacobi AB NewMatrix Kernel wrappers (Default non-threaded version)
// Kernel method for computing the local portion of C = (I-omega D^{-1} A)*B

template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class LocalOrdinalViewType>
void KernelWrappers2<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalOrdinalViewType>::jacobi_A_B_newmatrix_kernel_wrapper(Scalar omega,
                                                           const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & Dinv,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                           CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                           const LocalOrdinalViewType & targetMapToOrigRow,
                                                           const LocalOrdinalViewType & targetMapToImportRow,
                                                           const LocalOrdinalViewType & Bcol2Ccol,
                                                           const LocalOrdinalViewType & Icol2Ccol,
                                                           CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                           Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > Cimport,
                                                           const std::string& label,
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params) {

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  auto MM(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Nemwmatrix SerialCore")));
#endif

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_host_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  // Jacobi-specific
  typedef typename scalar_view_t::memory_space scalar_memory_space;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;

  typedef Map<LO,GO,NO>     map_type;
  size_t ST_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getLocalNumRows();
  size_t n = Ccolmap->getLocalNumElements();
  size_t b_max_nnz_per_row = Bview.origMatrix->getLocalMaxNumRowEntries();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS & Amat = Aview.origMatrix->getLocalMatrixHost();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrixHost();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;

  c_lno_view_t  Irowptr;
  lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixHost();
    Irowptr = lclB.graph.row_map;
    Icolind = lclB.graph.entries;
    Ivals   = lclB.values;
    b_max_nnz_per_row = std::max(b_max_nnz_per_row,Bview.importMatrix->getLocalMaxNumRowEntries());
  }

  // Jacobi-specific inner stuff
  auto Dvals = 
       Dinv.template getLocalView<scalar_memory_space>(Access::ReadOnly);

  // Teuchos::ArrayView::operator[].
  // The status array will contain the index into colind where this entry was last deposited.
  // c_status[i] < CSR_ip - not in the row yet.
  // c_status[i] >= CSR_ip, this is the entry where you can find the data
  // We start with this filled with INVALID's indicating that there are no entries yet.
  // Sadly, this complicates the code due to the fact that size_t's are unsigned.
  size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
  Array<size_t> c_status(n, ST_INVALID);

  // Classic csr assembly (low memory edition)
  //
  // mfh 27 Sep 2016: C_estimate_nnz does not promise an upper bound.
  // The method loops over rows of A, and may resize after processing
  // each row.  Chris Siefert says that this reflects experience in
  // ML; for the non-threaded case, ML found it faster to spend less
  // effort on estimation and risk an occasional reallocation.
  size_t CSR_alloc = std::max(C_estimate_nnz(*Aview.origMatrix, *Bview.origMatrix), n);
  lno_view_t Crowptr(Kokkos::ViewAllocateWithoutInitializing("Crowptr"),m+1);
  lno_nnz_view_t Ccolind(Kokkos::ViewAllocateWithoutInitializing("Ccolind"),CSR_alloc);
  scalar_view_t Cvals(Kokkos::ViewAllocateWithoutInitializing("Cvals"),CSR_alloc);
  size_t CSR_ip = 0, OLD_ip = 0;

  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
  // routine.  The routine computes
  //
  // C := (I - omega * D^{-1} * A) * (B_local + B_remote)).
  //
  // This corresponds to one sweep of (weighted) Jacobi.
  //
  // For column index Aik in row i of A, targetMapToOrigRow[Aik] tells
  // you whether the corresponding row of B belongs to B_local
  // ("orig") or B_remote ("Import").

  // For each row of A/C
  for (size_t i = 0; i < m; i++) {
    // mfh 27 Sep 2016: m is the number of rows in the input matrix A
    // on the calling process.
    Crowptr[i] = CSR_ip;
    SC minusOmegaDval = -omega*Dvals(i,0);

    // Entries of B
    for (size_t j = Browptr[i]; j < Browptr[i+1]; j++) {
      Scalar Bval = Bvals[j];
      if (Bval == SC_ZERO)
        continue;
      LO Bij = Bcolind[j];
      LO Cij = Bcol2Ccol[Bij];

      // Assume no repeated entries in B
      c_status[Cij]   = CSR_ip;
      Ccolind[CSR_ip] = Cij;
      Cvals[CSR_ip]   = Bvals[j];
      CSR_ip++;
    }

    // Entries of -omega * Dinv * A * B
    for (size_t k = Arowptr[i]; k < Arowptr[i+1]; k++) {
      LO Aik  = Acolind[k];
      const SC Aval = Avals[k];
      if (Aval == SC_ZERO)
        continue;

      if (targetMapToOrigRow[Aik] != LO_INVALID) {
        // Local matrix
        size_t Bk = static_cast<size_t> (targetMapToOrigRow[Aik]);

        for (size_t j = Browptr[Bk]; j < Browptr[Bk+1]; ++j) {
          LO Bkj = Bcolind[j];
          LO Cij = Bcol2Ccol[Bkj];

          if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip) {
            // New entry
            c_status[Cij]   = CSR_ip;
            Ccolind[CSR_ip] = Cij;
            Cvals[CSR_ip]   = minusOmegaDval* Aval * Bvals[j];
            CSR_ip++;

          } else {
            Cvals[c_status[Cij]] += minusOmegaDval* Aval * Bvals[j];
          }
        }

      } else {
        // Remote matrix
        size_t Ik = static_cast<size_t> (targetMapToImportRow[Aik]);
        for (size_t j = Irowptr[Ik]; j < Irowptr[Ik+1]; ++j) {
          LO Ikj = Icolind[j];
          LO Cij = Icol2Ccol[Ikj];

          if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip) {
            // New entry
            c_status[Cij]   = CSR_ip;
            Ccolind[CSR_ip] = Cij;
            Cvals[CSR_ip]   = minusOmegaDval* Aval * Ivals[j];
            CSR_ip++;
          } else {
            Cvals[c_status[Cij]] += minusOmegaDval* Aval * Ivals[j];
          }
        }
      }
    }

    // Resize for next pass if needed
   if (i+1 < m && CSR_ip + std::min(n,(Arowptr[i+2]-Arowptr[i+1]+1)*b_max_nnz_per_row) > CSR_alloc) {
     CSR_alloc *= 2;
     Kokkos::resize(Ccolind,CSR_alloc);
     Kokkos::resize(Cvals,CSR_alloc);
    }
    OLD_ip = CSR_ip;
  }
  Crowptr[m] = CSR_ip;

  // Downward resize
  Kokkos::resize(Ccolind,CSR_ip);
  Kokkos::resize(Cvals,CSR_ip);

  {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      auto MM2(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix Final Sort")));
#endif

      // Replace the column map
      //
      // mfh 27 Sep 2016: We do this because C was originally created
      // without a column Map.  Now we have its column Map.
      C.replaceColMap(Ccolmap);

      // Final sort & set of CRS arrays
      //
      // TODO (mfh 27 Sep 2016) Will the thread-parallel "local" sparse
      // matrix-matrix multiply routine sort the entries for us?
      // Final sort & set of CRS arrays
      if (params.is_null() || params->get("sort entries",true))
          Import_Util::sortCrsEntries(Crowptr,Ccolind, Cvals);
      C.setAllValues(Crowptr,Ccolind, Cvals);
  }
  {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      auto MM3(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Newmatrix ESFC")));
#endif

      // Final FillComplete
      //
      // mfh 27 Sep 2016: So-called "expert static fill complete" bypasses
      // Import (from domain Map to column Map) construction (which costs
      // lots of communication) by taking the previously constructed
      // Import object.  We should be able to do this without interfering
      // with the implementation of the local part of sparse matrix-matrix
      // multply above
      RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
      labelList->set("Timer Label",label);
      if(!params.is_null()) labelList->set("compute global constants",params->get("compute global constants",true));
      RCP<const Export<LO,GO,NO> > dummyExport;
      C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport,dummyExport,labelList);

  }
}


/*********************************************************************************************************/
// Kernel method for computing the local portion of C = (I-omega D^{-1} A)*B
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void jacobi_A_B_reuse(
  Scalar omega,
  const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
  const std::string& label,
  const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;

  typedef Import<LO,GO,NO>  import_type;
  typedef Map<LO,GO,NO>     map_type;

  // Kokkos typedefs
  typedef typename map_type::local_map_type local_map_type;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename NO::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
  typedef Kokkos::View<LO*, typename lno_view_t::array_layout, typename lno_view_t::device_type> lo_view_t;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  RCP<TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Reuse Cmap"))));
#endif
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // Grab all the maps
  RCP<const import_type> Cimport = C.getGraph()->getImporter();
  RCP<const map_type>    Ccolmap = C.getColMap();
  local_map_type Acolmap_local = Aview.colMap->getLocalMap();
  local_map_type Browmap_local = Bview.origMatrix->getRowMap()->getLocalMap();
  local_map_type Irowmap_local;  if(!Bview.importMatrix.is_null()) Irowmap_local = Bview.importMatrix->getRowMap()->getLocalMap();
  local_map_type Bcolmap_local = Bview.origMatrix->getColMap()->getLocalMap();
  local_map_type Icolmap_local;  if(!Bview.importMatrix.is_null()) Icolmap_local = Bview.importMatrix->getColMap()->getLocalMap();
  local_map_type Ccolmap_local = Ccolmap->getLocalMap();

  // Build the final importer / column map, hash table lookups for C
  lo_view_t Bcol2Ccol(Kokkos::ViewAllocateWithoutInitializing("Bcol2Ccol"),Bview.colMap->getLocalNumElements()), Icol2Ccol;
  {
    // Bcol2Col may not be trivial, as Ccolmap is compressed during fillComplete in newmatrix
    // So, column map of C may be a strict subset of the column map of B
    Kokkos::parallel_for(range_type(0,Bview.origMatrix->getColMap()->getLocalNumElements()),KOKKOS_LAMBDA(const LO i) {
        Bcol2Ccol(i) = Ccolmap_local.getLocalElement(Bcolmap_local.getGlobalElement(i));
      });

    if (!Bview.importMatrix.is_null()) {
      TEUCHOS_TEST_FOR_EXCEPTION(!Cimport->getSourceMap()->isSameAs(*Bview.origMatrix->getDomainMap()),
                                 std::runtime_error, "Tpetra::Jacobi: Import setUnion messed with the DomainMap in an unfortunate way");

      Kokkos::resize(Icol2Ccol,Bview.importMatrix->getColMap()->getLocalNumElements());
      Kokkos::parallel_for(range_type(0,Bview.importMatrix->getColMap()->getLocalNumElements()),KOKKOS_LAMBDA(const LO i) {
          Icol2Ccol(i) = Ccolmap_local.getLocalElement(Icolmap_local.getGlobalElement(i));
        });
    }
  }

  // Run through all the hash table lookups once and for all
  lo_view_t targetMapToOrigRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToOrigRow"),Aview.colMap->getLocalNumElements());
  lo_view_t targetMapToImportRow(Kokkos::ViewAllocateWithoutInitializing("targetMapToImportRow"),Aview.colMap->getLocalNumElements());
  Kokkos::parallel_for(range_type(Aview.colMap->getMinLocalIndex(), Aview.colMap->getMaxLocalIndex()+1),KOKKOS_LAMBDA(const LO i) {
      GO aidx = Acolmap_local.getGlobalElement(i);
      LO B_LID = Browmap_local.getLocalElement(aidx);
      if (B_LID != LO_INVALID) {
        targetMapToOrigRow(i)   = B_LID;
        targetMapToImportRow(i) = LO_INVALID;
      } else {
        LO I_LID = Irowmap_local.getLocalElement(aidx);
        targetMapToOrigRow(i)   = LO_INVALID;
        targetMapToImportRow(i) = I_LID;

      }
    });

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
#endif

  // Call the actual kernel.  We'll rely on partial template specialization to call the correct one ---
  // Either the straight-up Tpetra code (SerialNode) or the KokkosKernels one (other NGP node types)
  KernelWrappers2<Scalar,LocalOrdinal,GlobalOrdinal,Node,lo_view_t>::jacobi_A_B_reuse_kernel_wrapper(omega,Dinv,Aview,Bview,targetMapToOrigRow,targetMapToImportRow,Bcol2Ccol,Icol2Ccol,C,Cimport,label,params);
}



/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node,
         class LocalOrdinalViewType>
void KernelWrappers2<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalOrdinalViewType>::jacobi_A_B_reuse_kernel_wrapper(Scalar omega,
                                                                                               const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Dinv,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                                                                               CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                                                                               const LocalOrdinalViewType & targetMapToOrigRow,
                                                                                               const LocalOrdinalViewType & targetMapToImportRow,
                                                                                               const LocalOrdinalViewType & Bcol2Ccol,
                                                                                               const LocalOrdinalViewType & Icol2Ccol,
                                                                                               CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                                                                               Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > /* Cimport */,
                                                                                               const std::string& label,
                                                                                               const Teuchos::RCP<Teuchos::ParameterList>& /* params */) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  Teuchos::RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Reuse SerialCore"))));
  Teuchos::RCP<Teuchos::TimeMonitor> MM2;
#else
  (void)label;
#endif
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_host_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;
  typedef typename scalar_view_t::memory_space scalar_memory_space;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;
  typedef Map<LO,GO,NO>     map_type;
  const size_t ST_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  const SC SC_ZERO = Teuchos::ScalarTraits<Scalar>::zero();

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m = Aview.origMatrix->getLocalNumRows();
  size_t n = Ccolmap->getLocalNumElements();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS & Amat = Aview.origMatrix->getLocalMatrixHost();
  const KCRS & Bmat = Bview.origMatrix->getLocalMatrixHost();
  const KCRS & Cmat = C.getLocalMatrixHost();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map, Crowptr = Cmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries, Ccolind = Cmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  scalar_view_t Cvals = Cmat.values;

  c_lno_view_t  Irowptr;
  lno_nnz_view_t  Icolind;
  scalar_view_t  Ivals;
  if(!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixHost();
    Irowptr = lclB.graph.row_map;
    Icolind = lclB.graph.entries;
    Ivals   = lclB.values;
  }

  // Jacobi-specific inner stuff
  auto Dvals = 
       Dinv.template getLocalView<scalar_memory_space>(Access::ReadOnly);

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM2 = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Reuse SerialCore - Compare"))));
#endif

  // The status array will contain the index into colind where this entry was last deposited.
  //   c_status[i] <  CSR_ip - not in the row yet
  //   c_status[i] >= CSR_ip - this is the entry where you can find the data
  // We start with this filled with INVALID's indicating that there are no entries yet.
  // Sadly, this complicates the code due to the fact that size_t's are unsigned.
  std::vector<size_t> c_status(n, ST_INVALID);

  // For each row of A/C
  size_t CSR_ip = 0, OLD_ip = 0;
  for (size_t i = 0; i < m; i++) {

    // First fill the c_status array w/ locations where we're allowed to
    // generate nonzeros for this row
    OLD_ip = Crowptr[i];
    CSR_ip = Crowptr[i+1];
    for (size_t k = OLD_ip; k < CSR_ip; k++) {
      c_status[Ccolind[k]] = k;

      // Reset values in the row of C
      Cvals[k] = SC_ZERO;
    }

    SC minusOmegaDval = -omega*Dvals(i,0);

    // Entries of B
    for (size_t j = Browptr[i]; j < Browptr[i+1]; j++) {
      Scalar Bval = Bvals[j];
      if (Bval == SC_ZERO)
        continue;
      LO Bij = Bcolind[j];
      LO Cij = Bcol2Ccol[Bij];

      TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
        std::runtime_error, "Trying to insert a new entry into a static graph");

      Cvals[c_status[Cij]] = Bvals[j];
    }

    // Entries of -omega * Dinv * A * B
    for (size_t k = Arowptr[i]; k < Arowptr[i+1]; k++) {
      LO Aik  = Acolind[k];
      const SC Aval = Avals[k];
      if (Aval == SC_ZERO)
        continue;

      if (targetMapToOrigRow[Aik] != LO_INVALID) {
        // Local matrix
        size_t Bk = static_cast<size_t> (targetMapToOrigRow[Aik]);

        for (size_t j = Browptr[Bk]; j < Browptr[Bk+1]; ++j) {
          LO Bkj = Bcolind[j];
          LO Cij = Bcol2Ccol[Bkj];

          TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
            std::runtime_error, "Trying to insert a new entry into a static graph");

          Cvals[c_status[Cij]] += minusOmegaDval * Aval * Bvals[j];
        }

      } else {
        // Remote matrix
        size_t Ik = static_cast<size_t> (targetMapToImportRow[Aik]);
        for (size_t j = Irowptr[Ik]; j < Irowptr[Ik+1]; ++j) {
          LO Ikj = Icolind[j];
          LO Cij = Icol2Ccol[Ikj];

          TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
            std::runtime_error, "Trying to insert a new entry into a static graph");

          Cvals[c_status[Cij]] += minusOmegaDval * Aval * Ivals[j];
        }
      }
    }
  }

#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM2 = Teuchos::null;
  MM2 = rcp(new TimeMonitor (*TimeMonitor::getNewTimer(prefix_mmm + std::string("Jacobi Reuse ESFC"))));
#endif

  C.fillComplete(C.getDomainMap(), C.getRangeMap());
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM2 = Teuchos::null;
  MM =  Teuchos::null;
#endif

}



/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void import_and_extract_views(
  const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>&   A,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >   targetMap,
  CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>&   Aview,
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> > prototypeImporter,
  bool                                                          userAssertsThereAreNoRemotes,
  const std::string&                                            label,
  const Teuchos::RCP<Teuchos::ParameterList>&                   params)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;

  typedef Map<LO,GO,NO>             map_type;
  typedef Import<LO,GO,NO>          import_type;
  typedef CrsMatrix<SC,LO,GO,NO>    crs_matrix_type;

#ifdef HAVE_TPETRA_MMM_TIMINGS
  std::string prefix_mmm = std::string("TpetraExt ") + label + std::string(": ");
  using Teuchos::TimeMonitor;
  RCP<Teuchos::TimeMonitor> MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM I&X Alloc"))));
#endif
  // The goal of this method is to populate the 'Aview' struct with views of the
  // rows of A, including all rows that correspond to elements in 'targetMap'.
  //
  // If targetMap includes local elements that correspond to remotely-owned rows
  // of A, then those remotely-owned rows will be imported into
  // 'Aview.importMatrix', and views of them will be included in 'Aview'.
  Aview.deleteContents();

  Aview.origMatrix   = rcp(&A, false);
  // trigger creation of int-typed row pointer array for use in TPLs, but don't actually need it here
  Aview.origMatrix->getApplyHelper(); 
  Aview.origRowMap   = A.getRowMap();
  Aview.rowMap       = targetMap;
  Aview.colMap       = A.getColMap();
  Aview.domainMap    = A.getDomainMap();
  Aview.importColMap = null;
  RCP<const map_type> rowMap = A.getRowMap();
  const int numProcs = rowMap->getComm()->getSize();

  // Short circuit if the user swears there are no remotes (or if we're in serial)
  if (userAssertsThereAreNoRemotes || numProcs < 2)
    return;

  RCP<const import_type> importer;
  if (params != null && params->isParameter("importer")) {
    importer = params->get<RCP<const import_type> >("importer");

  } else {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM I&X RemoteMap"))));
#endif

    // Mark each row in targetMap as local or remote, and go ahead and get a view
    // for the local rows
    RCP<const map_type> remoteRowMap;
    size_t numRemote = 0;
    int mode = 0;
    if (!prototypeImporter.is_null() &&
        prototypeImporter->getSourceMap()->isSameAs(*rowMap)     &&
        prototypeImporter->getTargetMap()->isSameAs(*targetMap)) {
      // We have a valid prototype importer --- ask it for the remotes
#ifdef HAVE_TPETRA_MMM_TIMINGS
      TimeMonitor MM2 = *TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM I&X RemoteMap-Mode1"));
#endif
      ArrayView<const LO> remoteLIDs = prototypeImporter->getRemoteLIDs();
      numRemote = prototypeImporter->getNumRemoteIDs();

      Array<GO> remoteRows(numRemote);
      for (size_t i = 0; i < numRemote; i++)
        remoteRows[i] = targetMap->getGlobalElement(remoteLIDs[i]);

      remoteRowMap = rcp(new map_type(Teuchos::OrdinalTraits<global_size_t>::invalid(), remoteRows(),
                                      rowMap->getIndexBase(), rowMap->getComm()));
      mode = 1;

    } else if (prototypeImporter.is_null()) {
      // No prototype importer --- count the remotes the hard way
#ifdef HAVE_TPETRA_MMM_TIMINGS
      TimeMonitor MM2 = *TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM I&X RemoteMap-Mode2"));
#endif
      ArrayView<const GO> rows    = targetMap->getLocalElementList();
      size_t              numRows = targetMap->getLocalNumElements();

      Array<GO> remoteRows(numRows);
      for(size_t i = 0; i < numRows; ++i) {
        const LO mlid = rowMap->getLocalElement(rows[i]);

        if (mlid == Teuchos::OrdinalTraits<LO>::invalid())
          remoteRows[numRemote++] = rows[i];
      }
      remoteRows.resize(numRemote);
      remoteRowMap = rcp(new map_type(Teuchos::OrdinalTraits<global_size_t>::invalid(), remoteRows(),
                                      rowMap->getIndexBase(), rowMap->getComm()));
      mode = 2;

    } else {
      // PrototypeImporter is bad.  But if we're in serial that's OK.
      mode = 3;
    }

    if (numProcs < 2) {
      TEUCHOS_TEST_FOR_EXCEPTION(numRemote > 0, std::runtime_error,
            "MatrixMatrix::import_and_extract_views ERROR, numProcs < 2 but attempting to import remote matrix rows.");
      // If only one processor we don't need to import any remote rows, so return.
      return;
    }

    //
    // Now we will import the needed remote rows of A, if the global maximum
    // value of numRemote is greater than 0.
    //
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM I&X Collective-0"))));
#endif

    global_size_t globalMaxNumRemote = 0;
    Teuchos::reduceAll(*(rowMap->getComm()), Teuchos::REDUCE_MAX, (global_size_t)numRemote, Teuchos::outArg(globalMaxNumRemote) );

    if (globalMaxNumRemote > 0) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
      MM = Teuchos::null;
      MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM I&X Import-2"))));
#endif
      // Create an importer with target-map remoteRowMap and source-map rowMap.
      if (mode == 1)
        importer = prototypeImporter->createRemoteOnlyImport(remoteRowMap);
      else if (mode == 2)
        importer = rcp(new import_type(rowMap, remoteRowMap));
      else
        throw std::runtime_error("prototypeImporter->SourceMap() does not match A.getRowMap()!");
    }

    if (params != null)
      params->set("importer", importer);
  }

  if (importer != null) {
#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM I&X Import-3"))));
#endif

    // Now create a new matrix into which we can import the remote rows of A that we need.
    Teuchos::ParameterList labelList;
    labelList.set("Timer Label", label);
    auto & labelList_subList = labelList.sublist("matrixmatrix: kernel params",false);

    bool isMM = true;
    bool overrideAllreduce = false;
    int mm_optimization_core_count=::Tpetra::Details::Behavior::TAFC_OptimizationCoreCount();
    // Minor speedup tweak - avoid computing the global constants
    Teuchos::ParameterList params_sublist;
    if(!params.is_null()) {
        labelList.set("compute global constants", params->get("compute global constants",false));
        params_sublist = params->sublist("matrixmatrix: kernel params",false);
        mm_optimization_core_count = params_sublist.get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
        int mm_optimization_core_count2 = params->get("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
        if(mm_optimization_core_count2<mm_optimization_core_count) mm_optimization_core_count=mm_optimization_core_count2;
        isMM = params_sublist.get("isMatrixMatrix_TransferAndFillComplete",false);
        overrideAllreduce = params_sublist.get("MM_TAFC_OverrideAllreduceCheck",false);
    }
    labelList_subList.set("isMatrixMatrix_TransferAndFillComplete",isMM);
    labelList_subList.set("MM_TAFC_OptimizationCoreCount",mm_optimization_core_count);
    labelList_subList.set("MM_TAFC_OverrideAllreduceCheck",overrideAllreduce);

    Aview.importMatrix = Tpetra::importAndFillCompleteCrsMatrix<crs_matrix_type>(rcpFromRef(A), *importer,
                                    A.getDomainMap(), importer->getTargetMap(), rcpFromRef(labelList));
    // trigger creation of int-typed row pointer array for use in TPLs, but don't actually need it here
    Aview.importMatrix->getApplyHelper();

#if 0
    // Disabled code for dumping input matrices
    static int count=0;
    char str[80];
    sprintf(str,"import_matrix.%d.dat",count);
    Tpetra::MatrixMarket::Writer<crs_matrix_type>::writeSparseFile(str,Aview.importMatrix);
    count++;
#endif

#ifdef HAVE_TPETRA_MMM_STATISTICS
    printMultiplicationStatistics(importer, label + std::string(" I&X MMM"));
#endif


#ifdef HAVE_TPETRA_MMM_TIMINGS
    MM = Teuchos::null;
    MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix_mmm + std::string("MMM I&X Import-4"))));
#endif

    // Save the column map of the imported matrix, so that we can convert indices back to global for arithmetic later
    Aview.importColMap = Aview.importMatrix->getColMap();
#ifdef HAVE_TPETRA_MMM_TIMINGS
   MM = Teuchos::null;
#endif
  }
}

/*********************************************************************************************************/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
void import_and_extract_views(
  const BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& M,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >      targetMap,
  BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Mview,
  Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >   prototypeImporter,
  bool                                                             userAssertsThereAreNoRemotes)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::null;

  typedef Scalar            SC;
  typedef LocalOrdinal      LO;
  typedef GlobalOrdinal     GO;
  typedef Node              NO;

  typedef Map<LO,GO,NO>             map_type;
  typedef Import<LO,GO,NO>          import_type;
  typedef BlockCrsMatrix<SC,LO,GO,NO>    blockcrs_matrix_type;

  // The goal of this method is to populate the 'Mview' struct with views of the
  // rows of M, including all rows that correspond to elements in 'targetMap'.
  //
  // If targetMap includes local elements that correspond to remotely-owned rows
  // of M, then those remotely-owned rows will be imported into
  // 'Mview.importMatrix', and views of them will be included in 'Mview'.
  Mview.deleteContents();

  Mview.origMatrix   = rcp(&M, false);
  // trigger creation of int-typed row pointer array for use in TPLs, but don't actually need it here
  Mview.origMatrix->getApplyHelper();
  Mview.origRowMap   = M.getRowMap();
  Mview.rowMap       = targetMap;
  Mview.colMap       = M.getColMap();
  Mview.importColMap = null;
  RCP<const map_type> rowMap = M.getRowMap();
  const int numProcs = rowMap->getComm()->getSize();

  // Short circuit if the user swears there are no remotes (or if we're in serial)
  if (userAssertsThereAreNoRemotes || numProcs < 2) return;

  // Mark each row in targetMap as local or remote, and go ahead and get a view
  // for the local rows
  RCP<const map_type> remoteRowMap;
  size_t numRemote = 0;
  int mode = 0;
  if (!prototypeImporter.is_null() &&
      prototypeImporter->getSourceMap()->isSameAs(*rowMap)     &&
      prototypeImporter->getTargetMap()->isSameAs(*targetMap)) {

      // We have a valid prototype importer --- ask it for the remotes
      ArrayView<const LO> remoteLIDs = prototypeImporter->getRemoteLIDs();
      numRemote = prototypeImporter->getNumRemoteIDs();

      Array<GO> remoteRows(numRemote);
      for (size_t i = 0; i < numRemote; i++)
        remoteRows[i] = targetMap->getGlobalElement(remoteLIDs[i]);

      remoteRowMap = rcp(new map_type(Teuchos::OrdinalTraits<global_size_t>::invalid(), remoteRows(),
                                      rowMap->getIndexBase(), rowMap->getComm()));
      mode = 1;

    } else if (prototypeImporter.is_null()) {

      // No prototype importer --- count the remotes the hard way
      ArrayView<const GO> rows    = targetMap->getLocalElementList();
      size_t              numRows = targetMap->getLocalNumElements();

      Array<GO> remoteRows(numRows);
      for(size_t i = 0; i < numRows; ++i) {
        const LO mlid = rowMap->getLocalElement(rows[i]);

        if (mlid == Teuchos::OrdinalTraits<LO>::invalid())
          remoteRows[numRemote++] = rows[i];
      }
      remoteRows.resize(numRemote);
      remoteRowMap = rcp(new map_type(Teuchos::OrdinalTraits<global_size_t>::invalid(), remoteRows(),
                                      rowMap->getIndexBase(), rowMap->getComm()));
      mode = 2;

    } else {
      // PrototypeImporter is bad.  But if we're in serial that's OK.
      mode = 3;
    }

    if (numProcs < 2) {
      TEUCHOS_TEST_FOR_EXCEPTION(numRemote > 0, std::runtime_error,
            "MatrixMatrix::import_and_extract_views ERROR, numProcs < 2 but attempting to import remote matrix rows.");
      // If only one processor we don't need to import any remote rows, so return.
      return;
    }

    // Now we will import the needed remote rows of M, if the global maximum
    // value of numRemote is greater than 0.
    global_size_t globalMaxNumRemote = 0;
    Teuchos::reduceAll(*(rowMap->getComm()), Teuchos::REDUCE_MAX, (global_size_t)numRemote, Teuchos::outArg(globalMaxNumRemote) );

    RCP<const import_type> importer;

    if (globalMaxNumRemote > 0) {
      // Create an importer with target-map remoteRowMap and source-map rowMap.
      if (mode == 1)
        importer = prototypeImporter->createRemoteOnlyImport(remoteRowMap);
      else if (mode == 2)
        importer = rcp(new import_type(rowMap, remoteRowMap));
      else
        throw std::runtime_error("prototypeImporter->SourceMap() does not match M.getRowMap()!");
    }

  if (importer != null) {
    // Get import matrix
    // TODO: create the int-typed row-pointer here
    Mview.importMatrix = Tpetra::importAndFillCompleteBlockCrsMatrix<blockcrs_matrix_type>(rcpFromRef(M), *importer);
    // trigger creation of int-typed row pointer array for use in TPLs, but don't actually need it here
    Mview.importMatrix->getApplyHelper();
    // Save the column map of the imported matrix, so that we can convert indices
    // back to global for arithmetic later
    Mview.importColMap = Mview.importMatrix->getColMap();
  }
}

/*********************************************************************************************************/
 // This only merges matrices that look like B & Bimport, aka, they have no overlapping rows
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalOrdinalViewType>
const typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type
merge_matrices(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                    CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                    const LocalOrdinalViewType & Acol2Brow,
                    const LocalOrdinalViewType & Acol2Irow,
                    const LocalOrdinalViewType & Bcol2Ccol,
                    const LocalOrdinalViewType & Icol2Ccol,
                    const size_t mergedNodeNumCols) {

  using Teuchos::RCP;
  typedef typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;
  // Grab the  Kokkos::SparseCrsMatrices
  const KCRS & Ak = Aview.origMatrix->getLocalMatrixDevice();
  const KCRS & Bk = Bview.origMatrix->getLocalMatrixDevice();

  // We need to do this dance if either (a) We have Bimport or (b) We don't A's colMap is not the same as B's rowMap
  if(!Bview.importMatrix.is_null() || (Bview.importMatrix.is_null() && (&*Aview.origMatrix->getGraph()->getColMap() != &*Bview.origMatrix->getGraph()->getRowMap()))) {
    // We do have a Bimport
    // NOTE: We're going merge Borig and Bimport into a single matrix and reindex the columns *before* we multiply.
    // This option was chosen because we know we don't have any duplicate entries, so we can allocate once.
    RCP<const KCRS> Ik_;
    if(!Bview.importMatrix.is_null()) Ik_ = Teuchos::rcpFromRef<const KCRS>(Bview.importMatrix->getLocalMatrixDevice());
    const KCRS * Ik     = Bview.importMatrix.is_null() ? 0 : &*Ik_;
    KCRS Iks;
    if(Ik!=0) Iks = *Ik;
    size_t merge_numrows =  Ak.numCols();
    // The last entry of this at least, need to be initialized
    lno_view_t Mrowptr("Mrowptr", merge_numrows + 1);

    const LocalOrdinal LO_INVALID =Teuchos::OrdinalTraits<LocalOrdinal>::invalid();

    // Use a Kokkos::parallel_scan to build the rowptr
    typedef typename Node::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
    Kokkos::parallel_scan ("Tpetra_MatrixMatrix_merge_matrices_buildRowptr", range_type (0, merge_numrows),
      KOKKOS_LAMBDA(const size_t i, size_t& update, const bool final) {
        if(final) Mrowptr(i) = update;
        // Get the row count
        size_t ct=0;
        if(Acol2Brow(i)!=LO_INVALID)
          ct = Bk.graph.row_map(Acol2Brow(i)+1) - Bk.graph.row_map(Acol2Brow(i));
        else
          ct = Iks.graph.row_map(Acol2Irow(i)+1) - Iks.graph.row_map(Acol2Irow(i));
        update+=ct;

        if(final && i+1==merge_numrows)
          Mrowptr(i+1)=update;
      });

    // Allocate nnz
    size_t merge_nnz = ::Tpetra::Details::getEntryOnHost(Mrowptr,merge_numrows);
    lno_nnz_view_t Mcolind(Kokkos::ViewAllocateWithoutInitializing("Mcolind"),merge_nnz);
    scalar_view_t Mvalues(Kokkos::ViewAllocateWithoutInitializing("Mvals"),merge_nnz);

    // Use a Kokkos::parallel_for to fill the rowptr/colind arrays
    typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
    Kokkos::parallel_for ("Tpetra_MatrixMatrix_merg_matrices_buildColindValues", range_type (0, merge_numrows),KOKKOS_LAMBDA(const size_t i) {
        if(Acol2Brow(i)!=LO_INVALID) {
          size_t row   = Acol2Brow(i);
          size_t start = Bk.graph.row_map(row);
          for(size_t j= Mrowptr(i); j<Mrowptr(i+1); j++) {
            Mvalues(j) = Bk.values(j-Mrowptr(i)+start);
            Mcolind(j) = Bcol2Ccol(Bk.graph.entries(j-Mrowptr(i)+start));
          }
        }
        else {
          size_t row   = Acol2Irow(i);
          size_t start = Iks.graph.row_map(row);
          for(size_t j= Mrowptr(i); j<Mrowptr(i+1); j++) {
            Mvalues(j) = Iks.values(j-Mrowptr(i)+start);
            Mcolind(j) = Icol2Ccol(Iks.graph.entries(j-Mrowptr(i)+start));
          }
        }
      });

    KCRS newmat("CrsMatrix",merge_numrows,mergedNodeNumCols,merge_nnz,Mvalues,Mrowptr,Mcolind);
    return newmat;
  }
  else {
    // We don't have a Bimport (the easy case)
    return Bk;
  }
}//end merge_matrices

/*********************************************************************************************************/
 // This only merges matrices that look like B & Bimport, aka, they have no overlapping rows
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalOrdinalViewType>
const typename Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type
merge_matrices(BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
               BlockCrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
               const LocalOrdinalViewType & Acol2Brow,
               const LocalOrdinalViewType & Acol2Irow,
               const LocalOrdinalViewType & Bcol2Ccol,
               const LocalOrdinalViewType & Icol2Ccol,
               const size_t mergedNodeNumCols) 
{
  using Teuchos::RCP;
  typedef typename Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_device_type KBCRS;
  typedef typename KBCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KBCRS::values_type::non_const_type scalar_view_t;

  // Grab the  KokkosSparse::BsrMatrix
  const KBCRS & Ak = Aview.origMatrix->getLocalMatrixDevice();
  const KBCRS & Bk = Bview.origMatrix->getLocalMatrixDevice();

  // We need to do this dance if either (a) We have Bimport or (b) A's colMap is not the same as B's rowMap
  if(!Bview.importMatrix.is_null() ||
     (Bview.importMatrix.is_null() &&
     (&*Aview.origMatrix->getGraph()->getColMap() != &*Bview.origMatrix->getGraph()->getRowMap()))) {
    
    // We do have a Bimport
    // NOTE: We're going merge Borig and Bimport into a single matrix and reindex the columns *before* we multiply.
    // This option was chosen because we know we don't have any duplicate entries, so we can allocate once.
    RCP<const KBCRS> Ik_;
    if(!Bview.importMatrix.is_null()) Ik_ = Teuchos::rcpFromRef<const KBCRS>(Bview.importMatrix->getLocalMatrixDevice());
    const KBCRS * Ik = Bview.importMatrix.is_null() ? 0 : &*Ik_;
    KBCRS Iks;
    if(Ik!=0) Iks = *Ik;
    size_t merge_numrows =  Ak.numCols();
    
    // The last entry of this at least, need to be initialized
    lno_view_t Mrowptr("Mrowptr", merge_numrows + 1);

    const LocalOrdinal LO_INVALID =Teuchos::OrdinalTraits<LocalOrdinal>::invalid();

    // Use a Kokkos::parallel_scan to build the rowptr
    typedef typename Node::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
    Kokkos::parallel_scan ("Tpetra_MatrixMatrix_merge_matrices_buildRowptr", range_type (0, merge_numrows),
      KOKKOS_LAMBDA(const size_t i, size_t& update, const bool final) {
        if(final) Mrowptr(i) = update;
        // Get the row count
        size_t ct=0;
        if(Acol2Brow(i)!=LO_INVALID)
          ct = Bk.graph.row_map(Acol2Brow(i)+1) - Bk.graph.row_map(Acol2Brow(i));
        else
          ct = Iks.graph.row_map(Acol2Irow(i)+1) - Iks.graph.row_map(Acol2Irow(i));
        update+=ct;

        if(final && i+1==merge_numrows)
          Mrowptr(i+1)=update;
      });

    // Allocate nnz
    size_t merge_nnz = ::Tpetra::Details::getEntryOnHost(Mrowptr,merge_numrows);
    const int blocksize = Ak.blockDim();
    lno_nnz_view_t Mcolind(Kokkos::ViewAllocateWithoutInitializing("Mcolind"),merge_nnz);
    scalar_view_t Mvalues(Kokkos::ViewAllocateWithoutInitializing("Mvals"),merge_nnz*blocksize*blocksize);

    // Use a Kokkos::parallel_for to fill the rowptr/colind arrays
    typedef Kokkos::RangePolicy<execution_space, size_t> range_type;
    Kokkos::parallel_for ("Tpetra_MatrixMatrix_merg_matrices_buildColindValues", range_type (0, merge_numrows),KOKKOS_LAMBDA(const size_t i) {
        if(Acol2Brow(i)!=LO_INVALID) {
          size_t row   = Acol2Brow(i);
          size_t start = Bk.graph.row_map(row);
          for(size_t j= Mrowptr(i); j<Mrowptr(i+1); j++) {
            Mcolind(j) = Bcol2Ccol(Bk.graph.entries(j-Mrowptr(i)+start));

            for (int b=0; b<blocksize*blocksize; ++b) {
              const int val_indx = j*blocksize*blocksize + b;
              const int b_val_indx = (j-Mrowptr(i)+start)*blocksize*blocksize + b;
              Mvalues(val_indx) = Bk.values(b_val_indx);
            }
          }
        }
        else {
          size_t row   = Acol2Irow(i);
          size_t start = Iks.graph.row_map(row);
          for(size_t j= Mrowptr(i); j<Mrowptr(i+1); j++) {
            Mcolind(j) = Icol2Ccol(Iks.graph.entries(j-Mrowptr(i)+start));

            for (int b=0; b<blocksize*blocksize; ++b) {
              const int val_indx = j*blocksize*blocksize + b;
              const int b_val_indx = (j-Mrowptr(i)+start)*blocksize*blocksize + b;
              Mvalues(val_indx) = Iks.values(b_val_indx);
            }
          }
        }
      });

    // Build and return merged KokkosSparse matrix
    KBCRS newmat("CrsMatrix",merge_numrows,mergedNodeNumCols,merge_nnz,Mvalues,Mrowptr,Mcolind, blocksize);
    return newmat;
  }
  else {
    // We don't have a Bimport (the easy case)
    return Bk;
  }
}//end merge_matrices

/*********************************************************************************************************/
template<typename SC, typename LO, typename GO, typename NO>
void AddKernels<SC, LO, GO, NO>::
addSorted(
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::values_array& Avals,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::row_ptrs_array_const& Arowptrs,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::col_inds_array& Acolinds,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::impl_scalar_type scalarA,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::values_array& Bvals,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::row_ptrs_array_const& Browptrs,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::col_inds_array& Bcolinds,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::impl_scalar_type scalarB,
#if KOKKOSKERNELS_VERSION >= 40299
  GO numGlobalCols,
#endif
  typename MMdetails::AddKernels<SC, LO, GO, NO>::values_array& Cvals,
  typename MMdetails::AddKernels<SC, LO, GO, NO>::row_ptrs_array& Crowptrs,
  typename MMdetails::AddKernels<SC, LO, GO, NO>::col_inds_array& Ccolinds)
{
  using Teuchos::TimeMonitor;
  using AddKern = MMdetails::AddKernels<SC, LO, GO, NO>;
  TEUCHOS_TEST_FOR_EXCEPTION(Arowptrs.extent(0) != Browptrs.extent(0), std::runtime_error, "Can't add matrices with different numbers of rows.");
  auto nrows = Arowptrs.extent(0) - 1;
  Crowptrs = row_ptrs_array(Kokkos::ViewAllocateWithoutInitializing("C row ptrs"), nrows + 1);
  typename AddKern::KKH handle;
  handle.create_spadd_handle(true);
  auto addHandle = handle.get_spadd_handle();
#ifdef HAVE_TPETRA_MMM_TIMINGS
  auto MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("TpetraExt::MatrixMatrix::add() sorted symbolic")));
#endif
  KokkosSparse::Experimental::spadd_symbolic
    (&handle,
#if KOKKOSKERNELS_VERSION >= 40299
     nrows, numGlobalCols,
#endif
     Arowptrs, Acolinds, Browptrs, Bcolinds, Crowptrs);
  //KokkosKernels requires values to be zeroed
  Cvals = values_array("C values", addHandle->get_c_nnz());
  Ccolinds = col_inds_array(Kokkos::ViewAllocateWithoutInitializing("C colinds"), addHandle->get_c_nnz());
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("TpetraExt::MatrixMatrix::add() sorted numeric")));
#endif
  KokkosSparse::Experimental::spadd_numeric(&handle,
#if KOKKOSKERNELS_VERSION >= 40299
     nrows, numGlobalCols,
#endif
    Arowptrs, Acolinds, Avals, scalarA,
    Browptrs, Bcolinds, Bvals, scalarB,
    Crowptrs, Ccolinds, Cvals);
}

template<typename SC, typename LO, typename GO, typename NO>
void AddKernels<SC, LO, GO, NO>::
addUnsorted(
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::values_array& Avals,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::row_ptrs_array_const& Arowptrs,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::col_inds_array& Acolinds,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::impl_scalar_type scalarA,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::values_array& Bvals,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::row_ptrs_array_const& Browptrs,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::col_inds_array& Bcolinds,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::impl_scalar_type scalarB,
#if KOKKOSKERNELS_VERSION >= 40299
  GO numGlobalCols,
#else
  GO /* numGlobalCols */,
#endif
  typename MMdetails::AddKernels<SC, LO, GO, NO>::values_array& Cvals,
  typename MMdetails::AddKernels<SC, LO, GO, NO>::row_ptrs_array& Crowptrs,
  typename MMdetails::AddKernels<SC, LO, GO, NO>::col_inds_array& Ccolinds)
{
  using Teuchos::TimeMonitor;
  using AddKern = MMdetails::AddKernels<SC, LO, GO, NO>;
  TEUCHOS_TEST_FOR_EXCEPTION(Arowptrs.extent(0) != Browptrs.extent(0), std::runtime_error, "Can't add matrices with different numbers of rows.");
  auto nrows = Arowptrs.extent(0) - 1;
  Crowptrs = row_ptrs_array(Kokkos::ViewAllocateWithoutInitializing("C row ptrs"), nrows + 1);
  typedef MMdetails::AddKernels<SC, LO, GO, NO> AddKern;
  typename AddKern::KKH handle;
  handle.create_spadd_handle(false);
  auto addHandle = handle.get_spadd_handle();
#ifdef HAVE_TPETRA_MMM_TIMINGS
  auto MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("TpetraExt::MatrixMatrix::add() unsorted symbolic")));
#endif
  KokkosSparse::Experimental::spadd_symbolic
      (&handle,
#if KOKKOSKERNELS_VERSION >= 40299
       nrows, numGlobalCols,
#endif
       Arowptrs, Acolinds, Browptrs, Bcolinds, Crowptrs);
  //Cvals must be zeroed out
  Cvals = values_array("C values", addHandle->get_c_nnz());
  Ccolinds = col_inds_array(Kokkos::ViewAllocateWithoutInitializing("C colinds"), addHandle->get_c_nnz());
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("TpetraExt::MatrixMatrix::add() unsorted kernel: unsorted numeric")));
#endif
  KokkosSparse::Experimental::spadd_numeric(&handle,
#if KOKKOSKERNELS_VERSION >= 40299
    nrows, numGlobalCols,
#endif
    Arowptrs, Acolinds, Avals, scalarA,
    Browptrs, Bcolinds, Bvals, scalarB,
    Crowptrs, Ccolinds, Cvals);
}

template<typename GO,
         typename LocalIndicesType,
         typename GlobalIndicesType,
         typename ColMapType>
struct ConvertLocalToGlobalFunctor
{
  ConvertLocalToGlobalFunctor(
      const LocalIndicesType& colindsOrig_,
      const GlobalIndicesType& colindsConverted_,
      const ColMapType& colmap_) :
    colindsOrig (colindsOrig_),
    colindsConverted (colindsConverted_),
    colmap (colmap_)
  {}
  KOKKOS_INLINE_FUNCTION void
  operator() (const GO i) const
  {
    colindsConverted(i) = colmap.getGlobalElement(colindsOrig(i));
  }
  LocalIndicesType colindsOrig;
  GlobalIndicesType colindsConverted;
  ColMapType colmap;
};

template<typename SC, typename LO, typename GO, typename NO>
void MMdetails::AddKernels<SC, LO, GO, NO>::
convertToGlobalAndAdd(
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::KCRS& A,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::impl_scalar_type scalarA,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::KCRS& B,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::impl_scalar_type scalarB,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::local_map_type& AcolMap,
  const typename MMdetails::AddKernels<SC, LO, GO, NO>::local_map_type& BcolMap,
  typename MMdetails::AddKernels<SC, LO, GO, NO>::values_array& Cvals,
  typename MMdetails::AddKernels<SC, LO, GO, NO>::row_ptrs_array& Crowptrs,
  typename MMdetails::AddKernels<SC, LO, GO, NO>::global_col_inds_array& Ccolinds)
{
  using Teuchos::TimeMonitor;
  //Need to use a different KokkosKernelsHandle type than other versions,
  //since the ordinals are now GO
  using KKH_GO = KokkosKernels::Experimental::KokkosKernelsHandle<size_t, GO, impl_scalar_type,
              typename NO::execution_space, typename NO::memory_space, typename NO::memory_space>;

  const values_array Avals = A.values;
  const values_array Bvals = B.values;
  const col_inds_array Acolinds = A.graph.entries;
  const col_inds_array Bcolinds = B.graph.entries;
  auto Arowptrs = A.graph.row_map;
  auto Browptrs = B.graph.row_map;
  global_col_inds_array AcolindsConverted(Kokkos::ViewAllocateWithoutInitializing("A colinds (converted)"), Acolinds.extent(0));
  global_col_inds_array BcolindsConverted(Kokkos::ViewAllocateWithoutInitializing("B colinds (converted)"), Bcolinds.extent(0));
#ifdef HAVE_TPETRA_MMM_TIMINGS
  auto MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("TpetraExt::MatrixMatrix::add() diff col map kernel: " + std::string("column map conversion"))));
#endif
  ConvertLocalToGlobalFunctor<GO, col_inds_array, global_col_inds_array, local_map_type> convertA(Acolinds, AcolindsConverted, AcolMap);
  Kokkos::parallel_for("Tpetra_MatrixMatrix_convertColIndsA", range_type(0, Acolinds.extent(0)), convertA);
  ConvertLocalToGlobalFunctor<GO, col_inds_array, global_col_inds_array, local_map_type> convertB(Bcolinds, BcolindsConverted, BcolMap);
  Kokkos::parallel_for("Tpetra_MatrixMatrix_convertColIndsB", range_type(0, Bcolinds.extent(0)), convertB);
  KKH_GO handle;
  handle.create_spadd_handle(false);
  auto addHandle = handle.get_spadd_handle();
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("TpetraExt::MatrixMatrix::add() diff col map kernel: unsorted symbolic")));
#endif
  auto nrows = Arowptrs.extent(0) - 1;
  Crowptrs = row_ptrs_array(Kokkos::ViewAllocateWithoutInitializing("C row ptrs"), nrows + 1);
  KokkosSparse::Experimental::spadd_symbolic
    (&handle,
#if KOKKOSKERNELS_VERSION >= 40299
     nrows, A.numCols(),
#endif
     Arowptrs, AcolindsConverted, Browptrs, BcolindsConverted, Crowptrs);
  Cvals = values_array("C values", addHandle->get_c_nnz());
  Ccolinds = global_col_inds_array(Kokkos::ViewAllocateWithoutInitializing("C colinds"), addHandle->get_c_nnz());
#ifdef HAVE_TPETRA_MMM_TIMINGS
  MM = Teuchos::null;
  MM = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("TpetraExt::MatrixMatrix::add() diff col map kernel: unsorted numeric")));
#endif
  KokkosSparse::Experimental::spadd_numeric(&handle,
#if KOKKOSKERNELS_VERSION >= 40299
    nrows, A.numCols(),
#endif
    Arowptrs, AcolindsConverted, Avals, scalarA,
    Browptrs, BcolindsConverted, Bvals, scalarB,
    Crowptrs, Ccolinds, Cvals);
}


} //End namepsace MMdetails

} //End namespace Tpetra

/*********************************************************************************************************/
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//
namespace Tpetra {

#define TPETRA_MATRIXMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
template \
  void MatrixMatrix::Multiply( \
    const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
    bool transposeA, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& B, \
    bool transposeB, \
    CrsMatrix< SCALAR , LO , GO , NODE >& C, \
    bool call_FillComplete_on_result, \
    const std::string & label, \
    const Teuchos::RCP<Teuchos::ParameterList>& params); \
\
template \
  void MatrixMatrix::Multiply( \
    const Teuchos::RCP<const BlockCrsMatrix< SCALAR , LO , GO , NODE > >& A, \
    bool transposeA, \
    const Teuchos::RCP<const BlockCrsMatrix< SCALAR , LO , GO , NODE > >& B, \
    bool transposeB, \
    Teuchos::RCP<BlockCrsMatrix< SCALAR , LO , GO , NODE > >& C,	\
    const std::string & label); \
\
template \
  void MatrixMatrix::Jacobi( \
    SCALAR omega, \
    const Vector< SCALAR, LO, GO, NODE > & Dinv, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& B, \
    CrsMatrix< SCALAR , LO , GO , NODE >& C, \
    bool call_FillComplete_on_result, \
    const std::string & label, \
    const Teuchos::RCP<Teuchos::ParameterList>& params); \
\
  template \
  void MatrixMatrix::Add( \
    const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
    bool transposeA, \
    SCALAR scalarA, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& B, \
    bool transposeB, \
    SCALAR scalarB, \
    Teuchos::RCP<CrsMatrix< SCALAR , LO , GO , NODE > >& C); \
\
  template \
  void MatrixMatrix::Add( \
    const CrsMatrix< SCALAR , LO , GO , NODE >& A, \
    bool transposeA, \
    SCALAR scalarA, \
    const CrsMatrix< SCALAR , LO , GO , NODE >& B, \
    bool transposeB, \
    SCALAR scalarB, \
    const Teuchos::RCP<CrsMatrix< SCALAR , LO , GO , NODE > >& C); \
\
  template \
  void MatrixMatrix::Add( \
    const CrsMatrix<SCALAR, LO, GO, NODE>& A, \
    bool transposeA, \
    SCALAR scalarA, \
    CrsMatrix<SCALAR, LO, GO, NODE>& B, \
    SCALAR scalarB ); \
\
  template \
  Teuchos::RCP<CrsMatrix< SCALAR , LO , GO , NODE > > \
  MatrixMatrix::add<SCALAR, LO, GO, NODE> \
                    (const SCALAR & alpha, \
                     const bool transposeA, \
                     const CrsMatrix<SCALAR, LO, GO, NODE>& A, \
                     const SCALAR & beta, \
                     const bool transposeB, \
                     const CrsMatrix<SCALAR, LO, GO, NODE>& B, \
                     const Teuchos::RCP<const Map<LO, GO, NODE> >& domainMap, \
                     const Teuchos::RCP<const Map<LO, GO, NODE> >& rangeMap, \
                     const Teuchos::RCP<Teuchos::ParameterList>& params); \
\
  template \
  void \
  MatrixMatrix::add< SCALAR , LO, GO , NODE > \
                    (const SCALAR & alpha, \
                     const bool transposeA, \
                     const CrsMatrix< SCALAR , LO, GO , NODE >& A, \
                     const SCALAR& beta, \
                     const bool transposeB, \
                     const CrsMatrix< SCALAR , LO, GO , NODE >& B, \
                     CrsMatrix< SCALAR , LO, GO , NODE >& C, \
                     const Teuchos::RCP<const Map<LO, GO , NODE > >& domainMap, \
                     const Teuchos::RCP<const Map<LO, GO , NODE > >& rangeMap, \
                     const Teuchos::RCP<Teuchos::ParameterList>& params); \
\
  template struct MMdetails::AddKernels<SCALAR, LO, GO, NODE>; \
\
  template void MMdetails::import_and_extract_views<SCALAR, LO, GO, NODE>(const CrsMatrix<SCALAR, LO, GO, NODE>& M, \
                                                                          Teuchos::RCP<const Map<LO, GO, NODE> > targetMap, \
                                                                          CrsMatrixStruct<SCALAR, LO, GO, NODE>& Mview, \
                                                                          Teuchos::RCP<const Import<LO,GO, NODE> > prototypeImporter, \
                                                                          bool userAssertsThereAreNoRemotes, \
                                                                          const std::string& label, \
                                                                          const Teuchos::RCP<Teuchos::ParameterList>& params); \
\
  template void MMdetails::import_and_extract_views<SCALAR, LO, GO, NODE>(const BlockCrsMatrix<SCALAR, LO, GO, NODE>& M, \
                                                                          Teuchos::RCP<const Map<LO, GO, NODE> > targetMap, \
                                                                          BlockCrsMatrixStruct<SCALAR, LO, GO, NODE>& Mview, \
                                                                          Teuchos::RCP<const Import<LO,GO, NODE> > prototypeImporter, \
                                                                          bool userAssertsThereAreNoRemotes);
} //End namespace Tpetra

#endif // TPETRA_MATRIXMATRIX_DEF_HPP
