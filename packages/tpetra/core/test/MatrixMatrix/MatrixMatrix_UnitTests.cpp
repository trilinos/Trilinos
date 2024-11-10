// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "TpetraExt_MatrixMatrix.hpp"
#include "Teuchos_Assert.hpp"
#include "TpetraExt_TripleMatrixMultiply.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Import_Util.hpp"
#include <cmath>

namespace {
  bool verbose = false;
  std::string matnamesFile;

  using Tpetra::MatrixMarket::Reader;
  using Tpetra::CrsMatrix;
  using Tpetra::BlockCrsMatrix;
  using Tpetra::CrsMatrixMultiplyOp;
  using Tpetra::global_size_t;
  using Tpetra::Map;
  using Tpetra::Import;
  using Tpetra::RowMatrixTransposer;
  using Tpetra::Vector;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Teuchos::null;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::tuple;
  using std::endl;

TEUCHOS_STATIC_SETUP(){
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("matnames-file", &matnamesFile,
    "A file containing a list of matricies we'll import", true);
  clp.setOption("v", "not-verbose", &verbose,
    "Whether or not to use verbose output");
}

template<class SC, class LO, class GO,class NT>
RCP<CrsMatrix<SC, LO, GO, NT> >
getIdentityMatrix (Teuchos::FancyOStream& out,
                   const Tpetra::global_size_t globalNumRows,
                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> Matrix_t;
  typedef Tpetra::Map<LO, GO, NT> Map_t;

  Teuchos::OSTab tab0 (out);
  out << "getIdentityMatrix" << endl;
  Teuchos::OSTab tab1 (out);

  out << "Create row Map" << endl;
  RCP<const Map_t> identityRowMap =
    Tpetra::createUniformContigMapWithNode<LO, GO, NT> (globalNumRows, comm);

  out << "Create CrsMatrix" << endl;
  RCP<Matrix_t> identityMatrix =
    Tpetra::createCrsMatrix<SC, LO, GO, NT> (identityRowMap, 1);

  out << "Fill CrsMatrix" << endl;
  Teuchos::ArrayView<const GO> gblRows = identityRowMap->getLocalElementList ();
  for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
    Teuchos::Array<GO> col (1, *it);
    Teuchos::Array<SC> val (1, Teuchos::ScalarTraits<SC>::one ());
    identityMatrix->insertGlobalValues (*it, col (), val ());
  }

  out << "Call fillComplete" << endl;
  identityMatrix->fillComplete ();

  out << "Done!" << endl;
  return identityMatrix;
}

template<class SC, class LO, class GO,class NT>
RCP<CrsMatrix<SC, LO, GO, NT> >
getIdentityMatrixWithMap (Teuchos::FancyOStream& out,
                          Teuchos::RCP<const Tpetra::Map<LO,GO,NT> >& identityRowMap,
                          const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> Matrix_t;

  Teuchos::OSTab tab0 (out);
  out << "getIdentityMatrixWithMap" << endl;
  Teuchos::OSTab tab1 (out);

  out << "Create CrsMatrix" << endl;
  RCP<Matrix_t> identityMatrix =
    Tpetra::createCrsMatrix<SC, LO, GO, NT> (identityRowMap, 1);

  out << "Fill CrsMatrix" << endl;
  Teuchos::ArrayView<const GO> gblRows = identityRowMap->getLocalElementList ();
  for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
    Teuchos::Array<GO> col (1, *it);
    Teuchos::Array<SC> val (1, Teuchos::ScalarTraits<SC>::one ());
    identityMatrix->insertGlobalValues (*it, col (), val ());
  }

  out << "Call fillComplete" << endl;
  identityMatrix->fillComplete ();

  out << "Done!" << endl;
  return identityMatrix;
}

template<class Matrix_t>
RCP<Matrix_t>
copyMatrix(const RCP<Matrix_t> &A) {
  return rcp(new Matrix_t(*A,Teuchos::Copy));
}

typedef struct add_test_results_struct{
  double correctNorm;
  double computedNorm;
  double epsilon;
} add_test_results;

typedef struct mult_test_results_struct{
  double epsilon;
  double cNorm;
  double compNorm;
  bool   isImportValid;
} mult_test_results;


template<class Matrix_t>
add_test_results regular_add_test(
    const std::string& name,
    RCP<Matrix_t > A,
    RCP<Matrix_t > B,
    bool AT,
    bool BT,
    RCP<Matrix_t > C,
    RCP<const Comm<int> > comm)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;

  add_test_results toReturn;
  toReturn.correctNorm = C->getFrobeniusNorm ();

  RCP<const Map_t > rowmap = AT ? A->getDomainMap() : A->getRowMap();
  size_t estSize = A->getGlobalMaxNumRowEntries() + B->getGlobalMaxNumRowEntries();
         // estSize is upper bound for A, B; estimate only for AT, BT.
  RCP<Matrix_t> computedC = rcp( new Matrix_t(rowmap, estSize));

  SC one = Teuchos::ScalarTraits<SC>::one();
  Tpetra::MatrixMatrix::Add(*A, AT, one, *B, BT, one, computedC);
  computedC->fillComplete(A->getDomainMap(), A->getRangeMap());
  toReturn.computedNorm = computedC->getFrobeniusNorm ();
  toReturn.epsilon = fabs(toReturn.correctNorm - toReturn.computedNorm);

#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedC);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", C);
#endif


  return toReturn;
}

/// \brief Test the three-argument (A, B, C) version of CrsMatrix add,
///   where the output argument C is null on input.
///
/// \tparam Matrix_t A specialization of Tpetra::CrsMatrix.
template<class Matrix_t>
add_test_results
null_add_test_1 (const Matrix_t& A,
                 const Matrix_t& B,
                 const bool AT,
                 const bool BT,
                 const Matrix_t& C,
                 Teuchos::FancyOStream& out,
                 bool& success)
{
  typedef typename Matrix_t::scalar_type scalar_type;
  typedef typename Matrix_t::local_ordinal_type local_ordinal_type;
  typedef typename Matrix_t::global_ordinal_type global_ordinal_type;
  typedef typename Matrix_t::node_type NT;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, NT> map_type;
  typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, NT> export_type;
  const scalar_type one = STS::one ();

  RCP<const Comm<int> > comm = A.getMap ()->getComm ();
  const int myRank = comm->getRank ();

  out << "  Computing Frobenius norm of the expected result C" << endl;
  add_test_results toReturn;
  toReturn.correctNorm = C.getFrobeniusNorm ();

  out << "  Calling 3-arg add" << endl;
  RCP<const map_type> domainMap = BT ? B.getRangeMap () : B.getDomainMap ();
  RCP<const map_type> rangeMap = BT ? B.getDomainMap () : B.getRangeMap ();
  RCP<Matrix_t> C_computed;
  // for each MPI process to catch any exception message
  std::ostringstream errStrm;
  int lclSuccess = 1;
  int gblSuccess = 0; // output argument
  try {
    C_computed =
      Tpetra::MatrixMatrix::add (one, AT, A, one, BT, B,
                                 domainMap, rangeMap, Teuchos::null);
  }
  catch (std::exception& e) {
    errStrm << "Proc " << myRank << ": add threw an exception: "
      << e.what () << endl;
    lclSuccess = 0;
  }
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *comm);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    C_computed.is_null (), std::logic_error, "3-arg add returned null.");

  RCP<Matrix_t> C_exported;
  if (! C_computed->getRowMap ()->isSameAs (* (C.getRowMap ()))) {
    // Export C_computed to C's row Map, so we can compare the two.
    export_type exp (C_computed->getRowMap (), C.getRowMap ());
    C_exported =
      Tpetra::exportAndFillCompleteCrsMatrix<Matrix_t> (C_computed, exp,
                                                             C.getDomainMap (),
                                                             C.getRangeMap ());
  } else {
    C_exported = C_computed;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getRowMap ()->isSameAs (* (C.getRowMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different row Maps.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getDomainMap ()->isSameAs (* (C.getDomainMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different domain Maps.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getRangeMap ()->isSameAs (* (C.getRangeMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different range Maps.");

  toReturn.computedNorm = C_exported->getFrobeniusNorm ();
  toReturn.epsilon = STS::magnitude (toReturn.correctNorm - toReturn.computedNorm);
  return toReturn;
}

/// \brief Test the three-argument (A, B, C) version of CrsMatrix Add,
///   where the output argument C is null on input.
///
/// \tparam Matrix_t A specialization of Tpetra::CrsMatrix.
template<class Matrix_t>
add_test_results
null_add_test_2 (const Matrix_t& A,
                 const Matrix_t& B,
                 const bool AT,
                 const bool BT,
                 const Matrix_t& C,
                 Teuchos::FancyOStream& out,
                 bool& success)
{
  typedef typename Matrix_t::scalar_type scalar_type;
  typedef typename Matrix_t::local_ordinal_type local_ordinal_type;
  typedef typename Matrix_t::global_ordinal_type global_ordinal_type;
  typedef typename Matrix_t::node_type NT;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, NT> map_type;
  typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, NT> export_type;
  const scalar_type one = STS::one ();

  RCP<const Comm<int> > comm = A.getMap ()->getComm ();
  const int myRank = comm->getRank ();

  out << "  Computing Frobenius norm of the expected result C" << endl;
  add_test_results toReturn;
  toReturn.correctNorm = C.getFrobeniusNorm ();

  out << "  Calling 3-arg add" << endl;
  RCP<const map_type> domainMap = BT ? B.getRangeMap () : B.getDomainMap ();
  RCP<const map_type> rangeMap = BT ? B.getDomainMap () : B.getRangeMap ();
  RCP<Matrix_t> C_computed;
  // for each MPI process to catch any exception message
  std::ostringstream errStrm;
  int lclSuccess = 1;
  int gblSuccess = 0; // output argument
  try {
    Tpetra::MatrixMatrix::Add (A, AT, one, B, BT, one, C_computed);
  }
  catch (std::exception& e) {
    errStrm << "Proc " << myRank << ": add threw an exception: "
      << e.what () << endl;
    lclSuccess = 0;
  }
  // MatrixMatrix::Add, with C_computed null on input, should not call fillComplete.
  TEST_ASSERT(C_computed->isFillActive());
  // Call fillComplete now so that we can compare against C.
  C_computed->fillComplete(C.getDomainMap(), C.getRangeMap());
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *comm);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    C_computed.is_null (), std::logic_error, "3-arg add returned null.");

  RCP<Matrix_t> C_exported;
  if (! C_computed->getRowMap ()->isSameAs (* (C.getRowMap ()))) {
    // Export C_computed to C's row Map, so we can compare the two.
    export_type exp (C_computed->getRowMap (), C.getRowMap ());
    C_exported =
      Tpetra::exportAndFillCompleteCrsMatrix<Matrix_t> (C_computed, exp,
                                                             C.getDomainMap (),
                                                             C.getRangeMap ());
  } else {
    C_exported = C_computed;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getRowMap ()->isSameAs (* (C.getRowMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different row Maps.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getDomainMap ()->isSameAs (* (C.getDomainMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different domain Maps.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getRangeMap ()->isSameAs (* (C.getRangeMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different range Maps.");

  toReturn.computedNorm = C_exported->getFrobeniusNorm ();
  toReturn.epsilon = STS::magnitude (toReturn.correctNorm - toReturn.computedNorm);
  return toReturn;
}

template<class Matrix_t>
add_test_results add_into_test(
    RCP<Matrix_t > A,
    RCP<Matrix_t > B,
    bool AT,
    RCP<Matrix_t > C,
    RCP<const Comm<int> > comm)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;

  add_test_results toReturn;
  toReturn.correctNorm = C->getFrobeniusNorm ();

  RCP<const Map_t > rowmap =
    AT ? A->getDomainMap () : A->getRowMap ();
  RCP<Matrix_t> computedC = rcp (new Matrix_t (rowmap, 1));
  SC one = Teuchos::ScalarTraits<SC>::one();
  Tpetra::MatrixMatrix::Add (*A, AT, one, *B, one);
  B->fillComplete ();
  toReturn.computedNorm = B->getFrobeniusNorm ();
  toReturn.epsilon = fabs (toReturn.correctNorm - toReturn.computedNorm);

  return toReturn;
}


template<class Matrix_t>
add_test_results reuse_add_test(
    const std::string& name,
    RCP<Matrix_t > A,
    RCP<Matrix_t > B,
    bool AT,
    bool BT,
    RCP<Matrix_t > C,
    RCP<const Comm<int> > comm)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;

  add_test_results toReturn;
  toReturn.correctNorm = C->getFrobeniusNorm ();

  RCP<const Map_t > rowmap = AT ? A->getDomainMap() : A->getRowMap();
  size_t estSize = A->getGlobalMaxNumRowEntries() + B->getGlobalMaxNumRowEntries();
         // estSize is upper bound for A, B; estimate only for AT, BT.
  RCP<Matrix_t> computedC = rcp( new Matrix_t(rowmap, estSize));

  SC one = Teuchos::ScalarTraits<SC>::one();
  Tpetra::MatrixMatrix::Add(*A, AT, one, *B, BT, one, computedC);
  computedC->fillComplete(A->getDomainMap(), A->getRangeMap());

  // Call Add a second time
  Tpetra::MatrixMatrix::Add(*A, AT, one, *B, BT, one, computedC);

  TEUCHOS_ASSERT(A->getDomainMap()->isSameAs(*computedC->getDomainMap()));
  TEUCHOS_ASSERT(A->getRangeMap()->isSameAs(*computedC->getRangeMap()));

  toReturn.computedNorm = computedC->getFrobeniusNorm ();
  toReturn.epsilon = fabs(toReturn.correctNorm - toReturn.computedNorm);

#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedC);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", C);
#endif


  return toReturn;
}


template<class Matrix_t>
mult_test_results multiply_test_manualfc(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  bool AT,
  bool BT,
  RCP<Matrix_t> C,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;
  typedef Import<LO,GO,NT> Import_t;
  RCP<const Map_t> map = C->getRowMap();

  size_t maxNumEntriesPerRow = C->getGlobalMaxNumRowEntries();
  RCP<Matrix_t> computedC = rcp( new Matrix_t(map, maxNumEntriesPerRow));

  Tpetra::MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC, false);
  computedC->fillComplete(C->getDomainMap(), C->getRangeMap());

#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedC);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", C);
#endif
  SC one = Teuchos::ScalarTraits<SC>::one();

  RCP<Matrix_t> diffMatrix = Tpetra::createCrsMatrix<SC,LO,GO,NT>(C->getRowMap(),
                                                                  maxNumEntriesPerRow);
  Tpetra::MatrixMatrix::Add(*C, false, -one, *computedC, false, one, diffMatrix);
  diffMatrix->fillComplete(C->getDomainMap(), C->getRangeMap());

  mult_test_results results;
  results.cNorm    = C->getFrobeniusNorm ();
  results.compNorm = diffMatrix->getFrobeniusNorm ();
  results.epsilon  = results.compNorm/results.cNorm;

  // Check importer validity
  RCP<const Import_t> myImport = C->getGraph()->getImporter();
  if(!myImport.is_null()) results.isImportValid=Tpetra::Import_Util::checkImportValidity(*myImport);
  else results.isImportValid=true;

  return results;
}

template<class Matrix_t>
mult_test_results multiply_test_autofc(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  bool AT,
  bool BT,
  RCP<Matrix_t> C,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;
  typedef Import<LO,GO,NT> Import_t;
  RCP<const Map_t> map = C->getRowMap();

  RCP<Matrix_t> computedC = rcp( new Matrix_t(map, 0));
  SC one = Teuchos::ScalarTraits<SC>::one();

  Tpetra::MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC, true);

  if(!C->getDomainMap()->isSameAs(*computedC->getDomainMap())) throw std::runtime_error("Domain map mismatch");
  if(!C->getRangeMap()->isSameAs(*computedC->getRangeMap())) throw std::runtime_error("Range map mismatch");
  if(!C->getRowMap()->isSameAs(*computedC->getRowMap())) throw std::runtime_error("Row map mismatch");

#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedC);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", C);
#endif

#if 0
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  std::cout<<"*** C->colMap() ***"<<std::endl;
  C->getColMap()->describe(*fancy,Teuchos::VERB_EXTREME);
  std::cout<<"*** C->domainMap() ***"<<std::endl;
  C->getDomainMap()->describe(*fancy,Teuchos::VERB_EXTREME);

  if(C->getGraph()->getImporter().is_null()) std::cout<<"C->getImporter is null"<<std::endl;
  else {
    std::cout<<"*** C->getImporter()->getTargetMap() ***"<<std::endl;
    C->getGraph()->getImporter()->getTargetMap()->describe(*fancy,Teuchos::VERB_EXTREME);
  }

  std::cout<<"*** computedC->colMap() ***"<<std::endl;
  computedC->getColMap()->describe(*fancy,Teuchos::VERB_EXTREME);
  std::cout<<"*** computedC->domainMap() ***"<<std::endl;
  computedC->getDomainMap()->describe(*fancy,Teuchos::VERB_EXTREME);

  if(computedC->getGraph()->getImporter().is_null()) std::cout<<"computedC->getImporter is null"<<std::endl;
  else {
    std::cout<<"*** computedC->getImporter()->getTargetMap() ***"<<std::endl;
    computedC->getGraph()->getImporter()->getTargetMap()->describe(*fancy,Teuchos::VERB_EXTREME);
  }
#endif


  RCP<Matrix_t> diffMatrix =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(C->getRowMap(), C->getGlobalMaxNumRowEntries());
  Tpetra::MatrixMatrix::Add(*C, false, -one, *computedC, false, one, diffMatrix);
  diffMatrix->fillComplete(C->getDomainMap(), C->getRangeMap());

  mult_test_results results;
  results.cNorm    = C->getFrobeniusNorm ();
  results.compNorm = diffMatrix->getFrobeniusNorm ();
  results.epsilon  = results.compNorm/results.cNorm;

  if(results.epsilon > 1e-3) {
#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedC);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", C);
#endif
    if(!comm->getRank()) printf("ERROR: TEST %s FAILED\n",name.c_str());
  }

  // Check importer validity
  RCP<const Import_t> myImport = C->getGraph()->getImporter();
  if(!myImport.is_null()) results.isImportValid=Tpetra::Import_Util::checkImportValidity(*myImport);
  else results.isImportValid=true;


  return results;
}


template<class Matrix_t>
mult_test_results multiply_RAP_test_autofc(
  const std::string& name,
  RCP<Matrix_t> R,
  RCP<Matrix_t> A,
  RCP<Matrix_t> P,
  bool RT,
  bool AT,
  bool PT,
  RCP<Matrix_t> Ac,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;

  RCP<const Map_t> map = Ac->getRowMap();

  RCP<Matrix_t> computedAc = rcp( new Matrix_t(map, 0));
  SC one = Teuchos::ScalarTraits<SC>::one();

  Tpetra::TripleMatrixMultiply::MultiplyRAP(*R, RT, *A, AT, *P, PT, *computedAc, true);

  if(!Ac->getDomainMap()->isSameAs(*computedAc->getDomainMap())) throw std::runtime_error("Domain map mismatch");
  if(!Ac->getRangeMap()->isSameAs(*computedAc->getRangeMap())) throw std::runtime_error("Range map mismatch");
  if(!Ac->getRowMap()->isSameAs(*computedAc->getRowMap())) throw std::runtime_error("Row map mismatch");

#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedAc);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", Ac);
#endif

#if 0
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  std::cout<<"*** Ac->colMap() ***"<<std::endl;
  Ac->getColMap()->describe(*fancy,Teuchos::VERB_EXTREME);
  std::cout<<"*** Ac->domainMap() ***"<<std::endl;
  Ac->getDomainMap()->describe(*fancy,Teuchos::VERB_EXTREME);

  if(Ac->getGraph()->getImporter().is_null()) std::cout<<"Ac->getImporter is null"<<std::endl;
  else {
    std::cout<<"*** Ac->getImporter()->getTargetMap() ***"<<std::endl;
    Ac->getGraph()->getImporter()->getTargetMap()->describe(*fancy,Teuchos::VERB_EXTREME);
  }

  std::cout<<"*** computedAc->colMap() ***"<<std::endl;
  computedAc->getColMap()->describe(*fancy,Teuchos::VERB_EXTREME);
  std::cout<<"*** computedAc->domainMap() ***"<<std::endl;
  computedAc->getDomainMap()->describe(*fancy,Teuchos::VERB_EXTREME);

  if(computedAc->getGraph()->getImporter().is_null()) std::cout<<"computedAc->getImporter is null"<<std::endl;
  else {
    std::cout<<"*** computedAc->getImporter()->getTargetMap() ***"<<std::endl;
    computedAc->getGraph()->getImporter()->getTargetMap()->describe(*fancy,Teuchos::VERB_EXTREME);
  }
#endif


  RCP<Matrix_t> diffMatrix =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(Ac->getRowMap(),
                                         Ac->getGlobalMaxNumRowEntries());
  Tpetra::MatrixMatrix::Add(*Ac, false, -one, *computedAc, false, one, diffMatrix);
  diffMatrix->fillComplete(Ac->getDomainMap(), Ac->getRangeMap());

  mult_test_results results;
  results.cNorm    = Ac->getFrobeniusNorm ();
  results.compNorm = diffMatrix->getFrobeniusNorm ();
  results.epsilon  = results.compNorm/results.cNorm;

  if(results.epsilon > 1e-3) {
#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedAc);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", Ac);

#endif
  if(!comm->getRank()) printf("ERROR: TEST %s FAILED\n",name.c_str());
  }


  return results;
}

template<class Matrix_t>
mult_test_results multiply_RAP_reuse_test(
  const std::string& name,
  RCP<Matrix_t> R,
  RCP<Matrix_t> A,
  RCP<Matrix_t> P,
  bool RT,
  bool AT,
  bool PT,
  RCP<Matrix_t> Ac,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;

  RCP<const Map_t> map = Ac->getRowMap();

  RCP<Matrix_t> computedC1 = rcp( new Matrix_t(map, 0));
  SC one = Teuchos::ScalarTraits<SC>::one();

  // First cut
  Tpetra::TripleMatrixMultiply::MultiplyRAP(*R, RT, *A, AT, *P, PT, *computedC1);

  if(!Ac->getDomainMap()->isSameAs(*computedC1->getDomainMap())) throw std::runtime_error("Domain map mismatch");
  if(!Ac->getRangeMap()->isSameAs(*computedC1->getRangeMap())) throw std::runtime_error("Range map mismatch");
  if(!Ac->getRowMap()->isSameAs(*computedC1->getRowMap())) throw std::runtime_error("Row map mismatch");

  RCP<Matrix_t> computedC2 = rcp( new Matrix_t(computedC1->getCrsGraph()) );
  computedC2->fillComplete(computedC2->getDomainMap(), computedC2->getRangeMap());
  computedC2->resumeFill();
  Tpetra::TripleMatrixMultiply::MultiplyRAP(*R, RT, *A, AT, *P, PT, *computedC2);

  // Only check the second "reuse" matrix
  RCP<Matrix_t> diffMatrix =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(map, Ac->getGlobalMaxNumRowEntries());
  Tpetra::MatrixMatrix::Add(*Ac, false, -one, *computedC2, false, one, diffMatrix);
  diffMatrix->fillComplete(Ac->getDomainMap(), Ac->getRangeMap());

  mult_test_results results;
  results.cNorm    = Ac->getFrobeniusNorm ();
  results.compNorm = diffMatrix->getFrobeniusNorm ();
  results.epsilon  = results.compNorm/results.cNorm;

  if(results.epsilon > 1e-3) {
    if(!comm->getRank()) printf("ERROR: TEST %s FAILED\n",name.c_str());
  }


  return results;
}




template<class Matrix_t>
mult_test_results multiply_reuse_test(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  bool AT,
  bool BT,
  RCP<Matrix_t> C,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;
  typedef Import<LO,GO,NT> Import_t;
  typedef Vector<SC,LO,GO,NT> Vector_t;

  RCP<const Map_t> map = C->getRowMap();

  // Scaling vectors
  Teuchos::Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);

  RCP<Vector_t> leftScaling  = rcp( new Vector_t(C->getRangeMap()) );
  leftScaling->randomize();
  leftScaling->norm2(norms);
  leftScaling->scale(1.0/norms[0]);

  RCP<Vector_t> rightScaling = rcp( new Vector_t(C->getDomainMap()) );
  rightScaling->randomize();
  rightScaling->norm2(norms);
  rightScaling->scale(1.0/norms[0]);

  // As = leftScaling * op(A) =
  //   leftScaling * A, if AT=false
  //   A*leftScaling,   if AT=true
  RCP<Matrix_t> As = copyMatrix(A);
  if (AT == false) As->leftScale (*leftScaling);
  else             As->rightScale(*leftScaling);

  // Bs = op(B) * rightScaling =
  //   B * rightScaling, if BT=false
  //   rightScaling*B,   if BT=true
  RCP<Matrix_t> Bs = copyMatrix(B);
  if (BT == false) Bs->rightScale(*rightScaling);
  else             Bs->leftScale (*rightScaling);

  // computedC2 = op(As) * op(Bs)
  RCP<Matrix_t> computedC2 = rcp( new Matrix_t(C->getCrsGraph()) );
  computedC2->fillComplete(C->getDomainMap(), C->getRangeMap());

  computedC2->resumeFill();
  Tpetra::MatrixMatrix::Multiply(*As, AT, *Bs, BT, *computedC2, true/*call_FillCompleteOnResult*/);

  // computedC1 = leftScaling * (op(A)*op(B)) * rightScaling
  RCP<Matrix_t> computedC1 = rcp( new Matrix_t(map, computedC2->getGlobalMaxNumRowEntries()));
  Tpetra::MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC1, false/*call_FillCompleteOnResult*/);
  computedC1->fillComplete(C->getDomainMap(), C->getRangeMap());
  computedC1->leftScale (*leftScaling);
  computedC1->rightScale(*rightScaling);

  // diffMatrix = computedC2 - computedC1
  SC one = Teuchos::ScalarTraits<SC>::one();
  RCP<Matrix_t> diffMatrix =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(C->getRowMap(),
                                         computedC2->getGlobalMaxNumRowEntries());

  Kokkos::fence();
  Tpetra::MatrixMatrix::Add(*computedC1, false, -one, *computedC2, false, one, diffMatrix);
  diffMatrix->fillComplete(C->getDomainMap(), C->getRangeMap());

  mult_test_results results;
  results.cNorm    = C->getFrobeniusNorm ();
  results.compNorm = diffMatrix->getFrobeniusNorm ();
  results.epsilon  = results.compNorm/results.cNorm;

  // Check importer validity
  RCP<const Import_t> myImport = C->getGraph()->getImporter();
  if(!myImport.is_null()) results.isImportValid=Tpetra::Import_Util::checkImportValidity(*myImport);
  else results.isImportValid=true;

  return results;
}


template<class Matrix_t>
mult_test_results jacobi_test(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Vector<SC,LO,GO,NT> Vector_t;
  typedef Map<LO,GO,NT> Map_t;
  RCP<const Map_t> map = A->getRowMap();
  SC one = Teuchos::ScalarTraits<SC>::one();
  SC omega=Teuchos::ScalarTraits<SC>::one();
  Vector_t Dinv(B->getRowMap());
  Dinv.putScalar(1.0);


  // Jacobi version
  RCP<Matrix_t> C = rcp(new Matrix_t(B->getRowMap(),0));
  Tpetra::MatrixMatrix::Jacobi<SC, LO, GO, NT>(omega,Dinv,*A,*B,*C);

  // Multiply + Add version
  RCP<Matrix_t> C2;
  bool done=false;
#ifdef HAVE_TPETRA_INST_OPENMP
    if(std::is_same<NT,Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>::value) {
      Teuchos::ParameterList p;
      //p.set("openmp: jacobi algorithm","MSAK");
      p.set("openmp: jacobi algorithm","KK");
      C2 = rcp(new Matrix_t(B->getRowMap(),0));
      Tpetra::MatrixMatrix::Jacobi<SC, LO, GO, NT>(omega,Dinv,*A,*B,*C2,true,"jacobi_test_msak",rcp(&p,false));
      done=true;
    }
#endif
#ifdef HAVE_TPETRA_INST_CUDA
    if(std::is_same<NT,Tpetra::KokkosCompat::KokkosCudaWrapperNode>::value) {
      Teuchos::ParameterList p;
      //p.set("cuda: jacobi algorithm","MSAK");
      p.set("cuda: jacobi algorithm","KK");
      C2 = rcp(new Matrix_t(B->getRowMap(),0));
      Tpetra::MatrixMatrix::Jacobi<SC, LO, GO, NT>(omega,Dinv,*A,*B,*C2,true,"jacobi_test_msak",rcp(&p,false));
      done=true;
    }
#endif
    // Fallback
    if(!done) {
      Dinv.putScalar(omega);
      RCP<Matrix_t> AB = rcp(new Matrix_t(B->getRowMap(),0));

      Tpetra::MatrixMatrix::Multiply(*A,false,*B,false,*AB);
      AB->leftScale(Dinv);
      C2 = rcp(new Matrix_t(B->getRowMap(),
                            B->getGlobalMaxNumRowEntries() + AB->getGlobalMaxNumRowEntries())); // upper bound
      Tpetra::MatrixMatrix::Add(*AB,false,-one,*B,false,one,C2);
      if(!C2->isFillComplete()) C2->fillComplete(C->getDomainMap(),C->getRangeMap());
    }

  // Check the difference
  RCP<Matrix_t> C_check = rcp(new Matrix_t(B->getRowMap(), C->getGlobalMaxNumRowEntries()));
  Tpetra::MatrixMatrix::Add(*C, false, -one, *C2, false,one,C_check);
  C_check->fillComplete(B->getDomainMap(),B->getRangeMap());

  // Error Check
  double compNorm = C_check->getFrobeniusNorm();
  double cNorm = C->getFrobeniusNorm();
  mult_test_results results;
  results.epsilon  = compNorm/cNorm;
  results.cNorm    = cNorm;
  results.compNorm = compNorm;
  return results;
}


template<class Matrix_t>
mult_test_results jacobi_reuse_test(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Vector<SC,LO,GO,NT> Vector_t;
  typedef Map<LO,GO,NT> Map_t;

  RCP<const Map_t> map = A->getRowMap();

  // Scaling vectors
  Teuchos::Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
  RCP<Vector_t> rightScaling = rcp( new Vector_t(B->getDomainMap()) );
  rightScaling->randomize();
  rightScaling->norm2(norms);
  rightScaling->scale(1.0/norms[0]);
  SC one = Teuchos::ScalarTraits<SC>::one();

  SC omega = one;
  Vector_t Dinv(B->getRowMap());
  Dinv.putScalar(1.0);

  // Jacobi version
  RCP<Matrix_t> computedC1 = rcp(new Matrix_t(B->getRowMap(), 0));
  Tpetra::MatrixMatrix::Jacobi(omega, Dinv, *A, *B, *computedC1);
  computedC1->rightScale(*rightScaling);

  // Bs = B * rightScaling
  RCP<Matrix_t> Bs = copyMatrix(B);
  Bs->rightScale(*rightScaling);

  // computedC2 = (I - Dinv*A)*Bs
  RCP<Matrix_t> computedC2 = rcp( new Matrix_t(computedC1->getCrsGraph()) );
  Tpetra::MatrixMatrix::Jacobi(omega, Dinv, *A, *Bs, *computedC2);

  // diffMatrix = computedC2 - computedC1

  RCP<Matrix_t> diffMatrix =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(computedC1->getRowMap(),
                                         computedC1->getGlobalMaxNumRowEntries());
  Tpetra::MatrixMatrix::Add(*computedC1, false, -one, *computedC2, false, one, diffMatrix);
  diffMatrix->fillComplete(computedC1->getDomainMap(), computedC1->getRangeMap());

  mult_test_results results;
  results.cNorm    = computedC1->getFrobeniusNorm ();
  results.compNorm = diffMatrix->getFrobeniusNorm ();
  results.epsilon  = results.compNorm/results.cNorm;

  return results;
}


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMat, operations_test,SC,LO, GO, NT)  {
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  // NOTE: The matrix reader doesn't read real matrices into a complex data type, so we just swap down to MT here
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
  typedef Map<LO,GO,NT>                     map_type;
  typedef CrsMatrix<MT,LO,GO,NT> Matrix_t;
  const int myRank = comm->getRank ();
  //const int numProcs = comm->getSize();

  // Create an output stream that prints immediately, rather than
  // waiting until the test finishes.  This lets us get debug output
  // even if an unexpected exception or crash happens.  Unfortunately,
  // Teuchos::FancyOStream does not implement operator=, so we can't
  // replace 'out' temporarily.
  RCP<FancyOStream> newOutP = (myRank == 0) ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
    Teuchos::getFancyOStream (rcp (new Teuchos::oblackholestream ()));
  FancyOStream& newOut = *newOutP;

  newOut << "Tpetra sparse matrix-matrix multiply: operations_test" << endl;
  Teuchos::OSTab tab1 (newOut);

  newOut << "Get parameters from XML file" << endl;
  Teuchos::RCP<Teuchos::ParameterList> matrixSystems =
    Teuchos::getParametersFromXmlFile(matnamesFile);

  for (Teuchos::ParameterList::ConstIterator it = matrixSystems->begin();
       it != matrixSystems->end();
       ++it) {
    TEUCHOS_TEST_FOR_EXCEPTION(!it->second.isList(), std::runtime_error,
      "All top-level items in the list of matrix file names must be "
      "ParameterList instances.  " << endl <<
      "Bad tag's name: " << it->first <<
      "Type name: " << it->second.getAny().typeName() <<
      endl << endl);

    ParameterList currentSystem = matrixSystems->sublist (it->first);
    std::string name = currentSystem.name();
    std::string A_file = currentSystem.get<std::string> ("A");
    std::string A_domainmap_file = currentSystem.get<std::string> ("A_domainmap", "");
    std::string A_rangemap_file = currentSystem.get<std::string> ("A_rangemap", "");
    std::string A_rowmap_file = currentSystem.get<std::string> ("A_rowmap", "");
    std::string A_colmap_file = currentSystem.get<std::string> ("A_colmap", "");

    std::string B_file = currentSystem.get<std::string> ("B");
    std::string B_domainmap_file = currentSystem.get<std::string> ("B_domainmap", "");
    std::string B_rangemap_file = currentSystem.get<std::string> ("B_rangemap", "");
    std::string B_rowmap_file = currentSystem.get<std::string> ("B_rowmap", "");
    std::string B_colmap_file = currentSystem.get<std::string> ("B_colmap", "");

    std::string C_file = currentSystem.get<std::string> ("C");
    std::string C_domainmap_file = currentSystem.get<std::string> ("C_domainmap", "");
    std::string C_rangemap_file = currentSystem.get<std::string> ("C_rangemap", "");
    std::string C_rowmap_file = currentSystem.get<std::string> ("C_rowmap", "");
    std::string C_colmap_file = currentSystem.get<std::string> ("C_colmap", "");

    std::string D_file = currentSystem.get<std::string> ("D", "");
    std::string D_domainmap_file = currentSystem.get<std::string> ("D_domainmap", "");
    std::string D_rangemap_file = currentSystem.get<std::string> ("D_rangemap", "");
    std::string D_rowmap_file = currentSystem.get<std::string> ("D_rowmap", "");
    std::string D_colmap_file = currentSystem.get<std::string> ("D_colmap", "");

    const bool AT = currentSystem.get<bool> ("TransA");
    const bool BT = currentSystem.get<bool> ("TransB");
    const bool CT = currentSystem.get<bool> ("TransC", false);
    const int numProcs = currentSystem.get<int> ("numProcs", 0);

    // Don't run tests where the core count does not match the maps
    if(numProcs > 0 && comm->getSize() != numProcs) continue;


    double epsilon = currentSystem.get<double> ("epsilon",
                                       100. * Teuchos::ScalarTraits<SC>::eps());
    std::string op = currentSystem.get<std::string> ("op");

    RCP<Matrix_t> A, B, C, D;

    if (A_file != ""){
      if (A_domainmap_file == "" || A_rangemap_file == "" || A_rowmap_file == "" || A_colmap_file == "") {
        out << "Attempt to read sparse matrix file \"" << A_file
            << "\"" << endl;
        A = Reader<Matrix_t>::readSparseFile (A_file, comm);
      }
      else {
        RCP<const map_type> domainmap = Reader<Matrix_t>::readMapFile (A_domainmap_file, comm);
        RCP<const map_type> rangemap  = Reader<Matrix_t>::readMapFile (A_rangemap_file, comm);
        RCP<const map_type> rowmap    = Reader<Matrix_t>::readMapFile (A_rowmap_file, comm);
        RCP<const map_type> colmap    = Reader<Matrix_t>::readMapFile (A_colmap_file, comm);
        out << "Attempt to read sparse matrix file \""
            << A_file << "\"" << endl;
        A = Reader<Matrix_t>::readSparseFile (A_file, rowmap, colmap, domainmap, rangemap);
      }
    }
    if (B_domainmap_file == "" || B_rangemap_file == "" || B_rowmap_file == "" || B_colmap_file == "") {
      out << "Attempt to read sparse matrix file \"" << B_file
          << "\"" << endl;
      B = Reader<Matrix_t>::readSparseFile (B_file, comm);
    }
    else {
      RCP<const map_type> domainmap = Reader<Matrix_t>::readMapFile (B_domainmap_file, comm);
      RCP<const map_type> rangemap  = Reader<Matrix_t>::readMapFile (B_rangemap_file, comm);
      RCP<const map_type> rowmap    = Reader<Matrix_t>::readMapFile (B_rowmap_file, comm);
      RCP<const map_type> colmap    = Reader<Matrix_t>::readMapFile (B_colmap_file, comm);
      out << "Attempt to read sparse matrix file \"" << B_file
          << "\"" << endl;
      B = Reader<Matrix_t>::readSparseFile (B_file, rowmap, colmap, domainmap, rangemap);
    }
    if (C_domainmap_file == "" || C_rangemap_file == "" || C_rowmap_file == "" || C_colmap_file == "") {
      out << "Attempt to read sparse matrix file \"" << C_file
          << "\"" << endl;
      C = Reader<Matrix_t>::readSparseFile (C_file, comm);
    }
    else {
      RCP<const map_type> domainmap = Reader<Matrix_t>::readMapFile (C_domainmap_file, comm);
      RCP<const map_type> rangemap  = Reader<Matrix_t>::readMapFile (C_rangemap_file, comm);
      RCP<const map_type> rowmap    = Reader<Matrix_t>::readMapFile (C_rowmap_file, comm);
      RCP<const map_type> colmap    = Reader<Matrix_t>::readMapFile (C_colmap_file, comm);
      out << "Attempt to read sparse matrix file \"" << C_file
          << "\"" << endl;
      C = Reader<Matrix_t>::readSparseFile (C_file, rowmap, colmap, domainmap, rangemap);
    }
    if (D_file != "") {
      if (D_domainmap_file == "" || D_rangemap_file == "" || D_rowmap_file == "" || D_colmap_file == "") {
        out << "Attempt to read sparse matrix file \"" << D_file
            << "\"" << endl;
        D = Reader<Matrix_t>::readSparseFile (D_file, comm);
      }
      else {
        RCP<const map_type> domainmap = Reader<Matrix_t>::readMapFile (D_domainmap_file, comm);
        RCP<const map_type> rangemap  = Reader<Matrix_t>::readMapFile (D_rangemap_file, comm);
        RCP<const map_type> rowmap    = Reader<Matrix_t>::readMapFile (D_rowmap_file, comm);
        RCP<const map_type> colmap    = Reader<Matrix_t>::readMapFile (D_colmap_file, comm);
        out << "Attempt to read sparse matrix file \"" << D_file
            << "\"" << endl;
        D = Reader<Matrix_t>::readSparseFile (D_file, rowmap, colmap, domainmap, rangemap);
      }
    }

    if (A_file == "") {
      A = C;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(op != "multiply" && op != "add" && op != "RAP", std::runtime_error,
                               "Unrecognized Matrix Operation: " << op);

    if (op == "multiply") {
      if (verbose)
        newOut << "Running multiply test (manual FC) for " << currentSystem.name() << endl;

      mult_test_results results = multiply_test_manualfc(name, A, B, AT, BT, C, comm, newOut);

      if (verbose) {
        newOut << "Results:"     << endl;
        newOut << "\tEpsilon: "  << results.epsilon  << endl;
        newOut << "\tcNorm: "    << results.cNorm    << endl;
        newOut << "\tcompNorm: " << results.compNorm << endl;
        newOut << "\tisImportValid: " <<results.isImportValid << endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon);
      TEUCHOS_TEST_FOR_EXCEPTION(!results.isImportValid,std::logic_error,std::string("ManualFC: Import validity failed: ") + currentSystem.name());

      if (verbose)
        newOut << "Running multiply test (auto FC) for " << currentSystem.name() << endl;

      results = multiply_test_autofc(name, A, B, AT, BT, C, comm, newOut);

      if (verbose) {
        newOut << "Results:"     << endl;
        newOut << "\tEpsilon: "  << results.epsilon  << endl;
        newOut << "\tcNorm: "    << results.cNorm    << endl;
        newOut << "\tcompNorm: " << results.compNorm << endl;
        newOut << "\tisImportValid: " <<results.isImportValid << endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon);
      TEUCHOS_TEST_FOR_EXCEPTION(!results.isImportValid,std::logic_error,std::string("AutoFC: Import validity failed: ") + currentSystem.name());

      if (verbose)
        newOut << "Running multiply reuse test for " << currentSystem.name() << endl;

      results = multiply_reuse_test(name, A, B, AT, BT, C, comm, newOut);

      if (verbose) {
        newOut << "Results:"     << endl;
        newOut << "\tEpsilon: "  << results.epsilon  << endl;
        newOut << "\tcNorm: "    << results.cNorm    << endl;
        newOut << "\tcompNorm: " << results.compNorm << endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon);

      // Check if all diagonal entries are there, required for KokkosKernels Jacobi
      bool diagExists = true;
      auto rowMap = A->getRowMap();
      Tpetra::Vector<MT, LO, GO, NT> diags(rowMap);
      A->getLocalDiagCopy(diags);
      size_t diagLength = rowMap->getLocalNumElements();
      Teuchos::Array<MT> diagonal(diagLength);
      diags.get1dCopy(diagonal());

      for(size_t i = 0; i < diagLength; ++i) {
        if(diagonal[i] == Teuchos::ScalarTraits<SC>::zero()) {
          diagExists = false;
          break;
        }
      }

      // Do we try Jacobi?
      if (diagExists && AT == false && BT == false && A->getDomainMap()->isSameAs(*A->getRangeMap())) {
        if (verbose)
          newOut << "Running jacobi test for " << currentSystem.name() << endl;

        results = jacobi_test(name, A, B, comm, newOut);
        if (verbose) {
          newOut << "Results:"     << endl;
          newOut << "\tEpsilon: "  << results.epsilon  << endl;
          newOut << "\tcNorm: "    << results.cNorm    << endl;
          newOut << "\tcompNorm: " << results.compNorm << endl;
        }
        TEST_COMPARE(results.epsilon, <, epsilon)

        if (verbose)
          newOut << "Running jacobi reuse test for " << currentSystem.name() << endl;

        results = jacobi_reuse_test(name, A, B, comm, newOut);

        if (verbose) {
          newOut << "Results:"     << endl;
          newOut << "\tEpsilon: "  << results.epsilon  << endl;
          newOut << "\tcNorm: "    << results.cNorm    << endl;
          newOut << "\tcompNorm: " << results.compNorm << endl;
        }
        TEST_COMPARE(results.epsilon, <, epsilon)
      }
    }
    else if (op == "add") {
      if (verbose)
        newOut << "Running 3-argument add test (nonnull C on input) for "
               << currentSystem.name() << endl;

      add_test_results results = regular_add_test(name+"_add",A, B, AT, BT, C, comm);

      TEST_COMPARE(results.epsilon, <, epsilon);
      newOut << "Regular Add Test Results: " << endl;
      newOut << "\tCorrect Norm: " << results.correctNorm << endl;
      newOut << "\tComputed norm: " << results.computedNorm << endl;
      newOut << "\tEpsilon: " << results.epsilon << endl;

      if (verbose)
        newOut << "Running 3-argument add reuse test (nonnull C on input) for "
               << currentSystem.name() << endl;

      results = reuse_add_test(name+"_add",A, B, AT, BT, C, comm);

      TEST_COMPARE(results.epsilon, <, epsilon);
      newOut << "Reuse Add Test Results: " << endl;
      newOut << "\tCorrect Norm: " << results.correctNorm << endl;
      newOut << "\tComputed norm: " << results.computedNorm << endl;
      newOut << "\tEpsilon: " << results.epsilon << endl;

      // FIXME (mfh 09 May 2013) This test doesn't currently pass.  I
      // don't think anyone ever exercised the case where C is null on
      // input before.  I'm disabling this test for now until I have a
      // chance to fix that case.
      if (verbose)
        newOut << "Running 3-argument add test (null C on input) for "
               << currentSystem.name() << endl;

      TEUCHOS_TEST_FOR_EXCEPTION(A.is_null (), std::logic_error,
                                 "Before null_add_test_1: A is null");
      TEUCHOS_TEST_FOR_EXCEPTION(B.is_null (), std::logic_error,
                                 "Before null_add_test_1: B is null");
      TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
                                 "Before null_add_test_1: C is null");

      results = null_add_test_1<Matrix_t> (*A, *B, AT, BT, *C,
                                         newOut, success);

      TEST_COMPARE(results.epsilon, <, epsilon);
      newOut << "Null Add Test (1) Results: " << endl;
      newOut << "\tCorrect Norm: " << results.correctNorm << endl;
      newOut << "\tComputed norm: " << results.computedNorm << endl;
      newOut << "\tEpsilon: " << results.epsilon << endl;

      results = null_add_test_2<Matrix_t> (*A, *B, AT, BT, *C,
                                         newOut, success);

      TEST_COMPARE(results.epsilon, <, epsilon);
      newOut << "Null Add Test (2) Results: " << endl;
      newOut << "\tCorrect Norm: " << results.correctNorm << endl;
      newOut << "\tComputed norm: " << results.computedNorm << endl;
      newOut << "\tEpsilon: " << results.epsilon << endl;

      B = Reader<Matrix_t >::readSparseFile(B_file, comm, true, false, false);
      //declare E with enough entries to receive A + B
      size_t n = 0;
      for (size_t i = 0; i < B->getLocalNumRows(); i++) {
        if (n < B->getNumEntriesInLocalRow(i))
          n = B->getNumEntriesInLocalRow(i);
      }
      n += A->getLocalMaxNumRowEntries();

      RCP<const map_type> rm = B->getRowMap();
      RCP<Matrix_t> E = rcp (new Matrix_t(rm, n));
      auto one = Teuchos::ScalarTraits<MT>::one();
      Tpetra::MatrixMatrix::Add(*B, BT, one, *E, one);

      if (! BT) {
        if (verbose)
          newOut << "Running 2-argument add test for "
                 << currentSystem.name() << endl;

        results = add_into_test(A,E,AT,C,comm);
        TEST_COMPARE(results.epsilon, <, epsilon)
        newOut << "Add Into Test Results: " << endl;
        newOut << "\tCorrect Norm: " << results.correctNorm << endl;
        newOut << "\tComputed norm: " << results.computedNorm << endl;
        newOut << "\tEpsilon: " << results.epsilon << endl;
      }
      newOut << "Add tests complete" << endl;
    }
    else if (op == "RAP") {
      if (verbose)
        newOut << "Running multiply RAP test for " << currentSystem.name() << endl;

      mult_test_results results = multiply_RAP_test_autofc(name, A, B, C, AT, BT, CT, D, comm, newOut);

      if (verbose) {
        newOut << "Results:"     << endl;
        newOut << "\tEpsilon: "  << results.epsilon  << endl;
        newOut << "\tcNorm: "    << results.cNorm    << endl;
        newOut << "\tcompNorm: " << results.compNorm << endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon);

      if (verbose)
        newOut << "Running multiply RAP reuse test for " << currentSystem.name() << endl;

      results = multiply_RAP_reuse_test(name, A, B, C, AT, BT, CT, D, comm, newOut);

      if (verbose) {
        newOut << "Results:"     << endl;
        newOut << "\tEpsilon: "  << results.epsilon  << endl;
        newOut << "\tcNorm: "    << results.cNorm    << endl;
        newOut << "\tcompNorm: " << results.compNorm << endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon);
    }
  }

  const int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  newOut << "We made it through operations_test on all processes!" << endl;
  if (gblSuccess != 1) {
    newOut << "FAILED on at least one process!" << endl;
  }
  TEST_EQUALITY_CONST( gblSuccess, 1 );
}

/*
 * This test was created at the request of Chris Siefert to verify
 * that some inexplicable behaviour in MueLu was not due to a faulty
 * assumption in the Matrix Matrix Multiply Kernel.
 * KLN 15/06/2011
 */


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMat, range_row_test, SC, LO, GO, NT)  {
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  typedef Map<LO,GO,NT>         Map_t;
  typedef CrsMatrix<SC,LO,GO,NT> Matrix_t;

  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize();

  // Create an output stream that prints immediately, rather than
  // waiting until the test finishes.  This lets us get debug output
  // even if an unexpected exception or crash happens.  Unfortunately,
  // Teuchos::FancyOStream does not implement operator=, so we can't
  // replace 'out' temporarily.
  RCP<FancyOStream> newOutP = (myRank == 0) ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
    Teuchos::getFancyOStream (rcp (new Teuchos::oblackholestream ()));
  FancyOStream& newOut = *newOutP;

  newOut << "Tpetra sparse matrix-matrix multiply: range row test" << endl;
  Teuchos::OSTab tab1 (newOut);

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //THIS NUMBER MUST BE EVEN SO THAT WHEN I CALCULATE THE NUMBER
  //OF ROWS IN THE DOMAIN MAP I DON'T ENCOUNTER ANY
  //WEIRD RESULTS DUE TO INTEGER DIVISION
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  int numRowsPerProc = 4;
  int rank = comm->getRank();
  global_size_t globalNumRows = numRowsPerProc*numProcs;

  RCP<Matrix_t > dummy =
    getIdentityMatrix<SC,LO,GO,NT> (newOut, globalNumRows, comm);

//Create "B"
  Array<GO> myRows = tuple<GO>(
    rank*numRowsPerProc,
    rank*numRowsPerProc+1,
    rank*numRowsPerProc+2,
    rank*numRowsPerProc+3);
  Array<GO> rangeElements;
  if(rank == 0){
    rangeElements = tuple<GO>(
      (numProcs-1)*numRowsPerProc+1,
      (numProcs-1)*numRowsPerProc+2,
      (numProcs-1)*numRowsPerProc,
      (numProcs-1)*numRowsPerProc+3);
  }
  else{
    rangeElements = tuple<GO>(
      (rank-1)*numRowsPerProc+1,
      (rank-1)*numRowsPerProc+2,
      (rank-1)*numRowsPerProc,
      (rank-1)*numRowsPerProc+3);
  }

  newOut << "Create row, range, and domain Maps of B" << endl;
  RCP<const Map_t > bRowMap =
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(myRows, comm);
  RCP<const Map_t > bRangeMap =
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(rangeElements, comm);
  //We divide by 2 to make the matrix tall and "skinny"
  RCP<const Map_t > bDomainMap =
    Tpetra::createUniformContigMapWithNode<LO,GO,NT>(globalNumRows/2, comm);

  newOut << "Create identityMatrix" << endl;
  RCP<Matrix_t > identityMatrix =
    getIdentityMatrixWithMap<SC,LO,GO,NT> (newOut, bRowMap, comm);


  newOut << "Create bMatrix" << endl;
  RCP<Matrix_t > bMatrix =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(bRowMap, 1);

  newOut << "Fill bMatrix" << endl;
  {
    Teuchos::ArrayView<const GO> gblRows = bRowMap->getLocalElementList ();
    for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
      Array<GO> col(1,(*it)/2);
      Array<SC> val(1,3.0);
      bMatrix->insertGlobalValues(*it, col(), val());
    }
  }


  newOut << "Call fillComplete on bMatrix" << endl;
  bMatrix->fillComplete(bDomainMap, bRangeMap);

  newOut << "Regular I*P" << endl;
  mult_test_results results = multiply_test_manualfc(
    "Different Range and Row Maps",
    identityMatrix,
    bMatrix,
    false,
    false,
    bMatrix,
    comm,
    out);

  if(verbose){
    newOut << "Results:" << endl;
    newOut << "\tEpsilon: " << results.epsilon << endl;
    newOut << "\tcNorm: " << results.cNorm << endl;
    newOut << "\tcompNorm: " << results.compNorm << endl;
  }
  const double defaultEpsilon = 100. * Teuchos::ScalarTraits<SC>::eps();
  TEST_COMPARE(results.epsilon, <, defaultEpsilon);

  newOut << "Create identity2" << endl;
  RCP<Matrix_t > identity2 =
    getIdentityMatrix<SC,LO,GO,NT> (newOut, globalNumRows/2, comm);

  RCP<const Map_t > bTransRowMap =
    Tpetra::createUniformContigMapWithNode<LO,GO,NT>(globalNumRows/2,comm);

  Array<GO> bTransRangeElements;
  if(rank == 0){
    bTransRangeElements = tuple<GO>(
      (numProcs-1)*(numRowsPerProc/2)+1,
      (numProcs-1)*(numRowsPerProc/2));
  }
  else{
    bTransRangeElements = tuple<GO>(
      (rank-1)*(numRowsPerProc/2)+1,
      (rank-1)*(numRowsPerProc/2));
  }
  newOut << bTransRangeElements << endl;

  RCP<const Map_t > bTransRangeMap =
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(bTransRangeElements, comm);
  RCP<const Map_t > bTransDomainMap =
    Tpetra::createUniformContigMapWithNode<LO,GO,NT>(globalNumRows, comm);

  newOut << "Create and fill bTransTest" << endl;
  RCP<Matrix_t > bTransTest =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(bTransRowMap, 2);

  {
    Teuchos::ArrayView<const GO> gblRows = bRowMap->getLocalElementList ();
    for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
      Teuchos::Array<GO> col(1,*it);
      Teuchos::Array<SC> val(1,3.0);
      bTransTest->insertGlobalValues((*it)/2, col(), val());
    }
  }

  newOut << "Call fillComplete on bTransTest" << endl;
  bTransTest->fillComplete(bTransDomainMap, bTransRangeMap);

  newOut << "Create and fill bTrans" << endl;
  RCP<Matrix_t> bTrans =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(bTransRowMap,
                                         bTransTest->getGlobalMaxNumRowEntries());

  newOut << "Compute identity * transpose(bTrans)" << endl;
  Tpetra::MatrixMatrix::Multiply(*identity2,false,*bMatrix, true, *bTrans, false);

  newOut << "Call fillComplete on bTrans" << endl;
  bTrans->fillComplete(bTransDomainMap, bTransRangeMap);

  // FIXME (mfh 03 May 2016) I didn't write this output message.  I
  // don't know what it means, but I'm leaving it here in case it's
  // meaningful to someone.
  newOut << "Regular I*P^T" << endl;

  RCP<Matrix_t > bTransDiff =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(bTransRowMap, bTransTest->getGlobalMaxNumRowEntries());
  Tpetra::MatrixMatrix::Add<SC,LO,GO,NT>(*bTransTest, false, -1.0, *bTrans, false, 1.0,bTransDiff);

  newOut << "Call fillComplete on bTransDiff" << endl;
  bTransDiff->fillComplete(bTransDomainMap, bDomainMap);

  double diffNorm = bTransDiff->getFrobeniusNorm ();
  double realNorm = bTransTest->getFrobeniusNorm ();
  double calcEpsilon = diffNorm/realNorm;

  newOut << "B" << endl;

  if(verbose){
    out << "Results:" << endl;
    out << "\tEpsilon: " << calcEpsilon << endl;
    out << "\treal norm: " << realNorm << endl;
    out << "\tcompNorm: " << diffNorm << endl;
  }
  TEST_COMPARE(calcEpsilon, <, defaultEpsilon);

  const int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  newOut << "We made it through range_row_test on all processes!" << endl;
  if (gblSuccess != 1) {
    newOut << "FAILED on at least one process!" << endl;
  }
  TEST_EQUALITY_CONST( gblSuccess, 1 );
}

// Build upper diag matrix for test 
template<class BlockCrsMatrixType>
void build_A_matrix (const RCP<BlockCrsMatrixType>& A) {

  using LO = typename BlockCrsMatrixType::local_ordinal_type;
  using GO = typename BlockCrsMatrixType::global_ordinal_type;
  using Scalar = typename BlockCrsMatrixType::scalar_type;

  const typename BlockCrsMatrixType::map_type row_map = *(A->getRowMap());
  const typename BlockCrsMatrixType::map_type col_map = *(A->getColMap());

  int my_rank = row_map.getComm()->getRank();

  if(A->getBlockSize() != 2) {
    if (my_rank==0) std::cerr << "Error: A->getBlockSize != 2!" << std::endl;
    return;
  }

  const int blocksize = 2;

  for (LO localrow = row_map.getMinLocalIndex();
       localrow <= row_map.getMaxLocalIndex();
       ++localrow) {

    const GO globalrow = row_map.getGlobalElement(localrow);

    if (globalrow == 0) {

      LO local_col_indices[1];
      local_col_indices[0] = col_map.getLocalElement(globalrow);

      Scalar values[blocksize*blocksize];
      for (size_t b=0; b<blocksize*blocksize; ++b) {
        values[b] = 2*globalrow+1;
      }
      A->replaceLocalValues(localrow,
                            local_col_indices,
                            values,
                            1);

    } else if (globalrow == 1 || globalrow == 2) {

      LO local_col_indices[2];
      local_col_indices[0] = col_map.getLocalElement(globalrow-1);
      local_col_indices[1] = col_map.getLocalElement(globalrow);

      Scalar values[2*blocksize*blocksize];
      int local_indx = 0;
      for (GO globalcol=globalrow-1; globalcol<=globalrow; ++globalcol) {
        int start = local_indx*blocksize*blocksize;
        for (size_t b=0; b<blocksize*blocksize; ++b) {
          values[start+b] = (globalrow+1)+globalcol;
        }
        ++local_indx;
      }
      A->replaceLocalValues(localrow,
                            local_col_indices,
                            values,
                            2);

    } else { // globalrow == 3

      LO local_col_indices[1];
      local_col_indices[0] = col_map.getLocalElement(globalrow-1);

      Scalar values[blocksize*blocksize];
      for (size_t b=0; b<blocksize*blocksize; ++b) {
        values[b] = 2*globalrow;
      }
      A->replaceLocalValues(localrow,
                            local_col_indices,
                            values,
                            1);

    }
  }
}

// Build lower diag matrix for test   
template<class BlockCrsMatrixType>
void build_B_matrix (const RCP<BlockCrsMatrixType>& B) {

  using LO = typename BlockCrsMatrixType::local_ordinal_type;
  using GO = typename BlockCrsMatrixType::global_ordinal_type;
  using Scalar = typename BlockCrsMatrixType::scalar_type;

  const typename BlockCrsMatrixType::map_type row_map = *(B->getRowMap());
  const typename BlockCrsMatrixType::map_type col_map = *(B->getColMap());

  int my_rank = row_map.getComm()->getRank();

  if(B->getBlockSize() != 2) {
    if (my_rank==0) std::cerr << "Error: B->getBlockSize != 2!" << std::endl;
    return;
  }

  const int blocksize = 2;

  for (LO localrow = row_map.getMinLocalIndex();
       localrow <= row_map.getMaxLocalIndex();
       ++localrow) {

    const GO globalrow = row_map.getGlobalElement(localrow);

    LO local_col_indices[2];
    local_col_indices[0] = col_map.getLocalElement(globalrow);
    local_col_indices[1] = col_map.getLocalElement(globalrow+1);

    Scalar values[2*blocksize*blocksize];
    int local_indx = 0;
    for (GO globalcol=globalrow; globalcol<=globalrow+1; ++globalcol) {
      int start = local_indx*blocksize*blocksize;
      for (size_t b=0; b<blocksize*blocksize; ++b) {
        values[start+b] = (globalrow+1)+globalcol;
      }
      ++local_indx;
    }
    B->replaceLocalValues(localrow,
                          local_col_indices,
                          values,
                          2);
  }
}
// Build tri diag matrix for test 
template<class BlockCrsMatrixType>
void build_C_matrix (const RCP<BlockCrsMatrixType>& C) {

  using LO = typename BlockCrsMatrixType::local_ordinal_type;
  using GO = typename BlockCrsMatrixType::global_ordinal_type;
  using Scalar = typename BlockCrsMatrixType::scalar_type;

  const typename BlockCrsMatrixType::map_type row_map = *(C->getRowMap());
  const typename BlockCrsMatrixType::map_type col_map = *(C->getColMap());

  int my_rank = row_map.getComm()->getRank();

  if(C->getBlockSize() != 2) {
    if (my_rank==0) std::cerr << "Error: C->getBlockSize != 2!" << std::endl;
    return;
  }

  const int blocksize = 2;

  for (LO localrow = row_map.getMinLocalIndex();
       localrow <= row_map.getMaxLocalIndex();
       ++localrow) {

    const GO globalrow = row_map.getGlobalElement(localrow);

    if (globalrow == 0) {

      LO local_col_indices[2];
      local_col_indices[0] = col_map.getLocalElement(globalrow);
      local_col_indices[1] = col_map.getLocalElement(globalrow+1);

      Scalar values[2*blocksize*blocksize];
      int local_indx = 0;
      for (GO globalcol=globalrow; globalcol<=globalrow+1; ++globalcol) {
        int start = local_indx*blocksize*blocksize;
        for (size_t b=0; b<blocksize*blocksize; ++b) {
          values[start+b] = 2*(globalrow+1)*(globalcol+1);
        }
        ++local_indx;
      }
      C->replaceLocalValues(localrow,
                            local_col_indices,
                            values,
                            2);

    } else if (globalrow == 1 || globalrow == 2) {

      LO local_col_indices[3];
      local_col_indices[0] = col_map.getLocalElement(globalrow-1);
      local_col_indices[1] = col_map.getLocalElement(globalrow);
      local_col_indices[2] = col_map.getLocalElement(globalrow+1);

      Scalar values[3*blocksize*blocksize];
      int local_indx = 0;
      for (GO globalcol=globalrow-1; globalcol<=globalrow+1; ++globalcol) {
        Scalar vals_for_col[3];
        if (globalrow==1) {
          vals_for_col[0] = 4; vals_for_col[1] = 26; vals_for_col[2] = 24;
        } else if (globalrow==2) {
          vals_for_col[0] = 24; vals_for_col[1] = 82; vals_for_col[2] = 60;
        }
        int start = local_indx*blocksize*blocksize;
        for (size_t b=0; b<blocksize*blocksize; ++b) {
          values[start+b] = vals_for_col[local_indx];
        }
        ++local_indx;
      }
      C->replaceLocalValues(localrow,
                            local_col_indices,
                            values,
                            3);

    } else { // globalrow == 3

      LO local_col_indices[2];
      local_col_indices[0] = col_map.getLocalElement(globalrow-1);
      local_col_indices[1] = col_map.getLocalElement(globalrow);

      Scalar values[2*blocksize*blocksize];
      int local_indx = 0;
      for (GO globalcol=globalrow-1; globalcol<=globalrow; ++globalcol) {
        int start = local_indx*blocksize*blocksize;
        for (size_t b=0; b<blocksize*blocksize; ++b) {
          values[start+b] = 2*(globalrow+3)*(globalcol+3);
        }
        ++local_indx;
      }
      C->replaceLocalValues(localrow,
                            local_col_indices,
                            values,
                            2);

    }
  }
}

// Test that two matrices' rows have the same entries.
template<class BlockCrsMatrixType>
bool matrices_are_same(const RCP<BlockCrsMatrixType>& A1,
                       const RCP<BlockCrsMatrixType>& A2)
{
  // Loop through A1 and make sure each row has the same
  // entries as A2.  In the fully general case, the
  // redistribution may have added together values, resulting in
  // small rounding errors.  This is why we use an error tolerance
  // (with a little bit of wiggle room).

  int my_rank = A1->getRowMap()->getComm()->getRank();

  using LO = typename BlockCrsMatrixType::local_ordinal_type;
  using Scalar = typename BlockCrsMatrixType::scalar_type;
  using lids_type = typename BlockCrsMatrixType::local_inds_host_view_type;
  using vals_type = typename BlockCrsMatrixType::values_host_view_type;

  using ST = Teuchos::ScalarTraits<Scalar>;
  using magnitude_type = typename ST::magnitudeType;
  const magnitude_type tol =
    Teuchos::as<magnitude_type> (10) * Teuchos::ScalarTraits<magnitude_type>::eps ();

  const LO blocksize = A1->getBlockSize();
  // Verify the blocksizes are identical
  if (blocksize != A2->getBlockSize()) {
    if (my_rank==0) std::cout << "Error: Blocksizes are not the same!" << std::endl;
    return false;
  }

  lids_type A1LocalColInds;
  vals_type A1LocalRowVals;
  lids_type A2LocalColInds;
  vals_type A2LocalRowVals;
  for (LO localrow = A1->getRowMap()->getMinLocalIndex();
      localrow <= A1->getRowMap()->getMaxLocalIndex();
      ++localrow)
  {
    size_t A1NumEntries = A1->getNumEntriesInLocalRow (localrow);
    size_t A2NumEntries = A1->getNumEntriesInLocalRow (localrow);

    // Verify the same number of entries in each row
    if (A1NumEntries != A2NumEntries) {
      if (my_rank==0) std::cout << "Error: Matrices have different number of entries in at least one row!" << std::endl;
      return false;
    }

    A1->getLocalRowView (localrow, A1LocalColInds, A1LocalRowVals);
    A2->getLocalRowView (localrow, A2LocalColInds, A2LocalRowVals);

    // Verify the same number of values in each row
    if (A1LocalRowVals.extent(0) != A2LocalRowVals.extent(0)) {
      if (my_rank==0) std::cout << "Error: Matrices have different number of entries in at least one row!" << std::endl;
      return false;
    }

    // Col indices may be in different order. Store in a set and compare sets.
    std::set<LO> a1_inds;
    std::set<LO> a2_inds;
   typedef typename Array<Scalar>::size_type size_type;
    for (size_type k = 0; k < static_cast<size_type> (A1NumEntries); ++k) {
      a1_inds.insert(A1LocalColInds[k]);
      a2_inds.insert(A2LocalColInds[k]);
    }
    if(a1_inds!=a2_inds) {
      if (my_rank==0) std::cout << "Error: Matrices have different column indices!" << std::endl;
      return false;
    }

    // Loop over each local col entry of A1, find the corresponding col index of A2, and compare these value.
    for (size_type a1_k = 0; a1_k < static_cast<size_type> (A1NumEntries); ++a1_k) {
      LO a2_k=0;
      for (size_type i = 0; i < static_cast<size_type> (A2NumEntries); ++i) {
        if (A1LocalColInds[a1_k] == A2LocalColInds[i]) {
          a2_k = i;
          break;
        }
      }
      const int a1_start = a1_k*blocksize*blocksize;
      const int a2_start = a2_k*blocksize*blocksize;
      for (int b=0; b<blocksize*blocksize; ++b) {
        const magnitude_type rel_err = ST::magnitude(A1LocalRowVals[a1_start+b] - A2LocalRowVals[a2_start+b]);
        if(rel_err > tol) {
          if (my_rank==0) std::cout << "Error: Matrices have different values!" << std::endl;
          return false;
        }
      }
    }
  }

  return true;
}

/*
 * Test the matrix-matrix multiplication for a simple BlockCrsMatrix example.
 * The test is only run on 2 processors. The test is A*B=C, where A is 4x3 lower
 * diagonal matrix and  B is 3x4 upper diagonal matrix, each with up to 1 
 * off-diagonal entry. The result C is a tridiagonal matrix. Blocksize=2.
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMatMult, BlockCrsMult, SC,LO, GO, NT)  {
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

  // NOTE: The matrix reader doesn't read real matrices into a complex data type, so we just swap down to MT here
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
  typedef BlockCrsMatrix<MT,LO,GO,NT> BlockMatrix_t;
  typedef Tpetra::CrsGraph<LO, GO, NT> crs_graph_type;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize();

  // Create an output stream that prints immediately, rather than
  // waiting until the test finishes.  This lets us get debug output
  // even if an unexpected exception or crash happens.  Unfortunately,
  // Teuchos::FancyOStream does not implement operator=, so we can't
  // replace 'out' temporarily.
  RCP<FancyOStream> newOutP = (myRank == 0) ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
    Teuchos::getFancyOStream (rcp (new Teuchos::oblackholestream ()));
  FancyOStream& newOut = *newOutP;

  newOut << "Tpetra block sparse matrix-matrix multiply: operations_test" << endl;
  Teuchos::OSTab tab1 (newOut);

  // This test is 2 rank specific                                                             
  if (numProcs != 2) {
    return;
  }

  // Build A, B, C
  const GO indexBase = 0;
  LO num_local_elements_A = 2;
  LO num_local_elements_B = (myRank==0 ? 1 : 2);

  const int blocksize = 2;
  using GST = Tpetra::global_size_t;
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

  // Create row maps
  RCP<const map_type> row_map_A =
    rcp (new map_type (INVALID,
                       num_local_elements_A,
                       indexBase, comm));
  RCP<const map_type> row_map_B =
    rcp (new map_type (INVALID,
                       num_local_elements_B,
                       indexBase, comm));

  // Build graphs                                                        
  Teuchos::RCP<crs_graph_type> graph_A =
    Teuchos::rcp (new crs_graph_type (row_map_A, 2));
  Teuchos::RCP<crs_graph_type> graph_B =
    Teuchos::rcp (new crs_graph_type (row_map_B, 2));
  Teuchos::RCP<crs_graph_type> graph_C =
    Teuchos::rcp (new crs_graph_type (row_map_A, 3));
  {
    Array<GO> cols_A(2);
    for (LO localrow = row_map_A->getMinLocalIndex ();
         localrow <= row_map_A->getMaxLocalIndex (); ++localrow) {
      const GO globalrow = row_map_A->getGlobalElement(localrow);
      if (globalrow==0) {
        cols_A.resize(1);
        cols_A[0] = globalrow;
      } else if (globalrow==3) {
        cols_A.resize(1);
        cols_A[0] = globalrow-1;
      } else {
        cols_A.resize(2);
        cols_A[0] = globalrow-1; cols_A[1] = globalrow;
      }
      graph_A->insertGlobalIndices (globalrow, cols_A());
    }
    graph_A->fillComplete(row_map_B, row_map_A);
  }
  {
    Array<GO> cols_B(2);
    for (LO localrow = row_map_B->getMinLocalIndex ();
         localrow <= row_map_B->getMaxLocalIndex (); ++localrow) {
      const GO globalrow = row_map_B->getGlobalElement(localrow);
      cols_B.resize(2);
      cols_B[0] = globalrow; cols_B[1] = globalrow+1;
      graph_B->insertGlobalIndices (globalrow, cols_B());
    }
    graph_B->fillComplete(row_map_A, row_map_B);
  }
  {
    Array<GO> cols_C(2);
    for (LO localrow = row_map_A->getMinLocalIndex ();
         localrow <= row_map_A->getMaxLocalIndex (); ++localrow) {
      const GO globalrow = row_map_A->getGlobalElement(localrow);
      if (globalrow==0) {
        cols_C.resize(2);
        cols_C[0] = globalrow; cols_C[1] = globalrow+1;
      } else if (globalrow==3) {
        cols_C.resize(2);
        cols_C[0] = globalrow-1; cols_C[1] = globalrow;
      } else {
        cols_C.resize(3);
        cols_C[0] = globalrow-1; cols_C[1] = globalrow; cols_C[2] = globalrow+1; 
      }
      graph_C->insertGlobalIndices (globalrow, cols_C());
    }
    graph_C->fillComplete();
  }

  // Build matrices. C is used to test against the computed matrix
  RCP<BlockMatrix_t> A =
    rcp (new BlockMatrix_t (*graph_A, blocksize)); 
  RCP<BlockMatrix_t> B =
    rcp (new BlockMatrix_t (*graph_B, blocksize));
  RCP<BlockMatrix_t> C =
    rcp (new BlockMatrix_t (*graph_C, blocksize));

  build_A_matrix(A);  
  build_B_matrix(B);
  build_C_matrix(C);

  RCP<BlockMatrix_t> computedC;
  Tpetra::MatrixMatrix::Multiply(A.getConst(), false,
                                 B.getConst(), false, 
                                 computedC);

  // Test that C and the computed matrix are identical
  bool matrices_match = matrices_are_same<BlockMatrix_t>(computedC, C);
  TEST_ASSERT(matrices_match);

  const int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  newOut << "We made it through operations_test on all processes!" << endl;
  if (gblSuccess != 1) {
    newOut << "FAILED on at least one process!" << endl;
  }
  TEST_EQUALITY_CONST( gblSuccess, 1 );
}

/**
 * This test was written at the request of Chris Siefert
 * in order to verity that A^T * I produces correct results
 * when A's rowmap and rangemap are differnt.
 * KLN 23/06/2011
 */
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMat, ATI_range_row_test, SC, LO, GO, NT)  {
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  typedef Map<LO,GO,NT>          Map_t;
  typedef CrsMatrix<SC,LO,GO,NT> Matrix_t;

  // Create an output stream that prints immediately, rather than
  // waiting until the test finishes.  This lets us get debug output
  // even if an unexpected exception or crash happens.  Unfortunately,
  // Teuchos::FancyOStream does not implement operator=, so we can't
  // replace 'out' temporarily.
  const int myRank = comm->getRank ();
  RCP<FancyOStream> newOutP = (myRank == 0) ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
    Teuchos::getFancyOStream (rcp (new Teuchos::oblackholestream ()));
  FancyOStream& newOut = *newOutP;

  newOut << "Tpetra sparse matrix-matrix multiply: Test A^T * I, "
    "where A's row Map and range Map differ" << endl;
  Teuchos::OSTab tab1 (newOut);

  const int numProcs = comm->getSize();
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //THIS NUMBER MUST BE EVEN SO THAT WHEN I CALCULATE THE NUMBER
  //OF ROWS IN THE DOMAIN MAP I DON'T ENCOUNTER ANY
  //WEIRD RESULTS DUE TO INTEGER DIVISION
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  int numRowsPerProc = 4;
  const int rank = comm->getRank();
  global_size_t globalNumRows = numRowsPerProc*numProcs;

  newOut << "Create identity matrix" << endl;
  RCP<Matrix_t> identityMatrix =
    getIdentityMatrix<SC,LO,GO,NT> (newOut, globalNumRows, comm);

  newOut << "Create Maps for matrix aMat" << endl;
  Array<GO> aMyRows = tuple<GO>(
    rank*numRowsPerProc,
    rank*numRowsPerProc+1,
    rank*numRowsPerProc+2,
    rank*numRowsPerProc+3);
  RCP<const Map_t> aRowMap =
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(aMyRows, comm);
  RCP<const Map_t> aDomainMap =
    Tpetra::createUniformContigMapWithNode<LO,GO,NT>(globalNumRows/2, comm);
  Array<GO> aRangeElements;
  if(rank == 0){
    aRangeElements = tuple<GO>(
      (numProcs-1)*numRowsPerProc+1,
      (numProcs-1)*numRowsPerProc+2,
      (numProcs-1)*numRowsPerProc,
      (numProcs-1)*numRowsPerProc+3);
  }
  else{
    aRangeElements = tuple<GO>(
      (rank-1)*numRowsPerProc+1,
      (rank-1)*numRowsPerProc+2,
      (rank-1)*numRowsPerProc,
      (rank-1)*numRowsPerProc+3);
  }
  RCP<const Map<LO,GO,NT> > aRangeMap =
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(aRangeElements, comm);

  newOut << "Create matrix aMat" << endl;
  RCP<Matrix_t> aMat =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(aRowMap, 1);

  newOut << "Fill matrix aMat" << endl;
  {
    Teuchos::ArrayView<const GO> gblRows = aRowMap->getLocalElementList ();
    for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
      Teuchos::Array<GO> col(1,(*it)/2);
      Teuchos::Array<SC> val(1,3.0);
      aMat->insertGlobalValues(*it, col(), val());
    }
  }
  newOut << "Call fillComplete on matrix aMat" << endl;
  aMat->fillComplete(aDomainMap, aRangeMap);

  newOut << "Create RowMatrixTransposer with aMat" << endl;
  RowMatrixTransposer<SC,LO,GO,NT> transposer (aMat);

  newOut << "Use RowMatrixTransposer to create transpose of aMat" << endl;
  RCP<Matrix_t > knownAMat = transposer.createTranspose();

  // FIXME (mfh 03 May 2016) I'm not sure what this message means, so
  // I'll leave it.
  newOut << "Regular I*P" << endl;
  mult_test_results results = multiply_test_manualfc(
    "Different Range and Row Maps",
    aMat,
    identityMatrix,
    true,
    false,
    knownAMat,
    comm,
    newOut);
  if(verbose){
    newOut << "Results:" <<endl;
    newOut << "\tEpsilon: " << results.epsilon << endl;
    newOut << "\tcNorm: " << results.cNorm << endl;
    newOut << "\tcompNorm: " << results.compNorm << endl;
  }
  const double defaultEpsilon = 100. * Teuchos::ScalarTraits<SC>::eps();
  TEST_COMPARE(results.epsilon, <, defaultEpsilon);

  const int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  newOut << "We made it through ATI_range_row_test on all processes!" << endl;
  if (gblSuccess != 1) {
    newOut << "FAILED on at least one process!" << endl;
  }
  TEST_EQUALITY_CONST( gblSuccess, 1 );
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMat, threaded_add_sorted, SC, LO, GO, NT)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  using Teuchos::RCP;
  //First, make two local, random, sorted Kokkos sparse matrices
  size_t nrows = 1000;
  size_t nnzPerRow = 20;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using KCRS = typename crs_matrix_type::local_matrix_device_type;
  using ISC = typename crs_matrix_type::impl_scalar_type;
  using ValuesType = typename KCRS::values_type::non_const_type;
  using RowptrsType = typename KCRS::row_map_type::non_const_type;
  using ColindsType = typename KCRS::index_type::non_const_type;
  using Device = typename KCRS::device_type;
  using ExecSpace = typename Device::execution_space;
  //The 3 views are for A, B and C
  ValuesType valsCRS[3];
  RowptrsType rowptrsCRS[3];
  ColindsType colindsCRS[3];
  //Populate A and B
  typename ValuesType::HostMirror vals[3];
     vals[0] = typename ValuesType::HostMirror("vals0", nrows*nnzPerRow);
     vals[1] = typename ValuesType::HostMirror("vals1", nrows*nnzPerRow);
  typename RowptrsType::HostMirror rowptrs[3];
     rowptrs[0] = typename RowptrsType::HostMirror("rowptr0", nrows+1);
     rowptrs[1] = typename RowptrsType::HostMirror("rowptr1", nrows+1);
  typename ColindsType::HostMirror colinds[3]; 
     colinds[0] = typename ColindsType::HostMirror("colind0", nrows*nnzPerRow);
     colinds[1] = typename ColindsType::HostMirror("colind1", nrows*nnzPerRow);
  {
    //want consistent test results
    srand(12);
    for(LO m = 0; m < 2; m++)
    {
      for(size_t row = 0; row < nrows; row++)
      {
        //rows are sorted, so come up with nnzPerRow vals and cols and then sort them
        std::vector<ISC> rowVals(nnzPerRow);
        std::vector<LO> rowInds(nnzPerRow);
        for(size_t entry = 0; entry < nnzPerRow; entry++)
        {
          rowVals[entry] = ((double) rand()) / RAND_MAX;
          //don't allow repeats in col inds
          LO ind;
          do
          {
            ind = rand() % nrows;
          }
          while(std::find(rowInds.begin(), rowInds.end(), ind) != rowInds.end());
          rowInds[entry] = ind;
        }
        std::sort(rowInds.begin(), rowInds.end());
        //now insert new coordinates in big arrays
        for(size_t entry = 0; entry < nnzPerRow; entry++)
        {
          vals[m][nnzPerRow * row + entry] = rowVals[entry];
          colinds[m][nnzPerRow * row + entry] = rowInds[entry];
        }
      }
      //set rowptrs
      for(size_t row = 0; row <= nrows; row++)
      {
        rowptrs[m][row] = row * nnzPerRow;
      }
      valsCRS[m] = Kokkos::create_mirror_view_and_copy(
                                      typename NT::memory_space(), vals[m]);
      rowptrsCRS[m] = Kokkos::create_mirror_view_and_copy(
                                      typename NT::memory_space(), rowptrs[m]);
      colindsCRS[m] = Kokkos::create_mirror_view_and_copy(
                                      typename NT::memory_space(), colinds[m]);
    }
  }
  //now run the threaded addition on mats[0] and mats[1]
  ISC zero(0);
  ISC one(1);
  Tpetra::MMdetails::AddKernels<SC, LO, GO, NT>::addSorted(
                                valsCRS[0], rowptrsCRS[0], colindsCRS[0], one, 
                                valsCRS[1], rowptrsCRS[1], colindsCRS[1], one, 
#if KOKKOSKERNELS_VERSION >= 40299
                                nrows, // assumes square matrices
#endif
                                valsCRS[2], rowptrsCRS[2], colindsCRS[2]);

  ExecSpace().fence();

  vals[2] = Kokkos::create_mirror_view_and_copy(
                           typename ValuesType::HostMirror::memory_space(),
                           valsCRS[2]);
  rowptrs[2] = Kokkos::create_mirror_view_and_copy(
                              typename RowptrsType::HostMirror::memory_space(),
                              rowptrsCRS[2]);
  colinds[2] = Kokkos::create_mirror_view_and_copy(
                              typename ColindsType::HostMirror::memory_space(),
                              colindsCRS[2]);
  
  //the above function is an unfenced kernel launch, and the verification below relies on UVM, so fence here.
  //now scan through C's rows and entries to check they are correct
  TEST_ASSERT(rowptrs[0].extent(0) == rowptrs[2].extent(0));
  for(size_t i = 0; i < nrows; i++)
  {
    //also compute what C's row should be (as dense values)
    std::vector<ISC> correctVals(nrows, zero);
    std::vector<bool> correctEntries(nrows, false);
    for(size_t j = 0; j < nnzPerRow; j++)
    {
      int col1 = colinds[0](i * nnzPerRow + j);
      int col2 = colinds[1](i * nnzPerRow + j);
      correctVals[col1] += vals[0](i * nnzPerRow + j);
      correctEntries[col1] = true;
      correctVals[col2] += vals[1](i * nnzPerRow + j);
      correctEntries[col2] = true;
    }
    size_t actualNNZ = 0;
    for(size_t j = 0; j < nrows; j++)
    {
      if(correctEntries[j])
        actualNNZ++;
    }
    size_t Crowstart = rowptrs[2](i);
    size_t Crowlen = rowptrs[2](i + 1) - Crowstart;
    TEST_ASSERT(Crowlen == actualNNZ);
    for(size_t j = 0; j < Crowlen; j++)
    {
      size_t Coffset = Crowstart + j;
      if(j > 0)
      {
        //Check entries are sorted
        TEST_ASSERT(colinds[2](Coffset - 1) <= colinds[2](Coffset));
      }
      TEST_FLOATING_EQUALITY(vals[2](Coffset), correctVals[colinds[2](Coffset)], 1e-12);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMat, add_zero_rows, SC, LO, GO, NT)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  using Teuchos::RCP;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using MT = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  size_t nrows = 0;
  RCP<const map_type> emptyMap = rcp(new map_type(nrows, 0, comm));
  RCP<crs_matrix_type> A = rcp(new crs_matrix_type(emptyMap, 0));
  A->fillComplete(emptyMap, emptyMap);
  RCP<crs_matrix_type> B = rcp(new crs_matrix_type(emptyMap, 0));
  B->fillComplete(emptyMap, emptyMap);
  RCP<crs_matrix_type> C1 = rcp(new crs_matrix_type(emptyMap, 0));
  SC one = Teuchos::ScalarTraits<SC>::one();
  Tpetra::MatrixMatrix::add(one, false, *A, one, false, *B, *C1);
  RCP<crs_matrix_type> C2 = Tpetra::MatrixMatrix::add
    (one, false, *A, one, false, *B);
  MT magZero = Teuchos::ScalarTraits<MT>::zero();
  TEST_EQUALITY(C1->getFrobeniusNorm(), magZero);
  TEST_EQUALITY(C2->getFrobeniusNorm(), magZero);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMat, threaded_add_unsorted, SC, LO, GO, NT)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  using Teuchos::RCP;
  //First, make two local, random, unsorted Kokkos sparse matrices
  size_t nrows = 1000;
  size_t nnzPerRow = 20;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> crs_matrix_type;
  typedef typename crs_matrix_type::local_matrix_device_type KCRS;
  typedef typename crs_matrix_type::impl_scalar_type ISC;
  typedef typename KCRS::values_type::non_const_type ValuesType;
  typedef typename KCRS::row_map_type::non_const_type RowptrsType;
  typedef typename KCRS::index_type::non_const_type ColindsType;
  //The 3 views are for A, B and C
  ValuesType valsCRS[3];
  RowptrsType rowptrsCRS[3];
  ColindsType colindsCRS[3];
  //Populate A and B
  typename ValuesType::HostMirror vals[3];
     vals[0] = typename ValuesType::HostMirror("vals0", nrows*nnzPerRow);
     vals[1] = typename ValuesType::HostMirror("vals1", nrows*nnzPerRow);
     vals[2] = typename ValuesType::HostMirror("vals2", nrows*nnzPerRow);
  typename RowptrsType::HostMirror rowptrs[3];
     rowptrs[0] = typename RowptrsType::HostMirror("rowptr0", nrows+1);
     rowptrs[1] = typename RowptrsType::HostMirror("rowptr1", nrows+1);
     rowptrs[2] = typename RowptrsType::HostMirror("rowptr2", nrows+1);
  typename ColindsType::HostMirror colinds[3]; 
     colinds[0] = typename ColindsType::HostMirror("colind0", nrows*nnzPerRow);
     colinds[1] = typename ColindsType::HostMirror("colind1", nrows*nnzPerRow);
     colinds[2] = typename ColindsType::HostMirror("colind2", nrows*nnzPerRow);
  {
    //want consistent test results
    srand(12);
    for(LO m = 0; m < 2; m++)
    {
      for(size_t row = 0; row < nrows; row++)
      {
        //rows are sorted, so come up with nnzPerRow vals and cols and then sort them
        std::vector<ISC> rowVals(nnzPerRow);
        std::vector<LO> rowInds(nnzPerRow);
        for(size_t entry = 0; entry < nnzPerRow; entry++)
        {
          rowVals[entry] = ((double) rand()) / RAND_MAX;
          //don't allow repeats in col inds
          LO ind;
          do
          {
            ind = rand() % nrows;
          }
          while(std::find(rowInds.begin(), rowInds.end(), ind) != rowInds.end());
          rowInds[entry] = ind;
        }
        //now insert new coordinates in big arrays
        for(size_t entry = 0; entry < nnzPerRow; entry++)
        {
          vals[m][nnzPerRow * row + entry] = rowVals[entry];
          colinds[m][nnzPerRow * row + entry] = rowInds[entry];
        }
      }
      //set rowptrs
      for(size_t row = 0; row <= nrows; row++)
      {
        rowptrs[m][row] = row * nnzPerRow;
      }
      valsCRS[m] = Kokkos::create_mirror_view_and_copy(
                                      typename NT::memory_space(), vals[m]);
      rowptrsCRS[m] = Kokkos::create_mirror_view_and_copy(
                                      typename NT::memory_space(), rowptrs[m]);
      colindsCRS[m] = Kokkos::create_mirror_view_and_copy(
                                      typename NT::memory_space(), colinds[m]);
    }
  }
  //now run the threaded addition on mats[0] and mats[1]
  ISC zero(0);
  ISC one(1);
  Tpetra::MMdetails::AddKernels<SC, LO, GO, NT>::addUnsorted(
                               valsCRS[0], rowptrsCRS[0], colindsCRS[0], one, 
                               valsCRS[1], rowptrsCRS[1], colindsCRS[1], one, 
                               nrows, valsCRS[2], rowptrsCRS[2], colindsCRS[2]);

  //now scan through C's rows and entries to check they are correct
  TEST_ASSERT(rowptrsCRS[0].extent(0) == rowptrsCRS[2].extent(0));

  vals[2] = Kokkos::create_mirror_view_and_copy(
                           typename ValuesType::HostMirror::memory_space(),
                           valsCRS[2]);
  rowptrs[2] = Kokkos::create_mirror_view_and_copy(
                              typename RowptrsType::HostMirror::memory_space(),
                              rowptrsCRS[2]);
  colinds[2] = Kokkos::create_mirror_view_and_copy(
                              typename ColindsType::HostMirror::memory_space(),
                              colindsCRS[2]);

  for(size_t i = 0; i < nrows; i++)
  {
    //also compute what C's row should be (as dense values)
    std::vector<ISC> correctVals(nrows, zero);
    std::vector<bool> correctEntries(nrows, false);
    for(size_t j = 0; j < nnzPerRow; j++)
    {
      int col1 = colinds[0][i * nnzPerRow + j];
      int col2 = colinds[1][i * nnzPerRow + j];
      correctVals[col1] += vals[0](i * nnzPerRow + j);
      correctEntries[col1] = true;
      correctVals[col2] += vals[1](i * nnzPerRow + j);
      correctEntries[col2] = true;
    }
    size_t actualNNZ = 0;
    for(size_t j = 0; j < nrows; j++)
    {
      if(correctEntries[j])
        actualNNZ++;
    }
    size_t Crowstart = rowptrs[2](i);
    size_t Crowlen = rowptrs[2](i + 1) - Crowstart;
    TEST_ASSERT(Crowlen == actualNNZ);
    for(size_t j = 0; j < Crowlen; j++)
    {
      size_t Coffset = Crowstart + j;
      if(j > 0)
      {
        //Check entries are sorted
        TEST_ASSERT(colinds[2](Coffset - 1) <= colinds[2](Coffset));
      }
      TEST_FLOATING_EQUALITY(vals[2](Coffset), correctVals[colinds[2](Coffset)], 1e-12);
    }
  }
}

namespace AddTestUtils
{

//Local matrices of fill-complete matrices should always be sorted.
//Return true if sorted, false otherwise.
template<typename CrsMat>
bool checkLocallySorted(const CrsMat& A)
{
  using LO = typename CrsMat::local_ordinal_type;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  auto graph = A.getLocalMatrixHost().graph;
  LO numLocalRows = A.getLocalNumRows();
  int allSorted = 1;
  for(int i = 0; i < numLocalRows; i++)
  {
    for(size_t j = graph.row_map(i) + 1; j < graph.row_map(i+1); j++)
    {
      if(graph.entries(j - 1) > graph.entries(j))
      {
        allSorted = 0;
      }
    }
    if(!allSorted)
      break;
  }
  int globalAllSorted = 1;
  auto comm = A.getComm();
  reduceAll<int, int>(*comm, Teuchos::REDUCE_AND, allSorted, outArg(globalAllSorted));
  return globalAllSorted;
}

//Check that A+B=C, assuming all 3 have the smae row map
template<typename CrsMat>
bool verifySum(const CrsMat& A, const CrsMat& B, const CrsMat& C)
{
  using SC = typename CrsMat::scalar_type;
  using LO = typename CrsMat::local_ordinal_type;
  using GO = typename CrsMat::global_ordinal_type;
  using KAT = Kokkos::ArithTraits<SC>;
  using Teuchos::Array;
  typedef typename CrsMat::nonconst_global_inds_host_view_type gids_type;
  typedef typename CrsMat::nonconst_values_host_view_type vals_type;

  auto rowMap = A.getRowMap();
  LO numLocalRows = rowMap->getLocalNumElements();
  GO Amax = A.getGlobalMaxNumRowEntries();
  GO Bmax = B.getGlobalMaxNumRowEntries();
  GO Cmax = C.getGlobalMaxNumRowEntries();
  gids_type Ainds("Ainds",Amax);
  vals_type Avals("Avals",Amax);
  gids_type Binds("Binds",Bmax);
  vals_type Bvals("Bvals",Bmax);
  gids_type Cinds("Cinds",Cmax);
  vals_type Cvals("Cvals",Cmax);
  for(LO i = 0; i < numLocalRows; i++)
  {
    GO gid = rowMap->getGlobalElement(i);
    size_t Aentries;
    size_t Bentries;
    size_t Centries;
    A.getGlobalRowCopy(gid, Ainds, Avals, Aentries);
    B.getGlobalRowCopy(gid, Binds, Bvals, Bentries);
    C.getGlobalRowCopy(gid, Cinds, Cvals, Centries);
    Tpetra::sort2(Ainds, Aentries, Avals);
    Tpetra::sort2(Binds, Bentries, Bvals);
    Tpetra::sort2(Cinds, Centries, Cvals);
    //Now, scan through the row to make sure C's entries match
    size_t Ait = 0;
    size_t Bit = 0;
    for(size_t Cit = 0; Cit < Centries; Cit++)
    {
      GO col = Cinds[Cit];
      typename vals_type::value_type val = Cvals[Cit];
      typename vals_type::value_type goldVal = 0;
      if(Ait < Aentries && Ainds[Ait] == col)
        goldVal += Avals[Ait++];
      if(Bit < Bentries && Binds[Bit] == col)
        goldVal += Bvals[Bit++];
      //Any scalar magnitude should implicitly convert to double
      double err = KAT::abs(val - goldVal);
      if(err > 1e-13)
      {
        std::cerr << "On rank " << rowMap->getComm()->getRank() << ": global row " << gid << ", global col " << col << " has wrong value!" << std::endl;
        return false;
      }
    }
    //In looping over Centries, should have consumed all A and B entries too
    if(Ait != Aentries || Bit != Bentries)
    {
      std::cerr << "On rank " << rowMap->getComm()->getRank() << ": global row " << gid << " has too few entries!" << std::endl;
      return false;
    }
  }
  return true;
}

template<typename LO, typename GO, typename NT>
RCP<const Tpetra::Map<LO, GO, NT>> buildRandomColMap(
    RCP<const Tpetra::Map<LO, GO, NT>> domainMap,
    GO indexBase, GO minCol, GO maxCol, int seed, double proportion)
{
  using map_type = Tpetra::Map<LO, GO, NT>;
  using MemSpace = typename NT::memory_space;
  srand(seed * 21);
  std::vector<GO> present;
  for(GO i = minCol; i < maxCol; i++)
  {
    if(rand() < proportion * RAND_MAX)
      present.push_back(i);
  }
  Kokkos::View<GO*, Kokkos::HostSpace> globalIndsHost("Col GIDs", present.size());
  for(size_t i = 0; i < present.size(); i++)
    globalIndsHost(i) = present[i];
  Kokkos::View<GO*, MemSpace> globalInds = Kokkos::create_mirror_view_and_copy(MemSpace(), globalIndsHost);
  RCP<const map_type> zeroBased;
  Tpetra::Details::makeColMap(zeroBased, domainMap, globalInds);
  if(indexBase == 0)
    return zeroBased;
  //zeroBased always has 0 index base - build the real map with given indexBase
  auto indices = zeroBased->getMyGlobalIndices();
  typedef typename map_type::device_type device_type;
  auto device_indices = Kokkos::create_mirror_view_and_copy(device_type(), indices);
  return rcp(new map_type(zeroBased->getGlobalNumElements(), device_indices, indexBase, zeroBased->getComm()));
}

template<typename SC, typename LO, typename GO, typename NT>
RCP<Tpetra::CrsMatrix<SC, LO, GO, NT>> getTestMatrix(
    RCP<const Tpetra::Map<LO, GO, NT>> rowRangeMap,
    RCP<const Tpetra::Map<LO, GO, NT>> domainMap,
    RCP<const Tpetra::Map<LO, GO, NT>> colMap,
    int seed)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  const LO maxNnzPerRow = colMap->getLocalNumElements();
  auto mat = rcp(new Tpetra::CrsMatrix<SC, LO, GO, NT>(rowRangeMap, colMap, maxNnzPerRow));
  //get consistent results between runs on a given machine
  srand(seed);
  LO numLocalRows = mat->getLocalNumRows();
  for(LO i = 0; i < numLocalRows; i++)
  {
    int n = rand() % maxNnzPerRow;
    Teuchos::Array<SC> vals(n);
    Teuchos::Array<LO> inds(n);
    for(int j = 0; j < n; j++)
    {
      vals[j] = 10.0 * rand() / RAND_MAX;
      inds[j] = rand() % colMap->getLocalNumElements();
    }
    mat->insertLocalValues(i, inds(), vals());
  }
  mat->fillComplete(domainMap, rowRangeMap);
  return mat;
}

template<typename SC, typename LO, typename GO, typename NT>
RCP<Tpetra::CrsMatrix<SC, LO, GO, NT>> getUnsortedTestMatrix(
    RCP<const Tpetra::Map<LO, GO, NT>> rowRangeMap,
    RCP<const Tpetra::Map<LO, GO, NT>> domainMap,
    RCP<const Tpetra::Map<LO, GO, NT>> colMap,
    int seed)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using KCRS = typename crs_matrix_type::local_matrix_device_type;
  using size_type = typename KCRS::row_map_type::non_const_value_type;
  using lno_t =     typename KCRS::index_type::non_const_value_type;
  using kk_scalar_t =  typename KCRS::values_type::non_const_value_type;
  //Get consistent results between runs on a given machine
  srand(seed);
  lno_t numLocalRows = rowRangeMap->getLocalNumElements();
  lno_t numLocalCols = colMap->getLocalNumElements();
  Array<int> entriesPerRow(numLocalRows, 0);
  for(LO i = 0; i < numLocalRows; i++)
    entriesPerRow[i] = rand() % numLocalCols;
  size_type numEntries = 0;
  for(LO i = 0; i < numLocalRows; i++)
    numEntries += entriesPerRow[i];
  Array<lno_t> rowptrs(numLocalRows + 1);
  Array<lno_t> colinds(numEntries);
  Array<kk_scalar_t> values(numEntries);
  size_type accum = 0;
  for(lno_t i = 0; i <= numLocalRows; i++)
  {
    rowptrs[i] = accum;
    if(i == numLocalRows)
      break;
    accum += entriesPerRow[i];
  }
  for(lno_t i = 0; i < numLocalRows; i++)
  {
    Array<lno_t> unused(numLocalCols);
    for(lno_t j = 0; j < numLocalCols; j++)
      unused[j] = j;
    for(lno_t j = 0; j < entriesPerRow[i]; j++)
    {
      //Select a random column from the columns not used yet
      size_t index = rand() % unused.size();
      colinds[rowptrs[i] + j] = unused[index];
      values[rowptrs[i] + j] = 10.0 * rand() / RAND_MAX;
      //efficiently remove index element from unused, don't care about order
      unused[index] = unused.back();
      unused.pop_back();
    }
  }
  //Construct the local matrix
  KCRS lclMatrix(std::string("UnsortedLcl"), numLocalRows, numLocalCols, numEntries, values.data(), rowptrs.data(), colinds.data());
  auto matParams = rcp(new ParameterList);
  matParams->set("sorted", false);
  auto mat = rcp(new Tpetra::CrsMatrix<SC, LO, GO, NT>(rowRangeMap, colMap, lclMatrix, matParams));
  //mat will be returned as fill-complete
  TEUCHOS_TEST_FOR_EXCEPTION(!mat->isFillComplete(), std::logic_error,
      "Test matrix must be fill-complete");
  TEUCHOS_TEST_FOR_EXCEPTION(mat->getCrsGraph()->isSorted(), std::logic_error,
      "This test matrix isn't supposed to be locally sorted.");
  return mat;
}
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMatAdd, same_colmap, SC, LO, GO, NT)
{
  using namespace AddTestUtils;
  using exec_space = typename NT::execution_space;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using map_type = Tpetra::Map<LO, GO, NT>;
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int rank = comm->getRank();
  GO numGlobalRows = 1000;
  GO numGlobalCols = 1000;
  RCP<const map_type> rowDomainMap = rcp(new map_type(numGlobalRows, 0, comm));
  RCP<const map_type> colMap = buildRandomColMap<LO, GO, NT>(rowDomainMap, 0, 0, numGlobalCols, 234 + rank, 0.03);
  //Generate the first matrix without a given column map
  RCP<crs_matrix_type> A = getTestMatrix<SC, LO, GO, NT>(rowDomainMap, rowDomainMap, colMap, rank + 42);
  //Generate the second matrix using the first one's column map
  RCP<crs_matrix_type> B = getTestMatrix<SC, LO, GO, NT>(rowDomainMap, rowDomainMap, colMap, rank + 43);
  auto one = Teuchos::ScalarTraits<SC>::one();
  RCP<crs_matrix_type> C = Tpetra::MatrixMatrix::add(one, false, *A, one, false, *B);
  //Verify fillComplete and that local matrices are sorted
  exec_space().fence();
  TEST_ASSERT(C->isFillComplete());
  TEST_ASSERT(checkLocallySorted(*C));
  TEST_ASSERT(verifySum(*A, *B, *C));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMatAdd, different_col_maps, SC, LO, GO, NT)
{
  using namespace AddTestUtils;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using map_type = Tpetra::Map<LO, GO, NT>;
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  GO numGlobalRows = 1000;
  GO numGlobalCols = 1000;
  GO numLocalCols = 1000 / nprocs;
  //globally square; use same map for row/range/domain
  RCP<const map_type> rowDomainMap = rcp(new map_type(numGlobalRows, 0, comm));
  //Generate A/B column maps each covering a different uniform random subset (30 out of 1000) of all global cols
  GO colStart = numLocalCols * (rank - 1);
  if(colStart < 0)
    colStart = 0;
  GO colEnd = numLocalCols * (rank + 2);
  if(colEnd > numGlobalCols)
    colEnd = numGlobalCols;
  auto Acolmap = buildRandomColMap<LO, GO, NT>(rowDomainMap, 0, colStart, colEnd, 234 + rank, 0.03);
  auto Bcolmap = buildRandomColMap<LO, GO, NT>(rowDomainMap, 0, colStart, colEnd, 236 + rank, 0.03);
  auto A = getTestMatrix<SC, LO, GO, NT>(rowDomainMap, rowDomainMap, Acolmap, 123);
  auto B = getTestMatrix<SC, LO, GO, NT>(rowDomainMap, rowDomainMap, Bcolmap, 321);
  auto one = Teuchos::ScalarTraits<SC>::one();
  RCP<crs_matrix_type> C = Tpetra::MatrixMatrix::add(one, false, *A, one, false, *B);
  TEST_ASSERT(C->isFillComplete());
  TEST_ASSERT(checkLocallySorted(*C));
  TEST_ASSERT(verifySum(*A, *B, *C));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMatAdd, different_index_base, SC, LO, GO, NT)
{
  using namespace AddTestUtils;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using map_type = Tpetra::Map<LO, GO, NT>;
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  GO numGlobalRows = 333;
  GO numGlobalCols = 500;
  GO AindexBase = 100;
  GO numLocalCols = numGlobalCols / nprocs;
  RCP<const map_type> rowMap = rcp(new map_type(numGlobalRows, 0, comm));
  RCP<const map_type> domainMap = rcp(new map_type(numGlobalCols, 0, comm));
  //Generate A/B column maps each covering a different uniform random subset (30 out of 1000) of all global cols,
  //except A's column map has index base 100.
  GO AminCol = numLocalCols * (rank - 1);
  if(AminCol < AindexBase)
    AminCol = AindexBase;
  GO ABmaxCol = numLocalCols * (rank + 2);
  if(ABmaxCol > numGlobalCols)
    ABmaxCol = numGlobalCols;
  GO BminCol = numLocalCols * (rank - 1);
  if(BminCol < 0)
    BminCol = 0;
  auto Acolmap = buildRandomColMap<LO, GO, NT>(domainMap, AindexBase, AminCol, ABmaxCol, 234 + rank, 0.7);
  auto Bcolmap = buildRandomColMap<LO, GO, NT>(domainMap, 0, BminCol, ABmaxCol, 236 + rank, 0.7);
  auto A = getTestMatrix<SC, LO, GO, NT>(rowMap, domainMap, Acolmap, 123);
  auto B = getTestMatrix<SC, LO, GO, NT>(rowMap, domainMap, Bcolmap, 321);
  auto one = Teuchos::ScalarTraits<SC>::one();
  RCP<crs_matrix_type> C = Tpetra::MatrixMatrix::add(one, false, *A, one, false, *B);
  TEST_ASSERT(C->isFillComplete());
  TEST_ASSERT(checkLocallySorted(*C));
  TEST_ASSERT(verifySum(*A, *B, *C));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMatAdd, transposed_b, SC, LO, GO, NT)
{
  using namespace AddTestUtils;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using map_type = Tpetra::Map<LO, GO, NT>;
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  //To compute A + B^T, domain map of A must match range map of B and vice versa.
  GO numGlobalRows = 400;
  GO numGlobalCols = 400;
  GO numLocalCols = numGlobalCols / nprocs;
  //This is the row, domain and range map for both A and B
  RCP<const map_type> map = rcp(new map_type(numGlobalRows, 0, comm));
  auto colmap = buildRandomColMap<LO, GO, NT>(map, 0, numLocalCols * rank, numLocalCols * (rank + 1), 234 + rank, 1.0);
  auto A = getTestMatrix<SC, LO, GO, NT>(map, map, colmap, 123);
  auto B = getTestMatrix<SC, LO, GO, NT>(map, map, colmap, 321);
  auto one = Teuchos::ScalarTraits<SC>::one();
  RCP<crs_matrix_type> C = Tpetra::MatrixMatrix::add(one, false, *A, one, true, *B);
  TEST_ASSERT(C->isFillComplete());
  TEST_ASSERT(checkLocallySorted(*C));
  auto Btrans = Tpetra::RowMatrixTransposer<SC, LO, GO, NT>(B).createTranspose(Teuchos::null);
  TEST_ASSERT(verifySum(*A, *Btrans, *C));
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMatAdd, locally_unsorted, SC, LO, GO, NT)
{
  using namespace AddTestUtils;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using map_type = Tpetra::Map<LO, GO, NT>;
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  //To compute A + B^T, domain map of A must match range map of B and vice versa.
  GO numGlobalRows = 120;
  GO numGlobalCols = 120;
  GO numLocalCols = numGlobalCols / nprocs;
  //This is the row, domain and range map for both A and B
  RCP<const map_type> map = rcp(new map_type(numGlobalRows, 0, comm));
  auto colmap = buildRandomColMap<LO, GO, NT>(map, 0, numLocalCols * rank, numLocalCols * (rank + 1), 234 + rank, 1.0);
  //A and B must have the same column maps
  auto A = getTestMatrix<SC, LO, GO, NT>(map, map, colmap, 123);
  auto B = getUnsortedTestMatrix<SC, LO, GO, NT>(map, map, colmap, 321);
  auto one = Teuchos::ScalarTraits<SC>::one();
  RCP<crs_matrix_type> C = Tpetra::MatrixMatrix::add(one, false, *A, one, false, *B);
  TEST_ASSERT(C->isFillComplete());
  TEST_ASSERT(checkLocallySorted(*C));
  TEST_ASSERT(verifySum(*A, *B, *C));
}

/*
 * This test was added at the request of Chris Siefert
 * Turns out it fails because the original algorithm in
 * EpetraExt doesn't try to handle the case where
 * A has non-unique rows in it's row map. Perhaps I will
 * come back and address that some day.
 * KLN 28/06/2011
 */
/*TEUCHOS_UNIT_TEST(Tpetra_MatMat, Multiple_row_owners){
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  ParameterList defaultParameters;
  RCP<Matrix_t > A = Reader<Matrix_t >::readSparseFile("matrices/denserATa.mtx", comm);
  RCP<Matrix_t > C = Reader<Matrix_t >::readSparseFile("matrices/denserATc.mtx", comm, true, true);
  RowMatrixTransposer<double, int,int,NT> transposer(*A);
  RCP<Matrix_t> AT = transposer.createTranspose();

  RCP<Matrix_t> importedC =
    Tpetra::createCrsMatrix<double, int, int, NT>(AT->getRowMap());
  Tpetra::Import<int,int,NT> importer(C->getRowMap(), importedC->getRowMap());
  importedC->doImport(*C, importer, Tpetra::REPLACE);
  importedC->fillComplete(AT->getDomainMap(), AT->getRangeMap());

  mult_test_results results = multiply_test(
    "AT*A with multipl row owners",
    AT,
    A,
    false,
    false,
    importedC,
    comm,
    out);
  if(verbose){
    out << "Results:" <<endl;
    out << "\tEpsilon: " << results.epsilon << endl;
    out << "\tcNorm: " << results.cNorm << endl;
    out << "\tcompNorm: " << results.compNorm << endl;
  }
  TEST_COMPARE(results.epsilon, <, defaultEpsilon)

}*/


// These are the tests associated to UNIT_TEST_GROUP_SC_LO_GO_NO that work for all backends.
#define UNIT_TEST_GROUP_SC_LO_GO_NO_COMMON( SC, LO, GO, NT )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, range_row_test, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, ATI_range_row_test, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, threaded_add_sorted, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, threaded_add_unsorted, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, add_zero_rows, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMatAdd, same_colmap, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMatAdd, transposed_b, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMatAdd, locally_unsorted, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMatAdd, different_col_maps, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMatAdd, different_index_base, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMatMult, BlockCrsMult, SC, LO, GO, NT)

// FIXME_SYCL requires querying free device memory in KokkosKernels, see
// https://github.com/kokkos/kokkos-kernels/issues/1062.
// The SYCL specifications don't allow asking for that.
#if defined(HAVE_TPETRA_SYCL)
  #define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )  \
    UNIT_TEST_GROUP_SC_LO_GO_NO_COMMON( SC, LO, GO, NT )
#else
  #define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )                                 \
    TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, operations_test,SC, LO, GO, NT) \
    UNIT_TEST_GROUP_SC_LO_GO_NO_COMMON( SC, LO, GO, NT )
#endif

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO_NO )


  } //namespace Tpetra
