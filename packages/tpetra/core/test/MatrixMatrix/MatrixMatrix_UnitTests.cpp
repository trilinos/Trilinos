/*
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
*/
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "TpetraExt_MatrixMatrix.hpp"
#include "TpetraExt_TripleMatrixMultiply.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Import_Util.hpp"
#include <cmath>

namespace {
  static const double defaultEpsilon = 1e-10;
  bool verbose = false;
  std::string matnamesFile;

  using Tpetra::MatrixMarket::Reader;
  using Tpetra::CrsMatrix;
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
  Teuchos::ArrayView<const GO> gblRows = identityRowMap->getNodeElementList ();
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
  out << "getIdentityMatrix" << endl;
  Teuchos::OSTab tab1 (out);

  out << "Create CrsMatrix" << endl;
  RCP<Matrix_t> identityMatrix =
    Tpetra::createCrsMatrix<SC, LO, GO, NT> (identityRowMap, 1);

  out << "Fill CrsMatrix" << endl;
  Teuchos::ArrayView<const GO> gblRows = identityRowMap->getNodeElementList ();
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
  RCP<Matrix_t> computedC = rcp( new Matrix_t(rowmap, 1));

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
null_add_test (const Matrix_t& A,
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

  RCP<Matrix_t> computedC = rcp( new Matrix_t(map, 1));

  Tpetra::MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC, false);
  computedC->fillComplete(C->getDomainMap(), C->getRangeMap());

#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedC);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", C);
#endif
  SC one = Teuchos::ScalarTraits<SC>::one();

  RCP<Matrix_t> diffMatrix = Tpetra::createCrsMatrix<SC,LO,GO,NT>(C->getRowMap());
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


  RCP<Matrix_t> diffMatrix = Tpetra::createCrsMatrix<SC,LO,GO,NT>(C->getRowMap());
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


  RCP<Matrix_t> diffMatrix = Tpetra::createCrsMatrix<SC,LO,GO,NT>(Ac->getRowMap());
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

  // computedC1 = leftScaling * (op(A)*op(B)) * rightScaling
  RCP<Matrix_t> computedC1 = rcp( new Matrix_t(map, 1));
  Tpetra::MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC1, false/*call_FillCompleteOnResult*/);
  computedC1->fillComplete(C->getDomainMap(), C->getRangeMap());
  computedC1->leftScale (*leftScaling);
  computedC1->rightScale(*rightScaling);

  // NOTE (mfh 05 Jun 2016) This may even be null.  It exists at this
  // point only for the syntax.
  RCP<NT> node = map->getNode ();

  // As = leftScaling * op(A) =
  //   leftScaling * A, if AT=false
  //   A*leftScaling,   if AT=true
  RCP<Matrix_t> As = A->clone(node);
  if (AT == false) As->leftScale (*leftScaling);
  else             As->rightScale(*leftScaling);

  // Bs = op(B) * rightScaling =
  //   B * rightScaling, if BT=false
  //   rightScaling*B,   if BT=true
  RCP<Matrix_t> Bs = B->clone(node);
  if (BT == false) Bs->rightScale(*rightScaling);
  else             Bs->leftScale (*rightScaling);

  // computedC2 = op(As) * op(Bs)
  RCP<Matrix_t> computedC2 = rcp( new Matrix_t(C->getCrsGraph()) );
  computedC2->fillComplete(C->getDomainMap(), C->getRangeMap());

  computedC2->resumeFill();
  Tpetra::MatrixMatrix::Multiply(*As, AT, *Bs, BT, *computedC2, true/*call_FillCompleteOnResult*/);

  // diffMatrix = computedC2 - computedC1
  SC one = Teuchos::ScalarTraits<SC>::one();
  RCP<Matrix_t> diffMatrix = Tpetra::createCrsMatrix<SC,LO,GO,NT>(C->getRowMap());
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
  RCP<Matrix_t> C2 = rcp(new Matrix_t(B->getRowMap(),0));
  bool done=false;
#ifdef HAVE_TPETRA_INST_OPENMP
    if(std::is_same<NT,Kokkos::Compat::KokkosOpenMPWrapperNode>::value) {
      Teuchos::ParameterList p;
      p.set("openmp: jacobi algorithm","MSAK");
      Tpetra::MatrixMatrix::Jacobi<SC, LO, GO, NT>(omega,Dinv,*A,*B,*C2,true,"jacobi_test_msak",rcp(&p,false));
      done=true;
    }
#endif
#ifdef HAVE_TPETRA_INST_CUDA
    if(std::is_same<NT,Kokkos::Compat::KokkosCudaWrapperNode>::value) {
      Teuchos::ParameterList p;
      p.set("cuda: jacobi algorithm","MSAK");
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
      Tpetra::MatrixMatrix::Add(*AB,false,-one,*B,false,one,C2);
      if(!C2->isFillComplete()) C2->fillComplete(C->getDomainMap(),C->getRangeMap());
    }

  // Check the difference
  RCP<Matrix_t> C_check = rcp(new Matrix_t(B->getRowMap(),0));
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

  // NOTE (mfh 05 Jun 2016) This may even be null.  It exists at this
  // point only for the syntax.
  RCP<NT> node = map->getNode ();
  RCP<Matrix_t> Bs = B->clone(node);
  Bs->rightScale(*rightScaling);

  // computedC2 = (I - Dinv*A)*Bs
  RCP<Matrix_t> computedC2 = rcp( new Matrix_t(computedC1->getCrsGraph()) );
  Tpetra::MatrixMatrix::Jacobi(omega, Dinv, *A, *Bs, *computedC2);

  // diffMatrix = computedC2 - computedC1

  RCP<Matrix_t> diffMatrix = Tpetra::createCrsMatrix<SC,LO,GO,NT>(computedC1->getRowMap());
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


    double epsilon = currentSystem.get<double> ("epsilon", defaultEpsilon);
    std::string op = currentSystem.get<std::string> ("op");

    RCP<Matrix_t> A, B, C, D;

    if (A_file != ""){
      if (A_domainmap_file == "" || A_rangemap_file == "" || A_rowmap_file == "" || A_colmap_file == "")
        A = Reader<Matrix_t>::readSparseFile (A_file, comm);
      else {
        RCP<const map_type> domainmap = Reader<Matrix_t>::readMapFile (A_domainmap_file, comm);
        RCP<const map_type> rangemap  = Reader<Matrix_t>::readMapFile (A_rangemap_file, comm);
        RCP<const map_type> rowmap    = Reader<Matrix_t>::readMapFile (A_rowmap_file, comm);
        RCP<const map_type> colmap    = Reader<Matrix_t>::readMapFile (A_colmap_file, comm);
        A = Reader<Matrix_t>::readSparseFile (A_file, rowmap, colmap, domainmap, rangemap);
      }
    }
    if (B_domainmap_file == "" || B_rangemap_file == "" || B_rowmap_file == "" || B_colmap_file == "")
      B = Reader<Matrix_t>::readSparseFile (B_file, comm);
    else {
      RCP<const map_type> domainmap = Reader<Matrix_t>::readMapFile (B_domainmap_file, comm);
      RCP<const map_type> rangemap  = Reader<Matrix_t>::readMapFile (B_rangemap_file, comm);
      RCP<const map_type> rowmap    = Reader<Matrix_t>::readMapFile (B_rowmap_file, comm);
      RCP<const map_type> colmap    = Reader<Matrix_t>::readMapFile (B_colmap_file, comm);
      B = Reader<Matrix_t>::readSparseFile (B_file, rowmap, colmap, domainmap, rangemap);
    }
    if (C_domainmap_file == "" || C_rangemap_file == "" || C_rowmap_file == "" || C_colmap_file == "")
      C = Reader<Matrix_t>::readSparseFile (C_file, comm);
    else {
      RCP<const map_type> domainmap = Reader<Matrix_t>::readMapFile (C_domainmap_file, comm);
      RCP<const map_type> rangemap  = Reader<Matrix_t>::readMapFile (C_rangemap_file, comm);
      RCP<const map_type> rowmap    = Reader<Matrix_t>::readMapFile (C_rowmap_file, comm);
      RCP<const map_type> colmap    = Reader<Matrix_t>::readMapFile (C_colmap_file, comm);
      C = Reader<Matrix_t>::readSparseFile (C_file, rowmap, colmap, domainmap, rangemap);
    }
    if (D_file != "") {
      if (D_domainmap_file == "" || D_rangemap_file == "" || D_rowmap_file == "" || D_colmap_file == "")
        D = Reader<Matrix_t>::readSparseFile (D_file, comm);
      else {
        RCP<const map_type> domainmap = Reader<Matrix_t>::readMapFile (D_domainmap_file, comm);
        RCP<const map_type> rangemap  = Reader<Matrix_t>::readMapFile (D_rangemap_file, comm);
        RCP<const map_type> rowmap    = Reader<Matrix_t>::readMapFile (D_rowmap_file, comm);
        RCP<const map_type> colmap    = Reader<Matrix_t>::readMapFile (D_colmap_file, comm);
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


      // Do we try Jacobi?
      if (AT == false && BT == false && A->getDomainMap()->isSameAs(*A->getRangeMap())) {
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

      // FIXME (mfh 09 May 2013) This test doesn't currently pass.  I
      // don't think anyone ever exercised the case where C is null on
      // input before.  I'm disabling this test for now until I have a
      // chance to fix that case.
      if (verbose)
        newOut << "Running 3-argument add test (null C on input) for "
               << currentSystem.name() << endl;

      TEUCHOS_TEST_FOR_EXCEPTION(A.is_null (), std::logic_error,
                                 "Before null_add_test: A is null");
      TEUCHOS_TEST_FOR_EXCEPTION(B.is_null (), std::logic_error,
                                 "Before null_add_test: B is null");
      TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
                                 "Before null_add_test: C is null");

      results = null_add_test<Matrix_t> (*A, *B, AT, BT, *C,
                                         newOut, success);

      TEST_COMPARE(results.epsilon, <, epsilon);
      newOut << "Null Add Test Results: " << endl;
      newOut << "\tCorrect Norm: " << results.correctNorm << endl;
      newOut << "\tComputed norm: " << results.computedNorm << endl;
      newOut << "\tEpsilon: " << results.epsilon << endl;

      B = Reader<Matrix_t >::readSparseFile(B_file, comm, false);

      if (! BT) {
        if (verbose)
          newOut << "Running 2-argument add test for "
                 << currentSystem.name() << endl;

        results = add_into_test(A,B,AT,C,comm);
        TEST_COMPARE(results.epsilon, <, epsilon)
        newOut << "Add Into Test Results: " << endl;
        newOut << "\tCorrect Norm: " << results.correctNorm << endl;
        newOut << "\tComputed norm: " << results.computedNorm << endl;
        newOut << "\tEpsilon: " << results.epsilon << endl;
      }
    }
    else if (op == "RAP") {
      // if (verbose)
      //   newOut << "Running RAP multiply test (manual FC) for " << currentSystem.name() << endl;

      // mult_test_results results = multiply_RAP_test_manualfc(name, A, B, C, AT, BT, CT, D, comm, newOut);

      // if (verbose) {
      //   newOut << "Results:"     << endl;
      //   newOut << "\tEpsilon: "  << results.epsilon  << endl;
      //   newOut << "\tcNorm: "    << results.cNorm    << endl;
      //   newOut << "\tcompNorm: " << results.compNorm << endl;
      // }
      // TEST_COMPARE(results.epsilon, <, epsilon);

      if (verbose)
        newOut << "Running multiply test (auto FC) for " << currentSystem.name() << endl;

      mult_test_results results = multiply_RAP_test_autofc(name, A, B, C, AT, BT, CT, D, comm, newOut);

      if (verbose) {
        newOut << "Results:"     << endl;
        newOut << "\tEpsilon: "  << results.epsilon  << endl;
        newOut << "\tcNorm: "    << results.cNorm    << endl;
        newOut << "\tcompNorm: " << results.compNorm << endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon);

      // if (verbose)
      //   newOut << "Running multiply reuse test for " << currentSystem.name() << endl;

      // results = multiply_reuse_test(name, A, B, AT, BT, C, comm, newOut);

      // if (verbose) {
      //   newOut << "Results:"     << endl;
      //   newOut << "\tEpsilon: "  << results.epsilon  << endl;
      //   newOut << "\tcNorm: "    << results.cNorm    << endl;
      //   newOut << "\tcompNorm: " << results.compNorm << endl;
      // }
      // TEST_COMPARE(results.epsilon, <, epsilon);
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

  // This is just to fulfill syntax requirements.  It could even be null.
  RCP<NT> node = dummy->getNode ();

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
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(myRows, comm, node);
  RCP<const Map_t > bRangeMap =
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(rangeElements, comm, node);
  //We divide by 2 to make the matrix tall and "skinny"
  RCP<const Map_t > bDomainMap =
    Tpetra::createUniformContigMapWithNode<LO,GO,NT>(globalNumRows/2, comm, node);

  newOut << "Create identityMatrix" << endl;
  RCP<Matrix_t > identityMatrix =
    getIdentityMatrixWithMap<SC,LO,GO,NT> (newOut, bRowMap, comm);


  newOut << "Create bMatrix" << endl;
  RCP<Matrix_t > bMatrix =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(bRowMap, 1);

  newOut << "Fill bMatrix" << endl;
  {
    Teuchos::ArrayView<const GO> gblRows = bRowMap->getNodeElementList ();
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
  TEST_COMPARE(results.epsilon, <, defaultEpsilon);

  newOut << "Create identity2" << endl;
  RCP<Matrix_t > identity2 =
    getIdentityMatrix<SC,LO,GO,NT> (newOut, globalNumRows/2, comm);

  RCP<const Map_t > bTransRowMap =
    Tpetra::createUniformContigMapWithNode<LO,GO,NT>(globalNumRows/2,comm);

  newOut << "Create and fill bTrans" << endl;
  RCP<Matrix_t> bTrans = Tpetra::createCrsMatrix<SC,LO,GO,NT>(bTransRowMap, 1);
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
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(bTransRangeElements, comm, node);
  RCP<const Map_t > bTransDomainMap =
    Tpetra::createUniformContigMapWithNode<LO,GO,NT>(globalNumRows, comm, node);

  newOut << "Compute identity * transpose(bTrans)" << endl;
  Tpetra::MatrixMatrix::Multiply(*identity2,false,*bMatrix, true, *bTrans, false);

  newOut << "Call fillComplete on bTrans" << endl;
  bTrans->fillComplete(bTransDomainMap, bTransRangeMap);

  newOut << "Create and fill bTransTest" << endl;
  RCP<Matrix_t > bTransTest =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(bTransRowMap, 1);

  {
    Teuchos::ArrayView<const GO> gblRows = bRowMap->getNodeElementList ();
    for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
      Teuchos::Array<GO> col(1,*it);
      Teuchos::Array<SC> val(1,3.0);
      bTransTest->insertGlobalValues((*it)/2, col(), val());
    }
  }

  newOut << "Call fillComplete on bTransTest" << endl;
  bTransTest->fillComplete(bTransDomainMap, bTransRangeMap);

  // FIXME (mfh 03 May 2016) I didn't write this output message.  I
  // don't know what it means, but I'm leaving it here in case it's
  // meaningful to someone.
  newOut << "Regular I*P^T" << endl;

  RCP<Matrix_t > bTransDiff =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(bTransRowMap, 1);
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

  // NOTE (mfh 05 Jun 2016) This may even be null.  It exists at this
  // point only for the syntax.
  RCP<NT> node = identityMatrix->getNode ();

  newOut << "Create Maps for matrix aMat" << endl;
  Array<GO> aMyRows = tuple<GO>(
    rank*numRowsPerProc,
    rank*numRowsPerProc+1,
    rank*numRowsPerProc+2,
    rank*numRowsPerProc+3);
  RCP<const Map_t> aRowMap =
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(aMyRows, comm, node);
  RCP<const Map_t> aDomainMap =
    Tpetra::createUniformContigMapWithNode<LO,GO,NT>(globalNumRows/2, comm, node);
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
    Tpetra::createNonContigMapWithNode<LO,GO,NT>(aRangeElements, comm, node);

  newOut << "Create matrix aMat" << endl;
  RCP<Matrix_t> aMat =
    Tpetra::createCrsMatrix<SC,LO,GO,NT>(aRowMap, 1);

  newOut << "Fill matrix aMat" << endl;
  {
    Teuchos::ArrayView<const GO> gblRows = aRowMap->getNodeElementList ();
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
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> crs_matrix_type;
  typedef typename crs_matrix_type::local_matrix_type KCRS;
  typedef typename crs_matrix_type::impl_scalar_type ISC;
  typedef typename KCRS::values_type::non_const_type ValuesType;
  typedef typename KCRS::row_map_type::non_const_type RowptrsType;
  typedef typename KCRS::index_type::non_const_type ColindsType;
  //The 3 views are for A, B and C
  ValuesType valsCRS[3];
  RowptrsType rowptrsCRS[3];
  ColindsType colindsCRS[3];
  //Populate A and B
  {
    ISC* vals[2] = {new ISC[nrows * nnzPerRow], new ISC[nrows * nnzPerRow]};
    LO* rowptrs[2] = {new LO[nrows * nnzPerRow], new LO[nrows * nnzPerRow]};
    LO* colinds[2] = {new LO[nrows * nnzPerRow], new LO[nrows * nnzPerRow]};
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
      valsCRS[m] = ValuesType("Values", nrows * nnzPerRow);
      rowptrsCRS[m] = RowptrsType("RowPtrs", nrows + 1);
      colindsCRS[m] = ColindsType("ColInds", nrows * nnzPerRow);
      for(size_t i = 0; i < nrows + 1; i++)
      {
        rowptrsCRS[m](i) = rowptrs[m][i];
      }
      for(size_t i = 0; i < nrows * nnzPerRow; i++)
      {
        valsCRS[m](i) = vals[m][i];
        colindsCRS[m](i) = colinds[m][i];
      }
    }
    for(int i = 0; i < 2; i++)
    {
      delete[] vals[i];
      delete[] rowptrs[i];
      delete[] colinds[i];
    }
  }
  //now run the threaded addition on mats[0] and mats[1]
  ISC zero(0);
  ISC one(1);
  Tpetra::MatrixMatrix::AddDetails::AddKernels<SC, LO, GO, NT>::addSorted(valsCRS[0], rowptrsCRS[0], colindsCRS[0], one, valsCRS[1], rowptrsCRS[1], colindsCRS[1], one, valsCRS[2], rowptrsCRS[2], colindsCRS[2]);
  //now scan through C's rows and entries to check they are correct
  TEUCHOS_TEST_FOR_EXCEPTION(rowptrsCRS[0].extent(0) != rowptrsCRS[2].extent(0), std::logic_error,
      "Threaded addition of sorted Kokkos::CrsMatrix returned a matrix with the wrong number of rows.");
  for(size_t i = 0; i < nrows; i++)
  {
    //also compute what C's row should be (as dense values)
    std::vector<ISC> correctVals(nrows, zero);
    std::vector<bool> correctEntries(nrows, false);
    for(size_t j = 0; j < nnzPerRow; j++)
    {
      int col1 = colindsCRS[0](i * nnzPerRow + j);
      int col2 = colindsCRS[1](i * nnzPerRow + j);
      correctVals[col1] += valsCRS[0](i * nnzPerRow + j);
      correctEntries[col1] = true;
      correctVals[col2] += valsCRS[1](i * nnzPerRow + j);
      correctEntries[col2] = true;
    }
    size_t actualNNZ = 0;
    for(size_t j = 0; j < nrows; j++)
    {
      if(correctEntries[j])
        actualNNZ++;
    }
    size_t Crowstart = rowptrsCRS[2](i);
    size_t Crowlen = rowptrsCRS[2](i + 1) - Crowstart;
    TEUCHOS_TEST_FOR_EXCEPTION(Crowlen != actualNNZ, std::logic_error,
        std::string("Threaded addition of sorted Kokkos::CrsMatrix produced row ") + std::to_string(i) + " with the wrong number of entries.");
    for(size_t j = 0; j < Crowlen; j++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(valsCRS[2](Crowstart + j) != correctVals[colindsCRS[2](Crowstart + j)], std::logic_error,
          "Threaded addition of sorted Kokkos::CrsMatrix produced an incorrect value.");
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMat, threaded_add_unsorted, SC, LO, GO, NT)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  using Teuchos::RCP;
  //First, make two local, random, unsorted Kokkos sparse matrices
  size_t nrows = 1000;
  size_t nnzPerRow = 20;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> crs_matrix_type;
  typedef typename crs_matrix_type::local_matrix_type KCRS;
  typedef typename crs_matrix_type::impl_scalar_type ISC;
  typedef typename KCRS::values_type::non_const_type ValuesType;
  typedef typename KCRS::row_map_type::non_const_type RowptrsType;
  typedef typename KCRS::index_type::non_const_type ColindsType;
  //The 3 views are for A, B and C
  ValuesType valsCRS[3];
  RowptrsType rowptrsCRS[3];
  ColindsType colindsCRS[3];
  //Populate A and B
  {
    ISC* vals[2] = {new ISC[nrows * nnzPerRow], new ISC[nrows * nnzPerRow]};
    LO* rowptrs[2] = {new LO[nrows * nnzPerRow], new LO[nrows * nnzPerRow]};
    LO* colinds[2] = {new LO[nrows * nnzPerRow], new LO[nrows * nnzPerRow]};
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
      valsCRS[m] = ValuesType("Values", nrows * nnzPerRow);
      rowptrsCRS[m] = RowptrsType("RowPtrs", nrows + 1);
      colindsCRS[m] = ColindsType("ColInds", nrows * nnzPerRow);
      for(size_t i = 0; i < nrows + 1; i++)
      {
        rowptrsCRS[m](i) = rowptrs[m][i];
      }
      for(size_t i = 0; i < nrows * nnzPerRow; i++)
      {
        valsCRS[m](i) = vals[m][i];
        colindsCRS[m](i) = colinds[m][i];
      }
    }
    for(int i = 0; i < 2; i++)
    {
      delete[] vals[i];
      delete[] rowptrs[i];
      delete[] colinds[i];
    }
  }
  //now run the threaded addition on mats[0] and mats[1]
  ISC zero(0);
  ISC one(1);
  Tpetra::MatrixMatrix::AddDetails::AddKernels<SC, LO, GO, NT>::addUnsorted(valsCRS[0], rowptrsCRS[0], colindsCRS[0], one, valsCRS[1], rowptrsCRS[1], colindsCRS[1], one, nrows, valsCRS[2], rowptrsCRS[2], colindsCRS[2]);
  //now scan through C's rows and entries to check they are correct
  TEUCHOS_TEST_FOR_EXCEPTION(rowptrsCRS[0].extent(0) != rowptrsCRS[2].extent(0), std::logic_error,
      "Threaded addition of sorted Kokkos::CrsMatrix returned a matrix with the wrong number of rows.");
  for(size_t i = 0; i < nrows; i++)
  {
    //also compute what C's row should be (as dense values)
    std::vector<ISC> correctVals(nrows, zero);
    std::vector<bool> correctEntries(nrows, false);
    for(size_t j = 0; j < nnzPerRow; j++)
    {
      int col1 = colindsCRS[0][i * nnzPerRow + j];
      int col2 = colindsCRS[1][i * nnzPerRow + j];
      correctVals[col1] += valsCRS[0](i * nnzPerRow + j);
      correctEntries[col1] = true;
      correctVals[col2] += valsCRS[1](i * nnzPerRow + j);
      correctEntries[col2] = true;
    }
    size_t actualNNZ = 0;
    for(size_t j = 0; j < nrows; j++)
    {
      if(correctEntries[j])
        actualNNZ++;
    }
    size_t Crowstart = rowptrsCRS[2](i);
    size_t Crowlen = rowptrsCRS[2](i + 1) - Crowstart;
    TEUCHOS_TEST_FOR_EXCEPTION(Crowlen != actualNNZ, std::logic_error,
        std::string("Threaded addition of unsorted Kokkos::CrsMatrix produced row ") + std::to_string(i) + " with the wrong number of entries (is " + std::to_string(Crowlen) + ", should be " + std::to_string(actualNNZ) + ")");
    for(size_t j = 0; j < Crowlen; j++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(valsCRS[2](Crowstart + j) != correctVals[colindsCRS[2](Crowstart + j)], std::logic_error,
          "Threaded addition of unsorted Kokkos::CrsMatrix produced an incorrect value.");
    }
  }
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


#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, operations_test,SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, range_row_test, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, ATI_range_row_test, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, threaded_add_sorted, SC, LO, GO, NT) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMat, threaded_add_unsorted, SC, LO, GO, NT)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO_NO )


  } //namespace Tpetra
