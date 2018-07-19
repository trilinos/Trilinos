/*
//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
//@HEADER
*/

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "BelosConfigDefs.hpp"
#include "BelosMVOPTester.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosOutputManager.hpp"

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Tpetra::Map;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Tpetra::MultiVector;
  using Tpetra::CrsMatrix;
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Belos::OutputManager;
  using Belos::Warnings;
  using Teuchos::tuple;

  typedef Tpetra::Map<>::node_type Node;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return Tpetra::getDefaultComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  template<class Scalar, class O1, class O2>
  RCP<CrsMatrix<Scalar,O1,O2,Node> > constructDiagMatrix(const RCP<const Map<O1,O2,Node> > &map)
  {
    RCP<CrsMatrix<Scalar,O1,O2,Node> > op = rcp( new CrsMatrix<Scalar,O1,O2,Node>(map,1) );
    for (size_t i=0; i<map->getNodeNumElements(); ++i) {
      op->insertGlobalValues(map->getGlobalElement(i),tuple(map->getGlobalElement(i)), tuple(ScalarTraits<Scalar>::one()));
    }
    op->fillComplete();
    return op;
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, MVTestDist, O1, O2, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,O1,O2> MV;
    const O2 dim = 500;
    const Teuchos_Ordinal numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<Map<O1,O2,Node> > map = rcp( new Map<O1,O2,Node>(dim,0,comm) );
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestMultiVecTraits<Scalar,MV>(MyOM,mvec);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, MVTestLocal, O1, O2, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,O1,O2> MV;
    const O2 dim = 500;
    const Teuchos_Ordinal numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<Map<O1,O2,Node> > map = rcp(new Map<O1,O2,Node>(dim,0,comm,Tpetra::LocallyReplicated) );
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestMultiVecTraits<Scalar,MV>(MyOM,mvec);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, OPTestLocal, O1, O2, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,O1,O2> MV;
    typedef Tpetra::Operator<Scalar,O1,O2>    OP;
    const O2 dim = 500;
    const Teuchos_Ordinal numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<Map<O1,O2,Node> > map = rcp(new Map<O1,O2,Node>(dim,0,comm,Tpetra::LocallyReplicated) );
    // create a CrsMatrix
    RCP<OP> op = constructDiagMatrix<Scalar,O1,O2>(map);
    // create a multivector
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestOperatorTraits<Scalar,MV,OP>(MyOM,mvec,op);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, OPTestDist, O1, O2, Scalar )
  {
    typedef Tpetra::MultiVector<Scalar,O1,O2> MV;
    typedef Tpetra::Operator<Scalar,O1,O2>    OP;
    const O2 dim = 500;
    const Teuchos_Ordinal numVecs = 5;
    // Create an output manager to handle the I/O from the solver
    RCP<OutputManager<Scalar> > MyOM = rcp( new OutputManager<Scalar>(Warnings,rcp(&out,false)) );
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a uniform contiguous map
    RCP<Map<O1,O2,Node> > map = rcp( new Map<O1,O2,Node>(dim,0,comm) );
    // create a CrsMatrix
    RCP<OP> op = constructDiagMatrix<Scalar,O1,O2>(map);
    // create a multivector
    RCP<MV> mvec = rcp( new MV(map,numVecs,true) );
    bool res = Belos::TestOperatorTraits<Scalar,MV,OP>(MyOM,mvec,op);
    TEST_EQUALITY_CONST(res,true);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, Teuchos::outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( SCALAR, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, MVTestDist, LO, GO, SCALAR ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, MVTestLocal, LO, GO, SCALAR ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, OPTestDist, LO, GO, SCALAR ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, OPTestLocal, LO, GO, SCALAR )

#include "TpetraCore_ETIHelperMacros.h"

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )
}
