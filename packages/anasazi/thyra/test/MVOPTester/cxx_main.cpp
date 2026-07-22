// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test uses the MVOPTester.hpp functions to test the Anasazi adapters
//  to Tpetra and Thyra.
//

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMVOPTester.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicOutputManager.hpp"

#include "AnasaziThyraAdapter.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraLinearOp.hpp"

#include "Tpetra_Core.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char *argv[])
{
  bool ierr, gerr;
  gerr = true;

  Tpetra::ScopeGuard tpetraScope (&argc,&argv);
  {

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  int MyPID = comm->getRank();
  int NumProcs = comm->getSize();

  // number of global elements
  const int ROWS_PER_PROC = 50;
  int dim = ROWS_PER_PROC * NumProcs;
  const int blockSize = 5;

  bool verbose = false;
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
  }

  typedef double                                 ST;
  typedef Thyra::MultiVectorBase<ST>             MV;
  typedef Thyra::LinearOpBase<ST>                OP;
  typedef Tpetra::MultiVector<ST>               tMV;
  typedef tMV::global_ordinal_type               GO;
  typedef tMV::local_ordinal_type                LO;
  typedef tMV::node_type                       Node;

  // Create an output manager to handle the I/O from the solver
  Teuchos::RCP<Anasazi::OutputManager<ST> > MyOM = Teuchos::rcp( new Anasazi::BasicOutputManager<ST>() );
  if (verbose) {
    MyOM->setVerbosity( Anasazi::Warnings );
  }

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.
  Teuchos::RCP<const Tpetra::Map<LO,GO,Node> > Map = rcp (new Tpetra::Map<LO,GO,Node> (dim,ROWS_PER_PROC,0,comm));
  Teuchos::RCP<Tpetra::CrsMatrix<ST,LO,GO,Node> > A = rcp (new Tpetra::CrsMatrix<ST,LO,GO,Node> (Map, 4));
  LO base = MyPID * ROWS_PER_PROC;
  using Teuchos::tuple;
  if (MyPID != NumProcs-1) {
    for (LO i=0; i<ROWS_PER_PROC; ++i) {
      A->insertGlobalValues(static_cast<GO>(base+i  ), tuple<GO>(base+i  ), tuple<ST>( 2));
      A->insertGlobalValues(static_cast<GO>(base+i  ), tuple<GO>(base+i+1), tuple<ST>(-1));
      A->insertGlobalValues(static_cast<GO>(base+i+1), tuple<GO>(base+i  ), tuple<ST>(-1));
      A->insertGlobalValues(static_cast<GO>(base+i+1), tuple<GO>(base+i+1), tuple<ST>( 2));
    }
  }
  else {
    for (LO i=0; i<ROWS_PER_PROC-1; ++i) {
      A->insertGlobalValues(static_cast<GO>(base+i  ), tuple<GO>(base+i  ), tuple<ST>( 2));
      A->insertGlobalValues(static_cast<GO>(base+i  ), tuple<GO>(base+i+1), tuple<ST>(-1));
      A->insertGlobalValues(static_cast<GO>(base+i+1), tuple<GO>(base+i  ), tuple<ST>(-1));
      A->insertGlobalValues(static_cast<GO>(base+i+1), tuple<GO>(base+i+1), tuple<ST>( 2));
    }
  }
  A->fillComplete();
 
  // first, create a Thyra::VectorSpaceBase from an Tpetra::Map using the Tpetra-Thyra wrappers
  Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > domain_space = Thyra::createVectorSpace<ST>(A->getDomainMap());
  Teuchos::RCP<const Thyra::VectorSpaceBase<ST> > range_space = Thyra::createVectorSpace<ST>(A->getRangeMap());

  // then, create a Thyra::MultiVectorBase from the Thyra::VectorSpaceBase using Thyra creational functions
  Teuchos::RCP<Thyra::MultiVectorBase<ST> > thyra_ivec = Thyra::createMembers(domain_space,blockSize);

  // then, create a Thyra::LinearOpBase from the Tpetra::CrsMatrix using the Tpetra-Thyra wrappers
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > thyra_op = Thyra::tpetraLinearOp<ST,LO,GO,Node>(range_space, domain_space, A);

  // test the Thyra multivector adapter
  ierr = Anasazi::TestMultiVecTraits<ST,MV>(MyOM,thyra_ivec);
  gerr |= ierr;
  if (ierr) {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraAdapter PASSED TestMultiVecTraits()" << std::endl;
  }
  else {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraAdapter FAILED TestMultiVecTraits() ***" << std::endl << std::endl;
  }

  // test the Thyra operator adapter
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,thyra_ivec,thyra_op);
  gerr |= ierr;
  if (ierr) {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraAdapter PASSED TestOperatorTraits()" << std::endl;
  }
  else {
    MyOM->stream(Anasazi::Warnings) << "*** ThyraAdapter FAILED TestOperatorTraits() ***" << std::endl << std::endl;
  }

  if (gerr == false) {
    MyOM->print(Anasazi::Warnings,"End Result: TEST FAILED\n");
    return -1;
  }
  //
  // Default return value
  //
  MyOM->print(Anasazi::Warnings,"End Result: TEST PASSED\n");
  
  }  // end Tpetra::ScopeGuard 

  return 0;

}
