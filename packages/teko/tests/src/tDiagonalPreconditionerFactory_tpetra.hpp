// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tDiagonalPreconditionerFactory_tpetra_hpp__
#define __tDiagonalPreconditionerFactory_tpetra_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Teko_PreconditionerFactory.hpp"

#include <string>

#include "Test_Utils.hpp"

#include "Tpetra_CrsMatrix.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {
class DiagonalPreconditionerFactory;
class DiagonalPrecondState;
namespace Test {

class tDiagonalPreconditionerFactory_tpetra : public UnitTest {
 public:
  tDiagonalPreconditionerFactory_tpetra() : fact(0), pstate(0), block_starts(0), block_gids(0) {}
  virtual ~tDiagonalPreconditionerFactory_tpetra();

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream &stdstrm, std::ostream &failstrm, int &totalrun);
  virtual bool isParallel() const { return true; }

  void buildParameterList(int blocksize);

  bool test_createPrec(int verbosity, std::ostream &os, int blocksize);
  bool test_initializePrec(int verbosity, std::ostream &os);
  bool test_canApply(int verbosity, std::ostream &os);

 protected:
  double tolerance_;
  Teuchos::ParameterList List_;
  DiagonalPreconditionerFactory *fact;
  DiagonalPrecondState *pstate;
  LinearOp pop;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > F_;
  Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > tpetraF;
  GO *block_starts, *block_gids;
};

}  // end namespace Test
}  // end namespace Teko

#endif
