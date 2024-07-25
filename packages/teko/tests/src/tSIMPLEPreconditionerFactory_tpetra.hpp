// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tSIMPLEPreconditionerFactory_tpetra_hpp__
#define __tSIMPLEPreconditionerFactory_tpetra_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include <string>

#include "Test_Utils.hpp"
#include "Teko_InverseFactory.hpp"

#include "Teko_Utilities.hpp"

namespace Teko {
namespace Test {

class tSIMPLEPreconditionerFactory_tpetra : public UnitTest {
 public:
  virtual ~tSIMPLEPreconditionerFactory_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_createPrec(int verbosity, std::ostream& os);
  bool test_initializePrec(int verbosity, std::ostream& os, int use_blocking);
  bool test_uninitializePrec(int verbosity, std::ostream& os);
  bool test_isCompatable(int verbosity, std::ostream& os);
  bool test_iterativeSolves(int verbosity, std::ostream& os);
  bool test_hierarchicalSolves(int verbosity, std::ostream& os);

  // non-member tests
  bool test_result(int verbosity, std::ostream& os, int use_blocking);
  bool test_identity(int verbosity, std::ostream& os);
  bool test_diagonal(int verbosity, std::ostream& os, int use_blocking);

 protected:
  // some simple matrix subblocks
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > A_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > F_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > B_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > Bt_;
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > C_;

  Teuchos::RCP<InverseFactory> invF_;
  Teuchos::RCP<InverseFactory> invS_;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  // For block diagonal version
  Teuchos::ArrayRCP<GO> block_starts_;
  Teuchos::ArrayRCP<GO> block_gids_;

  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
