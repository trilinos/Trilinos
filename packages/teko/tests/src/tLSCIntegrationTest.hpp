// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tLSCIntegrationTest_hpp__
#define __tLSCIntegrationTest_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tLSCIntegrationTest : public UnitTest {
 public:
  virtual ~tLSCIntegrationTest() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_withmassStable(int verbosity, std::ostream& os);
  bool test_nomassStable(int verbosity, std::ostream& os);
  bool test_plConstruction(int verbosity, std::ostream& os);

 protected:
  void loadStableSystem();
  void solveList(Teuchos::ParameterList& paramList, int vcycles);

  double tolerance_;

  Teuchos::RCP<const Epetra_Map> velMap_;   // map of velocity space
  Teuchos::RCP<const Epetra_Map> prsMap_;   // map of pressure space
  Teuchos::RCP<const Epetra_Map> fullMap_;  // map of pressure space

  // stable discretizations matrices
  Teuchos::RCP<const Epetra_CrsMatrix> sF_;
  Teuchos::RCP<const Epetra_CrsMatrix> sB_;
  Teuchos::RCP<const Epetra_CrsMatrix> sBt_;
  Teuchos::RCP<const Epetra_CrsMatrix> sQu_;
  Teuchos::RCP<Epetra_Operator> sA_;

  // stable rhs and IFISS solution
  Teuchos::RCP<Epetra_Vector> rhs_;
  Teuchos::RCP<const Epetra_Vector> sExact_;
};

}  // namespace Test
}  // end namespace Teko

#endif
