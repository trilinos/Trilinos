// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tLSCIntegrationTest_tpetra_hpp__
#define __tLSCIntegrationTest_tpetra_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"

#include <string>

#include "Test_Utils.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {
namespace Test {

class tLSCIntegrationTest_tpetra : public UnitTest {
 public:
  virtual ~tLSCIntegrationTest_tpetra() {}

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

  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > velMap_;   // map of velocity space
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > prsMap_;   // map of pressure space
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > fullMap_;  // map of pressure space

  // stable discretizations matrices
  Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > sF_;
  Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > sB_;
  Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > sBt_;
  Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > sQu_;
  Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> > sA_;

  // stable rhs and IFISS solution
  Teuchos::RCP<Tpetra::Vector<ST, LO, GO, NT> > rhs_;
  Teuchos::RCP<const Tpetra::Vector<ST, LO, GO, NT> > sExact_;
};

}  // namespace Test
}  // end namespace Teko

#endif
