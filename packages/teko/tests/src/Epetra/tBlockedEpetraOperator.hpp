// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tBlockedEpetraOperator_hpp__
#define __tBlockedEpetraOperator_hpp__

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

class tBlockedEpetraOperator : public UnitTest {
 public:
  virtual ~tBlockedEpetraOperator() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_vector_constr(int verbosity, std::ostream& os);
  bool test_reorder(int verbosity, std::ostream& os, int total);

 protected:
  void buildBlockGIDs(std::vector<std::vector<int> >& blocks, const Epetra_Map& map) const;

  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
