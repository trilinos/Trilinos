// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tStridedEpetraOperator_hpp__
#define __tStridedEpetraOperator_hpp__

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

class tStridedEpetraOperator : public UnitTest {
 public:
  virtual ~tStridedEpetraOperator() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_numvars_constr(int verbosity, std::ostream& os);
  bool test_vector_constr(int verbosity, std::ostream& os);
  bool test_reorder(int verbosity, std::ostream& os, int total);

 protected:
  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
