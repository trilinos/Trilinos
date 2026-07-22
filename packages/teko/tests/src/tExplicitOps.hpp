// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tExplicitOps_hpp__
#define __tExplicitOps_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

// Teko includes
#include "Teko_Utilities.hpp"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tExplicitOps : public UnitTest {
 public:
  virtual ~tExplicitOps() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_mult_diagScaleMatProd(int verbosity, std::ostream& os);
  bool test_mult_diagScaling(int verbosity, std::ostream& os);
  bool test_mult_modScaleMatProd(int verbosity, std::ostream& os);

  bool test_add(int verbosity, std::ostream& os);
  bool test_add_mod(int verbosity, std::ostream& os);

 protected:
  double tolerance_;

  // matrix to invert
  Teko::ModifiableLinearOp F_;
  Teko::ModifiableLinearOp G_;
  Teko::LinearOp D_;
};

}  // namespace Test
}  // end namespace Teko

#endif
