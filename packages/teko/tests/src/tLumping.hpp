// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tLumping_hpp__
#define __tLumping_hpp__

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

class tLumping : public UnitTest {
 public:
  virtual ~tLumping() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_lumping(int verbosity, std::ostream& os);
  bool test_invLumping(int verbosity, std::ostream& os);

 protected:
  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
