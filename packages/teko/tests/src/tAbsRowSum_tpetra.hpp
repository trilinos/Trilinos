// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tAbsRowSum_tpetra_hpp__
#define __tAbsRowSum_tpetra_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"

// Teko includes
#include "Teko_Utilities.hpp"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tAbsRowSum_tpetra : public UnitTest {
 public:
  virtual ~tAbsRowSum_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_absRowSum(int verbosity, std::ostream& os);
  bool test_invAbsRowSum(int verbosity, std::ostream& os);

 protected:
  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
