// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tParallelInverse_tpetra_hpp__
#define __tParallelInverse_tpetra_hpp__

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

class tParallelInverse_tpetra : public UnitTest {
 public:
  virtual ~tParallelInverse_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_inverse(int verbosity, std::ostream& os);
  bool test_stridedInverse(int verbosity, std::ostream& os);

 protected:
  void loadMatrix();
  void loadStridedMatrix();

  ST tolerance_;

  // matrix to invert
  Teko::LinearOp F_;
};

}  // namespace Test
}  // end namespace Teko

#endif
