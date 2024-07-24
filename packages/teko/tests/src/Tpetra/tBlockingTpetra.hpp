// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tBlockingTpetra_hpp__
#define __tBlockingTpetra_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"

#include <string>

#include "Test_Utils.hpp"
#include "Teko_ConfigDefs.hpp"

namespace Teko {
namespace Test {

class tBlockingTpetra : public UnitTest {
 public:
  virtual ~tBlockingTpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_buildMaps(int verbosity, std::ostream& os);
  bool test_one2many(int verbosity, std::ostream& os);
  bool test_many2one(int verbosity, std::ostream& os);
  bool test_buildSubBlock(int verbosity, std::ostream& os);

 protected:
  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
