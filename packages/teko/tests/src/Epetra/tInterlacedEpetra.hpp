// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tInterlacedEpetra_hpp__
#define __tInterlacedEpetra_hpp__

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

class tInterlacedEpetra : public UnitTest {
 public:
  virtual ~tInterlacedEpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_buildSubMaps_num(int verbosity, std::ostream& os);
  bool test_buildSubMaps_vec(int verbosity, std::ostream& os);
  bool test_buildMaps(int verbosity, std::ostream& os);
  bool test_one2many(int verbosity, std::ostream& os);
  bool test_many2one(int verbosity, std::ostream& os);

 protected:
  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
