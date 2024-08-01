// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tEpetraThyraConverter_hpp__
#define __tEpetraThyraConverter_hpp__

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

class tEpetraThyraConverter : public UnitTest {
 public:
  virtual ~tEpetraThyraConverter() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_functionality(int verbosity, std::ostream& os);

  bool test_single_blockEpetraToThyra(int verbosity, std::ostream& os);
  bool test_single_blockThyraToEpetra(int verbosity, std::ostream& os);
  bool test_blockEpetraToThyra(int verbosity, std::ostream& os);
  bool test_blockThyraToEpetra(int verbosity, std::ostream& os);
};

}  // namespace Test
}  // end namespace Teko

#endif
