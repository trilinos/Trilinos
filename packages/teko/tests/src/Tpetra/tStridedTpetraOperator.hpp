// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tStridedTpetraOperator_hpp__
#define __tStridedTpetraOperator_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

#include <string>

#include "Teko_ConfigDefs.hpp"
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tStridedTpetraOperator : public UnitTest {
 public:
  virtual ~tStridedTpetraOperator() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_numvars_constr(int verbosity, std::ostream& os);
  bool test_vector_constr(int verbosity, std::ostream& os);
  bool test_reorder(int verbosity, std::ostream& os, int total);

 protected:
  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
