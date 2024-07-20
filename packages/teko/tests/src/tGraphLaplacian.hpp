// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tGraphLaplacian_hpp__
#define __tGraphLaplacian_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tGraphLaplacian : public UnitTest {
 public:
  virtual ~tGraphLaplacian() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_single_array(int verbosity, std::ostream& os);
  bool test_multi_array(int verbosity, std::ostream& os);

 protected:
  bool compareMatrix(const Epetra_CrsMatrix& gl, const std::string& name, int verbosity,
                     std::ostream& os) const;
  double tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
