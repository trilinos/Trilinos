// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tGraphLaplacian_tpetra_hpp__
#define __tGraphLaplacian_tpetra_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include "Tpetra_CrsMatrix.hpp"

#include <string>

#include "Test_Utils.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {
namespace Test {

class tGraphLaplacian_tpetra : public UnitTest {
 public:
  virtual ~tGraphLaplacian_tpetra() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return false; }

  bool test_single_array(int verbosity, std::ostream& os);
  bool test_multi_array(int verbosity, std::ostream& os);

 protected:
  bool compareMatrix(const Tpetra::CrsMatrix<ST, LO, GO, NT>& gl, const std::string& name,
                     int verbosity, std::ostream& os) const;
  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
