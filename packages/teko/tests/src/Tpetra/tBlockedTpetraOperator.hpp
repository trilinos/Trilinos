// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __tBlockedTpetraOperator_hpp__
#define __tBlockedTpetraOperator_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"

#include <string>

#include "Test_Utils.hpp"
#include "Teko_ConfigDefs.hpp"

namespace Teko {
namespace Test {

class tBlockedTpetraOperator : public UnitTest {
 public:
  virtual ~tBlockedTpetraOperator() {}

  virtual void initializeTest();
  virtual int runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm, int& totalrun);
  virtual bool isParallel() const { return true; }

  bool test_vector_constr(int verbosity, std::ostream& os);
  bool test_reorder(int verbosity, std::ostream& os, int total);

 protected:
  void buildBlockGIDs(std::vector<std::vector<GO> >& blocks,
                      const Tpetra::Map<LO, GO, NT>& map) const;

  ST tolerance_;
};

}  // namespace Test
}  // end namespace Teko

#endif
