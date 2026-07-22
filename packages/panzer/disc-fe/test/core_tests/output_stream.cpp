// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include "Panzer_OutputStream_DefaultImpl.hpp"
#include <sstream>

namespace panzer {

  class TestObject : public panzer::OutputStreamDefaultImpl {

  public:

    TestObject() {}

    void printValues() const
    {
      
      out() << "Test\n";
      
      if (doOutput(VERB_HIGH))
	out() << "Test2\n";

      // should not be written based on set verb level
      if (doOutput(VERB_LOW))
	out() << "Test4\n";

    }

  };
  

  TEUCHOS_UNIT_TEST(output, default_impl)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::FancyOStream;

    std::ostringstream os;

    RCP<FancyOStream> fos = 
      rcp(new FancyOStream(Teuchos::rcpFromRef(os)));

    TestObject t;
    
    t.setOStream(fos);
    t.setVerbosityLevel(VERB_MEDIUM);

    t.printValues();

    std::string gold_standard;

    gold_standard = "Test\nTest2\n";

    TEST_EQUALITY(gold_standard, os.str());

    TEST_ASSERT(nonnull(t.getOStream()));
    TEST_EQUALITY(t.getVerbosityLevel(), VERB_MEDIUM);

  }

}
