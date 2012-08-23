// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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
