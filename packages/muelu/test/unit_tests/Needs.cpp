// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"
#include "MueLu_Needs.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  using MueLu::Needs;

  //this macro declares the unit-test-class:
  TEUCHOS_UNIT_TEST(Needs, Constructor)
  {
    //we are now in a class method declared by the above macro, and
    //that method has these input arguments:
    //Teuchos::FancyOStream& out, bool& success
    out << "version: " << MueLu::Version() << std::endl;
    RCP<Needs> needs = rcp(new Needs() );
    TEST_EQUALITY(needs != Teuchos::null, true);
  }

  TEUCHOS_UNIT_TEST(Needs, NeedRequested)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    std::string aNeed = "knockNeed";
    TEST_EQUALITY(needs.IsRequested(aNeed, NULL), false);
    needs.Request(aNeed, NULL);
    TEST_EQUALITY(needs.IsRequested(aNeed, NULL), true);
  }

  TEUCHOS_UNIT_TEST(Needs, ValueIsAvailable)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    std::string aNeed = "knockNeed";
    TEST_EQUALITY(needs.IsAvailable(aNeed, NULL), false);
    needs.Request(aNeed, NULL);
    needs.Set(aNeed, 42, NULL);
    TEST_EQUALITY(needs.IsAvailable(aNeed, NULL), true);
  }

  TEUCHOS_UNIT_TEST(Needs, NumRequests_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    TEST_THROW( needs.NumRequests("nonExistentNeed", NULL), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, NumRequests)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    std::string aNeed = "knockNeed";
    needs.Request(aNeed, NULL);
    needs.Request(aNeed, NULL);
    TEST_EQUALITY(needs.NumRequests(aNeed, NULL), 2);
  }

  TEUCHOS_UNIT_TEST(Needs, Get_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    TEST_THROW( needs.Get<int>("nonExistentNeed", NULL), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, SetAndGet)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Request(aNeed, NULL);
    needs.Set(aNeed, trueValue, NULL);
    double expectedValue = 0;
    expectedValue = needs.Get<double>(aNeed, NULL);
    TEST_EQUALITY(trueValue, expectedValue);
  }

  TEUCHOS_UNIT_TEST(Needs, Release_Exception)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    TEST_THROW( needs.Release("nonExistentNeed", NULL), std::logic_error );
  }

  TEUCHOS_UNIT_TEST(Needs, Release_Without_Request)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Set(aNeed, trueValue, NULL);
    double expectedValue = 0;
    TEST_THROW( expectedValue = needs.Get<double>(aNeed, NULL), MueLu::Exceptions::RuntimeError );
    TEST_THROW( needs.Release(aNeed, NULL), MueLu::Exceptions::RuntimeError );

    //    RCP<MueLu::TentativePFactory> fac = rcp(new MueLu::TentativePFactory() );
    //    needs.SetData<double>("test", trueValue, &fac);
    //    TEST_THROW( expectedValue = needs.GetData<double>("test", &fac), MueLu::Exceptions::RuntimeError );
    //    TEST_THROW( needs.Release("test", &fac), MueLu::Exceptions::RuntimeError );

  }

  TEUCHOS_UNIT_TEST(Needs, Release)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    std::string aNeed = "knockNeed";
    double trueValue = 42;
    needs.Request(aNeed, NULL);         // TODO: write new test
    needs.Request(aNeed, NULL);
    needs.Set(aNeed, trueValue, NULL);
    double value = needs.Get<double>(aNeed, NULL);
    needs.Release(aNeed, NULL);
    TEST_EQUALITY(trueValue, value);
    TEST_EQUALITY(needs.NumRequests(aNeed, NULL), 1);
    value = needs.Get<double>(aNeed, NULL);
    needs.Release(aNeed, NULL);
    //try to get the need one too many times
    //JG TODO, disable for the moment    TEST_THROW( needs.Get(aNeed, value), std::logic_error );
  }

  class foobarClass {
  public:
    foobarClass() {}
    virtual ~foobarClass() {}
  };

  TEUCHOS_UNIT_TEST(Needs, nonPOD)
  {
    out << "version: " << MueLu::Version() << std::endl;
    Needs needs;
    RCP<foobarClass> trueValue = rcp(new foobarClass);
    needs.Request("foobar", NULL);
    needs.Set("foobar", trueValue, NULL);
    RCP<foobarClass> value;
    value = needs.Get<RCP<foobarClass> >("foobar", NULL);
    TEST_EQUALITY(trueValue, value);
  }

} // namespace MueLuTests
