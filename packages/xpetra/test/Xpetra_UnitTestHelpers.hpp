// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
#ifndef XPETRA_UNIT_TEST_HELPERS_HPP
#define XPETRA_UNIT_TEST_HELPERS_HPP

#include "Teuchos_UnitTestHelpers.hpp"

/** \brief Basic unit test creation macro for templated code on five template parameters. */
#define TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3, TYPE4, TYPE5) \
  template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5>         \
  class TEST_GROUP##_##TEST_NAME##_UnitTest : public Teuchos::UnitTestBase \
        { \
  public: \
    TEST_GROUP##_##TEST_NAME##_UnitTest( \
      const std::string& type1Name, \
      const std::string& type2Name, \
      const std::string& type3Name, \
      const std::string& type4Name, \
      const std::string& type5Name \
       ) \
      :Teuchos::UnitTestBase( \
         std::string(#TEST_GROUP)+"_"+type1Name+"_"+type2Name+"_"+type3Name+"_"+type4Name+"_"+type5Name, #TEST_NAME ) \
    {} \
    void runUnitTestImpl( Teuchos::FancyOStream &out, bool &success ) const; \
    virtual std::string unitTestFile() const { return __FILE__; } \
    virtual long int unitTestFileLineNumber() const { return __LINE__; } \
  }; \
  \
  template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5>         \
  void TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1,TYPE2,TYPE3,TYPE4, TYPE5>::runUnitTestImpl( \
    Teuchos::FancyOStream &out, bool &success ) const \

/** \brief Template instantiation for five templated types. */
#define TEUCHOS_UNIT_TEST_TEMPLATE_5_INSTANT(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3, TYPE4, TYPE5) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3, TYPE4, TYPE5 >; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3, TYPE4, TYPE5 >     \
  instance_##TEST_GROUP##_##TYPE1##_##TYPE2##_##TYPE3##_##TYPE4##_##TYPE5##_##TEST_NAME##_UnitTest(#TYPE1,#TYPE2,#TYPE3,#TYPE4,#TYPE5);

#endif  // XPETRA_UNIT_TEST_HELPERS_HPP
