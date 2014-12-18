/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "TpetraExt_TypeStack.hpp"

namespace {

using Teuchos::ParameterList;
using Teuchos::RCP;

template <class TS>
void recurseTypes(Teuchos::FancyOStream &os)
{
  using Teuchos::TypeNameTraits;
  //
  typedef typename TS::type type;
  os << TypeNameTraits<type>::name();
  if (TS::bottom) {
    os << std::endl;
  }
  else {
    os << " -> ";
    recurseTypes<typename TS::next>(os);
  }
}

TEUCHOS_UNIT_TEST( TypeStack, FullStack4 ) {
  // test the type stack
  TPETRAEXT_TYPESTACK4(FullStack,double,
                       float,
                       int,
                       short)
  TEST_EQUALITY_CONST( (int)FullStack::height, 4 );
  recurseTypes<FullStack>(out);
  // test that the macro is equivalent to the manual instantiation
  using Tpetra::Ext::TypeStack;
  typedef TypeStack< double,
          TypeStack< float,
          TypeStack< int ,
                     short  > > >
                     FullStackManual;
  const bool same = Teuchos::TypeTraits::is_same<FullStack,FullStackManual>::value;
  TEST_EQUALITY_CONST( same, true );
}

TEUCHOS_UNIT_TEST( TypeStack, FullStack3 ) {
  // test the type stack
  TPETRAEXT_TYPESTACK3(FullStack,double,
                       float,
                       int)
  TEST_EQUALITY_CONST( (int)FullStack::height, 3 );
  recurseTypes<FullStack>(out);
  // test that the macro is equivalent to the manual instantiation
  using Tpetra::Ext::TypeStack;
  typedef TypeStack< double,
          TypeStack< float,
                     int > >
                     FullStackManual;
  const bool same = Teuchos::TypeTraits::is_same<FullStack,FullStackManual>::value;
  TEST_EQUALITY_CONST( same, true );
}

TEUCHOS_UNIT_TEST( TypeStack, FullStack2 ) {
  // test the type stack
  TPETRAEXT_TYPESTACK2(FullStack,double,
                       float)
  TEST_EQUALITY_CONST( (int)FullStack::height, 2 );
  recurseTypes<FullStack>(out);
  // test that the macro is equivalent to the manual instantiation
  using Tpetra::Ext::TypeStack;
  typedef TypeStack< double,
                     float >
                     FullStackManual;
  const bool same = Teuchos::TypeTraits::is_same<FullStack,FullStackManual>::value;
  TEST_EQUALITY_CONST( same, true );
}

TEUCHOS_UNIT_TEST( TypeStack, FullStack1 ) {
  // test a boring type stack of one type
  // unlike above, this isn't actually a TypeStack object
  TPETRAEXT_TYPESTACK1(FullStack, int)
  TEST_EQUALITY_CONST( (int)FullStack::height, 1 );
  recurseTypes<FullStack>(out);
  // test that the macro is equivalent to the manual instantiation
  using Tpetra::Ext::TypeStackBottom;
  typedef TypeStackBottom<int>
                          FullStackManual;
  const bool same = Teuchos::TypeTraits::is_same<FullStack,FullStackManual>::value;
  TEST_EQUALITY_CONST( same, true );
}

struct TestDBInit {
  template <class T>
  RCP<ParameterList> initDB(ParameterList &params) {
    RCP<ParameterList> db = Teuchos::parameterList();
    db->set<std::string>("type",Teuchos::TypeNameTraits<T>::name());
    return db;
  }
};

TEUCHOS_UNIT_TEST( StackBuilder, DBBuilder ) {
  std::string xmlString(
    " <ParameterList>                   \n"
    "   <ParameterList name='child'>    \n"
    "     <ParameterList name='child'>  \n"
    "     </ParameterList>              \n"
    "   </ParameterList>                \n"
    " </ParameterList>                  \n"
  );
  ParameterList stackPL;
  Teuchos::updateParametersFromXmlString (xmlString, Teuchos::ptr (&stackPL));

  TPETRAEXT_TYPESTACK3( TestStack, double, float, int );
  TestDBInit testInit;
  RCP<ParameterList> stackDB = Tpetra::Ext::initStackDB<TestStack>(stackPL,testInit);
  TEST_EQUALITY_CONST( stackDB == Teuchos::null, false );
  out << *stackDB << std::endl;
}

}


