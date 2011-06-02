/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


namespace {


int n = 4;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "n", &n, "Number of elements in the vectors" );
}


template<class T>
std::vector<T> generatevector(const int n_in)
{
  using Teuchos::as;
  std::vector<T> a(n_in);
  for( int i = 0; i < n_in; ++i )
    a[i] = as<T>(i); // tests non-const operator[](i)
  return a;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, defaultConstruct, T )
{
  using std::vector;
  using Teuchos::as;
  vector<T> a2;
  TEST_EQUALITY_CONST( as<int>(a2.size()), 0 );
  TEST_EQUALITY_CONST( a2.empty(), true );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, sizedConstruct, T )
{
  using std::vector;
  using Teuchos::as;
  using Teuchos::getConst;
  typedef typename std::vector<T>::size_type size_type;
  vector<T> a(n);
  TEST_EQUALITY_CONST( a.empty(), false );
  TEST_EQUALITY( as<int>(a.size()), n );
  TEST_COMPARE( a.max_size(), >=, as<size_type>(n) );
  TEST_COMPARE( as<int>(a.capacity()), >=, n );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, operatorBracket, T )
{
  using std::vector;
  using Teuchos::as;
  out << "\nTest that a[i] == i ... ";
  vector<T> a = generatevector<T>(n);;
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEST_ARRAY_ELE_EQUALITY( a, i, as<T>(i) );
  }
  if (local_success) out << "passed\n";
  else success = false;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, constAt, T )
{
  using std::vector;
  using Teuchos::as;
  out << "\nTest that a.at(i) == i ...\n";
  vector<T> a = generatevector<T>(n);;
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEUCHOS_TEST_EQUALITY( a.at(i), as<T>(i), out, local_success );
  }
  if (local_success) out << "passed\n";
  else success = false;
}


//
// Instantiations
//


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, defaultConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, sizedConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, operatorBracket, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, constAt, T )


UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)


} // namespace
