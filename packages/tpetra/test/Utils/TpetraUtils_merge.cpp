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

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Util.hpp>
#include <functional>


TEUCHOS_UNIT_TEST( TpetraUtils, Merge2 )
{
  using Tpetra::merge2;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::OSTab;
  using Teuchos::toString;
  using std::endl;
  typedef Array<int>::size_type size_type;
  const size_type origNumEntries = 8;

  Array<int> ind (origNumEntries);
  ind[0] =  0;
  ind[1] =  1;
  ind[2] =  1;
  ind[3] =  3;
  ind[4] = -1;
  ind[5] = -1;
  ind[6] = -1;
  ind[7] =  0;
  Array<int> indCopy = ind; // deep copy

  Array<double> val (origNumEntries);
  val[0] =  42.0;
  val[1] =  -4.0;
  val[2] =  -3.0;
  val[3] =   1.5;
  val[4] =   1.0;
  val[5] =   2.0;
  val[6] =   3.0;
  val[7] = 100.0;
  Array<double> valCopy = val; // deep copy

  const int expNumEntries = 5;
  const int    indExp[] = { 0,    1,   3,  -1,     0  };
  const double valExp[] = {42.0, -7.0, 1.5, 6.0, 100.0};

  // Test merge2 with default merge policy (add).
  {
    Array<int>::iterator indEnd = ind.end ();
    Array<double>::iterator valEnd = val.end ();
    merge2 (indEnd, valEnd, ind.begin (), indEnd, val.begin (), valEnd);

    const size_type newIndLen = indEnd - ind.begin ();
    const size_type newValLen = valEnd - val.begin ();

    TEST_EQUALITY( newIndLen, expNumEntries );
    TEST_EQUALITY( newValLen, expNumEntries );

    const bool indEq = std::equal (ind.begin (), indEnd, indExp);
    const bool valEq = std::equal (val.begin (), valEnd, valExp);

    TEST_EQUALITY( indEq, true );
    TEST_EQUALITY( valEq, true );

    if (! valEq) {
      OSTab tab (out);
      out << "Input value range: " << toString (valCopy ()) << endl;
      out << "Expected output: "
          << toString (ArrayView<const double> ((const double*) valExp, expNumEntries))
          << endl;
      out << "Actual output: " << toString (val.view (0, newValLen)) << endl;
    }
  }

  ind = indCopy; // deep copy; restore original values
  val = valCopy; // deep copy; restore original values

  // Test merge2 with custom merge policy (also add).
  {
    Array<int>::iterator indEnd = ind.end ();
    Array<double>::iterator valEnd = val.end ();
    merge2 (indEnd, valEnd, ind.begin (), indEnd, val.begin (), valEnd, std::plus<double> ());

    const size_type newIndLen = indEnd - ind.begin ();
    const size_type newValLen = valEnd - val.begin ();

    TEST_EQUALITY( newIndLen, expNumEntries );
    TEST_EQUALITY( newValLen, expNumEntries );

    const bool indEq = std::equal (ind.begin (), indEnd, indExp);
    const bool valEq = std::equal (val.begin (), valEnd, valExp);

    TEST_EQUALITY( indEq, true );
    TEST_EQUALITY( valEq, true );

    if (! valEq) {
      OSTab tab (out);
      out << "Input value range: " << toString (valCopy ()) << endl;
      out << "Expected output: "
          << toString (ArrayView<const double> ((const double*) valExp, expNumEntries))
          << endl;
      out << "Actual output: " << toString (val.view (0, newValLen)) << endl;
    }
  }
}


