// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Util.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
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




