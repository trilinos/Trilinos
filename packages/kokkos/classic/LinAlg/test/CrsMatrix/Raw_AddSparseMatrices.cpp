//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include <Kokkos_ConfigDefs.hpp>
#include <Kokkos_Raw_addSparseMatrices_decl.hpp>
#include <Kokkos_Raw_addSparseMatrices_def.hpp>
#include <Teuchos_Array.hpp> // bounds checking is useful for debugging
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_toString.hpp>
#include <algorithm>

namespace {

/// \brief Test whether adding two sparse matrices works correctly.
///
/// Use linearity of sparse matrices' row sums (with respect to
/// multiplication by a dense vector on a right) to test sparse matrix
/// addition.
template<class OffsetType, class OrdinalType, class ScalarType>
typename Teuchos::ScalarTraits<ScalarType>::magnitudeType
testAddSparseMatrices (const Teuchos::ArrayView<const OffsetType>& ptrResult,
                       const Teuchos::ArrayView<const OrdinalType>& indResult,
                       const Teuchos::ArrayView<const ScalarType>& valResult,
                       const ScalarType alpha,
                       const Teuchos::ArrayView<const OffsetType> ptr1,
                       const Teuchos::ArrayView<const OrdinalType> ind1,
                       const Teuchos::ArrayView<const ScalarType> val1,
                       const ScalarType beta,
                       const Teuchos::ArrayView<const OffsetType> ptr2,
                       const Teuchos::ArrayView<const OrdinalType> ind2,
                       const Teuchos::ArrayView<const ScalarType> val2,
                       const OffsetType numRows,
                       const OffsetType numCols,
                       Teuchos::FancyOStream& out)
{
  using Teuchos::Array;
  using std::endl;
  typedef ScalarType ST;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef typename STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;

  out << "Testing sparse matrix add with alpha = " << alpha
      << " and beta = " << beta << endl;

  Array<ST> x (numCols);
  Array<ST> y (numRows);
  Array<ST> y_1 (numRows);
  Array<ST> y_2 (numRows);

  out << "Allocating vectors" << endl;

  // Filling the input vector with ones means that sparse mat-vec
  // computes row sums.  Sparse matrix addition must be linear with
  // respect to row sums.
  std::fill (x.begin (), x.end (), STS::one ());
  std::fill (y.begin (), y.end (), STS::zero ());
  std::fill (y_1.begin (), y_1.end (), STS::zero ());
  std::fill (y_2.begin (), y_2.end (), STS::zero ());

  out << "Doing sparse matrix-vector multiplies" << endl;

  for (OrdinalType i = 0; i < numRows; ++i) {
    for (OffsetType k = ptrResult[i]; k < ptrResult[i+1]; ++k) {
      y[i] += valResult[k] * x[indResult[k]];
    }
  }

  for (OrdinalType i = 0; i < numRows; ++i) {
    for (OffsetType k = ptr1[i]; k < ptr1[i+1]; ++k) {
      y_1[i] += val1[k] * x[ind1[k]];
    }
  }

  for (OrdinalType i = 0; i < numRows; ++i) {
    for (OffsetType k = ptr2[i]; k < ptr2[i+1]; ++k) {
      y_2[i] += val2[k] * x[ind2[k]];
    }
  }

  out << "Testing sums" << endl;

  // Now alpha*y_1 + beta*y_2 should equal y, up to a reasonable tolerance.
  MT maxDiff = STM::zero ();
  for (OrdinalType i = 0; i < numRows; ++i) {
    const MT curDiff = STS::magnitude (y[i] - (alpha * y_1[i] + beta * y_2[i]));
    maxDiff = std::max (maxDiff, curDiff);
  }

  return maxDiff;
}


template<class OffsetType, class OrdinalType, class ScalarType>
void
testMatrix1 (OffsetType*& ptr,
             OrdinalType*& ind,
             ScalarType*& val,
             const OrdinalType numRows,
             const OrdinalType numCols)
{
  using Teuchos::as;

  // Subdiagonal and diagonal filled on odd-numbered rows; no entries
  // on even-numbered rows.  Allocate extra space; doesn't matter.
  OffsetType* ptrOut = NULL;
  OrdinalType* indOut = NULL;
  ScalarType* valOut = NULL;
  try {
    ptrOut = new OffsetType [numRows+1];
    indOut = new OrdinalType [2*numRows];
    valOut = new ScalarType [2*numRows];
  } catch (...) {
    if (ptrOut != NULL) {
      delete [] ptrOut;
    }
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
    throw;
  }

  OffsetType count = 0;
  for (OrdinalType i = 0; i < numRows; ++i) {
    ptrOut[i] = count;
    if (i % 2 != 0 && i < numCols) {
      if (i > 0) {
        indOut[count] = i - 1;
        valOut[count] = as<ScalarType> (-17 + i / 2); // just some value
        ++count;
      }
      indOut[count] = i;
      valOut[count] = as<ScalarType> (i * 7); // just some value
      ++count;
    }
  }
  ptrOut[numRows] = count;

  ptr = ptrOut;
  ind = indOut;
  val = valOut;
}



template<class OffsetType, class OrdinalType, class ScalarType>
void
testMatrix2 (OffsetType*& ptr,
             OrdinalType*& ind,
             ScalarType*& val,
             const OrdinalType numRows,
             const OrdinalType numCols)
{
  using Teuchos::as;

  // Antidiagonal filled in every row.
  OffsetType* ptrOut = NULL;
  OrdinalType* indOut = NULL;
  ScalarType* valOut = NULL;
  try {
    ptrOut = new OffsetType [numRows+1];
    indOut = new OrdinalType [numRows];
    valOut = new ScalarType [numRows];
  } catch (...) {
    if (ptrOut != NULL) {
      delete [] ptrOut;
    }
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
    throw;
  }

  OffsetType count = 0;
  for (OrdinalType i = 0; i < numRows; ++i) {
    ptrOut[i] = count;
    if (i < numCols) {
      indOut[count] = (numRows - 1) - i; // antidiagonal entry
      valOut[count] = as<ScalarType> (-2 + 3 * i); // just some value
      ++count;
    }
  }
  ptrOut[numRows] = count;

  ptr = ptrOut;
  ind = indOut;
  val = valOut;
}

} // namespace (anonymous)


TEUCHOS_UNIT_TEST( KokkosRaw, AddSparseMatrices )
{
  using Teuchos::ArrayView;
  using Teuchos::toString;
  using std::endl;
  typedef ptrdiff_t OffsetType;
  typedef int OrdinalType;
  typedef double ScalarType;
  typedef ScalarType ST;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef Teuchos::ScalarTraits<ST>::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;

  const OrdinalType numRows = 42;
  const OrdinalType numCols = 50;

  out << "Constructing test matrix 1" << endl;
  OffsetType* ptr1 = NULL;
  OrdinalType* ind1 = NULL;
  ScalarType* val1 = NULL;
  try {
    testMatrix1<OffsetType, OrdinalType, ST> (ptr1, ind1, val1, numRows, numCols);
  } catch (...) {
    if (ptr1 != NULL) {
      delete [] ptr1;
    }
    if (ind1 != NULL) {
      delete [] ind1;
    }
    if (val1 != NULL) {
      delete [] val1;
    }
    throw;
  }

  ArrayView<const OffsetType> ptr1_view (const_cast<const OffsetType*> (ptr1), numRows+1);
  const OffsetType numEntries1 = ptr1_view[numRows];
  ArrayView<const OrdinalType> ind1_view (const_cast<const OrdinalType*> (ind1), numEntries1);
  ArrayView<const ScalarType> val1_view (const_cast<const ScalarType*> (val1), numEntries1);

  out << "Constructing test matrix 2" << endl;
  OffsetType* ptr2 = NULL;
  OrdinalType* ind2 = NULL;
  ScalarType* val2 = NULL;
  try {
    testMatrix2<OffsetType, OrdinalType, ST> (ptr2, ind2, val2, numRows, numCols);
  } catch (...) {
    if (ptr2 != NULL) {
      delete [] ptr2;
    }
    if (ind2 != NULL) {
      delete [] ind2;
    }
    if (val2 != NULL) {
      delete [] val2;
    }
    throw;
  }

  ArrayView<const OffsetType> ptr2_view (const_cast<const OffsetType*> (ptr2), numRows+1);
  const OffsetType numEntries2 = ptr2_view[numRows];
  ArrayView<const OrdinalType> ind2_view (const_cast<const OrdinalType*> (ind2), numEntries2);
  ArrayView<const ScalarType> val2_view (const_cast<const ScalarType*> (val2), numEntries2);

  out << "Test alpha = 1, beta = 1" << endl;
  {
    const ScalarType alpha = STS::one ();
    const ScalarType beta = STS::one ();
    OffsetType* ptrOut = NULL;
    OrdinalType* indOut = NULL;
    ScalarType* valOut = NULL;
    Kokkos::Raw::addSparseMatrices<OffsetType, OrdinalType, ScalarType> (
      ptrOut, indOut, valOut,
      alpha, ptr1, ind1, val1,
      beta, ptr2, ind2, val2,
      numRows);

    MT maxDiff = STM::zero ();
    try {
      TEUCHOS_TEST_FOR_EXCEPTION(numRows > 0 && ptrOut == NULL, std::logic_error,
        "addSparseMatrices failed to allocate ptrOut; it should have numRows+1 = "
        << numRows+1 << " entries.");
      ArrayView<const OffsetType> ptrOut_view (const_cast<const OffsetType*> (ptrOut), numRows+1);
      const OffsetType numEntriesOut = ptrOut_view[numRows];
      // out << "Output matrix has " << numRows << " rows and " << numEntriesOut
      //     << " entries." << endl << "ptrOut = " << toString (ptrOut_view) << endl;
      ArrayView<const OrdinalType> indOut_view (const_cast<const OrdinalType*> (indOut),
                                                numEntriesOut);
      // out << "indOut = " << toString (indOut_view) << endl;
      ArrayView<const ScalarType> valOut_view (const_cast<const ScalarType*> (valOut),
                                               numEntriesOut);
      // out << "valOut = " << toString (valOut_view) << endl;

      maxDiff = testAddSparseMatrices<OffsetType, OrdinalType, ST> (
        ptrOut_view, indOut_view, valOut_view,
        alpha, ptr1_view, ind1_view, val1_view,
        beta, ptr2_view, ind2_view, val2_view,
        numRows, numCols, out);
    } catch (...) {
      if (ptrOut != NULL) {
        delete [] ptrOut;
      }
      if (indOut != NULL) {
        delete [] indOut;
      }
      if (valOut != NULL) {
        delete [] valOut;
      }
      throw;
    }

    out << "maxDiff: " << maxDiff << endl;

    // The error bound is an upper bound for each row.  The sum of the
    // two test matrices contains at most three entries per row.
    // Thus, computing each entry of the vector resulting from
    // multiplying the sum by the vector of all ones (which is what
    // the test does) involves a sum of at most three nonzero values.
    // We multiply this by 10 to allow for a bit more error.
    const bool notTooBad = (maxDiff < Teuchos::as<ST> (30) * STS::eps ());
    TEST_EQUALITY_CONST( notTooBad, true );

    if (ptrOut != NULL) {
      delete [] ptrOut;
    }
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
  }

  out << "Test alpha = -3, beta = 4" << endl;
  {
    const ScalarType alpha = Teuchos::as<ST> (-3);
    const ScalarType beta = Teuchos::as<ST> (4);
    OffsetType* ptrOut = NULL;
    OrdinalType* indOut = NULL;
    ScalarType* valOut = NULL;
    Kokkos::Raw::addSparseMatrices<OffsetType, OrdinalType, ScalarType> (
      ptrOut, indOut, valOut,
      alpha, ptr1, ind1, val1,
      beta, ptr2, ind2, val2,
      numRows);
    MT maxDiff = STM::zero ();
    try {
      ArrayView<const OffsetType> ptr1_view (const_cast<const OffsetType*> (ptr1), numRows+1);
      const OffsetType numEntries1 = ptr1_view[numRows];
      ArrayView<const OrdinalType> ind1_view (const_cast<const OrdinalType*> (ind1), numEntries1);
      ArrayView<const ScalarType> val1_view (const_cast<const ScalarType*> (val1), numEntries1);

      ArrayView<const OffsetType> ptr2_view (const_cast<const OffsetType*> (ptr2), numRows+1);
      const OffsetType numEntries2 = ptr2_view[numRows];
      ArrayView<const OrdinalType> ind2_view (const_cast<const OrdinalType*> (ind2), numEntries2);
      ArrayView<const ScalarType> val2_view (const_cast<const ScalarType*> (val2), numEntries2);

      ArrayView<const OffsetType> ptrOut_view (const_cast<const OffsetType*> (ptrOut), numRows+1);
      const OffsetType numEntriesOut = ptrOut_view[numRows];
      ArrayView<const OrdinalType> indOut_view (const_cast<const OrdinalType*> (indOut), numEntriesOut);
      ArrayView<const ScalarType> valOut_view (const_cast<const ScalarType*> (valOut), numEntriesOut);

      maxDiff = testAddSparseMatrices<OffsetType, OrdinalType, ST> (
        ptrOut_view, indOut_view, valOut_view,
        alpha, ptr1_view, ind1_view, val1_view,
        beta, ptr2_view, ind2_view, val2_view,
        numRows, numCols, out);
    } catch (...) {
      if (ptrOut != NULL) {
        delete [] ptrOut;
      }
      if (indOut != NULL) {
        delete [] indOut;
      }
      if (valOut != NULL) {
        delete [] valOut;
      }
      throw;
    }

    out << "maxDiff: " << maxDiff << endl;

    // The error bound is an upper bound for each row.  The sum of the
    // two test matrices contains at most three entries per row.
    // Thus, computing each entry of the vector resulting from
    // multiplying the sum by the vector of all ones (which is what
    // the test does) involves a sum of at most three nonzero values.
    // Multiply this by 2 for alpha and beta not both one, and
    // multiply this by 10 to allow for a bit more error.
    const bool notTooBad = (maxDiff < Teuchos::as<ST> (60) * STS::eps ());
    TEST_EQUALITY_CONST( notTooBad, true );

    if (ptrOut != NULL) {
      delete [] ptrOut;
    }
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
  }


  out << "Test alpha = 3, beta = 0" << endl;
  {
    const ScalarType alpha = Teuchos::as<ST> (3);
    const ScalarType beta = Teuchos::as<ST> (0);
    OffsetType* ptrOut = NULL;
    OrdinalType* indOut = NULL;
    ScalarType* valOut = NULL;
    Kokkos::Raw::addSparseMatrices<OffsetType, OrdinalType, ScalarType> (
      ptrOut, indOut, valOut,
      alpha, ptr1, ind1, val1,
      beta, ptr2, ind2, val2,
      numRows);
    MT maxDiff = STM::zero ();
    try {
      ArrayView<const OffsetType> ptr1_view (const_cast<const OffsetType*> (ptr1), numRows+1);
      const OffsetType numEntries1 = ptr1_view[numRows];
      ArrayView<const OrdinalType> ind1_view (const_cast<const OrdinalType*> (ind1), numEntries1);
      ArrayView<const ScalarType> val1_view (const_cast<const ScalarType*> (val1), numEntries1);

      ArrayView<const OffsetType> ptr2_view (const_cast<const OffsetType*> (ptr2), numRows+1);
      const OffsetType numEntries2 = ptr2_view[numRows];
      ArrayView<const OrdinalType> ind2_view (const_cast<const OrdinalType*> (ind2), numEntries2);
      ArrayView<const ScalarType> val2_view (const_cast<const ScalarType*> (val2), numEntries2);

      ArrayView<const OffsetType> ptrOut_view (const_cast<const OffsetType*> (ptrOut), numRows+1);
      const OffsetType numEntriesOut = ptrOut_view[numRows];
      ArrayView<const OrdinalType> indOut_view (const_cast<const OrdinalType*> (indOut), numEntriesOut);
      ArrayView<const ScalarType> valOut_view (const_cast<const ScalarType*> (valOut), numEntriesOut);

      maxDiff = testAddSparseMatrices<OffsetType, OrdinalType, ST> (
        ptrOut_view, indOut_view, valOut_view,
        alpha, ptr1_view, ind1_view, val1_view,
        beta, ptr2_view, ind2_view, val2_view,
        numRows, numCols, out);
    } catch (...) {
      if (ptrOut != NULL) {
        delete [] ptrOut;
      }
      if (indOut != NULL) {
        delete [] indOut;
      }
      if (valOut != NULL) {
        delete [] valOut;
      }
      throw;
    }

    out << "maxDiff: " << maxDiff << endl;

    // The error bound is an upper bound for each row.  The sum of the
    // two test matrices contains at most three entries per row.
    // Thus, computing each entry of the vector resulting from
    // multiplying the sum by the vector of all ones (which is what
    // the test does) involves a sum of at most three nonzero values.
    // This is a generous upper bound, since beta = 0 in this case.
    // Multiply this by 10 to allow for a bit more error.
    const bool notTooBad = (maxDiff < Teuchos::as<ST> (30) * STS::eps ());
    TEST_EQUALITY_CONST( notTooBad, true );

    if (ptrOut != NULL) {
      delete [] ptrOut;
    }
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
  }

  out << "Test alpha = 0, beta = 42" << endl;
  {
    const ScalarType alpha = Teuchos::as<ST> (0);
    const ScalarType beta = Teuchos::as<ST> (42);
    OffsetType* ptrOut = NULL;
    OrdinalType* indOut = NULL;
    ScalarType* valOut = NULL;
    Kokkos::Raw::addSparseMatrices<OffsetType, OrdinalType, ScalarType> (
      ptrOut, indOut, valOut,
      alpha, ptr1, ind1, val1,
      beta, ptr2, ind2, val2,
      numRows);
    MT maxDiff = STM::zero ();
    try {
      ArrayView<const OffsetType> ptr1_view (const_cast<const OffsetType*> (ptr1), numRows+1);
      const OffsetType numEntries1 = ptr1_view[numRows];
      ArrayView<const OrdinalType> ind1_view (const_cast<const OrdinalType*> (ind1), numEntries1);
      ArrayView<const ScalarType> val1_view (const_cast<const ScalarType*> (val1), numEntries1);

      ArrayView<const OffsetType> ptr2_view (const_cast<const OffsetType*> (ptr2), numRows+1);
      const OffsetType numEntries2 = ptr2_view[numRows];
      ArrayView<const OrdinalType> ind2_view (const_cast<const OrdinalType*> (ind2), numEntries2);
      ArrayView<const ScalarType> val2_view (const_cast<const ScalarType*> (val2), numEntries2);

      ArrayView<const OffsetType> ptrOut_view (const_cast<const OffsetType*> (ptrOut), numRows+1);
      const OffsetType numEntriesOut = ptrOut_view[numRows];
      ArrayView<const OrdinalType> indOut_view (const_cast<const OrdinalType*> (indOut), numEntriesOut);
      ArrayView<const ScalarType> valOut_view (const_cast<const ScalarType*> (valOut), numEntriesOut);

      maxDiff = testAddSparseMatrices<OffsetType, OrdinalType, ST> (
        ptrOut_view, indOut_view, valOut_view,
        alpha, ptr1_view, ind1_view, val1_view,
        beta, ptr2_view, ind2_view, val2_view,
        numRows, numCols, out);
    } catch (...) {
      if (ptrOut != NULL) {
        delete [] ptrOut;
      }
      if (indOut != NULL) {
        delete [] indOut;
      }
      if (valOut != NULL) {
        delete [] valOut;
      }
      throw;
    }

    out << "maxDiff: " << maxDiff << endl;

    // The error bound is an upper bound for each row.  The sum of the
    // two test matrices contains at most three entries per row.
    // Thus, computing each entry of the vector resulting from
    // multiplying the sum by the vector of all ones (which is what
    // the test does) involves a sum of at most three nonzero values.
    // This is a generous upper bound, since alpha = 0 in this case.
    // Multiply this by 10 to allow for a bit more error.
    const bool notTooBad = (maxDiff < Teuchos::as<ST> (30) * STS::eps ());
    TEST_EQUALITY_CONST( notTooBad, true );

    if (ptrOut != NULL) {
      delete [] ptrOut;
    }
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
  }

  out << "Test alpha = 0, beta = 0" << endl;
  {
    const ScalarType alpha = Teuchos::as<ST> (0);
    const ScalarType beta = Teuchos::as<ST> (0);
    OffsetType* ptrOut = NULL;
    OrdinalType* indOut = NULL;
    ScalarType* valOut = NULL;
    Kokkos::Raw::addSparseMatrices<OffsetType, OrdinalType, ScalarType> (
      ptrOut, indOut, valOut,
      alpha, ptr1, ind1, val1,
      beta, ptr2, ind2, val2,
      numRows);
    MT maxDiff = STM::zero ();
    try {
      ArrayView<const OffsetType> ptr1_view (const_cast<const OffsetType*> (ptr1), numRows+1);
      const OffsetType numEntries1 = ptr1_view[numRows];
      ArrayView<const OrdinalType> ind1_view (const_cast<const OrdinalType*> (ind1), numEntries1);
      ArrayView<const ScalarType> val1_view (const_cast<const ScalarType*> (val1), numEntries1);

      ArrayView<const OffsetType> ptr2_view (const_cast<const OffsetType*> (ptr2), numRows+1);
      const OffsetType numEntries2 = ptr2_view[numRows];
      ArrayView<const OrdinalType> ind2_view (const_cast<const OrdinalType*> (ind2), numEntries2);
      ArrayView<const ScalarType> val2_view (const_cast<const ScalarType*> (val2), numEntries2);

      ArrayView<const OffsetType> ptrOut_view (const_cast<const OffsetType*> (ptrOut), numRows+1);
      const OffsetType numEntriesOut = ptrOut_view[numRows];
      ArrayView<const OrdinalType> indOut_view (const_cast<const OrdinalType*> (indOut), numEntriesOut);
      ArrayView<const ScalarType> valOut_view (const_cast<const ScalarType*> (valOut), numEntriesOut);

      maxDiff = testAddSparseMatrices<OffsetType, OrdinalType, ST> (
        ptrOut_view, indOut_view, valOut_view,
        alpha, ptr1_view, ind1_view, val1_view,
        beta, ptr2_view, ind2_view, val2_view,
        numRows, numCols, out);
    } catch (...) {
      if (ptrOut != NULL) {
        delete [] ptrOut;
      }
      if (indOut != NULL) {
        delete [] indOut;
      }
      if (valOut != NULL) {
        delete [] valOut;
      }
      throw;
    }

    out << "maxDiff: " << maxDiff << endl;

    // The error bound is an upper bound for each row.  The sum of the
    // two test matrices contains at most three entries per row.
    // Thus, computing each entry of the vector resulting from
    // multiplying the sum by the vector of all ones (which is what
    // the test does) involves a sum of at most three nonzero values.
    // This is a generous upper bound, since alpha = 0 and beta = 0 in
    // this case.  Multiply this by 10 to allow for a bit more error.
    const bool notTooBad = (maxDiff < Teuchos::as<ST> (30) * STS::eps ());
    TEST_EQUALITY_CONST( notTooBad, true );

    if (ptrOut != NULL) {
      delete [] ptrOut;
    }
    if (indOut != NULL) {
      delete [] indOut;
    }
    if (valOut != NULL) {
      delete [] valOut;
    }
  }

  out << "Freeing test matrices' memory" << endl;
  if (ptr1 != NULL) {
    delete [] ptr1;
  }
  if (ind1 != NULL) {
    delete [] ind1;
  }
  if (val1 != NULL) {
    delete [] val1;
  }
  if (ptr2 != NULL) {
    delete [] ptr2;
  }
  if (ind2 != NULL) {
    delete [] ind2;
  }
  if (val2 != NULL) {
    delete [] val2;
  }
}



