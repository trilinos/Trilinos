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

#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Sparse_impl_MKL.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <memory>

namespace { // (anonymous)

using Teuchos::TypeNameTraits;
using std::endl;
typedef ::Kokkos::View<int*, Kokkos::DefaultExecutionSpace>::HostMirror::device_type host_device_type;

template<class ScalarType, class LocalOrdinalType, class OffsetType>
struct TestMklSparseMatrixHandle {
  typedef typename Kokkos::Details::ArithTraits<ScalarType>::val_type val_type;
  typedef typename Kokkos::Details::ArithTraits<val_type>::mag_type mag_type;
  typedef KokkosSparse::CrsMatrix<val_type, LocalOrdinalType, host_device_type,
                                  void, OffsetType> matrix_type;
  typedef ::KokkosSparse::Impl::Mkl::WrappedTplMatrixHandle<matrix_type> handle_type;

  static void
  makeTestMatrixArrays (bool& success,
                        Teuchos::FancyOStream& out,
                        Kokkos::View<OffsetType*, host_device_type>& ptr,
                        Kokkos::View<LocalOrdinalType*, host_device_type>& ind,
                        Kokkos::View<val_type*, host_device_type>& val,
                        const LocalOrdinalType numRows,
                        const LocalOrdinalType numCols)
  {
    out << "Make test matrix:" << endl;
    Teuchos::OSTab tab0 (out);
    out << "ScalarType: " << TypeNameTraits<ScalarType>::name () << endl
        << "LocalOrdinalType: " << TypeNameTraits<LocalOrdinalType>::name () << endl
        << "OffsetType: " << TypeNameTraits<OffsetType>::name () << endl;

    const OffsetType numEnt = (numRows <= 2) ? (2*numRows) : (3*(numRows - 2) + 4);
    out << "numRows: " << numRows << ", numEnt: " << numEnt << endl;

    ptr = Kokkos::View<OffsetType*, host_device_type> ("ptr", numRows + 1);
    ind = Kokkos::View<LocalOrdinalType*, host_device_type> ("ind", numEnt);
    val = Kokkos::View<val_type*, host_device_type> ("val", numEnt);

    OffsetType curPos = 0;
    ptr[0] = curPos;
    for (LocalOrdinalType row = 0; row < numRows; ++row) {
      const LocalOrdinalType col0 = (row - 2) % numCols;
      const LocalOrdinalType col1 = row % numCols;
      const LocalOrdinalType col2 = (row + 2) % numCols;
      const val_type val0 = static_cast<val_type> (static_cast<mag_type> (col0));
      const val_type val1 = static_cast<val_type> (static_cast<mag_type> (col1));
      const val_type val2 = static_cast<val_type> (static_cast<mag_type> (col2));

      //out << " row: " << row << endl;

      if (row == 0) { // 2 entries
        //out << "   - curPos: " << curPos << endl;
        ind[curPos] = col1;
        val[curPos] = val1;
        ++curPos;
        //out << "   - curPos: " << curPos << endl;
        ind[curPos] = col2;
        val[curPos] = val2;
        ++curPos;
      }
      else if (row + 1 == numRows) { // 2 entries
        //out << "   - curPos: " << curPos << endl;
        ind[curPos] = col0;
        val[curPos] = val0;
        ++curPos;
        //out << "   - curPos: " << curPos << endl;
        ind[curPos] = col1;
        val[curPos] = val1;
        ++curPos;
      }
      else { // 3 entries
        //out << "   - curPos: " << curPos << endl;
        ind[curPos] = col0;
        val[curPos] = val0;
        ++curPos;
        //out << "   - curPos: " << curPos << endl;
        ind[curPos] = col1;
        val[curPos] = val1;
        ++curPos;
        //out << "   - curPos: " << curPos << endl;
        ind[curPos] = col2;
        val[curPos] = val2;
        ++curPos;
      }
      ptr[row+1] = curPos;
    }
    out << "Done!" << endl;
  }

  // Create a test matrix, and attempt to wrap it in an MKL TPL handle.
  static std::shared_ptr<handle_type>
  makeHandle (bool& success,
              Teuchos::FancyOStream& out,
              const LocalOrdinalType numRows,
              const LocalOrdinalType numCols)
  {
    Kokkos::View<OffsetType*, host_device_type> ptr;
    Kokkos::View<LocalOrdinalType*, host_device_type> ind;
    Kokkos::View<val_type*, host_device_type> val;

    makeTestMatrixArrays (success, out, ptr, ind, val, numRows, numCols);

    out << "Make KokkosSparse::CrsMatrix" << endl;
    matrix_type A ("A", numRows, numCols, val.dimension_0 (), val, ptr, ind);

    out << "Attempt to make MKL handle" << endl;
    std::shared_ptr<handle_type> handle;
#ifdef HAVE_TPETRAKERNELS_MKL
    TEST_NOTHROW( handle = std::shared_ptr<handle_type> (new handle_type (A, false)) );
    TEST_ASSERT( handle.get () != NULL );
#else
    TEST_THROW( handle = std::shared_ptr<handle_type> (new handle_type (A, false)), std::runtime_error );
    TEST_ASSERT( handle.get () == NULL );
#endif // HAVE_TPETRAKERNELS_MKL
    return handle;
  }
};

template<class ScalarType, class LocalOrdinalType, class OffsetType>
void
testMklSparseMatrixHandleOneCase (bool& success,
                                  Teuchos::FancyOStream& out,
                                  const LocalOrdinalType numRows,
                                  const LocalOrdinalType numCols)
{
  typedef TestMklSparseMatrixHandle<ScalarType, LocalOrdinalType, OffsetType> tester_type;
  (void) tester_type::makeHandle (success, out, numRows, numCols);
}

template<class LocalOrdinalType, class OffsetType>
void
testAllScalars (bool& success,
                Teuchos::FancyOStream& out,
                const LocalOrdinalType numRows,
                const LocalOrdinalType numCols)
{
  {
    typedef double scalar_type;
    out << "Test ScalarType=double" << endl;
    Teuchos::OSTab tab0 (out);
    testMklSparseMatrixHandleOneCase<scalar_type, LocalOrdinalType, OffsetType> (success, out, numRows, numCols);
  }
  {
    typedef float scalar_type;
    out << "Test ScalarType=float" << endl;
    Teuchos::OSTab tab0 (out);
    testMklSparseMatrixHandleOneCase<scalar_type, LocalOrdinalType, OffsetType> (success, out, numRows, numCols);
  }
  {
    typedef std::complex<double> scalar_type;
    out << "Test ScalarType=std::complex<double>" << endl;
    Teuchos::OSTab tab0 (out);
    testMklSparseMatrixHandleOneCase<scalar_type, LocalOrdinalType, OffsetType> (success, out, numRows, numCols);
  }
  {
    typedef std::complex<float> scalar_type;
    out << "Test ScalarType=std::complex<float>" << endl;
    Teuchos::OSTab tab0 (out);
    testMklSparseMatrixHandleOneCase<scalar_type, LocalOrdinalType, OffsetType> (success, out, numRows, numCols);
  }
}

template<class OffsetType>
void
testAllScalarsAndLocalOrdinals (bool& success,
                                Teuchos::FancyOStream& out)
{
  {
    typedef int LO;
    out << "Test LocalOrdinalType=int" << endl;
    Teuchos::OSTab tab0 (out);
    {
      const LO numRows = 30;
      const LO numCols = 15;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
    {
      const LO numRows = 1;
      const LO numCols = 3;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
    {
      const LO numRows = 2;
      const LO numCols = 3;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
  }

  {
    typedef long long LO;
    out << "Test LocalOrdinalType=long long" << endl;
    Teuchos::OSTab tab0 (out);
    {
      const LO numRows = 30;
      const LO numCols = 15;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
    {
      const LO numRows = 1;
      const LO numCols = 3;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
    {
      const LO numRows = 2;
      const LO numCols = 3;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
  }
}

void
testEverything (bool& success,
                Teuchos::FancyOStream& out)
{
  {
    typedef int offset_type;
    out << "Test OffsetType=int" << endl;
    Teuchos::OSTab tab0 (out);
    testAllScalarsAndLocalOrdinals<offset_type> (success, out);
  }
  {
    typedef size_t offset_type;
    out << "Test OffsetType=size_t" << endl;
    Teuchos::OSTab tab0 (out);
    testAllScalarsAndLocalOrdinals<offset_type> (success, out);
  }
  {
    typedef ptrdiff_t offset_type;
    out << "Test OffsetType=ptrdiff_t" << endl;
    Teuchos::OSTab tab0 (out);
    testAllScalarsAndLocalOrdinals<offset_type> (success, out);
  }
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using std::endl;

  Teuchos::RCP<Teuchos::FancyOStream> outPtr =
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  Teuchos::FancyOStream& out = *outPtr;

  out << "Call Kokkos::initialize" << endl;
  Kokkos::initialize (argc, argv);

  bool success = true;
  testEverything (success, out);

  out << "Call Kokkos::finalize" << endl;
  Kokkos::finalize ();

  if (success) {
    out << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}




