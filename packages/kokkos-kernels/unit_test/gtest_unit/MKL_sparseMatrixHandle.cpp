/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

//#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Sparse_impl_MKL.hpp"
//#include "Teuchos_TypeNameTraits.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <memory>

#include "KokkosKernels_Test_Macros.hpp"
namespace { // (anonymous)

//using Teuchos::TypeNameTraits;
using std::endl;
typedef ::Kokkos::View<int*>::HostMirror::execution_space host_execution_space;
typedef ::Kokkos::Device<host_execution_space, Kokkos::HostSpace> host_device_type;

template<class ScalarType, class LocalOrdinalType, class OffsetType>
struct TestMklSparseMatrixHandle {
  typedef typename Kokkos::Details::ArithTraits<ScalarType>::val_type val_type;
  typedef typename Kokkos::Details::ArithTraits<val_type>::mag_type mag_type;
  typedef KokkosSparse::CrsMatrix<val_type, LocalOrdinalType, host_device_type,
                                  void, OffsetType> matrix_type;
  typedef ::KokkosSparse::Impl::Mkl::WrappedTplMatrixHandle<matrix_type> handle_type;

  static void
  makeTestMatrixArrays (bool& success,
                        //Teuchos::FancyOStream& out,
                        std::ostream &out,
                        Kokkos::View<OffsetType*, host_device_type>& ptr,
                        Kokkos::View<LocalOrdinalType*, host_device_type>& ind,
                        Kokkos::View<val_type*, host_device_type>& val,
                        const LocalOrdinalType numRows,
                        const LocalOrdinalType numCols)
  {
    out << "Make test matrix:" << endl;
    //Teuchos::OSTab tab0 (out);
    out << "ScalarType: " << KokkosKernels::TypeNameTraits<ScalarType>::name () << endl
        << "LocalOrdinalType: " << KokkosKernels::TypeNameTraits<LocalOrdinalType>::name () << endl
        << "OffsetType: " << KokkosKernels::TypeNameTraits<OffsetType>::name () << endl;

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
              //Teuchos::FancyOStream& out,
			  std::ostream &out,
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
#ifdef HAVE_KOKKOSKERNELS_MKL
    //TEST_NOTHROW( handle = std::shared_ptr<handle_type> (new handle_type (A, false)) );
    KK_TEST_NOTHROW( handle = std::shared_ptr<handle_type> (new handle_type (A, false)), out, success );
    const bool l_result = handle.get () != NULL;
    //TEST_ASSERT( handle.get () != NULL );
    KK_TEST_ASSERT( l_result );
#else
    //TEST_THROW( handle = std::shared_ptr<handle_type> (new handle_type (A, false)), std::runtime_error );
    KK_TEST_THROW(  handle = std::shared_ptr<handle_type> (new handle_type (A, false)), std::runtime_error , out, success);
    //TEST_ASSERT( handle.get () == NULL );
    KK_TEST_ASSERT( handle.get () != NULL );
#endif // HAVE_KOKKOSKERNELS_MKL
    return handle;
  }
};

template<class ScalarType, class LocalOrdinalType, class OffsetType>
void
testMklSparseMatrixHandleOneCase (bool& success,
                                  //Teuchos::FancyOStream& out,
								  std::ostream &out,
                                  const LocalOrdinalType numRows,
                                  const LocalOrdinalType numCols)
{
  typedef TestMklSparseMatrixHandle<ScalarType, LocalOrdinalType, OffsetType> tester_type;
  (void) tester_type::makeHandle (success, out, numRows, numCols);
}

template<class LocalOrdinalType, class OffsetType>
void
testAllScalars (bool& success,
                //Teuchos::FancyOStream& out,
				std::ostream &out,
                const LocalOrdinalType numRows,
                const LocalOrdinalType numCols)
{
  {
    typedef double scalar_type;
    out << "Test ScalarType=double" << endl;
    //Teuchos::OSTab tab0 (out);
    testMklSparseMatrixHandleOneCase<scalar_type, LocalOrdinalType, OffsetType> (success, out, numRows, numCols);

    out << "Test conversion between our value_type and internal_value_type" << endl;

    typedef Kokkos::Details::ArithTraits<scalar_type>::val_type value_type;
    typedef ::KokkosSparse::Impl::Mkl::RawTplMatrixHandle<value_type> converter_type;
    static_assert (std::is_same<converter_type::value_type, value_type>::value,
                   "RawTplMatrixHandle<double>::value_type != double");
#ifdef HAVE_KOKKOSKERNELS_MKL
    static_assert (std::is_same<converter_type::internal_value_type, double>::value,
                   "RawTplMatrixHandle<double>::interval_value_type != double");
#endif // HAVE_KOKKOSKERNELS_MKL

    const value_type x_our (3.0);
    const auto x_mkl = converter_type::convertToInternalValue (x_our);
    static_assert (std::is_same<std::decay<decltype (x_mkl) >::type, converter_type::internal_value_type>::value,
                   "RawTplMatrixHandle<double>::convertToInternalValue returns the wrong type");
    const auto x_back = converter_type::convertFromInternalValue (x_mkl);
    static_assert (std::is_same<std::decay<decltype (x_back) >::type, converter_type::value_type>::value,
                   "RawTplMatrixHandle<double>::convertFromInternalValue returns the wrong type");
    //TEST_EQUALITY( x_back, x_our );
    EXPECT_TRUE( (x_back == x_our ));

  }
  {
    typedef float scalar_type;
    out << "Test ScalarType=float" << endl;
    //Teuchos::OSTab tab0 (out);
    testMklSparseMatrixHandleOneCase<scalar_type, LocalOrdinalType, OffsetType> (success, out, numRows, numCols);

    typedef Kokkos::Details::ArithTraits<scalar_type>::val_type value_type;
    typedef ::KokkosSparse::Impl::Mkl::RawTplMatrixHandle<value_type> converter_type;
    static_assert (std::is_same<converter_type::value_type, value_type>::value,
                   "RawTplMatrixHandle<float>::value_type != float");
#ifdef HAVE_KOKKOSKERNELS_MKL
    static_assert (std::is_same<converter_type::internal_value_type, float>::value,
                   "RawTplMatrixHandle<float>::interval_value_type != float");
#endif // HAVE_KOKKOSKERNELS_MKL

    const value_type x_our (3.0);
    const auto x_mkl = converter_type::convertToInternalValue (x_our);
    static_assert (std::is_same<std::decay<decltype (x_mkl) >::type, converter_type::internal_value_type>::value,
                   "RawTplMatrixHandle<float>::convertToInternalValue returns the wrong type");
    const auto x_back = converter_type::convertFromInternalValue (x_mkl);
    static_assert (std::is_same<std::decay<decltype (x_back) >::type, converter_type::value_type>::value,
                   "RawTplMatrixHandle<float>::convertFromInternalValue returns the wrong type");
    EXPECT_TRUE( (x_back == x_our ));
    //TEST_EQUALITY( x_back, x_our );
  }
  {
    typedef std::complex<double> scalar_type;
    out << "Test ScalarType=std::complex<double>" << endl;
    //Teuchos::OSTab tab0 (out);
    testMklSparseMatrixHandleOneCase<scalar_type, LocalOrdinalType, OffsetType> (success, out, numRows, numCols);

    out << "Test conversion between our value_type and internal_value_type" << endl;

    typedef Kokkos::Details::ArithTraits<scalar_type>::val_type value_type;
    typedef ::KokkosSparse::Impl::Mkl::RawTplMatrixHandle<value_type> converter_type;
    static_assert (std::is_same<converter_type::value_type, value_type>::value,
                   "RawTplMatrixHandle<Kokkos::complex<double> >::value_type != Kokkos::complex<double>");
#ifdef HAVE_KOKKOSKERNELS_MKL
    static_assert (std::is_same<converter_type::internal_value_type, MKL_Complex16>::value,
                   "RawTplMatrixHandle<Kokkos::complex<double> >::interval_value_type != MKL_Complex16");
#endif // HAVE_KOKKOSKERNELS_MKL

    const value_type x_our (3.0, -4.0);
    const auto x_mkl = converter_type::convertToInternalValue (x_our);
    static_assert (std::is_same<std::decay<decltype (x_mkl) >::type, converter_type::internal_value_type>::value,
                   "RawTplMatrixHandle<Kokkos::complex<double> >::convertToInternalValue returns the wrong type");
    const auto x_back = converter_type::convertFromInternalValue (x_mkl);
    static_assert (std::is_same<std::decay<decltype (x_back) >::type, converter_type::value_type>::value,
                   "RawTplMatrixHandle<Kokkos::complex<double> >::convertFromInternalValue returns the wrong type");
    EXPECT_TRUE( (x_back == x_our ));
    //TEST_EQUALITY( x_back, x_our );
  }
  {
    typedef std::complex<float> scalar_type;
    out << "Test ScalarType=std::complex<float>" << endl;
    //Teuchos::OSTab tab0 (out);
    testMklSparseMatrixHandleOneCase<scalar_type, LocalOrdinalType, OffsetType> (success, out, numRows, numCols);

    out << "Test conversion between our value_type and internal_value_type" << endl;

    typedef Kokkos::Details::ArithTraits<scalar_type>::val_type value_type;
    typedef ::KokkosSparse::Impl::Mkl::RawTplMatrixHandle<value_type> converter_type;
    static_assert (std::is_same<converter_type::value_type, value_type>::value,
                   "RawTplMatrixHandle<Kokkos::complex<float> >::value_type != Kokkos::complex<float>");
#ifdef HAVE_KOKKOSKERNELS_MKL
    static_assert (std::is_same<converter_type::internal_value_type, MKL_Complex8>::value,
                   "RawTplMatrixHandle<Kokkos::complex<float> >::interval_value_type != MKL_Complex8");
#endif // HAVE_KOKKOSKERNELS_MKL

    const value_type x_our (3.0, -4.0);
    const auto x_mkl = converter_type::convertToInternalValue (x_our);
    static_assert (std::is_same<std::decay<decltype (x_mkl) >::type, converter_type::internal_value_type>::value,
                   "RawTplMatrixHandle<Kokkos::complex<float> >::convertToInternalValue returns the wrong type");
    const auto x_back = converter_type::convertFromInternalValue (x_mkl);
    static_assert (std::is_same<std::decay<decltype (x_back) >::type, converter_type::value_type>::value,
                   "RawTplMatrixHandle<Kokkos::complex<float> >::convertFromInternalValue returns the wrong type");
    //TEST_EQUALITY( x_back, x_our );
    EXPECT_TRUE( (x_back == x_our ));

  }
}

template<class OffsetType>
void
testAllScalarsAndLocalOrdinals (bool& success,std::ostream &out)
                                //Teuchos::FancyOStream& out)
{
  {
    typedef int LO;
    out << "Test LocalOrdinalType=int" << endl;
    //Teuchos::OSTab tab0 (out);
    {
      const LO numRows = 30;
      const LO numCols = 15;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      //Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
    {
      const LO numRows = 1;
      const LO numCols = 3;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      //Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
    {
      const LO numRows = 2;
      const LO numCols = 3;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      //Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
  }

  {
    typedef long long LO;
    out << "Test LocalOrdinalType=long long" << endl;
    //Teuchos::OSTab tab0 (out);
    {
      const LO numRows = 30;
      const LO numCols = 15;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      //Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
    {
      const LO numRows = 1;
      const LO numCols = 3;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      //Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
    {
      const LO numRows = 2;
      const LO numCols = 3;
      out << "Test numRows=" << numRows << ", numCols=" << numCols << endl;
      //Teuchos::OSTab tab1 (out);
      testAllScalars<LO, OffsetType> (success, out, numRows, numCols);
    }
  }
}

void
testEverything (bool& success, std::ostream &out)
                //Teuchos::FancyOStream& out)
{
  {
    typedef int offset_type;
    out << "Test OffsetType=int" << endl;
    //Teuchos::OSTab tab0 (out);
    testAllScalarsAndLocalOrdinals<offset_type> (success, out);
  }
  {
    typedef size_t offset_type;
    out << "Test OffsetType=size_t" << endl;
    //Teuchos::OSTab tab0 (out);
    testAllScalarsAndLocalOrdinals<offset_type> (success, out);
  }
  {
    typedef ptrdiff_t offset_type;
    out << "Test OffsetType=ptrdiff_t" << endl;
    //Teuchos::OSTab tab0 (out);
    testAllScalarsAndLocalOrdinals<offset_type> (success, out);
  }
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using std::endl;
  std::ostream &out = std::cout;

  /*
  Teuchos::RCP<Teuchos::FancyOStream> outPtr =
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  Teuchos::FancyOStream& out = *outPtr;
  */

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




