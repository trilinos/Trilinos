/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Teuchos_UnitTestHarness.hpp>
#include <KokkosCompat_TMM.hpp>
#include <Teuchos_Array.hpp>

// CMakeList.txt is not set up for compiling Cuda
// so choose an appropriate non-cuda device:
#if defined( KOKKOS_ENABLE_CUDA )
typedef Kokkos::HostSpace::execution_space TestDevice ;
#else
typedef Kokkos::DefaultExecutionSpace TestDevice ;
#endif

namespace {

  // Test functor, that fills the given 1-D vector in parallel.  The
  // functor is like APL's iota, except that it starts with an initial
  // constant instead of with zero.  That makes it more "interesting"
  // so that tests are more likely to catch bugs.
  template<class Device>
  class FillFunctor {
  public:
    // Functors need a execution_space typedef.
    typedef Device execution_space;

    // Constructor accepts a View of a 1-D array.
    FillFunctor (const Kokkos::View<double*, Device>& vec) : vec_ (vec) {}

    // Initialize the array.
    KOKKOS_INLINE_FUNCTION
    void operator () (int i) const {
      vec_[i] = 42.0 + static_cast<double> (i);
    }

  private:
    Kokkos::View<double*, Device> vec_;
  };

  // Custom deallocator for Teuchos::ArrayRCP.  It doesn't actually
  // deallocate the ArrayRCP's data (which belongs to the View passed
  // into Deallocator's constructor, not to the ArrayRCP).  Instead,
  // it just keeps a reference to the View.  That way, the data won't
  // go away until the View's reference count goes to zero.
  template<class Scalar, class ViewType>
  class Deallocator {
  public:
    // Constructor takes the View that owns the memory.
    Deallocator (const ViewType& view) : view_ (view) {}

    // "Deallocation function" doesn't actually deallocate its input
    // pointer; the View is responsible for deallocation of its
    // memory.
    void free (Scalar*) {}
  private:
    ViewType view_; // View that owns the memory
  };

} // namespace (anonymous)


// Just test whether Teuchos memory management objects and Kokkos
// Array Views can coexist in the same program.  This test does not
// have the Teuchos and Kokkos objects interact with each other.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, NoInteraction ) {
  TestDevice::initialize();
  typedef Teuchos::Array<double>::size_type size_type;

  const size_type numElts = 10;
  Teuchos::Array<double> x (numElts);
  for (size_type k = 0; k < numElts; ++k) {
    x[k] = 42.0 + static_cast<double> (k);
  }
  Teuchos::ArrayView<double> xView = x.view (3, 5); // view of [3, 4, 5, 6, 7]
  TEST_EQUALITY_CONST( xView.size(), 5 );
  for (size_type k = 0; k < xView.size (); ++k) {
    TEST_EQUALITY( xView[k], x[k+3] );
  }

  typedef Kokkos::View<double*, TestDevice> ka_view_type;
  ka_view_type y ("y", numElts);
  Kokkos::parallel_for (y.extent (0), FillFunctor<TestDevice> (y));
  TestDevice::finalize();
}


// Get a Teuchos::ArrayView of a Kokkos::View, and make sure that
// it points to the same data.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, ArrayViewOfView ) {
  TestDevice::initialize();
  typedef Teuchos::Array<double>::size_type size_type;
  typedef Kokkos::View<double*, Kokkos::LayoutLeft, TestDevice> ka_view_type;

  const size_type numElts = 10;
  ka_view_type y ("y", numElts);
  Kokkos::parallel_for (y.extent (0), FillFunctor<TestDevice> (y));

  double* const y_raw = y.data ();
  const size_type y_size = static_cast<size_type> (y.extent (0));

  Teuchos::ArrayView<double> y_view (y_raw, y_size);
  TEST_EQUALITY_CONST( y_view.size(), y_size );
  for (size_type k = 0; k < y_size; ++k) {
    TEST_EQUALITY( y_view[k], y[k] );
  }
  TestDevice::finalize();
}


// Get a Kokkos::View of a Teuchos::ArrayView, and make sure that
// it points to the same data.  Thanks to Christian Trott for
// implementing the necessary functionality (View constructor for
// certain View specializations, that takes a raw pointer and
// dimensions) in Kokkos Array.
//
// This example will be useful for implementing the
// Tpetra::MultiVector methods get1dCopy and get2dCopy.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, ViewOfArrayView ) {
  TestDevice::initialize();
  typedef Teuchos::Array<double>::size_type size_type;
  typedef Kokkos::View<double*, Kokkos::LayoutLeft, TestDevice, Kokkos::MemoryUnmanaged> ka_view_type;
  typedef Kokkos::View<const double*, Kokkos::LayoutLeft, TestDevice, Kokkos::MemoryUnmanaged> ka_const_view_type;

  const size_type numElts = 10;
  Teuchos::Array<double> x (numElts);
  for (size_type k = 0; k < numElts; ++k) {
    x[k] = 42.0 + static_cast<double> (k);
  }
  // You can make an (unmanaged) View of a raw array with left or
  // right (Fortran or C) layout on the HostSpace device, just by passing
  // the array and stride(s) into the constructor.  If you want the
  // dimensions to differ from the strides, you'll have to create a
  // subview with the desired dimensions.
  ka_view_type x_view (x.getRawPtr (), x.size ());

  TEST_EQUALITY( x.size(), static_cast<size_type>(x_view.extent(0)) );
  for (size_type k = 0; k < x.size (); ++k) {
    TEST_EQUALITY( x_view[k], x[k] );
  }

  // This ensures that conversions from double* to const double* work correctly.
  // x.getRawPtr() returns double*.
  ka_const_view_type x_view_const ( (const double*) x.getRawPtr (), x.size ());

  TEST_EQUALITY( x.size(), static_cast<size_type>(x_view_const.extent(0)) );
  for (size_type k = 0; k < x.size (); ++k) {
    TEST_EQUALITY( x_view_const[k], x[k] );
  }

  TestDevice::finalize();
}


// Create a 2-D View (LayoutLeft, in HostSpace memory), and then create an
// ArrayRCP which owns that View (using a custom destructor).  This
// will be useful for implementing Tpetra::MultiVector's getData,
// getDataNonConst, get1dView, and get1dViewNonConst methods.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, ArrayRCP1D_of_2DView ) {
  TestDevice::initialize();

  typedef Kokkos::View<double**, Kokkos::LayoutLeft, TestDevice> ka_view_type;

  const size_t numRows = 75;
  const size_t numCols = 5;
  const size_t stride = 100;
  const size_t ZERO = static_cast<size_t> (0);

  ka_view_type X ("X", stride, numCols);
  ka_view_type X_view = Kokkos::subview (X, std::make_pair (ZERO, numRows), std::make_pair (ZERO, numCols));
  TEST_EQUALITY(X_view.extent(0), numRows);
  TEST_EQUALITY(X_view.extent(1), numCols);

  // Test that the strides of X_view are correct, for int.  Kokkos
  // Array templates the stride() method on the integer type of the
  // input.  I just want to make sure that this works for different
  // integer types; I'll test int here and size_t.
  {
    int strides[3] = { 0 , 0 , 0 };
    X_view.stride (strides);
    TEST_EQUALITY_CONST(strides[0], 1); // stride between X_view(i,j) and X_view(i+1,j)
    // The stride must be at least as given, but can be greater (due to possible padding for alignment).
    TEST_ASSERT(static_cast<size_t>(strides[1]) >= stride); // stride between X_view(i,j) and X_view(i,j+1)
  }

  // Test that the strides of X_view are correct, for size_t.
  {
    size_t strides[3] = { 0 , 0 , 0 };
    X_view.stride (strides);
    TEST_EQUALITY_CONST(strides[0], static_cast<size_t>(1)); // stride between X_view(i,j) and X_view(i+1,j)
    // The stride must be at least as given, but can be greater (due to possible padding for alignment).
    TEST_ASSERT(strides[1] >= stride); // stride between X_view(i,j) and X_view(i,j+1)
  }

  // Create a nonowning "view" of X's data.  The ArrayRCP has a custom
  // destructor which simply holds a reference to the View and doesn't
  // actually deallocate memory.  That way, the ArrayRCP can use the
  // View's raw pointer, but still defers to the View for memory
  // management.

  Teuchos::ArrayRCP<double> Y_values (X.data (), 0, stride*numCols, Deallocator<double, ka_view_type> (X), true);
  TEST_EQUALITY(Y_values.getRawPtr(), X.data());
  TEST_EQUALITY(Y_values.getRawPtr(), X_view.data());
  TEST_EQUALITY(Y_values.size(), stride*numCols);

  TestDevice::finalize();
}


// Create a 2-D View (LayoutLeft, in HostSpace memory), and then create an
// ArrayRCP<ArrayRCP<double> >, each entry of which owns the
// corresponding column of that View (using a custom destructor).
// This will be useful for implementing Tpetra::MultiVector's
// get2dView and get2dViewNonConst methods.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, ArrayRCP2D_of_2DView ) {
  // View<double*, LayoutLeft, ...> and View<double*, LayoutRight,
  // ...> always have unit stride.  Furthermore, subview(X,j) returns
  // a view of a column for LayoutLeft, and a view of a row for
  // LayoutRight.  Thus, to support both layouts, each column of the
  // multivector still needs to be a View<double**, ...>, and we have
  // to use the subview() overload that takes two ranges of row and
  // column indices.
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, TestDevice> ka_view_type;
  TestDevice::initialize();

  const size_t numRows = 75;
  const size_t numCols = 5;
  const size_t stride = 100;
  const size_t ZERO = static_cast<size_t> (0);

  ka_view_type X ("X", stride, numCols);
  ka_view_type X_view = Kokkos::subview (X, std::make_pair (ZERO, numRows), std::make_pair (ZERO, numCols));
  TEST_EQUALITY( & X(0,0) , & X_view(0,0) );
  TEST_EQUALITY(X_view.extent(0), numRows);
  TEST_EQUALITY(X_view.extent(1), numCols);

  // Test that the strides of X_view are correct, for size_t.
  {
    size_t strides[3]; // stride[rank] is the maximum span
    strides[0] = static_cast<size_t> (0);
    strides[1] = static_cast<size_t> (0);
    strides[2] = static_cast<size_t> (0);
    X_view.stride (strides);
    TEST_EQUALITY_CONST(strides[0], static_cast<size_t>(1)); // stride between X_view(i,j) and X_view(i+1,j)
    // The stride must be at least as given, but can be greater (due to possible padding for alignment).
    TEST_ASSERT(strides[1] >= stride); // stride between X_view(i,j) and X_view(i,j+1)
  }

  // Make a 2-D "view" (array of arrays) of X_view.  This is how we
  // will implement Tpetra::MultiVector methods like get2dView.
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > Y_2D (X_view.extent (1));
  for (size_t j = 0; j < static_cast<size_t> (X_view.extent (1)); ++j) {
    ka_view_type X_j = Kokkos::subview (X_view, std::make_pair (ZERO, numRows), std::make_pair (j, j+1));
    TEST_EQUALITY( & X_view(0,j) , & X_j(0,0) );
    TEST_EQUALITY(static_cast<size_t>(X_j.extent(0)), numRows);
    TEST_EQUALITY_CONST(X_j.extent(1), 1);

    // Test that the strides of X_j are correct.
    {
      size_t strides[3]; // stride[rank] is the maximum span
      strides[0] = static_cast<size_t> (0);
      strides[1] = static_cast<size_t> (0);
      strides[2] = static_cast<size_t> (0);
      X_j.stride (strides);
      TEST_EQUALITY_CONST(strides[0], static_cast<size_t>(1)); // stride between X_j(i,j) and X_j(i+1,j)
      // Stride between X_j(i,j) and X_j(i,j+1), even though X_j only
      // has one column.  The stride must be at least as given, but
      // can be greater (due to possible padding for alignment).
      TEST_ASSERT(strides[1] >= stride);
    }

    // Create a nonowning "view" of X_j's data.  The ArrayRCP has a
    // custom destructor which simply holds a reference to the View
    // X_j, and doesn't actually deallocate memory.  That way, the
    // ArrayRCP can use the View's raw pointer, but still defers to
    // the View for memory management.
    Teuchos::ArrayRCP<double> Y_j (X_j.data (), 0, numRows, Deallocator<double, ka_view_type> (X_j), true);
    Y_2D[j] = Y_j;
  }
  TestDevice::finalize();
}


