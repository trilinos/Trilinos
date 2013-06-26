/*
//@HEADER
// ************************************************************************
//
//                             KokkosArray
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

namespace {

  // Test functor, that fills the given 1-D vector in parallel.  The
  // functor is like APL's iota, except that it starts with an initial
  // constant instead of with zero.  That makes it more "interesting"
  // so that tests are more likely to catch bugs.
  template<class Device>
  class FillFunctor {
  public:
    // Functors need a device_type typedef.
    typedef Device device_type;
    
    // Constructor accepts a View of a 1-D array.
    FillFunctor (const KokkosArray::View<double*, Device>& vec) : vec_ (vec) {}

    // Initialize the array.
    KOKKOSARRAY_INLINE_FUNCTION
    void operator () (int i) const {
      vec_[i] = 42.0 + static_cast<double> (i);
    }

  private:
    KokkosArray::View<double*, Device> vec_;
  };

} // namespace (anonymous)


// Just test whether Teuchos memory management objects and Kokkos
// Array Views can coexist in the same program.  This test does not
// have the Teuchos and Kokkos objects interact with each other.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkosArray, NoInteraction ) {
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

  typedef KokkosArray::View<double*, KokkosArray::Host> ka_view_type;
  ka_view_type y ("y", numElts);
  KokkosArray::parallel_for (y.dimension_0 (), FillFunctor<KokkosArray::Host> (y));
}


// First test of interaction between Teuchos and Kokkos objects:
//
// Get a Teuchos::ArrayView of a KokkosArray::View, 
// and make sure that it points to the same data.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkosArray, ArrayViewOfView ) {
  typedef Teuchos::Array<double>::size_type size_type;
  typedef KokkosArray::View<double*, KokkosArray::LayoutLeft, KokkosArray::Host> ka_view_type;

  const size_type numElts = 10;
  ka_view_type y ("y", numElts);
  KokkosArray::parallel_for (y.dimension_0 (), FillFunctor<KokkosArray::Host> (y));

  // It's possible to get the View's raw pointer because we know its
  // layout.  Not every kind of View necessarily implements the
  // ptr_on_device() method, but certainly Views in Host memory with
  // left or right (Fortran or C) layout implement this method.
  double* const y_raw = y.ptr_on_device ();
  const size_type y_size = static_cast<size_type> (y.dimension_0 ());

  Teuchos::ArrayView<double> y_view (y_raw, y_size);
  TEST_EQUALITY_CONST( y_view.size(), y_size );
  for (size_type k = 0; k < y_size; ++k) {
    TEST_EQUALITY( y_view[k], y[k] );
  }
}


// Second test of interaction between Teuchos and Kokkos objects:
//
// Get a KokkosArray::View of a Teuchos::ArrayView, and make sure that
// it points to the same data.  Thanks to Christian Trott for
// implementing the necessary functionality (View constructor for
// certain View specializations, that takes a raw pointer and
// dimensions) in Kokkos Array.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkosArray, ViewOfArrayView ) {
  typedef Teuchos::Array<double>::size_type size_type;
  typedef KokkosArray::View<double*, KokkosArray::LayoutLeft, KokkosArray::Host, KokkosArray::MemoryUnmanaged> ka_view_type;
  typedef KokkosArray::View<const double*, KokkosArray::LayoutLeft, KokkosArray::Host, KokkosArray::MemoryUnmanaged> ka_const_view_type;

  const size_type numElts = 10;
  Teuchos::Array<double> x (numElts);
  for (size_type k = 0; k < numElts; ++k) {
    x[k] = 42.0 + static_cast<double> (k);
  }
  // You can make an (unmanaged) View of a raw array with left or
  // right (Fortran or C) layout on the Host device, just by passing
  // the array and stride(s) into the constructor.  If you want the
  // dimensions to differ from the strides, you'll have to create a
  // subview with the desired dimensions.
  ka_view_type x_view (x.getRawPtr (), x.size ());

  TEST_EQUALITY( x.size(), static_cast<size_type>(x_view.dimension_0()) );
  for (size_type k = 0; k < x.size (); ++k) {
    TEST_EQUALITY( x_view[k], x[k] );
  }

  // This ensures that conversions from double* to const double* work correctly.
  // x.getRawPtr() returns double*.
  ka_const_view_type x_view_const (x.getRawPtr (), x.size ());

  TEST_EQUALITY( x.size(), static_cast<size_type>(x_view_const.dimension_0()) );
  for (size_type k = 0; k < x.size (); ++k) {
    TEST_EQUALITY( x_view_const[k], x[k] );
  }
}


