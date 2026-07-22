//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <KokkosCompat_TMM.hpp>
#include <Teuchos_Array.hpp>

namespace TestTeuchosKokkosCompat {

  // Test functor, that fills the given 1-D vector in parallel.  The
  // functor is like APL's iota, except that it starts with an initial
  // constant instead of with zero.  That makes it more "interesting"
  // so that tests are more likely to catch bugs.
  template<class Device>
  class FillFunctor {
  public:
    FillFunctor (const Kokkos::View<double*, Device>& vec) : vec_ (vec) {}
    KOKKOS_INLINE_FUNCTION void operator () (const int i) const {
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

} // namespace TestTeuchosKokkosCompat

// Just test whether Teuchos memory management objects and Kokkos
// Array Views can coexist in the same program.  This test does not
// have the Teuchos and Kokkos objects interact with each other.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, NoInteraction ) {
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

  typedef Kokkos::View<double*> ka_view_type;
  ka_view_type y ("y", numElts);

  using execution_space = ka_view_type::execution_space;
  Kokkos::RangePolicy<execution_space, int> range (0, y.extent (0));
  using functor_type = TestTeuchosKokkosCompat::FillFunctor<execution_space>;
  Kokkos::parallel_for ("FillFunctor", range, functor_type (y));
}


// Get a Teuchos::ArrayView of a Kokkos::View, and make sure that
// it points to the same data.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, ArrayViewOfView ) {
  typedef Teuchos::Array<double>::size_type size_type;
  typedef Kokkos::View<double*> ka_view_type;

  const size_type numElts = 10;
  ka_view_type y ("y", numElts);

  using execution_space = ka_view_type::execution_space;
  Kokkos::RangePolicy<execution_space, int> range (0, y.extent (0));
  using functor_type = TestTeuchosKokkosCompat::FillFunctor<execution_space>;
  Kokkos::parallel_for ("FillFunctor", range, functor_type (y));

  auto y_host = Kokkos::create_mirror_view (y);
  Kokkos::deep_copy (y_host, y);

  Teuchos::ArrayView<double> y_view (y_host.data (), y_host.extent (0));
  TEST_EQUALITY_CONST( size_t (y_view.size ()), size_t (y_host.extent (0)) );
  if (success) {
    for (size_type k = 0; k < size_type (y_host.extent (0)); ++k) {
      TEST_EQUALITY( y_view[k], y_host[k] );
    }
  }
}

template<class DataType, class LayoutType = Kokkos::LayoutLeft>
using ManagedHostView =
  Kokkos::View<DataType, LayoutType,
               Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace> >;

template<class DataType, class LayoutType = Kokkos::LayoutLeft>
using UnmanagedHostView =
  Kokkos::View<DataType, LayoutType,
               Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>,
               Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

// Get a Kokkos::View of a Teuchos::ArrayView, and make sure that it
// points to the same data.  Thanks to Christian Trott for
// implementing the necessary functionality (View constructor for
// certain View specializations, that takes a raw pointer and
// dimensions) in Kokkos Array.
//
// This example will be useful for implementing the
// Tpetra::MultiVector methods get1dCopy and get2dCopy.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, ViewOfArrayView ) {
  using size_type = Teuchos::Array<double>::size_type;

  const size_type numElts = 10;
  Teuchos::Array<double> x (numElts);
  for (size_type k = 0; k < numElts; ++k) {
    x[k] = 42.0 + static_cast<double> (k);
  }
  // You can make an (unmanaged) View of a raw array with left or
  // right (Fortran or C) layout on the HostSpace device, just by
  // passing the array and stride(s) into the constructor.  If you
  // want the dimensions to differ from the strides, you'll have to
  // create a subview with the desired dimensions.
  UnmanagedHostView<double*> x_view (x.getRawPtr (), x.size ());

  TEST_EQUALITY( x.size(), static_cast<size_type> (x_view.extent (0)) );
  for (size_type k = 0; k < x.size (); ++k) {
    TEST_EQUALITY( x_view(k), x[k] );
  }

  UnmanagedHostView<const double*> x_view_const (x.getRawPtr (), x.size ());

  TEST_EQUALITY( x.size (), static_cast<size_type> (x_view_const.extent (0)) );
  for (size_type k = 0; k < x.size (); ++k) {
    TEST_EQUALITY( x_view_const(k), x[k] );
  }
}

// Create a 2-D View (LayoutLeft, in HostSpace memory), and then create an
// ArrayRCP which owns that View (using a custom destructor).  This
// will be useful for implementing Tpetra::MultiVector's getData,
// getDataNonConst, get1dView, and get1dViewNonConst methods.
TEUCHOS_UNIT_TEST( LinkTeuchosAndKokkos, ArrayRCP1D_of_2DView ) {
  using ka_view_type = ManagedHostView<double**>;

  const size_t numRows = 75;
  const size_t numCols = 5;
  const size_t stride = 100;
  const size_t ZERO = static_cast<size_t> (0);

  ka_view_type X ("X", stride, numCols);
  ka_view_type X_view = Kokkos::subview (X, std::make_pair (ZERO, numRows),
                                         std::make_pair (ZERO, numCols));
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

  using dealloc_type = TestTeuchosKokkosCompat::Deallocator<double, ka_view_type>;
  Teuchos::ArrayRCP<double> Y_values (X.data (), 0, stride*numCols,
                                      dealloc_type (X), true);
  TEST_EQUALITY(Y_values.getRawPtr(), X.data());
  TEST_EQUALITY(Y_values.getRawPtr(), X_view.data());
  TEST_EQUALITY(Y_values.size(), stride*numCols);
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
  using ka_view_type = ManagedHostView<double**>;

  const size_t numRows = 75;
  const size_t numCols = 5;
  const size_t stride = 100;
  const size_t ZERO = static_cast<size_t> (0);

  ka_view_type X ("X", stride, numCols);
  ka_view_type X_view = Kokkos::subview (X, std::make_pair (ZERO, numRows),
                                         std::make_pair (ZERO, numCols));
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
    ka_view_type X_j = Kokkos::subview (X_view, std::make_pair (ZERO, numRows),
                                        std::make_pair (j, j+1));
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
    using dealloc_type = TestTeuchosKokkosCompat::Deallocator<double, ka_view_type>;
    Teuchos::ArrayRCP<double> Y_j (X_j.data (), 0, numRows,
                                   dealloc_type (X_j), true);
    Y_2D[j] = Y_j;
  }
}


int
main (int argc, char* argv[])
{
  // Initialize MPI (if enabled) before initializing Kokkos.  This
  // lets MPI control things like pinning processes to sockets.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
