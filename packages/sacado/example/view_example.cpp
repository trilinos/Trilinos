// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Sacado.hpp"

//
// A simple example showing how to use Sacado with Kokkos
//
// The code would be simpler using lambda's instead of functors, but that
// can be problematic for Cuda and for versions of g++ that claim to, but
// don't really support C++11
//

// This code relies on the Kokkos view specializations being enabled, which
// is the default, but can be disabled by the user
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

// Functor to compute matrix-vector product c = A*b using Kokkos
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
class MatVec {
  const ViewTypeA A;
  const ViewTypeB b;
  const ViewTypeC c;
  const size_t m, n;

public:

  MatVec(const ViewTypeA& A_, const ViewTypeB& b_, const ViewTypeC& c_) :
    A(A_), b(b_), c(c_),
    m(A.extent(0)), n(A.extent(1))
  {
    typedef typename ViewTypeC::execution_space execution_space;
    Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,m), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t i) const {
    typedef typename ViewTypeC::value_type scalar_type;

    scalar_type t = 0.0;
    for (size_t j=0; j<n; ++j) {
      scalar_type bb = b(j); // fix for Intel 17.0.1 with optimization
      t += A(i,j)*bb;
    }
    c(i) = t;
  }
};

// Function to run mat-vec functor above
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void run_mat_vec(const ViewTypeA& A, const ViewTypeB& b,const ViewTypeC& c)
{
  MatVec<ViewTypeA,ViewTypeB,ViewTypeC> f(A,b,c);
}

// Functor to compute the derivative of the matrix-vector product
// c = A*b using Kokkos.
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
class MatVecDeriv {
  const ViewTypeA A;
  const ViewTypeB b;
  const ViewTypeC c;
  const size_t m, n, p;

public:

  MatVecDeriv(const ViewTypeA& A_, const ViewTypeB& b_, const ViewTypeC& c_) :
    A(A_), b(b_), c(c_),
    m(A.extent(0)), n(A.extent(1)), p(A.extent(2)-1)
  {
    typedef typename ViewTypeC::execution_space execution_space;
    Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0,m), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t i) const {
    typedef typename ViewTypeC::value_type scalar_type;

    // Derivative portion
    for (size_t k=0; k<p; ++k) {
      scalar_type t = 0.0;
      for (size_t j=0; j<n; ++j)
        t += A(i,j,k)*b(j,p) + A(i,j,p)*b(j,k);
      c(i,k) = t;
    }

    // Value portion
    scalar_type t = 0.0;
    for (size_t j=0; j<n; ++j)
      t += A(i,j,p)*b(j,p);
    c(i,p) = t;
  }
};

// Function to run mat-vec derivative functor above
template <typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
void run_mat_vec_deriv(const ViewTypeA& A, const ViewTypeB& b,
                       const ViewTypeC& c)
{
  MatVecDeriv<ViewTypeA,ViewTypeB,ViewTypeC> f(A,b,c);
}

#endif

int main(int argc, char* argv[]) {
  int ret = 0;

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

  Kokkos::initialize(argc, argv);
  {

    const size_t m = 10; // Number of rows
    const size_t n = 2;  // Number of columns
    const size_t p = 1;  // Derivative dimension

    // Allocate Kokkos view's for matrix, input vector, and output vector
    // Derivative dimension must be specified as the last constructor argument
    typedef Sacado::Fad::SFad<double,p> FadType;
    Kokkos::View<FadType**> A("A",m,n,p+1);
    Kokkos::View<FadType*>  b("b",n,p+1);
    Kokkos::View<FadType*>  c("c",m,p+1);

    // Initialize values
    Kokkos::deep_copy( A, FadType(p, 0, 2.0) );
    Kokkos::deep_copy( b, FadType(p, 0, 3.0) );
    Kokkos::deep_copy( c, 0.0 );

    // Run mat-vec
    run_mat_vec(A,b,c);

    // Print result
    std::cout << "\nc = A*b:  Differentiated using Sacado:" << std::endl;
    auto h_c = Kokkos::create_mirror_view(c);
    Kokkos::deep_copy(h_c, c);
    for (size_t i=0; i<m; ++i)
      std::cout << "\tc(" << i << ") = " << h_c(i) << std::endl;

    // Now compute derivative analytically.  Any Sacado view can be flattened
    // into a standard view of one higher rank, with the extra dimension equal
    // to p+1
    Kokkos::View<FadType*> c2("c",m,p+1);
    Kokkos::View<double***> A_flat = A;
    Kokkos::View<double**>  b_flat = b;
    Kokkos::View<double**>  c_flat = c2;
    run_mat_vec_deriv(A_flat, b_flat, c_flat);

    // Print result
    std::cout << "\nc = A*b:  Differentiated analytically:" << std::endl;
    auto h_c2 = Kokkos::create_mirror_view(c2);
    Kokkos::deep_copy(h_c2, c2);
    for (size_t i=0; i<m; ++i)
      std::cout << "\tc(" << i << ") = " << h_c2(i) << std::endl;

    // Compute the error
    double err = 0.0;
    for (size_t i=0; i<m; ++i) {
      for (size_t k=0; k<p; ++k)
        err += std::abs(h_c(i).fastAccessDx(k)-h_c2(i).fastAccessDx(k));
      err += std::abs(h_c(i).val()-h_c2(i).val());
    }

    double tol = 1.0e-14;
    if (err < tol) {
      ret = 0;
      std::cout << "\nExample passed!" << std::endl;
    }
    else {
      ret = 1;
      std::cout <<"\nSomething is wrong, example failed!" << std::endl;
    }

  }
  Kokkos::finalize();

#endif

  return ret;
}
