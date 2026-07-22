// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Sacado.hpp"

#if !defined(__CUDACC__)
#include "Kokkos_DynRankView_Fad.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Teuchos_Assert.hpp"

// Example to demonstrate the use of Kokkos::ViewFactory for constructing
// view's of Fad's without explicitly referencing the sacado dimension

// An example function that takes two views of various ranks and scalar types
// and produces a third allocated using the ViewFactory
template <class View1, class View2>
Kokkos::View< typename Kokkos::ViewFactory<View1,View2>::value_type*,
              typename View1::device_type >
my_func(const View1& v1, const View2& v2)
{
  typedef Kokkos::ViewFactory<View1,View2> ViewFac;
  typedef typename ViewFac::value_type value_type;
  typedef typename View1::device_type device_type;
  typedef typename View1::size_type size_type;
  typedef Kokkos::View<value_type*,device_type> ViewTmp;

  static_assert( unsigned(View1::rank) == 2, "" );
  static_assert( unsigned(View2::rank) == 1, "" );

  const size_type m = v1.extent(0);
  const size_type n = v1.extent(1);
  ViewTmp vtmp = ViewFac::template create_view<ViewTmp>(v1,v2,"tmp",m);

  Kokkos::parallel_for(m, KOKKOS_LAMBDA (const size_type i) {
    value_type t = 0;
    for (size_type j=0; j<n; ++j)
      t += v1(i,j)*v2(j);
    vtmp(i) = t;
  });

  return vtmp;
}

#if defined(HAVE_SACADO_KOKKOS)
// An example function that takes two dynamic-rank views of various ranks and
// scalar types and produces a third allocated using the ViewFactory
template <class View1, class View2>
Kokkos::DynRankView< typename Kokkos::ViewFactory<View1,View2>::value_type,
                     typename View1::device_type >
my_func_dynamic(const View1& v1, const View2& v2)
{
  typedef Kokkos::ViewFactory<View1,View2> ViewFac;
  typedef typename ViewFac::value_type value_type;
  typedef typename View1::device_type device_type;
  typedef typename View1::size_type size_type;
  typedef Kokkos::DynRankView<value_type,device_type> ViewTmp;

  TEUCHOS_ASSERT( v1.rank() == 2 );
  TEUCHOS_ASSERT( v2.rank() == 1 );

  const size_type m = v1.extent(0);
  const size_type n = v1.extent(1);
  ViewTmp vtmp = ViewFac::template create_view<ViewTmp>(v1,v2,"tmp",m);

  Kokkos::parallel_for(m, KOKKOS_LAMBDA (const size_type i) {
    value_type t = 0;
    for (size_type j=0; j<n; ++j)
      t += v1(i,j)*v2(j);
    vtmp(i) = t;
  });

  return vtmp;
}
#endif

#endif

int main(int argc, char* argv[]) {

#if !defined(__CUDACC__)
  typedef Sacado::Fad::DFad<double> FadType;
  typedef Kokkos::DefaultExecutionSpace execution_space;

  Kokkos::initialize(argc, argv);

  const size_t m = 10;
  const size_t n = 2;
  const size_t deriv_dim = 1;

  // First use the static-rank view
  {

    // Calculation with double's
    Kokkos::View<double**,execution_space> v1("v1",m,n);
    Kokkos::View<double* ,execution_space> v2("v2",n);

    Kokkos::deep_copy(v1, 2.0);
    Kokkos::deep_copy(v2, 3.0);

    auto v3 = my_func(v1,v2);

    std::cout << "v3 = " << std::endl;
    for (size_t i=0; i<m; ++i) {
      std::cout << "\t" << v3(i) << std::endl;
    }

    // Calculation with Fad's (initialize each component of v2_fad with a
    // Fad object with value == 3.0, one derivative == 1
    Kokkos::View<FadType*,execution_space> v2_fad("v2_fad",n,deriv_dim+1);
    Kokkos::deep_copy(v2_fad, FadType(deriv_dim, 0, 3.0));

    auto v3_fad = my_func(v1,v2_fad);

    std::cout << "v3_fad = " << std::endl;
    for (size_t i=0; i<m; ++i) {
      std::cout << "\t" << v3_fad(i) << std::endl;
    }

  }

#if defined(HAVE_SACADO_KOKKOS)
  // Now use the dynamic-rank view
  {

    // Calculation with double's
    Kokkos::DynRankView<double,execution_space> v1("v1",m,n);
    Kokkos::DynRankView<double,execution_space> v2("v2",n);

    Kokkos::deep_copy(v1, 2.0);
    Kokkos::deep_copy(v2, 3.0);

    auto v3 = my_func_dynamic(v1,v2);

    std::cout << "v3 = " << std::endl;
    for (size_t i=0; i<m; ++i) {
      std::cout << "\t" << v3(i) << std::endl;
    }

    // Calculation with Fad's (initialize each component of v2_fad with a
    // Fad object with value == 3.0, one derivative == 1
    Kokkos::DynRankView<FadType,execution_space> v2_fad("v2_fad",n,deriv_dim+1);
    Kokkos::deep_copy(v2_fad, FadType(deriv_dim, 0, 3.0));

    auto v3_fad = my_func_dynamic(v1,v2_fad);

    std::cout << "v3_fad = " << std::endl;
    for (size_t i=0; i<m; ++i) {
      std::cout << "\t" << v3_fad(i) << std::endl;
    }

  }
#endif

  Kokkos::finalize();
#endif

  return 0;
}
