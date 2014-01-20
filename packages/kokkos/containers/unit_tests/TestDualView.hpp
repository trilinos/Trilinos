#ifndef KOKKOS_TEST_DUALVIEW_HPP
#define KOKKOS_TEST_DUALVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <impl/Kokkos_Timer.hpp>

namespace Test {

namespace Impl {

  template <typename Scalar, class Device>
  struct test_dualview_combinations
  {
    typedef test_dualview_combinations<Scalar,Device> self_type;

    typedef Scalar scalar_type;
    typedef Device device_type;

    Scalar reference;
    Scalar result;

    template <typename ViewType>
    Scalar run_me(unsigned int n,unsigned int m){
      if(n<10) n = 10;
      if(m<3) m = 3;
      ViewType a("A",n,m);

      Kokkos::Impl::ViewFill<typename ViewType::t_dev>(a.d_view,1);
      a.template modify<typename ViewType::device_type>();
      a.template sync<typename ViewType::host_mirror_device_type>();

      a.h_view(5,1) = 3;
      a.h_view(6,1) = 4;
      a.h_view(7,2) = 5;
      a.template modify<typename ViewType::host_mirror_device_type>();
      ViewType b = Kokkos::subview<ViewType>(a,std::pair<unsigned int, unsigned int>(6,9),std::pair<unsigned int, unsigned int>(0,1));
      a.template sync<typename ViewType::device_type>();
      b.template modify<typename ViewType::device_type>();
      Kokkos::Impl::ViewFill<typename ViewType::t_dev>(b.d_view,2);
      a.template sync<typename ViewType::host_mirror_device_type>();
      Scalar count = 0;
      for(unsigned int i = 0; i<a.d_view.dimension_0(); i++)
        for(unsigned int j = 0; j<a.d_view.dimension_1(); j++)
          count += a.h_view(i,j);
      return count -  a.d_view.dimension_0()*a.d_view.dimension_1()-2-4-3*2;
    }


    test_dualview_combinations(unsigned int size)
    {
      result = run_me< Kokkos::DualView<Scalar**,Kokkos::LayoutLeft,Device> >(size,3);
    }

   };

} // namespace Impl




template <typename Scalar, typename Device>
void test_dualview_combinations(unsigned int size)
{
  Impl::test_dualview_combinations<Scalar,Device> test(size);
  ASSERT_EQ( test.result,0);

}


} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP
