#ifndef KOKKOS_TEST_DUALVIEW_HPP
#define KOKKOS_TEST_DUALVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Atomic.hpp>
#include <cmath>

namespace Test {

namespace Impl{

struct RandomProperties {
  uint64_t count;
  double mean;
  double variance;
  double covariance;
  RandomProperties() {
    count = 0;
    mean = 0.0;
    variance = 0.0;
    covariance = 0.0;
  }
  RandomProperties& operator+=(const RandomProperties& add) {
    count      += add.count;
    mean       += add.mean;
    variance   += add.variance;
    covariance += add.covariance;
    return *this;
  }
  volatile RandomProperties& operator+=(const volatile RandomProperties& add) volatile {
    count      += add.count;
    mean       += add.mean;
    variance   += add.variance;
    covariance += add.covariance;
    return *this;
  }
};

template<class GeneratorPool, class Scalar>
struct test_random_functor {
  typedef typename GeneratorPool::generator_type rnd_type;

  typedef RandomProperties value_type;
  typedef typename GeneratorPool::device_type device_type;

  GeneratorPool rand_pool;
  const double mean;
  typedef Kokkos::View<int[1024000],typename GeneratorPool::device_type> type_1d;
  type_1d density_1d;
  typedef Kokkos::View<int[128][128][128],typename GeneratorPool::device_type> type_3d;
  type_3d density_3d;

  test_random_functor(GeneratorPool rand_pool_,type_1d d1d, type_3d d3d):
      rand_pool(rand_pool_),mean(0.5*Kokkos::rand<rnd_type,Scalar>::max()),
      density_1d(d1d),density_3d(d3d) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, RandomProperties& prop) const {
    rnd_type rand_gen = rand_pool.get_state();
    for(int k = 0;k<1024;k++) {
      const Scalar tmp = Kokkos::rand<rnd_type,Scalar>::draw(rand_gen);
      prop.count++;
      prop.mean += tmp;
      prop.variance += (tmp-mean)*(tmp-mean);
      const Scalar tmp2 = Kokkos::rand<rnd_type,Scalar>::draw(rand_gen);
      prop.count++;
      prop.mean += tmp2;
      prop.variance += (tmp2-mean)*(tmp2-mean);
      prop.covariance += (tmp-mean)*(tmp2-mean);
      const Scalar tmp3 = Kokkos::rand<rnd_type,Scalar>::draw(rand_gen);
      prop.count++;
      prop.mean += tmp3;
      prop.variance += (tmp3-mean)*(tmp3-mean);
      prop.covariance += (tmp2-mean)*(tmp3-mean);

      Kokkos::atomic_fetch_add(&density_1d(int(1024000.0*tmp/Kokkos::rand<rnd_type,Scalar>::max())),1);
      Kokkos::atomic_fetch_add(&density_1d(int(1024000.0*tmp2/Kokkos::rand<rnd_type,Scalar>::max())),1);
      Kokkos::atomic_fetch_add(&density_1d(int(1024000.0*tmp3/Kokkos::rand<rnd_type,Scalar>::max())),1);
      Kokkos::atomic_fetch_add(&density_3d(int(128.0*tmp/Kokkos::rand<rnd_type,Scalar>::max()),
                 int(128.0*tmp2/Kokkos::rand<rnd_type,Scalar>::max()),
                 int(128.0*tmp3/Kokkos::rand<rnd_type,Scalar>::max())),1);
    }
    rand_pool.free_state(rand_gen);
  }
};

template <class RandomGenerator,class Scalar>
struct test_random_scalar {
  typedef typename RandomGenerator::generator_type rnd_type;

  int pass_mean,pass_var,pass_covar;
  test_random_scalar(
      typename test_random_functor<RandomGenerator,int>::type_1d& density_1d,
      typename test_random_functor<RandomGenerator,int>::type_3d& density_3d,
      RandomGenerator& pool, unsigned int num_draws) {
    RandomProperties result;

    Kokkos::parallel_reduce(num_draws/1024,test_random_functor<RandomGenerator,Scalar>(pool,density_1d,density_3d),result);

    //printf("Result: %lf %lf %lf\n",result.mean/num_draws/3,result.variance/num_draws/3,result.covariance/num_draws/2);
    double tolerance = 1.6*sqrt(1.0/num_draws);
    double mean_expect = 0.5*Kokkos::rand<rnd_type,Scalar>::max();
    double variance_expect = 1.0/3.0*mean_expect*mean_expect;
    double mean_eps = mean_expect/(result.mean/num_draws/3)-1.0;
    double variance_eps = variance_expect/(result.variance/num_draws/3)-1.0;
    double covariance_eps = result.covariance/num_draws/2/variance_expect;
    pass_mean  = ((-tolerance < mean_eps) &&
                  ( tolerance > mean_eps)) ? 1:0;
    pass_var   = ((-tolerance < variance_eps) &&
                  ( tolerance > variance_eps)) ? 1:0;
    pass_covar = ((-1.4*tolerance < covariance_eps) &&
                  ( 1.4*tolerance > covariance_eps)) ? 1:0;
    printf("Pass: %i %i %e %e %e || %e\n",pass_mean,pass_var,mean_eps,variance_eps,covariance_eps,tolerance);
  }
};

template <class RandomGenerator>
void test_random(unsigned int num_draws)
{
  typedef typename RandomGenerator::generator_type rnd_type;
  typename test_random_functor<RandomGenerator,int>::type_1d density_1d("D1d");
  typename test_random_functor<RandomGenerator,int>::type_3d density_3d("D3d");

  RandomGenerator pool(31891);
  test_random_scalar<RandomGenerator,int> test_int(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_int.pass_mean,1);
  ASSERT_EQ( test_int.pass_var,1);
  ASSERT_EQ( test_int.pass_covar,1);
  test_random_scalar<RandomGenerator,unsigned int> test_uint(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_uint.pass_mean,1);
  ASSERT_EQ( test_uint.pass_var,1);
  ASSERT_EQ( test_uint.pass_covar,1);
  test_random_scalar<RandomGenerator,int64_t> test_int64(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_int64.pass_mean,1);
  ASSERT_EQ( test_int64.pass_var,1);
  ASSERT_EQ( test_int64.pass_covar,1);
  test_random_scalar<RandomGenerator,uint64_t> test_uint64(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_uint64.pass_mean,1);
  ASSERT_EQ( test_uint64.pass_var,1);
  ASSERT_EQ( test_uint64.pass_covar,1);
  test_random_scalar<RandomGenerator,float> test_float(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_float.pass_mean,1);
  ASSERT_EQ( test_float.pass_var,1);
  ASSERT_EQ( test_float.pass_covar,1);
  test_random_scalar<RandomGenerator,double> test_double(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_double.pass_mean,1);
  ASSERT_EQ( test_double.pass_var,1);
  ASSERT_EQ( test_double.pass_covar,1);

}
}

} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP
