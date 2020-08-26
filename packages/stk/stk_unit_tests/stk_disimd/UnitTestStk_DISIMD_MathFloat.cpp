#include <gtest/gtest.h>

#ifndef STK_KOKKOS_SIMD
#define STK_KOKKOS_SIMD
#endif

#include <stk_simd/Simd.hpp>
#include <cmath>
#include <algorithm>

#include "DISimdFixture.hpp"

// Math ToolKit unit tests

TEST(StkDISimd, SimdAddSubtractMultDivideFloat) 
{
  int N = 400000;
  double t0; // timing variable

  std::vector<float> x(N);
  std::vector<float> y(N);
 
  std::vector<float> out1(N);
  std::vector<float> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.5)/RAND_MAX;
    y[n] = 26*(rand()-0.4)/RAND_MAX;
  }
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::simd::load(&y[n]);
    const stk::simd::Float c = ( b*(a-b) )/a + 1;
    stk::simd::store(&out1[n],c);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "DI-SIMD ADD,SUB,MUL,DIV took " << t0 << " seconds" <<  std::endl;
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    const float b = y[n];
    out2[n] = ( b*(a-b) )/a + 1;
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real ADD,SUB,MUL,DIV took " << t0 << " seconds" <<  std::endl;

  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_EQ( maxerr, 0.0 );
}

TEST(StkDISimd, SimdMiscSelfAddSubEtcFloat) 
{
  int N = 10000;

  std::vector<float> x(N);
  std::vector<float> y(N);
 
  std::vector<float> out1(N);
  std::vector<float> out2(N);

  float maxerr = 0.0;

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.5)/RAND_MAX;
    y[n] = 26*(rand()-0.4)/RAND_MAX;
  }
  
  // *= -a and traits

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::simd::load(&y[n]);
    stk::simd::Float c = 3.0f;
    c *= -b + a + stk::Traits<stk::simd::Float>::TWO/5;
    stk::simd::store(&out1[n],c);
  }
  
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    const float b = y[n];
    out2[n] = 3.0f;
    out2[n] *= -b + a + stk::Traits<float>::TWO/5;
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_EQ( maxerr, 0.0 );
  
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::simd::load(&y[n]);
    stk::simd::Float c = 3.2f;
    c /= -b+a+5.6f*c;
    stk::simd::store(&out1[n],c);
  }
  
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    const float b = y[n];
    out2[n] = 3.2f;
    out2[n] /= -b+a+5.6f*out2[n];
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_EQ( maxerr, 0.0 );

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::simd::load(&y[n]);
    stk::simd::Float c = stk::Traits<stk::simd::Float>::THIRD;
    c += -(b/(1.0f-a)+(c+5.2f));
    stk::simd::store(&out1[n],c);
  }
  
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    const float b = y[n];
    float c = 1.0/3.0;
    c += -(b/(1.0f-a)+(c+5.2f));
    out2[n] = c;
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_EQ( maxerr, 0.0);

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::simd::load(&y[n]);
    stk::simd::Float c = stk::simd::Float(-0.5);
    c -= 5.2f+(b/(a-5.4f)+3.5f/c);
    stk::simd::store(&out1[n],c);
  }
  
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    const float b = y[n];
    out2[n] = -0.5f;
    out2[n] -= 5.2f+(b/(a-5.4f)+3.5f/out2[n]);
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_EQ( maxerr, 0.0 );

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::simd::load(&y[n]);
    const stk::simd::Float c = -0.3f;
    const stk::simd::Float d = -c + 1.4f*a/b;
    stk::simd::store(&out1[n],d);
  }
  
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    const float b = y[n];
    const float c = -0.3f;
    out2[n] = -c + 1.4f*a/b;
  }
 
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_EQ( maxerr, 0.0 );
}

TEST(StkDISimd, Simd_fmaddFloat)
{

  int N = 400000;
  double t0; // timing variable

  std::vector<float> x(N);
  std::vector<float> y(N);
  std::vector<float> z(N);

  std::vector<float> out1(N);
  std::vector<float> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 2.1*(rand()-0.5)/RAND_MAX;
    y[n] = 6.0*(rand()-0.4)/RAND_MAX;
    z[n] = 6.0*(rand()-0.4)/RAND_MAX;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::simd::load(&y[n]);
    const stk::simd::Float c = stk::simd::load(&z[n]);
    stk::simd::Float d = stk::math::fmadd(a,b,c);
    stk::simd::store(&out1[n],d);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD ADD,SUB,MUL,DIV took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
#if defined(__INTEL_COMPILER)
#pragma novector
#endif
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    const float b = y[n];
    const float c = z[n];
    out2[n] = a*b+c;
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real ADD,SUB,MUL,DIV took " << t0 << " seconds" <<  std::endl;

  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 1e-5 );

}

TEST(StkDISimd, SimdSqrtFloat) 
{

  int N = 400000;
  double t0; // timing variable

  std::vector<float> x(N);
 
  std::vector<float> out1(N);
  std::vector<float> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 21.0f*(rand()+1)/RAND_MAX;
  }
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::math::sqrt(a);
    stk::simd::store(&out1[n],b);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "DI-SIMD SQRT took " << t0 << " seconds" <<  std::endl;
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    out2[n] = std::sqrt(a);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real SQRT took " << t0 << " seconds" <<  std::endl;

  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
}

TEST(StkDISimd, SimdLogFloat) 
{

  int N = 400000;
  double t0; // timing variable

  std::vector<float> x(N);
 
  std::vector<float> out1(N);
  std::vector<float> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 21.0f*(rand()+1.0f)/RAND_MAX;
  }
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::math::log(a);
    stk::simd::store(&out1[n],b);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "DI-SIMD Log took " << t0 << " seconds" <<  std::endl;
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    out2[n] = std::log(a);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real Log took " << t0 << " seconds" <<  std::endl;

  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 1.0e-6 );
}

TEST(StkDISimd, SimdExpFloat) 
{

  int N = 400000;
  double t0; // timing variable

  std::vector<float> x(N);
 
  std::vector<float> out1(N);
  std::vector<float> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 5.4f*(rand()/RAND_MAX-0.5);
  }
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::math::exp(a);
    stk::simd::store(&out1[n],b);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "DI-SIMD Exp took " << t0 << " seconds" <<  std::endl;
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    out2[n] = std::exp(a);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real Exp took " << t0 << " seconds" <<  std::endl;

  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
}

TEST(StkDISimd, SimdPowAFloat) 
{

  int N = 400000;
  double t0; // timing variable

  std::vector<float> x(N);
  std::vector<float> y(N/stk::simd::nfloats);
 
  std::vector<float> out1(N);
  std::vector<float> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 5.4*(rand()/RAND_MAX-0.5);
  }
  for (int n=0; n < N/stk::simd::nfloats; ++n) {
    y[n] = 3.2*(rand()/RAND_MAX-0.5);
  }
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const float b = y[n/stk::simd::nfloats];
    const stk::simd::Float c = stk::math::pow(a,b);
    stk::simd::store(&out1[n],c);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "DI-SIMD Exp took " << t0 << " seconds" <<  std::endl;
  
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N/stk::simd::nfloats; ++n) {
    float exp = y[n];
    for (int i=0; i < stk::simd::nfloats; ++i) {
      const float a = x[stk::simd::nfloats*n+i];
      out2[stk::simd::nfloats*n+i] = std::pow(a,exp);
    }
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real Exp took " << t0 << " seconds" <<  std::endl;

  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
}

TEST(StkDISimd, SimdPowBFloat)
{

  int N = 400000;
  double t0; // timing variable

  std::vector<float> x(N);
  std::vector<float> y(N);

  std::vector<float> out1(N);
  std::vector<float> out2(N);

  for (int n=0; n < N; ++n) {
    x[n] = 5.4*(rand()/RAND_MAX-0.5);
  }
  for (int n=0; n < N; ++n) {
    y[n] = 3.2*(rand()/RAND_MAX-0.5);
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    const stk::simd::Float b = stk::simd::load(&y[n]);
    const stk::simd::Float c = stk::math::pow(a,b);
    stk::simd::store(&out1[n],c);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD Exp took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    const float a = x[n];
    const float b = y[n];
    out2[n] = std::pow(a,b);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real Exp took " << t0 << " seconds" <<  std::endl;

  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

}

TEST(StkDISimd, SimdCbrtFloat)
{

  int N = 800000;
  double t0; // timing variable

  std::vector<float> x(N);

  std::vector<float> out1(N);
  std::vector<float> out2(N);


  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.5)/RAND_MAX;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float d = stk::math::cbrt(xl);
    stk::simd::store(&out1[n],d);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "SIMD cbrt took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    out2[n] = std::cbrt(x[n]);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Real cbrt took " << t0 << " seconds" <<  std::endl;

  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 2.5e-7 );
}

TEST(StkDISimd, SimdTimeLoadStoreDataLayoutFloat)
{

  double t0;

  const int N = 20000;
  const int sz = 32;

  std::vector<float,non_std::AlignedAllocator<float,64> > x(sz*N);
  std::vector<float,non_std::AlignedAllocator<float,64> > y(sz*N);
  std::vector<float,non_std::AlignedAllocator<float,64> > z(sz*N);
  std::vector<float,non_std::AlignedAllocator<float,64> > w(sz*N);

  for (int n=0; n < N; ++n) {
    for (int i=0; i < sz; ++i) {
      x[sz*n+i] = 1+std::sqrt(n+n*(i+1));
    }
  }

  float* X = x.data();
  float* Y = y.data();

  stk::simd::Float a[sz];

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats, X+=sz*stk::simd::nfloats, Y+=sz*stk::simd::nfloats) {
    // load the slow way
    for (int i=0; i < sz; ++i) a[i] = stk::simd::load( X+i, sz );
    // store the slow way
    for (int i=0; i < sz; ++i) stk::simd::store( Y+i, a[i]+stk::simd::Float(1.0), sz );
  }

  t0 += stk::get_time_in_seconds();
  std::cout << "Method 1: Offset load/store, took " << t0 << " seconds" <<  std::endl;

  // reorder the arrays...
  const int ssz = sz*stk::simd::nfloats;
  float tmp[ssz];

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        tmp[sz*j+i] = x[sz*n+sz*j+i];
      }
    }
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        x[sz*n+stk::simd::nfloats*i+j] = tmp[sz*j+i];
      }
    }
  }

  X = x.data();
  Y = z.data();

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats, X+=sz*stk::simd::nfloats, Y+=sz*stk::simd::nfloats) {

    for (int i=0; i < sz; ++i) {
      a[i] = stk::simd::load( X+stk::simd::nfloats*i );
    }

    for (int i=0; i < sz; ++i) {
      stk::simd::store( Y+stk::simd::nfloats*i, a[i]+stk::simd::Float(1.0) );
    }

  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Method 2: Reordered load, took " << t0 << " seconds" << std::endl;

  X = x.data();
  Y = w.data();

  stk::simd::Float* XX = stk::simd::simd_ptr_cast<float>(X);
  stk::simd::Float* YY = stk::simd::simd_ptr_cast<float>(Y);

  t0 = -stk::get_time_in_seconds();

  for (int n=0,i=0; n < N; n+=stk::simd::nfloats) {
    for (int j=0; j < sz; ++j,++i) {
      YY[i] = XX[i] + stk::simd::Float(1.0);
    }

  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Method 3: Load/Store in place, took " << t0 << " seconds" << std::endl;

  //reorder back! (y and w)
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        tmp[sz*j+i] = z[sz*n+stk::simd::nfloats*i+j];
      }
    }
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        z[sz*n+sz*j+i] = tmp[sz*j+i];
      }
    }
  }

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        tmp[sz*j+i] = w[sz*n+stk::simd::nfloats*i+j];
      }
    }
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        w[sz*n+sz*j+i] = tmp[sz*j+i];
      }
    }
  }

  // figure out error
  float maxerr = 0.0;
  for (int n=0; n < sz*N; ++n) {
    float err = stk::math::abs(y[n]-z[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  //printf("maxerror = %g\n",maxerr);

  ASSERT_NEAR( maxerr, 0.0, 1.0e-16 );

  // figure out error
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(y[n]-w[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  //printf("maxerror = %g\n",maxerr);

  ASSERT_NEAR( maxerr, 0.0, 1.0e-16 );

}

TEST(StkDISimd, SimdTimeLoadStoreInnerProductFloat)
{

  double t0;

  const int N = 20000;
  const int sz = 32;

  std::vector<float,non_std::AlignedAllocator<float,64> > x(sz*N);
  std::vector<float,non_std::AlignedAllocator<float,64> > y(N);
  std::vector<float,non_std::AlignedAllocator<float,64> > z(N);
  std::vector<float,non_std::AlignedAllocator<float,64> > w(N);

  for (int n=0; n < N; ++n) {
    for (int i=0; i < sz; ++i) {
      x[sz*n+i] = 1+std::sqrt(n+n*(i+1));
    }
  }

  float* X = x.data();
  float* Y = y.data();

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats, X+=sz*stk::simd::nfloats, Y+=stk::simd::nfloats) {
    // load the slow way
    stk::simd::Float dot = stk::simd::Float(0.0);
    for (int i=0; i < sz; ++i) {
      stk::simd::Float tmp = stk::simd::load( X+i, sz );
      dot += tmp*tmp;
    }
    stk::simd::store( Y, dot );
  }

  t0 += stk::get_time_in_seconds();
  std::cout << "Method 1: Offset load/store, took " << t0 << " seconds" <<  std::endl;

  // reorder the arrays...
  const int ssz = sz*stk::simd::nfloats;
  float tmp[ssz];

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        tmp[sz*j+i] = x[sz*n+sz*j+i];
      }
    }
    for (int i=0; i < sz; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        x[sz*n+stk::simd::nfloats*i+j] = tmp[sz*j+i];
      }
    }
  }

  X = x.data();
  Y = z.data();

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats, X+=sz*stk::simd::nfloats, Y+=stk::simd::nfloats) {
    stk::simd::Float dot = stk::simd::Float(0.0);
    for (int i=0; i < sz; ++i) {
      const stk::simd::Float tmpFloats = stk::simd::load( X+stk::simd::nfloats*i );
      dot += tmpFloats*tmpFloats;
    }

    stk::simd::store( Y, dot );

  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Method 2: Reordered load, took " << t0 << " seconds" << std::endl;

  X = x.data();
  Y = w.data();

  const stk::simd::Float* XX = stk::simd::simd_ptr_cast<float>(X);
  stk::simd::Float* YY = stk::simd::simd_ptr_cast<float>(Y);

  t0 = -stk::get_time_in_seconds();

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    stk::simd::Float dot = stk::simd::Float(0.0);
    for (int i=0; i < sz; ++i) {
      const stk::simd::Float tmpFloats = XX[i];
      dot += tmpFloats*tmpFloats;
    }
    XX+=sz;
    *(YY++) = dot;
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Method 3: Load/Store in place, took " << t0 << " seconds" << std::endl;

  //reorder back! (y and w)
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < 1; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        tmp[1*j+i] = z[1*n+stk::simd::nfloats*i+j];
      }
    }
    for (int i=0; i < 1; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        z[1*n+1*j+i] = tmp[1*j+i];
      }
    }
  }

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < 1; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        tmp[1*j+i] = w[1*n+stk::simd::nfloats*i+j];
      }
    }
    for (int i=0; i < 1; ++i) {
      for (int j=0; j < stk::simd::nfloats; ++j) {
        w[1*n+1*j+i] = tmp[1*j+i];
      }
    }
  }

  // figure out error
  float maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(y[n]-z[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  //printf("maxerror = %g\n",maxerr);

  ASSERT_NEAR( maxerr, 0.0, 1.0e-16 );

  // figure out error
  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(y[n]-w[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  //printf("maxerror = %g\n",maxerr);

  ASSERT_NEAR( maxerr, 0.0, 1.0e-16 );

}


class SimdFloatBool : public ::testing::Test
{
public:
  stk::simd::Float zeros_ones()
  {
    for (int i=0; i<stk::simd::nfloats; i++) {
      m_data[i] = i%2;
    }
    return stk::simd::load(m_data);
  }

private:
  float m_data[stk::simd::nfloats];
};

TEST_F(SimdFloatBool, allFalse)
{
  stk::simd::Boolf simdBool(false);

  for (int i=0; i<stk::simd::nfloats; i++) {
    EXPECT_FALSE(simdBool[i]);
  }
}

TEST_F(SimdFloatBool, allTrue)
{
  stk::simd::Boolf simdTrue(true);

  for (int i=0; i<stk::simd::nfloats; i++) {
    EXPECT_TRUE(simdTrue[i]);
  }
}

TEST_F(SimdFloatBool, someTrue)
{
  stk::simd::Float half(0.5);
  stk::simd::Float zeroOne = zeros_ones();

  stk::simd::Boolf simdBool = (zeroOne < half);

  for (int i=0; i<stk::simd::nfloats; i++) {
    if (zeroOne[i] < half[i]) {
      EXPECT_TRUE(simdBool[i]);
    }
    else {
      EXPECT_FALSE(simdBool[i]);
    }
  }
}


TEST(StkDISimd, SimdIfThenBoolFloat)
{

  typedef stk::Traits<float>::bool_type float_bool;
  typedef stk::Traits<stk::simd::Float>::bool_type Floats_bool;

  ASSERT_TRUE( stk::simd::are_all( stk::Traits<stk::simd::Float>::TRUE_VAL ) );
  ASSERT_TRUE( stk::Traits<float>::TRUE_VAL );

  ASSERT_FALSE( stk::simd::are_all( stk::Traits<stk::simd::Float>::FALSE_VAL ) );
  ASSERT_FALSE( stk::Traits<float>::FALSE_VAL );

  const int N = 2000;
  const float a = 5.1;
  const float b = -3.2;
  float maxerr = 0.0;

  std::vector<float> x(N);
  std::vector<float> y(N);
  std::vector<float> x2(N);
  std::vector<float> y2(N);
  std::vector<float> z1(N);
  std::vector<float> z2(N);
  std::vector<float> z3(N);

  for (int n=0; n < N; ++n) {
    x[n] = (rand()-0.5)/RAND_MAX;
    y[n] = rand()/RAND_MAX;
    x2[n] = (rand()-0.1)/RAND_MAX;
    y2[n] = (rand()+0.5)/RAND_MAX;
  }

  y[10] = x[10] = 5;
  y[33] = x[33] = 6.4;
  y[101]= x[101]= -3;

  // less than

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = x[n+i] < y[n+i] ? a : b;
      float_bool tmp = x[n+i] < y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);

    stk::simd::Float zl = stk::math::if_then_else(xl < yl, stk::simd::Float(a), stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // less than equal

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = x[n+i] <= y[n+i] ? a : b;
      float_bool tmp = x[n+i] <= y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);

    stk::simd::Float zl = stk::math::if_then_else(xl <= yl, stk::simd::Float(a), stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // equal

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = x[n+i] == y[n+i] ? a : b;
      float_bool tmp = x[n+i] == y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);

    stk::simd::Float zl = stk::math::if_then_else(xl == yl, stk::simd::Float(a), stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );


  // greater than equal

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = x[n+i] >= y[n+i] ? a : b;
      float_bool tmp = x[n+i] >= y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);

    stk::simd::Float zl = stk::math::if_then_else(xl >= yl, stk::simd::Float(a), stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // greater than

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = x[n+i] > y[n+i] ? a : b;
      float_bool tmp = x[n+i] > y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);

    stk::simd::Float zl = stk::math::if_then_else(xl > yl, stk::simd::Float(a), stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // not equal

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = x[n+i] != y[n+i] ? a : b;
      float_bool tmp = x[n+i] != y[n+i];
      z3[n+i] = stk::math::if_then_else(tmp,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);

    stk::simd::Float zl = stk::math::if_then_else(xl != yl, stk::simd::Float(a), stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // if then zero

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = x[n+i] < y[n+i] ? a : 0;
      float_bool tmp = x[n+i] < y[n+i];
      z3[n+i] = stk::math::if_then_else_zero(tmp,a);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);

    stk::simd::Float zl = stk::math::if_then_else_zero(xl < yl, stk::simd::Float(a));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // if ! then

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = !(x[n+i] < y[n+i]) ? a : b;
      float_bool tmp = x[n+i] < y[n+i];
      z3[n+i] = stk::math::if_then_else(!tmp,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    Floats_bool tmp = xl < yl;
    stk::simd::Float zl = stk::math::if_then_else(!tmp, stk::simd::Float(a),stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // &&

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      float_bool tmp  = x[n+i] < y[n+i];
      float_bool tmp2 = x2[n+i] > y2[n+i];
      z2[n+i] = (tmp&tmp2) ? a : b;
      z3[n+i] = stk::math::if_then_else(tmp&&tmp2,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    const stk::simd::Float xl2 = stk::simd::load(&x2[n]);
    const stk::simd::Float yl2 = stk::simd::load(&y2[n]);
    Floats_bool tmp = xl < yl;
    Floats_bool tmp2 = xl2 > yl2;
    stk::simd::Float zl = stk::math::if_then_else(tmp&&tmp2, stk::simd::Float(a),stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // ||

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      float_bool tmp  = x[n+i] < y[n+i];
      float_bool tmp2 = x2[n+i] > y2[n+i];
      z2[n+i] = (tmp|tmp2) ? a : b;
      z3[n+i] = stk::math::if_then_else(tmp||tmp2,a,b);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    const stk::simd::Float xl2 = stk::simd::load(&x2[n]);
    const stk::simd::Float yl2 = stk::simd::load(&y2[n]);
    Floats_bool tmp = xl < yl;
    Floats_bool tmp2 = xl2 > yl2;
    stk::simd::Float zl = stk::math::if_then_else(tmp||tmp2,
                                                   stk::simd::Float(a),
                                                   stk::simd::Float(b));
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // any

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    bool anyl=false;
    for (int i=0; i < stk::simd::nfloats; ++i) {
      float_bool tmp  = x[n+i] < y[n+i];
      anyl = tmp|anyl;
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    Floats_bool tmp = xl < yl;
    bool anyl_simd = stk::simd::are_any(tmp,stk::simd::nfloats);

    ASSERT_TRUE(anyl_simd==anyl);

  }

  // all

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    bool alll=true;
    for (int i=0; i < stk::simd::nfloats; ++i) {
      float_bool tmp  = x[n+i] < y[n+i];
      alll = tmp&alll;
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    Floats_bool tmp = xl < yl;
    bool alll_simd = stk::simd::are_all(tmp,stk::simd::nfloats);

    ASSERT_TRUE(alll_simd==alll);

  }

#if defined(StkSimd_SIMD) // these don't make sense for non-simd


  for (int n=0; n < N; n+=stk::simd::nfloats) {
    bool anyl=false;
    for (int i=0; i < stk::simd::nfloats-1; ++i) {
      float_bool tmp  = x[n+i] < y[n+i];
      anyl = tmp|anyl;
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    Floats_bool tmp = xl < yl;
    bool anyl_simd = stk::simd::are_any(tmp,stk::simd::nfloats-1);

    ASSERT_TRUE(anyl_simd==anyl);

  }

  // all (partial)

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    bool alll=true;
    for (int i=0; i < stk::simd::nfloats-1; ++i) {
      float_bool tmp  = x[n+i] < y[n+i];
      alll = tmp&alll;
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    Floats_bool tmp = xl < yl;
    bool alll_simd = stk::simd::are_all(tmp,stk::simd::nfloats-1);

    ASSERT_TRUE(alll_simd==alll);

  }

#endif

}

TEST(StkDISimd, SimdSpecialFunctionsFloat)
{

  typedef stk::Traits<float>::bool_type float_bool;
  typedef stk::Traits<stk::simd::Float>::bool_type Floats_bool;

  const int N = 2000;
  float maxerr = 0.0;

  std::vector<float> x(N);
  std::vector<float> y(N);
  std::vector<float> z1(N);
  std::vector<float> z2(N);
  std::vector<float> z3(N);

  for (int n=0; n < N; ++n) {
    x[n] = (rand()-0.5)/RAND_MAX;
    y[n] = (rand()-0.5)/RAND_MAX;
  }

  // abs

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = std::abs(x[n+i]);
      z3[n+i] = stk::math::abs(x[n+i]);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    stk::simd::Float zl = stk::math::abs(xl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // multiplysign

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = x[n+i]*( (0.0 <= y[n+i]) - (y[n+1] < 0.0) );
      z3[n+i] = stk::math::multiplysign(x[n+i],y[n+i]);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    stk::simd::Float zl = stk::math::multiplysign(xl,yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // copysign

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = std::abs(x[n+i])*( (0.0 <= y[n+i]) - (y[n+1] < 0.0) );
      z3[n+i] = stk::math::copysign(x[n+i],y[n+i]);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    stk::simd::Float zl = stk::math::copysign(xl,yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // max

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = stk::math::max(x[n+i],y[n+i]);
      z3[n+i] = stk::math::max(x[n+i],y[n+i]);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    stk::simd::Float zl = stk::math::max(xl,yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

  // min

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    for (int i=0; i < stk::simd::nfloats; ++i) {
      z2[n+i] = stk::math::min(x[n+i],y[n+i]);
      z3[n+i] = stk::math::min(x[n+i],y[n+i]);
    }
    const stk::simd::Float xl = stk::simd::load(&x[n]);
    const stk::simd::Float yl = stk::simd::load(&y[n]);
    stk::simd::Float zl = stk::math::min(xl,yl);
    stk::simd::store(&z1[n],zl);
  }

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(z1[n]-z2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(z1[n]-z3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

}


TEST(StkDISimd, SimdTimeSet1VsConstFloats)
{
  int N = 1000000;
  double t0; // timing variable
  float maxerr;

  std::vector<float> x(N);

  std::vector<float> out1(N);
  std::vector<float> out2(N);
  std::vector<float> out3(N);
  std::vector<float> out4(N);

  const float three = 3.0;
  const float two = 2.0;

  for (int n=0; n < N; ++n) {
    x[n] = 21*(rand()-0.5)/RAND_MAX;
  }

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    stk::simd::store(&out1[n],stk::Traits<stk::simd::Float>::TWO+(a+stk::Traits<stk::simd::Float>::THREE));
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Adding a const Floats took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    stk::simd::store(&out2[n],stk::simd::Float(2.0)+(a+stk::simd::Float(3.0)));
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Adding a local const Floats took " << t0 << " seconds" <<  std::endl;

  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    const stk::simd::Float a = stk::simd::load(&x[n]);
    stk::simd::store(&out3[n],2.0f+(a+3.0f));
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Adding a 3.0 (with load1_pd) took " << t0 << " seconds" <<  std::endl;
  t0 = -stk::get_time_in_seconds();
  for (int n=0; n < N; ++n) {
    out4[n] = two+(x[n] + three);
  }
  t0 += stk::get_time_in_seconds();
  std::cout << "Non simd took " << t0 << " seconds" <<  std::endl;

  maxerr = 0.0;
  for (int n=0; n < N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(out1[n]-out3[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(out1[n]-out4[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );

}

template <typename REAL_TYPE> 
void NegateVec(REAL_TYPE * const in, REAL_TYPE * const out) {
  REAL_TYPE tmp[12];
  tmp[9] = -in[0];
  tmp[10] = -in[1];
  tmp[11] = -in[2];
  
  out[0] = tmp[9];
  out[1] = tmp[10];
  out[2] = tmp[11];
}

template <typename REAL_TYPE> static void NegateVec2(REAL_TYPE* const in, REAL_TYPE* const out) {
  static const REAL_TYPE ZERO(0.0);
  out[0] = ZERO - in[0];
  out[1] = ZERO - in[1];
  out[2] = ZERO - in[2];
}


TEST(StkDISimd, NegatingAVectorFloat)
{

  int N = 8000;

  std::vector<float> x(3*N);
  
  std::vector<float> out1(3*N);
  std::vector<float> out2(3*N);
  std::vector<float> out3(3*N);

  for (int n=0; n < 3*N; ++n) {
    x[n] = 21.0*(rand()-0.5)/RAND_MAX;
  }
 
  for (int n=0; n < N; ++n) {
    NegateVec(&x[3*n],&out1[3*n]);
  }

  stk::simd::Float a[3];
  stk::simd::Float b[3];
  for (int n=0; n < N; n+=stk::simd::nfloats) {
    stk::simd::load_array<3>(a,&x[3*n]);
    NegateVec(a,b);
    stk::simd::store_array<3>(&out2[3*n],b);
  }

  for (int n=0; n < N; n+=stk::simd::nfloats) {
    stk::simd::load_array<3>(a,&x[3*n]);
    NegateVec2(a,b);
    stk::simd::store_array<3>(&out3[3*n],b);
  }

  float maxerr = 0.0;
  for (int n=0; n < 3*N; ++n) {
    float err = stk::math::abs(out1[n]-out2[n]);
    maxerr = stk::math::max(err,maxerr);
    err = stk::math::abs(out1[n]-out3[n]);
    maxerr = stk::math::max(err,maxerr);
  }

  ASSERT_NEAR( maxerr, 0.0, 0.0 );
  
}

TEST(StkDISimd, SimdIsnanFloat)
{

  const int N = stk::simd::nfloats;

  std::vector<float> x(N);
  for (int i=0; i < stk::simd::nfloats; ++i) {
    x[i] = i;
  }

  stk::simd::Float X = stk::simd::load(x.data());
  stk::simd::Float Y = stk::simd::Float(1.0)/X;

  stk::simd::Boolf IsNaN = stk::math::isnan(Y);

  ASSERT_FALSE( stk::simd::are_any(IsNaN) );
  ASSERT_FALSE( stk::simd::are_all(IsNaN) );

  Y*= stk::simd::Float(0.0);

  IsNaN = stk::math::isnan(Y);

  ASSERT_TRUE( stk::simd::are_any(IsNaN) );
  ASSERT_TRUE( !stk::simd::are_all(IsNaN) || (stk::simd::nfloats==1) );
  
}

TEST(StkDISimd, SimdChooseFloat)
{
  std::cout<<"Using simd_abi = '"<<stk::unit_test_simd::get_simd_abi_native_name()<<"'"<<std::endl;

  SIMD_NAMESPACE::simd<float, SIMD_NAMESPACE::simd_abi::native> one(1.0f);
  SIMD_NAMESPACE::simd<float, SIMD_NAMESPACE::simd_abi::native> zero(0.0f);

  std::vector<float> maskVec(stk::simd::nfloats, 0.0f);
  maskVec[0] = 1.0f;
  SIMD_NAMESPACE::simd<float, SIMD_NAMESPACE::simd_abi::native> maskFloat(maskVec.data(), SIMD_NAMESPACE::element_aligned_tag());
  //SIMD_NAMESPACE::simd_mask<float, SIMD_NAMESPACE::simd_abi::native> mask(maskFloat.get());
  auto mask = (maskFloat == one);

  SIMD_NAMESPACE::simd<float, SIMD_NAMESPACE::simd_abi::native> masked = SIMD_NAMESPACE::choose(mask, one, zero);
  std::vector<float> maskedVec(stk::simd::nfloats);
  masked.copy_to(maskedVec.data(), SIMD_NAMESPACE::element_aligned_tag());

  ASSERT_EQ(maskVec, maskedVec);
}
