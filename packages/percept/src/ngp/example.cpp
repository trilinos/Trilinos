#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include <iostream>

typedef Kokkos::View<double*>  ViewVectorType;
typedef Kokkos::View<double**> ViewMatrixType;

struct TestKernel {
  int S = 4194304, M = 1024, N = 4096, nrepeat = 100;

  ViewVectorType y;
  ViewVectorType x;
  ViewMatrixType A;

  TestKernel() 
  {
    S = pow( 2, 22);
    M = 1024;
    N = S/N;
    nrepeat = 100;

    y = ViewVectorType("y", N);
    x = ViewVectorType("x", M);
    A = ViewMatrixType("A", N, M);

    printf("  Initialize\n");
    initialize ();
  }

  void initialize () 
  {
    // Initialize y vector
    Kokkos::parallel_for(N, KOKKOS_LAMBDA (int i){
        y( i ) = 1; 
      });
    
    printf("  first pfor\n");
    // Initialize x vector 
    Kokkos::parallel_for(M, KOKKOS_LAMBDA (int i){
        x( i ) = 1; 
      });
    
    // Initialize A matrix, note 2D indexing
    Kokkos::parallel_for(N, KOKKOS_LAMBDA (int j){
        for ( int i = 0 ; i < M ; ++i ) {
          A( j , i ) = 1; 
        }
      });    
  }

  void run () 
  {
    const double solution = (double)N *(double)M;

    // Timer products
    struct timeval begin,end;
    
    gettimeofday(&begin,NULL);
    
    for ( int repeat = 0; repeat < nrepeat; repeat++) {
      
      //Application: <y,Ax> = y^T*A*x
      double result = 0;
      Kokkos::parallel_reduce(N, KOKKOS_LAMBDA ( int j, double &update ) {
          double temp2 = 0;
          for ( int i = 0 ; i < M ; ++i ) {
            temp2 += A( j , i ) * x( i );
          }
          update += y( j ) * temp2;
        }, result );
      
      EXPECT_NEAR(result, solution, 1.e-12);
    }

    gettimeofday(&end,NULL);
    
    // Calculate time
    double time = 1.0*(end.tv_sec-begin.tv_sec) +
      1.0e-6*(end.tv_usec-begin.tv_usec);
    
    // Calculate bandwidth.
    // Each matrix A row (each of length M) is read once.
    // The x vector (of length M) is read N times.
    // The y vector (of length N) is read once.
    // double Gbytes = 1.0e-9 * double(sizeof(double) * ( 2 * M * N + N ));
    double Gbytes = 1.0e-9 * double(sizeof(double) * ( M + M * N + N )) ;
    
    // Print results (problem size, time and bandwidth in GB/s)
    printf("  M( %d ) N( %d ) nrepeat ( %d ) problem( %g MB ) time( %g s ) bandwidth( %g GB/s )\n",
           M , N, nrepeat, Gbytes * 1000, time, Gbytes * nrepeat / time );
  }
};

TEST(test,foo)
{
  TestKernel().run();
}

void parallel_reduce_and_modify_array()
{
  const int N=3;
  ViewVectorType x = ViewVectorType("x", N);

  Kokkos::parallel_for(N, KOKKOS_LAMBDA (int i){
      x(i) = i; 
    });

  int sum=0;

  Kokkos::parallel_reduce(N, KOKKOS_LAMBDA (int i, int &local_sum) {
      local_sum += x(i); 
      x(i) *= 2;
    }, sum);

  EXPECT_EQ(sum, 0+1+2);
}

TEST(test,parallel_reduce_and_modify_array)
{
  parallel_reduce_and_modify_array();
}

KOKKOS_INLINE_FUNCTION
double metric(const double x, bool& valid)
{
  double val = 0.0;
  valid = true;
  for (unsigned i=0; i<100; i++)
    val += std::sin(M_PI*i*x);
  if (val<0) valid = false;
  return val;
}

void test_metric() {

  int nrepeat = 500;

  unsigned M = 10000;
  ViewVectorType x = ViewVectorType("x", M);

  Kokkos::parallel_for(M, KOKKOS_LAMBDA (int i){
      x(i) = (double)i/(double)M;
    });

  // Timer products
  struct timeval begin,end;
  
  gettimeofday(&begin,NULL);
  
  for ( int repeat = 0; repeat < nrepeat; repeat++) {
    
    double total_metric = 0;
    Kokkos::parallel_reduce(M, KOKKOS_LAMBDA (int i, double &local_metric){
        bool valid;
        local_metric += metric(x(i),valid);
      }, total_metric);
    
    if (repeat==nrepeat-1)
      std::cout << "total_metric = " << total_metric << std::endl;
  }

  gettimeofday(&end,NULL);
  
  // Calculate time
  double time = 1.0*(end.tv_sec-begin.tv_sec) +
    1.0e-6*(end.tv_usec-begin.tv_usec);
  
  // Calculate bandwidth.
  // Each matrix A row (each of length M) is read once.
  // The x vector (of length M) is read N times.
  // The y vector (of length N) is read once.
  // double Gbytes = 1.0e-9 * double(sizeof(double) * ( 2 * M * N + N ));
  double Gbytes = 1.0e-9 * double(sizeof(double) * ( M )) ;
  
  // Print results (problem size, time and bandwidth in GB/s)
  printf("  M( %d ) nrepeat ( %d ) problem( %g MB ) time( %g s ) bandwidth( %g GB/s )\n",
         M , nrepeat, Gbytes * 1000, time, Gbytes * nrepeat / time );
}

TEST(foo, metric)
{
  test_metric();
}
