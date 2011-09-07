#include <iostream>
#include <cstdlib>

namespace test{

  void test_Host(int b, int e, int r);
  void test_TPI(int b, int e, int r, int t);
  void test_Cuda(int b, int e, int r);

}

int main( int argc , char ** argv ){

  int threads, runs;
  int beg = 22;
  int end = 23;

  if (argc > 1)
    threads = atoi(argv[1]);
  else
    threads = 4;

  if (argc > 2)
    runs = atoi(argv[2]);
  else
    runs = 1;

/************************************************************/
/*  argument 1 specifies the number of threads used by     	*/
/*  TPI and OpenMP; agrument 2 sets the number of runs    	*/
/*  to be averaged for the output timings          			*/
/*                             				 				*/
/*  Ex: to run with 8 threads and average 10 runs(LINUX):  	*/  
/*  ./fe_test.exe 8 10                    					*/  
/************************************************************/
  
  std::cout << std::endl << "\tStarting..." << std::endl;

  test::test_Host(beg, end, runs);
  test::test_TPI (beg, end, runs, threads);
//  test::test_Cuda(beg, end, runs);  

  return 0;

}

