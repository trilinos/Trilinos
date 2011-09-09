#include <iostream>
#include <cstdlib>

namespace Test{

	void test_Host(int beg, int end, int r);
	void test_TPI (int beg, int end, int r, int t);
	void test_TBB(int beg, int end, int r, int t);
	void test_Cuda(int beg, int end, int r);
}

int main(int argc, char ** argv){

	if ( argc != 5) 
	{
		int beg = 5;
		int end = 6;
		int runs = 1;
		int threads = 4;

		Test::test_Host(beg, end, runs);
		Test::test_TPI (beg, end, runs, threads);
		Test::test_TBB (beg, end, runs, threads);
		Test::test_Cuda(beg , end, runs);
		
	}
	else {
		int beg = atoi(argv[1]);
		int end = atoi(argv[2]);
		int runs = atoi(argv[3]);
		int threads = atoi(argv[4]);

		Test::test_Host(beg, end, runs);
		Test::test_TPI (beg, end, runs, threads);
		Test::test_TBB (beg, end, runs, threads);
		Test::test_Cuda(beg , end, runs);
	}
	return 0;

}
