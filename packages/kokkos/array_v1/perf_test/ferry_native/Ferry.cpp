#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include "offload.h"
#include "omp.h"
#ifdef __MIC__
#pragma offload_attribute(push, target(mic))
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#include <cilk/cilk.h>
#pragma offload_attribute(pop)
using namespace tbb;

#pragma warning( disable: 588)
#endif


#ifdef __MIC__
#include <micvec.h>
#endif

#define TYPE double

#define FERRY __declspec(target(mic))
struct 
FERRY
HexFill {

	TYPE * coord;
	HexFill(TYPE* arg_coords) : coord(arg_coords) {}

	void operator() (int ielem) const 
	{
	#ifdef __MIC__
		coord[(ielem * 0) + 0] = 0.;
		coord[(ielem * 1) + 0] = 0.;
		coord[(ielem * 2) + 0] = 0.;
		
		coord[(ielem * 0) + 1] = 1.;
		coord[(ielem * 1) + 1] = 0.;
		coord[(ielem * 2) + 1] = 0.;
		
		coord[(ielem * 0) + 2] = 1.;
		coord[(ielem * 1) + 2] = 1.;
		coord[(ielem * 2) + 2] = 0.;
		
		coord[(ielem * 0) + 3] = 0.;
		coord[(ielem * 1) + 3] = 1.;
		coord[(ielem * 2) + 3] = 0.;
		
		coord[(ielem * 0) + 4] = 0.;
		coord[(ielem * 1) + 4] = 0.;
		coord[(ielem * 2) + 4] = 1.;
		
		coord[(ielem * 0) + 5] = 1.;
		coord[(ielem * 1) + 5] = 0.;
		coord[(ielem * 2) + 5] = 1.;
		
		coord[(ielem * 0) + 6] = 1.;
		coord[(ielem * 1) + 6] = 1.;
		coord[(ielem * 2) + 6] = 1.;
		
		coord[(ielem * 0) + 7] = 0.;
		coord[(ielem * 1) + 7] = 1.;
		coord[(ielem * 2) + 7] = 1.;	
	#endif
	}
};

struct
FERRY
HexGrad {
	TYPE * coord ;
	TYPE * grad ;
	
	HexGrad(TYPE* arg_coord , TYPE* arg_grad) : coord(arg_coord) , grad(arg_grad) { }

	void operator()(int ielem) const {
	#ifdef __MIC__
#if 1
	TYPE x[8], y[8], z[8];
		
	for(int j = 0; j < 8 ; j++)
	{
		x[j] = coord[(ielem * 0) + j];
		y[j] = coord[(ielem * 1) + j];
		z[j] = coord[(ielem * 2) + j];
	}

	TYPE R42=(z[3] - z[1]);
    TYPE R52=(z[4] - z[1]);
    TYPE R54=(z[4] - z[3]);

    TYPE R63=(z[5] - z[2]);
    TYPE R83=(z[7] - z[2]);
    TYPE R86=(z[7] - z[5]);

    TYPE R31=(z[2] - z[0]);
    TYPE R61=(z[5] - z[0]);
    TYPE R74=(z[6] - z[3]);

    TYPE R72=(z[6] - z[1]);
    TYPE R75=(z[6] - z[4]);
    TYPE R81=(z[7] - z[0]);

    TYPE t1=(R63 + R54);
    TYPE t2=(R61 + R74);
    TYPE t3=(R72 + R81);

    TYPE t4 =(R86 + R42);
    TYPE t5 =(R83 + R52);
    TYPE t6 =(R75 + R31);

	int count = ielem;
	grad[count++] = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5)  + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54); 
    grad[count++] = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1)  - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61);
    grad[count++] = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2)  - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72);
    grad[count++] = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3)  + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83);
    grad[count++] = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2)  - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61);
    grad[count++] = (y[6] *  t5) - (y[4] *  t3)  - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
    grad[count++] = (y[7] *  t1) - (y[5] *  t5)  - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
    grad[count++] = (y[4] *  t2) - (y[6] *  t1)  + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54);

	R42=(x[3] - x[1]);
    R52=(x[4] - x[1]);
    R54=(x[4] - x[3]);

    R63=(x[5] - x[2]);
    R83=(x[7] - x[2]);
    R86=(x[7] - x[5]);

    R31=(x[2] - x[0]);
    R61=(x[5] - x[0]);
    R74=(x[6] - x[3]);

    R72=(x[6] - x[1]);
    R75=(x[6] - x[4]);
    R81=(x[7] - x[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

    grad[count++] = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
    grad[count++] = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
    grad[count++] = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
    grad[count++] = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
    grad[count++] = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
    grad[count++] = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
    grad[count++] = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
    grad[count++] = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);

    R42=(y[3] - y[1]);
    R52=(y[4] - y[1]);
    R54=(y[4] - y[3]);

    R63=(y[5] - y[2]);
    R83=(y[7] - y[2]);
    R86=(y[7] - y[5]);

    R31=(y[2] - y[0]);
    R61=(y[5] - y[0]);
    R74=(y[6] - y[3]);

    R72=(y[6] - y[1]);
    R75=(y[6] - y[4]);
    R81=(y[7] - y[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

    grad[count++] = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54);
    grad[count++] = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61);
    grad[count++] = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72);
    grad[count++] = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83);
    grad[count++] = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61);
    grad[count++] = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
    grad[count++] = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
    grad[count++] = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54);
#endif
#endif
	}
};

FERRY
void test( TYPE * coord, int size)
{
	#ifdef __MIC__
	HexFill hex(coord);

#if 0
	parallel_for(blocked_range<int>(0,size), [=] (const blocked_range<int> & r) {
		for(int i = r.begin() ; i != r.end(); ++i) 
			hex(i);
		} , auto_partitioner() );		
#endif
#if 0
	cilk_for(int i = 0 ; i < size ; i++) {
		hex(i);
	}
#endif
#if 0
	#pragma omp parallel for
	for(int i = 0 ; i < size ; i++) {
		hex(i);
	}
#endif
#if 1
	for(int i = 0 ; i < size ; i++) {
		hex(i);
	}
#endif
	#endif 
}

FERRY
void Grad(TYPE * coord, TYPE * grad , int size)
{
#ifdef __MIC__
	HexGrad hex(coord,grad);
#if 0
	parallel_for(blocked_range<int>(0,size), [=] (const blocked_range<int> & r) {
		for(int i = r.begin() ; i != r.end(); ++i) 
			hex(i);
		} , auto_partitioner() );		
#endif
#if 0
	cilk_for(int i = 0 ; i < size ; i++) {
		hex(i);
	}
#endif
#if 0
	#pragma omp parallel for
	for(int i = 0 ; i < size ; i++) {
		hex(i);
	}
#endif
#if 1
	#pragma simd
	for(int i = 0 ; i < size ; i++) {
		hex(i);
	}
#endif
#endif
}


int main(int argc , char* argv[])
{
	timeval start, stop, result;
	for(int i = 1024 ; i < 300000 ; i = i << 1 )	{

		//1D approach
		size_t data_size = sizeof(TYPE)*i*24;
		TYPE* coord = new TYPE[i];
		TYPE* grad  = new TYPE[i];
		//omp_set_num_threads_target(TARGET_MIC,0,64);
		#pragma offload target(mic) in(i) inout(coord : length(i))
		{
			test(coord, i);
		}

		gettimeofday(&start,NULL);
		#pragma offload target(mic) in(i) in(coord : length(i)) out(grad : length(i)) 
		{
			Grad(coord, grad,i);
		}
		gettimeofday(&stop,NULL);
		timersub(&stop, &start, &result);
		double time = (result.tv_sec + (result.tv_usec/1000000.0));
		std::cout<<"Compute Time: "<<time<<" "<< i<<std::endl;
		delete coord;
		delete grad;
	
	}
	std::cout<<"Finished"<<std::endl;
	return 0;
}

