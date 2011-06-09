#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <iomanip>

__global__ void HexGrad (double * f, double * g){

//////////////////////////////////////////////////////////////////////////
//	Original kernel														//
//																		//
//	this kernel computes x, y, and z gradients for a finite element		//
//	routine																//
//////////////////////////////////////////////////////////////////////////

	unsigned int index = blockIdx.x*(blockDim.x) + threadIdx.x;
	unsigned int e_count = gridDim.x * blockDim.x;

	double x[8], y[8], z[8];

	x[0] = f[index]; index += e_count;
	x[1] = f[index]; index += e_count;
	x[2] = f[index]; index += e_count;
	x[3] = f[index]; index += e_count;
	x[4] = f[index]; index += e_count;
	x[5] = f[index]; index += e_count;
	x[6] = f[index]; index += e_count;
	x[7] = f[index]; index += e_count;

	y[0] = f[index]; index += e_count;
	y[1] = f[index]; index += e_count;
	y[2] = f[index]; index += e_count;
	y[3] = f[index]; index += e_count;
	y[4] = f[index]; index += e_count;
	y[5] = f[index]; index += e_count;
	y[6] = f[index]; index += e_count;
	y[7] = f[index]; index += e_count;

	z[0] = f[index]; index += e_count;
	z[1] = f[index]; index += e_count;
	z[2] = f[index]; index += e_count;
	z[3] = f[index]; index += e_count;
	z[4] = f[index]; index += e_count;
	z[5] = f[index]; index += e_count;
	z[6] = f[index]; index += e_count;
	z[7] = f[index]; 

	// z difference vectors
    double R42=(z[3] - z[1]);
    double R52=(z[4] - z[1]);
    double R54=(z[4] - z[3]);

    double R63=(z[5] - z[2]);
    double R83=(z[7] - z[2]);
    double R86=(z[7] - z[5]);

    double R31=(z[2] - z[0]);
    double R61=(z[5] - z[0]);
    double R74=(z[6] - z[3]);

    double R72=(z[6] - z[1]);
    double R75=(z[6] - z[4]);
    double R81=(z[7] - z[0]);

    double t1=(R63 + R54);
    double t2=(R61 + R74);
    double t3=(R72 + R81);

    double t4 =(R86 + R42);
    double t5 =(R83 + R52);
    double t6 =(R75 + R31);

    g[index] = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5)  + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54); index -= e_count;
    g[index] = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1)  - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61); index -= e_count;
    g[index] = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2)  - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72); index -= e_count;
    g[index] = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3)  + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83); index -= e_count;
    g[index] = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2)  - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61); index -= e_count;
    g[index] = (y[6] *  t5) - (y[4] *  t3)  - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72); index -= e_count;
    g[index] = (y[7] *  t1) - (y[5] *  t5)  - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83); index -= e_count;
    g[index] = (y[4] *  t2) - (y[6] *  t1)  + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54); index -= e_count;

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

    g[index]  = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54); index -= e_count;
    g[index]  = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61); index -= e_count;
    g[index]  = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72); index -= e_count;
    g[index]  = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83); index -= e_count;
    g[index]  = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61); index -= e_count;
    g[index]  = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72); index -= e_count;
    g[index]  = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83); index -= e_count;
    g[index]  = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54); index -= e_count;

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

    g[index] = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54); index -= e_count;
    g[index] = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61); index -= e_count;
    g[index] = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72); index -= e_count;
    g[index] = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83); index -= e_count;
    g[index] = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61); index -= e_count;
    g[index] = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72); index -= e_count;
    g[index] = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83); index -= e_count;
    g[index] = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54);

}

__global__ void HexGrad2 (double * f, double * g){

//////////////////////////////////////////////////////////////////////////
//	first optimization of the original kernel							//
//																		//
//	this kernel features staggered reads, which improves performance	//
//	by only performing memory reads immediately before the data is		//
//	needed																//
//////////////////////////////////////////////////////////////////////////


	unsigned int index = blockIdx.x*(blockDim.x) + threadIdx.x;
	unsigned int e_count = gridDim.x * blockDim.x;

	double x[8], y[8], z[8];

	z[0] = f[index]; index += e_count;
	z[1] = f[index]; index += e_count;
	z[2] = f[index]; index += e_count;
	z[3] = f[index]; index += e_count;
	z[4] = f[index]; index += e_count;
	z[5] = f[index]; index += e_count;
	z[6] = f[index]; index += e_count;
	z[7] = f[index]; index += e_count;

	// z difference vectors
    double R42=(z[3] - z[1]);
    double R52=(z[4] - z[1]);
    double R54=(z[4] - z[3]);

    double R63=(z[5] - z[2]);
    double R83=(z[7] - z[2]);
    double R86=(z[7] - z[5]);

    double R31=(z[2] - z[0]);
    double R61=(z[5] - z[0]);
    double R74=(z[6] - z[3]);

    double R72=(z[6] - z[1]);
    double R75=(z[6] - z[4]);
    double R81=(z[7] - z[0]);

    double t1=(R63 + R54);
    double t2=(R61 + R74);
    double t3=(R72 + R81);

    double t4 =(R86 + R42);
    double t5 =(R83 + R52);
    double t6 =(R75 + R31);

	y[0] = f[index]; index += e_count;
	y[1] = f[index]; index += e_count;
	y[2] = f[index]; index += e_count;
	y[3] = f[index]; index += e_count;
	y[4] = f[index]; index += e_count;
	y[5] = f[index]; index += e_count;
	y[6] = f[index]; index += e_count;
	y[7] = f[index]; index += e_count;

    g[index] = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5)  + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54); index += e_count;
    g[index] = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1)  - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61); index += e_count;
    g[index] = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2)  - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72); index += e_count;
    g[index] = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3)  + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83); index += e_count;
    g[index] = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2)  - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61); index += e_count;
    g[index] = (y[6] *  t5) - (y[4] *  t3)  - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72); index += e_count;
    g[index] = (y[7] *  t1) - (y[5] *  t5)  - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83); index += e_count;
    g[index] = (y[4] *  t2) - (y[6] *  t1)  + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54); 

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

	x[7] = f[index]; index -= e_count;
	x[6] = f[index]; index -= e_count;
	x[5] = f[index]; index -= e_count;
	x[4] = f[index]; index -= e_count;
	x[3] = f[index]; index -= e_count;
	x[2] = f[index]; index -= e_count;
	x[1] = f[index]; index -= e_count;
	x[0] = f[index]; index -= e_count;

    g[index] = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54); index -= e_count;
    g[index] = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61); index -= e_count;
    g[index] = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72); index -= e_count;
    g[index] = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83); index -= e_count;
    g[index] = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61); index -= e_count;
    g[index] = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72); index -= e_count;
    g[index] = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83); index -= e_count;
    g[index] = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54); index -= e_count;

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

    g[index]  = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54); index -= e_count;
    g[index]  = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61); index -= e_count;
    g[index]  = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72); index -= e_count;
    g[index]  = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83); index -= e_count;
    g[index]  = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61); index -= e_count;
    g[index]  = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72); index -= e_count;
    g[index]  = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83); index -= e_count;
    g[index]  = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);



}

__global__ void HexGrad3 (double * f, double * g){

//////////////////////////////////////////////////////////////////////////
//	second optimization of the original kernel							//
//																		//
//	this kernel features truly staggered reads, which improves 			//
//	performance by only performing memory reads immediately before 		//
//	the data is needed													//
//////////////////////////////////////////////////////////////////////////


	unsigned int index = blockIdx.x*(blockDim.x) + threadIdx.x;
	unsigned int e_count = gridDim.x * blockDim.x;

	double x[8], y[8], z[8];

	// z difference vectors

	z[0] = f[index]; index += 2 * e_count;
	z[2] = f[index]; index += 3 * e_count;
	z[5] = f[index]; index += 2 * e_count;

    double R31=(z[2] - z[0]);
    double R61=(z[5] - z[0]);
	double R63=(z[5] - z[2]);

	z[7] = f[index]; index -= 3 * e_count;

    double R83=(z[7] - z[2]);
    double R86=(z[7] - z[5]);
    double R81=(z[7] - z[0]);

	z[4] = f[index]; index -= e_count;
	z[3] = f[index]; index -= 2 * e_count;
	z[1] = f[index]; index += 5 * e_count;

    double R42=(z[3] - z[1]);
	double t4 =(R86 + R42);
    double R52=(z[4] - z[1]);
	double t5 =(R83 + R52);
    double R54=(z[4] - z[3]);
   	double t1=(R63 + R54);

	z[6] = f[index]; index += 2 * e_count;

    double R72=(z[6] - z[1]);
	double t3=(R72 + R81);
    double R75=(z[6] - z[4]);
	double t6 =(R75 + R31);
    double R74=(z[6] - z[3]);
	double t2=(R61 + R74);

	y[0] = f[index]; index += e_count;
	y[1] = f[index]; index += e_count;
	y[2] = f[index]; index += e_count;
	y[3] = f[index]; index += e_count;
	y[4] = f[index]; index += e_count;
	y[5] = f[index]; index += e_count;
	y[6] = f[index]; index += e_count;
	y[7] = f[index]; index += e_count;

    g[index] = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5)  + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54); index += e_count;
    g[index] = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1)  - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61); index += e_count;
    g[index] = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2)  - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72); index += e_count;
    g[index] = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3)  + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83); index += e_count;
    g[index] = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2)  - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61); index += e_count;
    g[index] = (y[6] *  t5) - (y[4] *  t3)  - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72); index += e_count;
    g[index] = (y[7] *  t1) - (y[5] *  t5)  - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83); index += e_count;
    g[index] = (y[4] *  t2) - (y[6] *  t1)  + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54); 

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

	x[7] = f[index]; index -= e_count;
	x[6] = f[index]; index -= e_count;
	x[5] = f[index]; index -= e_count;
	x[4] = f[index]; index -= e_count;
	x[3] = f[index]; index -= e_count;
	x[2] = f[index]; index -= e_count;
	x[1] = f[index]; index -= e_count;
	x[0] = f[index]; index -= e_count;

    g[index] = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54); index -= e_count;
    g[index] = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61); index -= e_count;
    g[index] = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72); index -= e_count;
    g[index] = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83); index -= e_count;
    g[index] = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61); index -= e_count;
    g[index] = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72); index -= e_count;
    g[index] = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83); index -= e_count;
    g[index] = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54); index -= e_count;

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

    g[index]  = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54); index -= e_count;
    g[index]  = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61); index -= e_count;
    g[index]  = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72); index -= e_count;
    g[index]  = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83); index -= e_count;
    g[index]  = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61); index -= e_count;
    g[index]  = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72); index -= e_count;
    g[index]  = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83); index -= e_count;
    g[index]  = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);

}

__global__ void HexGrad4 (double * f, double * g){

//////////////////////////////////////////////////////////////////////////
//	third optimization of the original kernel							//
//																		//
//	this kernel features staggered reads and a reduced memory			//
//	footprint, by using only 1 local array to store coordinate info		//
//////////////////////////////////////////////////////////////////////////

	unsigned int index = blockIdx.x*(blockDim.x) + threadIdx.x;
	unsigned int e_count = gridDim.x * blockDim.x;

	double a[8];

	// Z
	a[0] = f[index]; index += e_count;
	a[1] = f[index]; index += e_count;
	a[2] = f[index]; index += e_count;
	a[3] = f[index]; index += e_count;
	a[4] = f[index]; index += e_count;
	a[5] = f[index]; index += e_count;
	a[6] = f[index]; index += e_count;
	a[7] = f[index]; index += e_count;
 
	// z difference vectors
    double R42=(a[3] - a[1]);
    double R52=(a[4] - a[1]);
    double R54=(a[4] - a[3]);

    double R63=(a[5] - a[2]);
    double R83=(a[7] - a[2]);
    double R86=(a[7] - a[5]);

    double R31=(a[2] - a[0]);
    double R61=(a[5] - a[0]);
    double R74=(a[6] - a[3]);

    double R72=(a[6] - a[1]);
    double R75=(a[6] - a[4]);
    double R81=(a[7] - a[0]);

    double t1=(R63 + R54);
    double t2=(R61 + R74);
    double t3=(R72 + R81);

    double t4 =(R86 + R42);
    double t5 =(R83 + R52);
    double t6 =(R75 + R31);

	// Y

	a[0] = f[index]; index += e_count;
	a[1] = f[index]; index += e_count;
	a[2] = f[index]; index += e_count;
	a[3] = f[index]; index += e_count;
	a[4] = f[index]; index += e_count;
	a[5] = f[index]; index += e_count;
	a[6] = f[index]; index += e_count;
	a[7] = f[index]; index += e_count;

	// X grad
    g[index] = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54); index += e_count;
    g[index] = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61); index += e_count;
    g[index] = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72); index += e_count;
    g[index] = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83); index += e_count;
    g[index] = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61); index += e_count;
    g[index] = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72); index += e_count;
    g[index] = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83); index += e_count;
    g[index] = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54); index += e_count;

    R42=(a[3] - a[1]);
    R52=(a[4] - a[1]);
    R54=(a[4] - a[3]);

    R63=(a[5] - a[2]);
    R83=(a[7] - a[2]);
    R86=(a[7] - a[5]);

    R31=(a[2] - a[0]);
    R61=(a[5] - a[0]);
    R74=(a[6] - a[3]);

    R72=(a[6] - a[1]);
    R75=(a[6] - a[4]);
    R81=(a[7] - a[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

	// X

	a[7] = f[index]; index -= e_count;
	a[6] = f[index]; index -= e_count;
	a[5] = f[index]; index -= e_count;
	a[4] = f[index]; index -= e_count;
	a[3] = f[index]; index -= e_count;
	a[2] = f[index]; index -= e_count;
	a[1] = f[index]; index -= e_count;
	a[0] = f[index]; index -= e_count;

	index -= 8 * e_count;

	// Z grad
    g[index] = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54); index -= e_count;
    g[index] = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83); index -= e_count;
    g[index] = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72); index -= e_count;
    g[index] = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61); index -= e_count;
    g[index] = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83); index -= e_count;
    g[index] = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72); index -= e_count;
    g[index] = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61); index -= e_count;
    g[index] = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54); index -= e_count;
 
    R42=(a[3] - a[1]);
    R52=(a[4] - a[1]);
    R54=(a[4] - a[3]);

    R63=(a[5] - a[2]);
    R83=(a[7] - a[2]);
    R86=(a[7] - a[5]);

    R31=(a[2] - a[0]);
    R61=(a[5] - a[0]);
    R74=(a[6] - a[3]);

    R72=(a[6] - a[1]);
    R75=(a[6] - a[4]);
    R81=(a[7] - a[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

	//  Z

	a[0] = f[index]; index += e_count;
	a[1] = f[index]; index += e_count;
	a[2] = f[index]; index += e_count;
	a[3] = f[index]; index += e_count;
	a[4] = f[index]; index += e_count;
	a[5] = f[index]; index += e_count;
	a[6] = f[index]; index += e_count;
	a[7] = f[index]; index += e_count;

	//	Y grad

    g[index]  = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54); index += e_count;
    g[index]  = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61); index += e_count;
    g[index]  = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72); index += e_count;
    g[index]  = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83); index += e_count;
    g[index]  = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61); index += e_count;
    g[index]  = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72); index += e_count;
    g[index]  = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83); index += e_count;
    g[index]  = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54); index += e_count;

}

__global__ void HexGrad5 (double * f, double * g){

//////////////////////////////////////////////////////////////////////////
//	fourth optimization of the original kernel							//
//																		//
//	this kernel features truly staggered reads and a reduced memory		//
//	footprint, by using only 1 local array to store coordinate info		//
//////////////////////////////////////////////////////////////////////////

	unsigned int index = blockIdx.x*(blockDim.x) + threadIdx.x;
	unsigned int e_count = gridDim.x * blockDim.x;

	double a[8];

	a[0] = f[index]; index += 2 * e_count;
	a[2] = f[index]; index += 3 * e_count;
	a[5] = f[index]; index += 2 * e_count;

    double R31=(a[2] - a[0]);
    double R61=(a[5] - a[0]);
	double R63=(a[5] - a[2]);

	a[7] = f[index]; index -= 3 * e_count;

    double R83=(a[7] - a[2]);
    double R86=(a[7] - a[5]);
    double R81=(a[7] - a[0]);

	a[4] = f[index]; index -= e_count;
	a[3] = f[index]; index -= 2 * e_count;
	a[1] = f[index]; index += 5 * e_count;

    double R42=(a[3] - a[1]);
	double t4 =(R86 + R42);
    double R52=(a[4] - a[1]);
	double t5 =(R83 + R52);
    double R54=(a[4] - a[3]);
   	double t1=(R63 + R54);

	a[6] = f[index]; index += 2 * e_count;

    double R72=(a[6] - a[1]);
	double t3=(R72 + R81);
    double R75=(a[6] - a[4]);
	double t6 =(R75 + R31);
    double R74=(a[6] - a[3]);
	double t2=(R61 + R74);

	// Y

	a[0] = f[index]; index += e_count;
	a[1] = f[index]; index += e_count;
	a[2] = f[index]; index += e_count;
	a[3] = f[index]; index += e_count;
	a[4] = f[index]; index += e_count;
	a[5] = f[index]; index += e_count;
	a[6] = f[index]; index += e_count;
	a[7] = f[index]; index += e_count;

	// X grad
    g[index] = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54); index += e_count;
    g[index] = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61); index += e_count;
    g[index] = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72); index += e_count;
    g[index] = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83); index += e_count;
    g[index] = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61); index += e_count;
    g[index] = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72); index += e_count;
    g[index] = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83); index += e_count;
    g[index] = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54); index += e_count;

    R42=(a[3] - a[1]);
    R52=(a[4] - a[1]);
    R54=(a[4] - a[3]);

    R63=(a[5] - a[2]);
    R83=(a[7] - a[2]);
    R86=(a[7] - a[5]);

    R31=(a[2] - a[0]);
    R61=(a[5] - a[0]);
    R74=(a[6] - a[3]);

    R72=(a[6] - a[1]);
    R75=(a[6] - a[4]);
    R81=(a[7] - a[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

	// X

	a[7] = f[index]; index -= e_count;
	a[6] = f[index]; index -= e_count;
	a[5] = f[index]; index -= e_count;
	a[4] = f[index]; index -= e_count;
	a[3] = f[index]; index -= e_count;
	a[2] = f[index]; index -= e_count;
	a[1] = f[index]; index -= e_count;
	a[0] = f[index]; index -= e_count;

	index -= 8 * e_count;

	// Z grad
    g[index] = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54); index -= e_count;
    g[index] = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83); index -= e_count;
    g[index] = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72); index -= e_count;
    g[index] = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61); index -= e_count;
    g[index] = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83); index -= e_count;
    g[index] = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72); index -= e_count;
    g[index] = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61); index -= e_count;
    g[index] = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54); index -= e_count;
 
    R42=(a[3] - a[1]);
    R52=(a[4] - a[1]);
    R54=(a[4] - a[3]);

    R63=(a[5] - a[2]);
    R83=(a[7] - a[2]);
    R86=(a[7] - a[5]);

    R31=(a[2] - a[0]);
    R61=(a[5] - a[0]);
    R74=(a[6] - a[3]);

    R72=(a[6] - a[1]);
    R75=(a[6] - a[4]);
    R81=(a[7] - a[0]);

    t1=(R63 + R54);
    t2=(R61 + R74);
    t3=(R72 + R81);

    t4 =(R86 + R42);
    t5 =(R83 + R52);
    t6 =(R75 + R31);

	//  Z

	a[0] = f[index]; index += e_count;
	a[1] = f[index]; index += e_count;
	a[2] = f[index]; index += e_count;
	a[3] = f[index]; index += e_count;
	a[4] = f[index]; index += e_count;
	a[5] = f[index]; index += e_count;
	a[6] = f[index]; index += e_count;
	a[7] = f[index]; index += e_count;

	//	Y grad

    g[index]  = (a[1] *  t1) - (a[2] * R42) - (a[3] *  t5)  + (a[4] *  t4) + (a[5] * R52) - (a[7] * R54); index += e_count;
    g[index]  = (a[2] *  t2) + (a[3] * R31) - (a[0] *  t1)  - (a[5] *  t6) + (a[6] * R63) - (a[4] * R61); index += e_count;
    g[index]  = (a[3] *  t3) + (a[0] * R42) - (a[1] *  t2)  - (a[6] *  t4) + (a[7] * R74) - (a[5] * R72); index += e_count;
    g[index]  = (a[0] *  t5) - (a[1] * R31) - (a[2] *  t3)  + (a[7] *  t6) + (a[4] * R81) - (a[6] * R83); index += e_count;
    g[index]  = (a[5] *  t3) + (a[6] * R86) - (a[7] *  t2)  - (a[0] *  t4) - (a[3] * R81) + (a[1] * R61); index += e_count;
    g[index]  = (a[6] *  t5) - (a[4] *  t3)  - (a[7] * R75) + (a[1] *  t6) - (a[0] * R52) + (a[2] * R72); index += e_count;
    g[index]  = (a[7] *  t1) - (a[5] *  t5)  - (a[4] * R86) + (a[2] *  t4) - (a[1] * R63) + (a[3] * R83); index += e_count;
    g[index]  = (a[4] *  t2) - (a[6] *  t1)  + (a[5] * R75) - (a[3] *  t6) - (a[2] * R74) + (a[0] * R54); index += e_count;

}

__global__ void HexGrad6 (	double * f, 
							double * g){


//////////////////////////////////////////////////////////////////////////
//	fifth optimization of the original kernel							//
//																		//
//	this kernel further reduces the memory footprint by eliminating 	//
//	stored intermediates. Instead, this kernel opts to calculate		//
//	all needed values on the fly.										//
//////////////////////////////////////////////////////////////////////////

	unsigned int index = blockIdx.x*(blockDim.x) + threadIdx.x;
	unsigned int e_count = gridDim.x * blockDim.x;

	double a[8], b[8];

	// Z
	a[0] = f[index]; index += e_count;
	a[1] = f[index]; index += e_count;
	a[2] = f[index]; index += e_count;
	a[3] = f[index]; index += e_count;
	a[4] = f[index]; index += e_count;
	a[5] = f[index]; index += e_count;
	a[6] = f[index]; index += e_count;
	a[7] = f[index]; index += e_count;

	// Y
	b[0] = f[index]; index += e_count;
	b[1] = f[index]; index += e_count;
	b[2] = f[index]; index += e_count;
	b[3] = f[index]; index += e_count;
	b[4] = f[index]; index += e_count;
	b[5] = f[index]; index += e_count;
	b[6] = f[index]; index += e_count;
	b[7] = f[index]; index += e_count;
 
	// X grad
    g[index] = (b[1] *  ((a[5] - a[2]) + (a[4] - a[3]))) - (b[2] * (a[3] - a[1])) - (b[3] *  ((a[7] - a[2]) + (a[4] - a[1]))) \
			 + (b[4] *  ((a[7] - a[5]) + (a[3] - a[1]))) + (b[5] * (a[4] - a[1])) - (b[7] * (a[4] - a[3])); index += e_count;
    g[index] = (b[2] *  ((a[5] - a[0]) + (a[6] - a[3]))) + (b[3] * (a[2] - a[0])) - (b[0] *  ((a[5] - a[2]) + (a[4] - a[3]))) \
			 - (b[5] *  ((a[6] - a[4]) + (a[2] - a[0]))) + (b[6] * (a[5] - a[2])) - (b[4] * (a[5] - a[0])); index += e_count;
    g[index] = (b[3] *  ((a[6] - a[1]) + (a[7] - a[0]))) + (b[0] * (a[3] - a[1])) - (b[1] *  ((a[5] - a[0]) + (a[6] - a[3]))) \
			 - (b[6] *  ((a[7] - a[5]) + (a[3] - a[1]))) + (b[7] * (a[6] - a[3])) - (b[5] * (a[6] - a[1])); index += e_count;
    g[index] = (b[0] *  ((a[7] - a[2]) + (a[4] - a[1]))) - (b[1] * (a[2] - a[0])) - (b[2] *  ((a[6] - a[1]) + (a[7] - a[0]))) \
			 + (b[7] *  ((a[6] - a[4]) + (a[2] - a[0]))) + (b[4] * (a[7] - a[0])) - (b[6] * (a[7] - a[2])); index += e_count;
    g[index] = (b[5] *  ((a[6] - a[1]) + (a[7] - a[0]))) + (b[6] * (a[7] - a[5])) - (b[7] *  ((a[5] - a[0]) + (a[6] - a[3]))) \
			 - (b[0] *  ((a[7] - a[5]) + (a[3] - a[1]))) - (b[3] * (a[7] - a[0])) + (b[1] * (a[5] - a[0])); index += e_count;
    g[index] = (b[6] *  ((a[7] - a[2]) + (a[4] - a[1]))) - (b[4] * ((a[6] - a[1]) + (b[7] - a[0])))  - (a[7] * (a[6] - a[4]))\
			 + (b[1] *  ((a[6] - a[4]) + (a[2] - a[0]))) - (b[0] * (a[4] - a[1])) + (b[2] * (a[6] - a[1])); index += e_count;
    g[index] = (b[7] *  ((a[5] - a[2]) + (a[4] - a[3]))) - (b[5] * ((a[7] - a[2]) + (b[4] - a[1])))  - (a[4] * (a[7] - a[5]))\
			 + (b[2] *  ((a[7] - a[5]) + (a[3] - a[1]))) - (b[1] * (a[5] - a[2])) + (b[3] * (a[7] - a[2])); index += e_count;
    g[index] = (b[4] *  ((a[5] - a[0]) + (a[6] - a[3]))) - (b[6] * ((a[5] - a[2]) + (b[4] - a[3])))  + (a[5] * (a[6] - a[4]))\
			 - (b[3] *  ((a[6] - a[4]) + (a[2] - a[0]))) - (b[2] * (a[6] - a[3])) + (b[0] * (a[4] - a[3])); index += e_count;

	// X

	a[7] = f[index]; index -= e_count;
	a[6] = f[index]; index -= e_count;
	a[5] = f[index]; index -= e_count;
	a[4] = f[index]; index -= e_count;
	a[3] = f[index]; index -= e_count;
	a[2] = f[index]; index -= e_count;
	a[1] = f[index]; index -= e_count;
	a[0] = f[index]; index -= e_count;

	index -= 8 * e_count;

	// Z grad
    g[index] = (a[4] *  ((b[5] - b[0]) + (b[6] - b[3]))) - (a[6] *  ((b[5] - b[2]) + (b[4] - b[3]))) + (a[5] * (b[6] - b[4])) \
			 - (a[3] *  ((b[6] - b[4]) + (b[2] - b[0]))) - (a[2] * (b[6] - b[3])) + (a[0] * (b[4] - b[3])); index -= e_count;
    g[index] = (a[7] *  ((b[5] - b[2]) + (b[4] - b[3]))) - (a[5] *  ((b[7] - b[2]) + (b[4] - b[1])))  - (a[4] * (b[7] - b[5]))\
			 + (a[2] *  ((b[7] - b[5]) + (b[3] - b[1]))) - (a[1] * (b[5] - b[2])) + (a[3] * (b[7] - b[2])); index -= e_count;
    g[index] = (a[6] *  ((b[7] - b[2]) + (b[4] - b[1]))) - (a[4] *  ((b[6] - b[1]) + (b[7] - b[0])))  - (a[7] * (b[6] - b[4]))\
			 + (a[1] *  ((b[6] - b[4]) + (b[2] - b[0]))) - (a[0] * (b[4] - b[1])) + (a[2] * (b[6] - b[1])); index -= e_count;
    g[index] = (a[5] *  ((b[6] - b[1]) + (b[7] - b[0]))) + (a[6] * (b[7] - b[5])) - (a[7] *  ((b[5] - b[0]) + (b[6] - b[3]))) \
			 - (a[0] *  ((b[7] - b[5]) + (b[3] - b[1]))) - (a[3] * (b[7] - b[0])) + (a[1] * (b[5] - b[0])); index -= e_count;
    g[index] = (a[0] *  ((b[7] - b[2]) + (b[4] - b[1]))) - (a[1] * (b[2] - b[0])) - (a[2] *  ((b[6] - b[1]) + (b[7] - b[0]))) \
			 + (a[7] *  ((b[6] - b[4]) + (b[2] - b[0]))) + (a[4] * (b[7] - b[0])) - (a[6] * (b[7] - b[2])); index -= e_count;
    g[index] = (a[3] *  ((b[6] - b[1]) + (b[7] - b[0]))) + (a[0] * (b[3] - b[1])) - (a[1] *  ((b[5] - b[0]) + (b[6] - b[3]))) \
			 - (a[6] *  ((b[7] - b[5]) + (b[3] - b[1]))) + (a[7] * (b[6] - b[3])) - (a[5] * (b[6] - b[1])); index -= e_count;
    g[index] = (a[2] *  ((b[5] - b[0]) + (b[6] - b[3]))) + (a[3] * (b[2] - b[0])) - (a[0] *  ((b[5] - b[2]) + (b[4] - b[3]))) \
			 - (a[5] *  ((b[6] - b[4]) + (b[2] - b[0]))) + (a[6] * (b[5] - b[2])) - (a[4] * (b[5] - b[0])); index -= e_count;
    g[index] = (a[1] *  ((b[5] - b[2]) + (b[4] - b[3]))) - (a[2] * (b[3] - b[1])) - (a[3] *  ((b[7] - b[2]) + (b[4] - b[1]))) \
			 + (a[4] *  ((b[7] - b[5]) + (b[3] - b[1]))) + (a[5] * (b[4] - b[1])) - (a[7] * (b[4] - b[3])); index -= e_count;
 
	// Z

	b[0] = f[index]; index += e_count;
	b[1] = f[index]; index += e_count;
	b[2] = f[index]; index += e_count;
	b[3] = f[index]; index += e_count;
	b[4] = f[index]; index += e_count;
	b[5] = f[index]; index += e_count;
	b[6] = f[index]; index += e_count;
	b[7] = f[index]; index += e_count;

	//	Y grad

    g[index] = (a[1] *  ((a[5] - a[2]) + (a[4] - a[3]))) - (a[2] * (a[3] - a[1])) - (a[3] *  ((a[7] - a[2]) + (a[4] - a[1]))) \
			 + (a[4] *  ((a[7] - a[5]) + (a[3] - a[1]))) + (a[5] * (a[4] - a[1])) - (a[7] * (a[4] - a[3])); index += e_count;
    g[index] = (a[2] *  ((a[5] - a[0]) + (a[6] - a[3]))) + (a[3] * (a[2] - a[0])) - (a[0] *  ((a[5] - a[2]) + (a[4] - a[3]))) \
			 - (a[5] *  ((a[6] - a[4]) + (a[2] - a[0]))) + (a[6] * (a[5] - a[2])) - (a[4] * (a[5] - a[0])); index += e_count;
    g[index] = (a[3] *  ((a[6] - a[1]) + (a[7] - a[0]))) + (a[0] * (a[3] - a[1])) - (a[1] *  ((a[5] - a[0]) + (a[6] - a[3]))) \
			 - (a[6] *  ((a[7] - a[5]) + (a[3] - a[1]))) + (a[7] * (a[6] - a[3])) - (a[5] * (a[6] - a[1])); index += e_count;
    g[index] = (a[0] *  ((a[7] - a[2]) + (a[4] - a[1]))) - (a[1] * (a[2] - a[0])) - (a[2] *  ((a[6] - a[1]) + (a[7] - a[0]))) \
			 + (a[7] *  ((a[6] - a[4]) + (a[2] - a[0]))) + (a[4] * (a[7] - a[0])) - (a[6] * (a[7] - a[2])); index += e_count;
    g[index] = (a[5] *  ((a[6] - a[1]) + (a[7] - a[0]))) + (a[6] * (a[7] - a[5])) - (a[7] *  ((a[5] - a[0]) + (a[6] - a[3]))) \
			 - (a[0] *  ((a[7] - a[5]) + (a[3] - a[1]))) - (a[3] * (a[7] - a[0])) + (a[1] * (a[5] - a[0])); index += e_count;
    g[index] = (a[6] *  ((a[7] - a[2]) + (a[4] - a[1]))) - (a[4] * ((a[6] - a[1]) + (a[7] - a[0])))  - (a[7] * (a[6] - a[4]))\
			 + (a[1] *  ((a[6] - a[4]) + (a[2] - a[0]))) - (a[0] * (a[4] - a[1])) + (a[2] * (a[6] - a[1])); index += e_count;
    g[index] = (a[7] *  ((a[5] - a[2]) + (a[4] - a[3]))) - (a[5] * ((a[7] - a[2]) + (a[4] - a[1])))  - (a[4] * (a[7] - a[5]))\
			 + (a[2] *  ((a[7] - a[5]) + (a[3] - a[1]))) - (a[1] * (a[5] - a[2])) + (a[3] * (a[7] - a[2])); index += e_count;
    g[index] = (a[4] *  ((a[5] - a[0]) + (a[6] - a[3]))) - (a[6] * ((a[5] - a[2]) + (a[4] - a[3])))  + (a[5] * (a[6] - a[4]))\
			 - (a[3] *  ((a[6] - a[4]) + (a[2] - a[0]))) - (a[2] * (a[6] - a[3])) + (a[0] * (a[4] - a[3])); index += e_count;


}

__global__ void scale (	double * input, 
							double * scalar_array, 
							const int location, 
							const int length){

//	scale a vector of size "length" by a scalar found at scalar_array[location]
//	
//	note: each thread services multiple components within the vector.

	unsigned int i = blockIdx.x*(blockDim.x) + threadIdx.x;
	unsigned int gridSize = blockDim.x*gridDim.x;

	double scalar = scalar_array[location];

	while(i < length){

		input[i] = input[i] * scalar;
		i += gridSize;	
	
	}

}

template <bool add>
__global__ void scale_subtract (	double * minuend, 
									double * subtrahend, 
									double * scalar_array, 
									const int location, 
									const int length){

//	scale a vector (subtrahend) of size "length" by a scalar found at scalar_array[location],
//	then add/subtract that vector from minuend
//	note: each thread services multiple components within the vector.

	unsigned int i = blockIdx.x*(blockDim.x) + threadIdx.x;
	unsigned int gridSize = blockDim.x*gridDim.x;
	double scalar = scalar_array[location];

	while(i < length){

		if(add)
			minuend[i] += scalar * subtrahend[i];
		if(!add)
			minuend[i] -= scalar * subtrahend[i];

		i += gridSize;	
	
	}

}

template <unsigned int blockSize>
__global__ void reduction (	double * input, 
								double * scratch, 
								const int length){

//	Brent's theorem optimized parallel reduction mechanism,
//	courtesy of NVidia GPU toolkit

    extern __shared__ double sdata[];

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;
    sdata[tid] = 0.0f;

    // reduce vector into blocks and shared memory
    while (i < length)
    {
        sdata[tid] += input[i] * input[i] + input[i + blockSize] * input[i + blockSize];  
        i += gridSize;
    } 
    __syncthreads();

    // do reduction of shared memory
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    
    if (tid < 32)
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; __syncthreads(); }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; __syncthreads(); }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; __syncthreads(); }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; __syncthreads(); }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; __syncthreads(); }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; __syncthreads(); }
    }
    
    // write result for this block to global mem 
    if (tid == 0) scratch[blockIdx.x] = sdata[0];

	// note:this kernel does not reduce the vectors to a single quantity,
	// 		but instead writes N partial sums to scratch, where N is given 
	//		by the number of blocks.

}

template <unsigned int blockSize>
__global__ void dot_product (	double * a, 
								double * b, 	
								double * scratch, 
								const int length){

//	Brent's theorem optimized parallel reduction mechanism,
//	courtesy of NVidia GPU toolkit, modified to compute a vector
//	dot product, instead of an arithmetic sum

    extern __shared__ double sdata[];

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;
    sdata[tid] = 0;

    // reduce vector into blocks and shared memory

    while (i < length)
    {
        sdata[tid] += a[i] * b[i] + a[i + blockSize] * b[i + blockSize];  
        i += gridSize;
    } 
    __syncthreads();

    // do reduction of shared memory
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    
    if (tid < 32)
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; __syncthreads(); }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; __syncthreads(); }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; __syncthreads(); }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; __syncthreads(); }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; __syncthreads(); }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; __syncthreads(); }
    }
    
    // write result for this block to global mem 
    if (tid == 0) scratch[blockIdx.x] = sdata[0];


	// note:this kernel does not reduce the vectors to a single quantity,
	// 		but instead writes N partial sums to scratch, where N is given 
	//		by the number of blocks.

}

template <bool root>
__global__ void block_reduction (double * scratch, double * R, const int location){

//	this is a kernel that aims to complete the reduction started in the 
//	dot_product and reduction kernels. 


    extern __shared__ double sdata[];

	uint tid = threadIdx.x;

//	read in scratch data to shared memory space
	sdata[tid] = scratch[tid];

	for (int offset = blockDim.x / 2; offset > 0; offset /= 2){

//		perform reduction in shared memory, to a single value
		if(tid < offset){
			sdata[tid] += sdata[tid + offset];
		}

		__syncthreads();

	}


//	finally, return the final calculated value
	if (tid == 0){
		if(root){
			R[location] = sdata[0];
			scratch[0] = 1.0f / sqrt(sdata[0]);
		}		
		if(!root)
			scratch[0] = sdata[0];
	}

}

double serial_reduction(double * array, int n){

//	host-side implementation of a serial reduction
	double temp = 0;

	for (int i = 0; i < n; i++)
		temp += array[i];
	
	std::cout << temp << std::endl;

	return temp;
}

int CheckCudaError(cudaError_t const err){

    if(err)
        printf("CUDA Error: %s\nExiting...\n", cudaGetErrorString(err));
    return err;

}

void GpuProperties(){

    int device;
    cudaDeviceProp prop;
    cudaError_t deviceError, propertiesError;
    
    deviceError = cudaGetDevice(&device);
    propertiesError = cudaGetDeviceProperties(&prop, device);
    
    if (deviceError || propertiesError)
    {
        printf("Error detecting CUDA-capable GPU\n");
    }
    else
    {
        printf("\n\tGPU Information:\n");
        printf("\tDevice: %s\n", prop.name);
        printf("\tCompute Capability: %d.%d\n", prop.major, prop.minor);
        printf("\tMemory (MB): %d\n\n", prop.totalGlobalMem >> 20);
    }
}

void test_gram_schmidt(){

	cudaError_t err = cudaSuccess;
	timeval start,stop,result;

	std::cout << "\tGram Schmidt algorithm: \n";
	std::cout << "Work Size     Time(s)     Time/Work\n\n";

//////////////////////////////////////////////////
//	reduction and future gram schmit kernels	//
//	i = number of vectors						//
//	j = vector length							//
//////////////////////////////////////////////////


	for(int i = 32; i < 33; i += 3){
		for(int j = 1024; j < 9000000; j *= 2){
	
			int threads = 256;
			int blocks = 64;

			double *Q[i], *R[i], * Q_D[i], * R_D[i];
			double *zeroes, *scratch;
					
			int Q_vector_size = sizeof(double) * (j + threads);
			int R_vector_size = sizeof(double) * (i);
			int S_vector_size = sizeof(double) * blocks;

			int shared_size = sizeof(double) * threads;

		//	init Q
			for (int a = 0; a < i; a++){

				double * temp = (double*)(malloc(Q_vector_size));
				
				for(int b = 0; b < j + threads; b++){


					temp[b] = (rand() % 1000) / (250.0 * M_PI);	
					
					if(a == 0)
						temp[b] = 1337;
					if(a == 1)
						temp[b] = 2;
					
				//	To use the Brent's Theorem optimized
				// 	reduction kernel, add a thread-sized 
				// 	buffer to the end of each vector.

					if (b >= j)
						temp[b] = 0.0;

				}

				Q[a] = temp;

			}

		//	init R
			for (int a = 0; a < i; a++){

				double * temp = (double*)(malloc(R_vector_size));

				for(int b = 0; b < i; b++){

					temp[b] = (rand() % 1000) / (250.0 * M_PI);
	
				}
					
				R[a] = temp;
				
			}

		//	init zeroes/scratch information

			zeroes = (double*)malloc(S_vector_size);

			for (int a = 0; a < blocks; a++){
				zeroes[a] = 0;
			}


		//	write multivectors to device		
			for(int a = 0; a < i; a++){
				err = cudaMalloc((void**)&Q_D[a], Q_vector_size);
				CheckCudaError(err);
				err = cudaMalloc((void**)&R_D[a], R_vector_size);
				CheckCudaError(err);
			}

			err = cudaMalloc((void**)&scratch, S_vector_size);
			CheckCudaError(err);

		//	copy information into device memory
			for(int a = 0; a < i; a++){
				err = cudaMemcpy(Q_D[a], Q[a], Q_vector_size, cudaMemcpyHostToDevice);
				CheckCudaError(err);
				err = cudaMemcpy(R_D[a], R[a], R_vector_size, cudaMemcpyHostToDevice);
				CheckCudaError(err);
			}

			err = cudaMemcpy(scratch, zeroes, S_vector_size, cudaMemcpyHostToDevice);
			CheckCudaError(err);


			gettimeofday(&start, NULL);

		//	perform modified Gram-Schmidt Orthanormalization
			for(int m = 0; m < i; ++m){

      		// 	Reduction   : tmp = dot( Q(:,j) , Q(:,j) );
      		// 	PostProcess : tmp = sqrt( tmp ); R(j,j) = tmp ; tmp = 1 / tmp ;
				dot_product<256><<<blocks, threads, shared_size>>>(Q_D[m], Q_D[m], scratch, j);
				block_reduction<true><<<1, blocks, shared_size>>>(scratch, R_D[m], m);

			// 	Q(:,j) *= ( 1 / R(j,j) ); => Q(:,j) *= tmp ;
				scale <<<blocks, threads>>>(Q_D[m], scratch, 0, j);


				for(int n = m + 1; n < i; ++n){

      	 		// 	Reduction   : R(j,k) = dot( Q(:,j) , Q(:,k) );
        		// 	PostProcess : tmp = - R(j,k);
        			dot_product<256><<<blocks, threads, shared_size>>>(Q_D[m], Q_D[n], scratch, j);
					block_reduction<false><<<1, blocks, shared_size>>>(scratch, R_D[m], n);

        		// 	Q(:,k) -= R(j,k) * Q(:,j); => Q(:,k) += tmp * Q(:,j)
					scale_subtract<true><<<blocks, threads>>>(Q_D[n], Q_D[m], R_D[m], n, j);

				}

			}

			cudaThreadSynchronize();

			gettimeofday(&stop, NULL);
			timersub(&stop, &start, &result);
			double time = (result.tv_sec + result.tv_usec/1000000.0);
			
			std::cout << std::setw(9) <<  j << ",  " << std::setw(9) <<  time << ",  " << std::setw(9) << time / j << std::endl;


		//	cleanup				
			for (int a = 0; a < i; a++){

				free(Q[a]);
				free(R[a]);
				cudaFree(Q_D[a]);
				cudaFree(R_D[a]);

			}	

			free(zeroes);
			cudaFree(scratch);

		}
	}

	std::cout << std::endl;


}

void test_hex_grad(){

	timeval start,stop,result;
	cudaError_t err = cudaSuccess;

	std::cout << "\tHexahedral Gradient: \n";
	std::cout << "Work Size     Time(s)     Time/Work\n\n";


	for(int i = 1024; i < 9000000; i*=2){

	//	work allocation parameters=-
		int threads = 256;
		int blocks = i / threads;

	//	sizes and data structures to be passed into kernel
	//	3 components * 8 coordinates = 24 doubles / element
		uint data_size = sizeof(double) * 24 * i;
		
		double * f = (double*)malloc(data_size);
		double * fD, * gD;

	//	allocate memory on the device
		err = cudaMalloc((void**)&fD, data_size);
		CheckCudaError(err);

		err = cudaMalloc((void**)&gD, data_size);
		CheckCudaError(err);

	//	copy dummy variables to the device
		err = cudaMemcpy(fD, f, data_size, cudaMemcpyHostToDevice);
		CheckCudaError(err);

		gettimeofday(&start, NULL);
		
		HexGrad <<<blocks, threads>>>(fD, gD);
	//	HexGrad2<<<blocks, threads>>>(fD, gD);
	//	HexGrad3<<<blocks, threads>>>(fD, gD);
	//	HexGrad4<<<blocks, threads>>>(fD, gD);
	//	HexGrad5<<<blocks, threads>>>(fD, gD);
	//	HexGrad6<<<blocks, threads>>>(fD, gD);

		cudaThreadSynchronize();

		gettimeofday(&stop, NULL);
		timersub(&stop, &start, &result);
		double time = (result.tv_sec + result.tv_usec/1000000.0);

		err = cudaMemcpy(f, fD, data_size, cudaMemcpyDeviceToHost);
		CheckCudaError(err);

	//	display work size, time, and time/work information

		std::cout << std::setw(9) << i << ",  " << std::setw(9) <<  time << ",  " << std::setw(9) << time / i << std::endl;

	//	clean up
		free (f);
		cudaFree(fD);
		cudaFree(gD);
	}

}

int main(){

//	print out device information
	GpuProperties();

	test_gram_schmidt();
	test_hex_grad();

    return 0;
 }
