/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <sys/time.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>
//
//  Defined numerical constants
//
#define _M * stride
#define _N * blockDim.x
//
//
//
const double ONE3RD =  (1.0 /  3.0);
const double ONE12TH = (1.0 / 12.0);
//
//  Indexes into a 4 by 3 matrix of hourglass constants
//
const int HG_X1 = 0;
const int HG_Y1 = 1;
const int HG_Z1 = 2;
const int HG_X2 = 3;
const int HG_Y2 = 4;
const int HG_Z2 = 5;
const int HG_X3 = 6;
const int HG_Y3 = 7;
const int HG_Z3 = 8;
const int HG_X4 = 9;
const int HG_Y4 = 10;
const int HG_Z4 = 11;
//
//  Indexes into a 3 by 3 symmetric tensor stored as a length 6 vector
//
const int K_S_XX = 0;
const int K_S_YY = 1;
const int K_S_ZZ = 2;
const int K_S_XY = 3;
const int K_S_YZ = 4;
const int K_S_ZX = 5;
const int K_S_YX = 3;
const int K_S_ZY = 4;
const int K_S_XZ = 5;
//
//  Indexes into a full 3 by 3 tensor stored as a length 9 vector
//
const int K_F_XX = 0;
const int K_F_YY = 1;
const int K_F_ZZ = 2;
const int K_F_XY = 3;
const int K_F_YZ = 4;
const int K_F_ZX = 5;
const int K_F_YX = 6;
const int K_F_ZY = 7;
const int K_F_XZ = 8;
//
//  Indexes into a 3 by 3 skew symmetric tensor stored as a length 9 vector
//
const int K_V_XY = 0;
const int K_V_YZ = 1;
const int K_V_ZX = 2;


/********************************************************************/

int CheckCudaError(cudaError_t const err){

    if(err)
        printf("CUDA Error: %s\nExiting...\n", cudaGetErrorString(err));
    return err;

}

/********************************************************************/

__global__ void cuda_fill(double * position, double * velocity, double * rotation){

	int stride = blockDim.x;
	int index_24 = blockIdx.x * blockDim.x * 24 + threadIdx.x;
	int index_9  = blockIdx.x * blockDim.x *  9 + threadIdx.x;

	//	unit cube in first octant, first node at origin

	//	X coords
		position[index_24 +  0 _M] = 0.0;
		position[index_24 +  1 _M] = 1.0;
		position[index_24 +  2 _M] = 1.0;
		position[index_24 +  3 _M] = 0.0;
		position[index_24 +  4 _M] = 0.0;
		position[index_24 +  5 _M] = 1.0;
		position[index_24 +  6 _M] = 1.0;
		position[index_24 +  7 _M] = 0.0;

	//	Y coords
		position[index_24 +  8 _M] = 0.0;
		position[index_24 +  9 _M] = 0.0;
		position[index_24 + 10 _M] = 1.0;
		position[index_24 + 11 _M] = 1.0;
		position[index_24 + 12 _M] = 0.0;
		position[index_24 + 13 _M] = 0.0;
		position[index_24 + 14 _M] = 1.0;
		position[index_24 + 15 _M] = 1.0;

	//	Z coords
		position[index_24 + 16 _M] = 0.0;
		position[index_24 + 17 _M] = 0.0;
		position[index_24 + 18 _M] = 0.0;
		position[index_24 + 19 _M] = 0.0;
		position[index_24 + 20 _M] = 1.0;
		position[index_24 + 21 _M] = 1.0;
		position[index_24 + 22 _M] = 1.0;
		position[index_24 + 23 _M] = 1.0;


	//	first node with velocity in the direction of
	//	the cube's centroid, all other nodes static

	//	X vel
		velocity[index_24 +  0 _M] = 0.25;	
		velocity[index_24 +  1 _M] = 0.0;	
		velocity[index_24 +  2 _M] = 0.0;	
		velocity[index_24 +  3 _M] = 0.0;	
		velocity[index_24 +  4 _M] = 0.0;	
		velocity[index_24 +  5 _M] = 0.0;	
		velocity[index_24 +  6 _M] = 0.0;	
		velocity[index_24 +  7 _M] = 0.0;	

	//	Y vel
		velocity[index_24 +  8 _M] = 0.25;	
		velocity[index_24 +  9 _M] = 0.0;	
		velocity[index_24 + 10 _M] = 0.0;	
		velocity[index_24 + 11 _M] = 0.0;	
		velocity[index_24 + 12 _M] = 0.0;	
		velocity[index_24 + 13 _M] = 0.0;	
		velocity[index_24 + 14 _M] = 0.0;	
		velocity[index_24 + 15 _M] = 0.0;

	//	Z vel
		velocity[index_24 + 16 _M] = 0.25;	
		velocity[index_24 + 17 _M] = 0.0;	
		velocity[index_24 + 18 _M] = 0.0;	
		velocity[index_24 + 19 _M] = 0.0;	
		velocity[index_24 + 20 _M] = 0.0;	
		velocity[index_24 + 21 _M] = 0.0;	
		velocity[index_24 + 22 _M] = 0.0;	
		velocity[index_24 + 23 _M] = 0.0;



	//	No rotations
		rotation[index_9 + 0 _M] = 0.0;	
		rotation[index_9 + 1 _M] = 0.0;	
		rotation[index_9 + 2 _M] = 0.0;	
		rotation[index_9 + 3 _M] = 0.0;	
		rotation[index_9 + 4 _M] = 0.0;	
		rotation[index_9 + 5 _M] = 0.0;	
		rotation[index_9 + 6 _M] = 0.0;	
		rotation[index_9 + 7 _M] = 0.0;
		rotation[index_9 + 8 _M] = 0.0;


}



__global__ void cuda_comp_grad(double * pos, double * vel, double * m_pos, double * grad, double * m_vol, double dt){

	int disp_index = blockIdx.x * blockDim.x * 24 + threadIdx.x;

	double * m_pos_ptr = &m_pos[disp_index];
	double * pos_ptr = &pos[disp_index];
	double * vel_ptr = &vel[disp_index];
	double * grad_ptr = &grad[disp_index];
	
	double x[8], y[8], z[8];
	double delta_t = -0.5 * dt;

	m_pos_ptr[ 0 _N] = x[0] = vel_ptr[ 0 _N] * delta_t + pos_ptr[ 0 _N];
	m_pos_ptr[ 1 _N] = x[1] = vel_ptr[ 1 _N] * delta_t + pos_ptr[ 1 _N];
	m_pos_ptr[ 2 _N] = x[2] = vel_ptr[ 2 _N] * delta_t + pos_ptr[ 2 _N];
	m_pos_ptr[ 3 _N] = x[3] = vel_ptr[ 3 _N] * delta_t + pos_ptr[ 3 _N];
	m_pos_ptr[ 4 _N] = x[4] = vel_ptr[ 4 _N] * delta_t + pos_ptr[ 4 _N];
	m_pos_ptr[ 5 _N] = x[5] = vel_ptr[ 5 _N] * delta_t + pos_ptr[ 5 _N];
	m_pos_ptr[ 6 _N] = x[6] = vel_ptr[ 6 _N] * delta_t + pos_ptr[ 6 _N];
	m_pos_ptr[ 7 _N] = x[7] = vel_ptr[ 7 _N] * delta_t + pos_ptr[ 7 _N];

	m_pos_ptr[ 8 _N] = y[0] = vel_ptr[ 8 _N] * delta_t + pos_ptr[ 8 _N];
	m_pos_ptr[ 9 _N] = y[1] = vel_ptr[ 9 _N] * delta_t + pos_ptr[ 9 _N];
	m_pos_ptr[10 _N] = y[2] = vel_ptr[10 _N] * delta_t + pos_ptr[10 _N];
	m_pos_ptr[11 _N] = y[3] = vel_ptr[11 _N] * delta_t + pos_ptr[11 _N];
	m_pos_ptr[12 _N] = y[4] = vel_ptr[12 _N] * delta_t + pos_ptr[12 _N];
	m_pos_ptr[13 _N] = y[5] = vel_ptr[13 _N] * delta_t + pos_ptr[13 _N];
	m_pos_ptr[14 _N] = y[6] = vel_ptr[14 _N] * delta_t + pos_ptr[14 _N];
	m_pos_ptr[15 _N] = y[7] = vel_ptr[15 _N] * delta_t + pos_ptr[15 _N];

  	double t1 =((y[5] - y[2]) + (y[4] - y[3]));
  	double t2 =((y[5] - y[0]) + (y[6] - y[3]));
  	double t3 =((y[6] - y[1]) + (y[7] - y[0]));
  	double t4 =((y[7] - y[5]) + (y[3] - y[1]));
  	double t5 =((y[7] - y[2]) + (y[4] - y[1]));
  	double t6 =((y[6] - y[4]) + (y[2] - y[0]));

//	Z grad

    grad_ptr[16 _N] = (x[1] *  t1) - (x[2] * (y[3] - y[1])) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * (y[4] - y[1])) - (x[7] * (y[4] - y[3]));
    grad_ptr[17 _N] = (x[2] *  t2) + (x[3] * (y[2] - y[0])) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * (y[5] - y[2])) - (x[4] * (y[5] - y[0]));
    grad_ptr[18 _N] = (x[3] *  t3) + (x[0] * (y[3] - y[1])) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * (y[6] - y[3])) - (x[5] * (y[6] - y[1]));
    grad_ptr[19 _N] = (x[0] *  t5) - (x[1] * (y[2] - y[0])) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * (y[7] - y[0])) - (x[6] * (y[7] - y[2])); 
    grad_ptr[20 _N] = (x[5] *  t3) + (x[6] * (y[7] - y[5])) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * (y[7] - y[0])) + (x[1] * (y[5] - y[0]));
    grad_ptr[21 _N] = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * (y[6] - y[4])) + (x[1] *  t6) - (x[0] * (y[4] - y[1])) + (x[2] * (y[6] - y[1])); 
    grad_ptr[22 _N] = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * (y[7] - y[5])) + (x[2] *  t4) - (x[1] * (y[5] - y[2])) + (x[3] * (y[7] - y[2])); 
    grad_ptr[23 _N] = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * (y[6] - y[4])) - (x[3] *  t6) - (x[2] * (y[6] - y[3])) + (x[0] * (y[4] - y[3]));

	m_pos_ptr[16 _N] = z[0] = vel_ptr[16 _N] * delta_t + pos_ptr[16 _N];
	m_pos_ptr[17 _N] = z[1] = vel_ptr[17 _N] * delta_t + pos_ptr[17 _N];
	m_pos_ptr[18 _N] = z[2] = vel_ptr[18 _N] * delta_t + pos_ptr[18 _N];
	m_pos_ptr[19 _N] = z[3] = vel_ptr[19 _N] * delta_t + pos_ptr[19 _N];
	m_pos_ptr[20 _N] = z[4] = vel_ptr[20 _N] * delta_t + pos_ptr[20 _N];
	m_pos_ptr[21 _N] = z[5] = vel_ptr[21 _N] * delta_t + pos_ptr[21 _N];
	m_pos_ptr[22 _N] = z[6] = vel_ptr[22 _N] * delta_t + pos_ptr[22 _N];
	m_pos_ptr[23 _N] = z[7] = vel_ptr[23 _N] * delta_t + pos_ptr[23 _N];

  	t1 =((x[5] - x[2]) + (x[4] - x[3]));
  	t2 =((x[5] - x[0]) + (x[6] - x[3]));
  	t3 =((x[6] - x[1]) + (x[7] - x[0]));
  	t4 =((x[7] - x[5]) + (x[3] - x[1]));
  	t5 =((x[7] - x[2]) + (x[4] - x[1]));
  	t6 =((x[6] - x[4]) + (x[2] - x[0]));

//	Y grad

    grad_ptr[ 8 _N] = (z[1] *  t1) - (z[2] * (x[3] - x[1])) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * (x[4] - x[1])) - (z[7] * (x[4] - x[3])); 
    grad_ptr[ 9 _N] = (z[2] *  t2) + (z[3] * (x[2] - x[0])) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * (x[5] - x[2])) - (z[4] * (x[5] - x[0])); 
    grad_ptr[10 _N] = (z[3] *  t3) + (z[0] * (x[3] - x[1])) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * (x[6] - x[3])) - (z[5] * (x[6] - x[1]));
    grad_ptr[11 _N] = (z[0] *  t5) - (z[1] * (x[2] - x[0])) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * (x[7] - x[0])) - (z[6] * (x[7] - x[2])); 
    grad_ptr[12 _N] = (z[5] *  t3) + (z[6] * (x[7] - x[5])) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * (x[7] - x[0])) + (z[1] * (x[5] - x[0])); 
    grad_ptr[13 _N] = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * (x[6] - x[4])) + (z[1] *  t6) - (z[0] * (x[4] - x[1])) + (z[2] * (x[6] - x[1])); 
	grad_ptr[14 _N] = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * (x[7] - x[5])) + (z[2] *  t4) - (z[1] * (x[5] - x[2])) + (z[3] * (x[7] - x[2])); 
    grad_ptr[15 _N] = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * (x[6] - x[4])) - (z[3] *  t6) - (z[2] * (x[6] - x[3])) + (z[0] * (x[4] - x[3])); 


  	t1 =((z[5] - z[2]) + (z[4] - z[3]));
  	t2 =((z[5] - z[0]) + (z[6] - z[3]));
  	t3 =((z[6] - z[1]) + (z[7] - z[0]));
  	t4 =((z[7] - z[5]) + (z[3] - z[1]));
  	t5 =((z[7] - z[2]) + (z[4] - z[1]));
  	t6 =((z[6] - z[4]) + (z[2] - z[0]));

//	X grad

    grad_ptr[0 _N] = (y[1] *  t1) - (y[2] * (z[3] - z[1])) - (y[3] *  t5)  + (y[4] *  t4) + (y[5] * (z[4] - z[1])) - (y[7] * (z[4] - z[3]));
    grad_ptr[1 _N] = (y[2] *  t2) + (y[3] * (z[2] - z[0])) - (y[0] *  t1)  - (y[5] *  t6) + (y[6] * (z[5] - z[2])) - (y[4] * (z[5] - z[0])); 
    grad_ptr[2 _N] = (y[3] *  t3) + (y[0] * (z[3] - z[1])) - (y[1] *  t2)  - (y[6] *  t4) + (y[7] * (z[6] - z[3])) - (y[5] * (z[6] - z[1]));
    grad_ptr[3 _N] = (y[0] *  t5) - (y[1] * (z[2] - z[0])) - (y[2] *  t3)  + (y[7] *  t6) + (y[4] * (z[7] - z[0])) - (y[6] * (z[7] - z[2])); 
    grad_ptr[4 _N] = (y[5] *  t3) + (y[6] * (z[7] - z[5])) - (y[7] *  t2)  - (y[0] *  t4) - (y[3] * (z[7] - z[0])) + (y[1] * (z[5] - z[0]));
    grad_ptr[5 _N] = (y[6] *  t5) - (y[4] *  t3)  - (y[7] * (z[6] - z[4])) + (y[1] *  t6) - (y[0] * (z[4] - z[1])) + (y[2] * (z[6] - z[1]));
    grad_ptr[6 _N] = (y[7] *  t1) - (y[5] *  t5)  - (y[4] * (z[7] - z[5])) + (y[2] *  t4) - (y[1] * (z[5] - z[2])) + (y[3] * (z[7] - z[2])); 
    grad_ptr[7 _N] = (y[4] *  t2) - (y[6] *  t1)  + (y[5] * (z[6] - z[4])) - (y[3] *  t6) - (y[2] * (z[6] - z[3])) + (y[0] * (z[4] - z[3]));


	delta_t  = x[0] * grad_ptr[0 _N]; 
	delta_t += x[1] * grad_ptr[1 _N]; 
	delta_t += x[2] * grad_ptr[2 _N]; 
	delta_t += x[3] * grad_ptr[3 _N];
	delta_t += x[4] * grad_ptr[4 _N]; 
	delta_t += x[5] * grad_ptr[5 _N];
	delta_t += x[6] * grad_ptr[6 _N]; 
	delta_t += x[7] * grad_ptr[7 _N];
	delta_t /= 12.0;


}

__global__ void cuda_comp_grad_simple(double * pos, double * grad, double * m_vol, double dt){

	int disp_index = blockIdx.x * blockDim.x * 24 + threadIdx.x;

	double * pos_ptr = &pos[disp_index];
	double * grad_ptr = &grad[disp_index];
	
	double x[8], y[8], z[8];
	double delta_t = -0.5 * dt;

	x[0] = pos_ptr[ 0 _N];
	x[1] = pos_ptr[ 1 _N];
	x[2] = pos_ptr[ 2 _N];
	x[3] = pos_ptr[ 3 _N];
	x[4] = pos_ptr[ 4 _N];
	x[5] = pos_ptr[ 5 _N];
	x[6] = pos_ptr[ 6 _N];
	x[7] = pos_ptr[ 7 _N];

	y[0] = pos_ptr[ 8 _N];
	y[1] = pos_ptr[ 9 _N];
	y[2] = pos_ptr[10 _N];
	y[3] = pos_ptr[11 _N];
	y[4] = pos_ptr[12 _N];
	y[5] = pos_ptr[13 _N];
	y[6] = pos_ptr[14 _N];
	y[7] = pos_ptr[15 _N];

  	double t1 =((y[5] - y[2]) + (y[4] - y[3]));
  	double t2 =((y[5] - y[0]) + (y[6] - y[3]));
  	double t3 =((y[6] - y[1]) + (y[7] - y[0]));
  	double t4 =((y[7] - y[5]) + (y[3] - y[1]));
  	double t5 =((y[7] - y[2]) + (y[4] - y[1]));
  	double t6 =((y[6] - y[4]) + (y[2] - y[0]));

//	Z grad

    grad_ptr[16 _N] = (x[1] *  t1) - (x[2] * (y[3] - y[1])) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * (y[4] - y[1])) - (x[7] * (y[4] - y[3]));
    grad_ptr[17 _N] = (x[2] *  t2) + (x[3] * (y[2] - y[0])) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * (y[5] - y[2])) - (x[4] * (y[5] - y[0]));
    grad_ptr[18 _N] = (x[3] *  t3) + (x[0] * (y[3] - y[1])) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * (y[6] - y[3])) - (x[5] * (y[6] - y[1]));
    grad_ptr[19 _N] = (x[0] *  t5) - (x[1] * (y[2] - y[0])) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * (y[7] - y[0])) - (x[6] * (y[7] - y[2])); 
    grad_ptr[20 _N] = (x[5] *  t3) + (x[6] * (y[7] - y[5])) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * (y[7] - y[0])) + (x[1] * (y[5] - y[0]));
    grad_ptr[21 _N] = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * (y[6] - y[4])) + (x[1] *  t6) - (x[0] * (y[4] - y[1])) + (x[2] * (y[6] - y[1])); 
    grad_ptr[22 _N] = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * (y[7] - y[5])) + (x[2] *  t4) - (x[1] * (y[5] - y[2])) + (x[3] * (y[7] - y[2])); 
    grad_ptr[23 _N] = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * (y[6] - y[4])) - (x[3] *  t6) - (x[2] * (y[6] - y[3])) + (x[0] * (y[4] - y[3]));

	z[0] = pos_ptr[16 _N];
	z[1] = pos_ptr[17 _N];
	z[2] = pos_ptr[18 _N];
	z[3] = pos_ptr[19 _N];
	z[4] = pos_ptr[20 _N];
	z[5] = pos_ptr[21 _N];
	z[6] = pos_ptr[22 _N];
	z[7] = pos_ptr[23 _N];

  	t1 =((x[5] - x[2]) + (x[4] - x[3]));
  	t2 =((x[5] - x[0]) + (x[6] - x[3]));
  	t3 =((x[6] - x[1]) + (x[7] - x[0]));
  	t4 =((x[7] - x[5]) + (x[3] - x[1]));
  	t5 =((x[7] - x[2]) + (x[4] - x[1]));
  	t6 =((x[6] - x[4]) + (x[2] - x[0]));

//	Y grad

    grad_ptr[ 8 _N] = (z[1] *  t1) - (z[2] * (x[3] - x[1])) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * (x[4] - x[1])) - (z[7] * (x[4] - x[3])); 
    grad_ptr[ 9 _N] = (z[2] *  t2) + (z[3] * (x[2] - x[0])) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * (x[5] - x[2])) - (z[4] * (x[5] - x[0])); 
    grad_ptr[10 _N] = (z[3] *  t3) + (z[0] * (x[3] - x[1])) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * (x[6] - x[3])) - (z[5] * (x[6] - x[1]));
    grad_ptr[11 _N] = (z[0] *  t5) - (z[1] * (x[2] - x[0])) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * (x[7] - x[0])) - (z[6] * (x[7] - x[2])); 
    grad_ptr[12 _N] = (z[5] *  t3) + (z[6] * (x[7] - x[5])) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * (x[7] - x[0])) + (z[1] * (x[5] - x[0])); 
    grad_ptr[13 _N] = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * (x[6] - x[4])) + (z[1] *  t6) - (z[0] * (x[4] - x[1])) + (z[2] * (x[6] - x[1])); 
	grad_ptr[14 _N] = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * (x[7] - x[5])) + (z[2] *  t4) - (z[1] * (x[5] - x[2])) + (z[3] * (x[7] - x[2])); 
    grad_ptr[15 _N] = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * (x[6] - x[4])) - (z[3] *  t6) - (z[2] * (x[6] - x[3])) + (z[0] * (x[4] - x[3])); 


  	t1 =((z[5] - z[2]) + (z[4] - z[3]));
  	t2 =((z[5] - z[0]) + (z[6] - z[3]));
  	t3 =((z[6] - z[1]) + (z[7] - z[0]));
  	t4 =((z[7] - z[5]) + (z[3] - z[1]));
  	t5 =((z[7] - z[2]) + (z[4] - z[1]));
  	t6 =((z[6] - z[4]) + (z[2] - z[0]));

//	X grad

    grad_ptr[0 _N] = (y[1] *  t1) - (y[2] * (z[3] - z[1])) - (y[3] *  t5)  + (y[4] *  t4) + (y[5] * (z[4] - z[1])) - (y[7] * (z[4] - z[3]));
    grad_ptr[1 _N] = (y[2] *  t2) + (y[3] * (z[2] - z[0])) - (y[0] *  t1)  - (y[5] *  t6) + (y[6] * (z[5] - z[2])) - (y[4] * (z[5] - z[0])); 
    grad_ptr[2 _N] = (y[3] *  t3) + (y[0] * (z[3] - z[1])) - (y[1] *  t2)  - (y[6] *  t4) + (y[7] * (z[6] - z[3])) - (y[5] * (z[6] - z[1]));
    grad_ptr[3 _N] = (y[0] *  t5) - (y[1] * (z[2] - z[0])) - (y[2] *  t3)  + (y[7] *  t6) + (y[4] * (z[7] - z[0])) - (y[6] * (z[7] - z[2])); 
    grad_ptr[4 _N] = (y[5] *  t3) + (y[6] * (z[7] - z[5])) - (y[7] *  t2)  - (y[0] *  t4) - (y[3] * (z[7] - z[0])) + (y[1] * (z[5] - z[0]));
    grad_ptr[5 _N] = (y[6] *  t5) - (y[4] *  t3)  - (y[7] * (z[6] - z[4])) + (y[1] *  t6) - (y[0] * (z[4] - z[1])) + (y[2] * (z[6] - z[1]));
    grad_ptr[6 _N] = (y[7] *  t1) - (y[5] *  t5)  - (y[4] * (z[7] - z[5])) + (y[2] *  t4) - (y[1] * (z[5] - z[2])) + (y[3] * (z[7] - z[2])); 
    grad_ptr[7 _N] = (y[4] *  t2) - (y[6] *  t1)  + (y[5] * (z[6] - z[4])) - (y[3] *  t6) - (y[2] * (z[6] - z[3])) + (y[0] * (z[4] - z[3]));


	delta_t  = x[0] * grad_ptr[0 _N]; 
	delta_t += x[1] * grad_ptr[1 _N]; 
	delta_t += x[2] * grad_ptr[2 _N]; 
	delta_t += x[3] * grad_ptr[3 _N];
	delta_t += x[4] * grad_ptr[4 _N]; 
	delta_t += x[5] * grad_ptr[5 _N];
	delta_t += x[6] * grad_ptr[6 _N]; 
	delta_t += x[7] * grad_ptr[7 _N];
	delta_t /= 12.0;

	m_vol[blockIdx.x * blockDim.x + threadIdx.x] = delta_t;

}

__global__ void cuda_v_grad(double * vel, double * grad, double * m_vol, double * v_grad){

	int stride =   blockDim.x;
	int index  =   blockIdx.x * blockDim.x * 24 + threadIdx.x;
	int v_index =  blockIdx.x * blockDim.x *  9 + threadIdx.x;

	double volinv12th = 1.0 / (m_vol[blockIdx.x * blockDim.x + threadIdx.x] * 12.0);

	double x[8][2], y[8][2], z[8][2];

	x[0][0] = vel[index + 0 * stride];	 x[0][1] = grad[v_index + 0 * stride];
	x[1][0] = vel[index + 1 * stride];	 x[1][1] = grad[v_index + 1 * stride];
	x[2][0] = vel[index + 2 * stride];	 x[2][1] = grad[v_index + 2 * stride];
	x[3][0] = vel[index + 3 * stride];	 x[3][1] = grad[v_index + 3 * stride];
	x[4][0] = vel[index + 4 * stride];	 x[4][1] = grad[v_index + 4 * stride];
	x[5][0] = vel[index + 5 * stride];	 x[5][1] = grad[v_index + 5 * stride];
	x[6][0] = vel[index + 6 * stride];	 x[6][1] = grad[v_index + 6 * stride];
	x[7][0] = vel[index + 7 * stride];	 x[7][1] = grad[v_index + 7 * stride];

	y[0][0] = vel[index +  8 * stride];	 y[0][1] = grad[v_index +  8 * stride];
	y[1][0] = vel[index +  9 * stride];	 y[1][1] = grad[v_index +  9 * stride];
	y[2][0] = vel[index + 10 * stride];	 y[2][1] = grad[v_index + 10 * stride];
	y[3][0] = vel[index + 11 * stride];	 y[3][1] = grad[v_index + 11 * stride];
	y[4][0] = vel[index + 12 * stride];	 y[4][1] = grad[v_index + 12 * stride];
	y[5][0] = vel[index + 13 * stride];	 y[5][1] = grad[v_index + 13 * stride];
	y[6][0] = vel[index + 14 * stride];	 y[6][1] = grad[v_index + 14 * stride];
	y[7][0] = vel[index + 15 * stride];	 y[7][1] = grad[v_index + 15 * stride];

	z[0][0] = vel[index + 16 * stride];	 z[0][1] = grad[v_index + 16 * stride];
	z[1][0] = vel[index + 17 * stride];	 z[1][1] = grad[v_index + 17 * stride];
	z[2][0] = vel[index + 18 * stride];	 z[2][1] = grad[v_index + 18 * stride];
	z[3][0] = vel[index + 19 * stride];	 z[3][1] = grad[v_index + 19 * stride];
	z[4][0] = vel[index + 20 * stride];	 z[4][1] = grad[v_index + 20 * stride];
	z[5][0] = vel[index + 21 * stride];	 z[5][1] = grad[v_index + 21 * stride];
	z[6][0] = vel[index + 22 * stride];	 z[6][1] = grad[v_index + 22 * stride];
	z[7][0] = vel[index + 23 * stride];	 z[7][1] = grad[v_index + 23 * stride];


	v_grad[v_index + K_F_XX * stride] =	   (x[0][1] * x[0][0] + \
											x[1][1] * x[1][0] + \
											x[2][1] * x[2][0] + \
											x[3][1] * x[3][0] + \
											x[4][1] * x[4][0] + \
											x[5][1] * x[5][0] + \
											x[6][1] * x[6][0] + \
											x[7][1] * x[7][0]) * volinv12th;

	v_grad[v_index + K_F_XY * stride] =	   (y[0][1] * x[0][0] + \
											y[1][1] * x[1][0] + \
											y[2][1] * x[2][0] + \
											y[3][1] * x[3][0] + \
											y[4][1] * x[4][0] + \
											y[5][1] * x[5][0] + \
											y[6][1] * x[6][0] + \
											y[7][1] * x[7][0]) * volinv12th;

	v_grad[v_index + K_F_XZ * stride] =	   (z[0][1] * x[0][0] + \
											z[1][1] * x[1][0] + \
											z[2][1] * x[2][0] + \
											z[3][1] * x[3][0] + \
											z[4][1] * x[4][0] + \
											z[5][1] * x[5][0] + \
											z[6][1] * x[6][0] + \
											z[7][1] * x[7][0]) * volinv12th;



	v_grad[v_index + K_F_YX * stride] =	   (x[0][1] * y[0][0] + \
											x[1][1] * y[1][0] + \
											x[2][1] * y[2][0] + \
											x[3][1] * y[3][0] + \
											x[4][1] * y[4][0] + \
											x[5][1] * y[5][0] + \
											x[6][1] * y[6][0] + \
											x[7][1] * y[7][0]) * volinv12th;

	v_grad[v_index + K_F_YY * stride] =	   (y[0][1] * y[0][0] + \
											y[1][1] * y[1][0] + \
											y[2][1] * y[2][0] + \
											y[3][1] * y[3][0] + \
											y[4][1] * y[4][0] + \
											y[5][1] * y[5][0] + \
											y[6][1] * y[6][0] + \
											y[7][1] * y[7][0]) * volinv12th;

	v_grad[v_index + K_F_YZ * stride] =	   (z[0][1] * y[0][0] + \
											z[1][1] * y[1][0] + \
											z[2][1] * y[2][0] + \
											z[3][1] * y[3][0] + \
											z[4][1] * y[4][0] + \
											z[5][1] * y[5][0] + \
											z[6][1] * y[6][0] + \
											z[7][1] * y[7][0]) * volinv12th;



	v_grad[v_index + K_F_ZX * stride] =	   (x[0][1] * z[0][0] + \
											x[1][1] * z[1][0] + \
											x[2][1] * z[2][0] + \
											x[3][1] * z[3][0] + \
											x[4][1] * z[4][0] + \
											x[5][1] * z[5][0] + \
											x[6][1] * z[6][0] + \
											x[7][1] * z[7][0]) * volinv12th;

	v_grad[v_index + K_F_ZY * stride] =	   (y[0][1] * z[0][0] + \
											y[1][1] * z[1][0] + \
											y[2][1] * z[2][0] + \
											y[3][1] * z[3][0] + \
											y[4][1] * z[4][0] + \
											y[5][1] * z[5][0] + \
											y[6][1] * z[6][0] + \
											y[7][1] * z[7][0]) * volinv12th;

	v_grad[v_index + K_F_ZZ * stride] =	   (z[0][1] * z[0][0] + \
											z[1][1] * z[1][0] + \
											z[2][1] * z[2][0] + \
											z[3][1] * z[3][0] + \
											z[4][1] * z[4][0] + \
											z[5][1] * z[5][0] + \
											z[6][1] * z[6][0] + \
											z[7][1] * z[7][0]) * volinv12th;


}

__global__ void cuda_additive_decomp(double * grad, double * stretch, double * vort){

	int stride =   blockDim.x;
	int index  =   blockIdx.x * blockDim.x * 24 + threadIdx.x;
	int s_index =  blockIdx.x * blockDim.x *  6 + threadIdx.x;
	int v_index =  blockIdx.x * blockDim.x *  3 + threadIdx.x;

	double * g = &grad[index];
	double * s = &stretch[s_index];
	double * v = &vort[v_index];

	s[K_S_XX * stride] = g[K_F_XX * stride];
	s[K_S_YY * stride] = g[K_F_YY * stride];
	s[K_S_ZZ * stride] = g[K_F_ZZ * stride];
  	s[K_S_XY * stride] = 0.5*(g[K_F_XY * stride] + g[K_F_YX * stride]);
  	s[K_S_YZ * stride] = 0.5*(g[K_F_YZ * stride] + g[K_F_ZY * stride]);
  	s[K_S_ZX * stride] = 0.5*(g[K_F_ZX * stride] + g[K_F_XZ * stride]);

  	v[K_V_XY * stride] = 0.5*(g[K_F_XY * stride] - g[K_F_YX * stride]);
  	v[K_V_YZ * stride] = 0.5*(g[K_F_YZ * stride] - g[K_F_ZY * stride]);
  	v[K_V_ZX * stride] = 0.5*(g[K_F_ZX * stride] - g[K_F_XZ * stride]);

}

__global__ void cuda_polar_decomp(double dt, double * stretch1, double * stretch2, double * rot_old, double * vort, double *rot_new){

	int stride =   blockDim.x;
	int r_index =  blockDim.x * blockIdx.x * 9 + threadIdx.x;
	int s_index =  blockDim.x * blockIdx.x * 6 + threadIdx.x;
	int v_index =  blockDim.x * blockIdx.x * 3 + threadIdx.x;

	double * s1 = &stretch1[s_index];
	double * s2 = &stretch2[s_index];
	double * r_o = &rot_old[r_index];
	double * r_n = &rot_new[r_index];
	double * v   =    &vort[v_index];

//	stretching = s1
//	stretch    = s2

	double z1 = 	s1[stride * K_S_XY] * s2[stride * K_S_ZX] - s1[stride * K_S_ZX] * s2[stride * K_S_XY] + s1[stride * K_S_YY] * s2[stride * K_S_YZ] -
    				s1[stride * K_S_YZ] * s2[stride * K_S_YY] + s1[stride * K_S_YZ] * s2[stride * K_S_ZZ] - s1[stride * K_S_ZZ] * s2[stride * K_S_YZ];
  	double z2 = 	s1[stride * K_S_ZX] * s2[stride * K_S_XX] - s1[stride * K_S_XX] * s2[stride * K_S_ZX] + s1[stride * K_S_YZ] * s2[stride * K_S_XY] -
    				s1[stride * K_S_XY] * s2[stride * K_S_YZ] + s1[stride * K_S_ZZ] * s2[stride * K_S_ZX] - s1[stride * K_S_ZX] * s2[stride * K_S_ZZ];
  	double z3 = 	s1[stride * K_S_XX] * s2[stride * K_S_XY] - s1[stride * K_S_XY] * s2[stride * K_S_XX] + s1[stride * K_S_XY] * s2[stride * K_S_YY] -
    				s1[stride * K_S_YY] * s2[stride * K_S_XY] + s1[stride * K_S_ZX] * s2[stride * K_S_YZ] - s1[stride * K_S_YZ] * s2[stride * K_S_ZX];

	//
  	// forward elimination
  	//
  	double a1inv = 1.0 / (s2[stride * K_S_YY] + s2[stride * K_S_ZZ]);

  	double a4BYa1 = -s2[stride * K_S_XY] * a1inv;
  	double a2inv = 1.0 / (s2[stride * K_S_ZZ] + s2[stride * K_S_XX] + s2[stride * K_S_XY] * a4BYa1);

  	double a5 =  -s2[stride * K_S_YZ] + s2[stride * K_S_ZX] * a4BYa1;
  	z2 -= z1 * a4BYa1;
  	double a6BYa1 = -s2[stride * K_S_ZX] * a1inv;
  	double a5BYa2 = a5 * a2inv;
  	z3 -= z1 * a6BYa1 - z2 * a5BYa2;
  	//
  	// backward substitution -
  	//
  	z3 /= (s2[stride * K_S_XX] + s2[stride * K_S_YY] + s2[stride * K_S_ZX] * a6BYa1 + a5 * a5BYa2);
  	z2 = (z2 - a5 * z3) * a2inv;
  	z1 = (z1*a1inv - a6BYa1 * z3 -a4BYa1 * z2);
	
	//
  	// calculate rotation rates - recall that spin_rate is an asymmetric tensor,
  	// so compute spin rate vector as dual of spin rate tensor,
  	// i.e   w_i = e_ijk * spin_rate_jk
  	//
  	z1 += v[stride * K_V_YZ];
  	z2 += v[stride * K_V_ZX];
  	z3 += v[stride * K_V_XY];
  	//
  	// update rotation tensor:
  	//   1) premultiply old rotation tensor to get right-hand side.
  	//
  	const double dt_half = 0.5 * dt;
  	double r_XX = r_o[stride * K_F_XX] + dt_half*( z3 * r_o[stride * K_F_YX] - z2 * r_o[stride * K_F_ZX] );
  	double r_YX = r_o[stride * K_F_YX] + dt_half*( z1 * r_o[stride * K_F_ZX] - z3 * r_o[stride * K_F_XX] );
  	double r_ZX = r_o[stride * K_F_ZX] + dt_half*( z2 * r_o[stride * K_F_XX] - z1 * r_o[stride * K_F_YX] );
  	double r_XY = r_o[stride * K_F_XY] + dt_half*( z3 * r_o[stride * K_F_YY] - z2 * r_o[stride * K_F_ZY] );
  	double r_YY = r_o[stride * K_F_YY] + dt_half*( z1 * r_o[stride * K_F_ZY] - z3 * r_o[stride * K_F_XY] );
  	double r_ZY = r_o[stride * K_F_ZY] + dt_half*( z2 * r_o[stride * K_F_XY] - z1 * r_o[stride * K_F_YY] );
  	double r_XZ = r_o[stride * K_F_XZ] + dt_half*( z3 * r_o[stride * K_F_YZ] - z2 * r_o[stride * K_F_ZZ] );
  	double r_YZ = r_o[stride * K_F_YZ] + dt_half*( z1 * r_o[stride * K_F_ZZ] - z3 * r_o[stride * K_F_XZ] );
  	double r_ZZ = r_o[stride * K_F_ZZ] + dt_half*( z2 * r_o[stride * K_F_XZ] - z1 * r_o[stride * K_F_YZ] );
  	//
  	//    2) solve for new rotation tensor via gauss elimination.
  	// forward elimination -
  	//
  	double a12 = - dt_half * z3;
  	double a13 =   dt_half * z2;
  	double b32 = - dt_half * z1;
  	double a22inv = 1.0 / (1.0 + a12 * a12);
	
  	double a13a12 = a13*a12;
  	double a23 = b32 + a13a12;
  	r_YX += r_XX * a12;
  	r_YY += r_XY * a12;
  	r_YZ += r_XZ * a12;
  	b32 = (b32 - a13a12) * a22inv;
  	r_ZX += r_XX * a13 + r_YX * b32;
  	r_ZY += r_XY * a13 + r_YY * b32;
  	r_ZZ += r_XZ * a13 + r_YZ * b32;
  	//
  	// backward substitution -
  	//
  	double a33inv = 1.0 / (1.0 + a13 * a13 + a23 * b32);

  	r_n[stride * K_F_ZX]  = r_ZX * a33inv;
  	r_n[stride * K_F_ZY]  = r_ZY * a33inv;
  	r_n[stride * K_F_ZZ]  = r_ZZ * a33inv;
  	r_n[stride * K_F_YX]  = ( r_YX - r_n[stride * K_F_ZX] * a23 ) * a22inv;
  	r_n[stride * K_F_YY]  = ( r_YY - r_n[stride * K_F_ZY] * a23 ) * a22inv;
  	r_n[stride * K_F_YZ]  = ( r_YZ - r_n[stride * K_F_ZZ] * a23 ) * a22inv;
  	r_n[stride * K_F_XX]  = r_XX - r_n[stride * K_F_ZX] * a13 - r_n[stride * K_F_YX] * a12;
  	r_n[stride * K_F_XY]  = r_XY - r_n[stride * K_F_ZY] * a13 - r_n[stride * K_F_YY] * a12;
  	r_n[stride * K_F_XZ]  = r_XZ - r_n[stride * K_F_ZZ] * a13 - r_n[stride * K_F_YZ] * a12;
  	//
  	// update stretch tensor in the new configuration -
  	//
  	double a1 = s1[stride * K_S_XY] + v[stride * K_V_XY];
  	double a2 = s1[stride * K_S_YZ] + v[stride * K_V_YZ];
  	double a3 = s1[stride * K_S_ZX] + v[stride * K_V_ZX];
  	double b1 = s1[stride * K_S_ZX] - v[stride * K_V_ZX];
  	double b2 = s1[stride * K_S_XY] - v[stride * K_V_XY];
  	double b3 = s1[stride * K_S_YZ] - v[stride * K_V_YZ];
	
  	double s_XX = s2[stride * K_S_XX];
  	double s_YY = s2[stride * K_S_YY];
  	double s_ZZ = s2[stride * K_S_ZZ];
  	double s_XY = s2[stride * K_S_XY];
  	double s_YZ = s2[stride * K_S_YZ];
  	double s_ZX = s2[stride * K_S_ZX];
	
  	s2[stride * K_S_XX] += dt * (s1[stride * K_S_XX] * s_XX + ( a1 + z3 ) * s_XY + ( b1 - z2 ) * s_ZX);
  	s2[stride * K_S_YY] += dt * (s1[stride * K_S_YY] * s_YY + ( a2 + z1 ) * s_YZ + ( b2 - z3 ) * s_XY);
  	s2[stride * K_S_ZZ] += dt * (s1[stride * K_S_ZZ] * s_ZZ + ( a3 + z2 ) * s_ZX + ( b3 - z1 ) * s_YZ);
  	s2[stride * K_S_XY] += dt * (s1[stride * K_S_XX] * s_XY + ( a1 )      * s_YY + ( b1      ) * s_YZ - z3 * s_XX + z1 * s_ZX);
  	s2[stride * K_S_YZ] += dt * (s1[stride * K_S_YY] * s_YZ + ( a2 )      * s_ZZ + ( b2      ) * s_ZX - z1 * s_YY + z2 * s_XY);
  	s2[stride * K_S_ZX] += dt * (s1[stride * K_S_ZZ] * s_ZX + ( a3 )      * s_XX + ( b3      ) * s_XY - z2 * s_ZZ + z3 * s_YZ);

}

__global__ void cuda_rotate_tensor_forward(double * rot_new, double * stretch, double * rot_str){

	int stride =   blockDim.x;
	int r_index =  blockDim.x * blockIdx.x * 9 + threadIdx.x;
	int s_index =  blockDim.x * blockIdx.x * 6 + threadIdx.x;

	double * s = &stretch[s_index];
	double * r =  &rot_new[r_index];
	double * r_s = &rot_str[s_index];

	double t[9];

  	t[0] = s[stride * K_S_XX]*r[stride * K_F_XX] + s[stride * K_S_XY]*r[stride * K_F_YX] + s[stride * K_S_XZ]*r[stride * K_F_ZX];
  	t[1] = s[stride * K_S_YX]*r[stride * K_F_XX] + s[stride * K_S_YY]*r[stride * K_F_YX] + s[stride * K_S_YZ]*r[stride * K_F_ZX];
  	t[2] = s[stride * K_S_ZX]*r[stride * K_F_XX] + s[stride * K_S_ZY]*r[stride * K_F_YX] + s[stride * K_S_ZZ]*r[stride * K_F_ZX];
  	t[3] = s[stride * K_S_XX]*r[stride * K_F_XY] + s[stride * K_S_XY]*r[stride * K_F_YY] + s[stride * K_S_XZ]*r[stride * K_F_ZY];
  	t[4] = s[stride * K_S_YX]*r[stride * K_F_XY] + s[stride * K_S_YY]*r[stride * K_F_YY] + s[stride * K_S_YZ]*r[stride * K_F_ZY];
  	t[5] = s[stride * K_S_ZX]*r[stride * K_F_XY] + s[stride * K_S_ZY]*r[stride * K_F_YY] + s[stride * K_S_ZZ]*r[stride * K_F_ZY];
  	t[6] = s[stride * K_S_XX]*r[stride * K_F_XZ] + s[stride * K_S_XY]*r[stride * K_F_YZ] + s[stride * K_S_XZ]*r[stride * K_F_ZZ];
  	t[7] = s[stride * K_S_YX]*r[stride * K_F_XZ] + s[stride * K_S_YY]*r[stride * K_F_YZ] + s[stride * K_S_YZ]*r[stride * K_F_ZZ];
  	t[8] = s[stride * K_S_ZX]*r[stride * K_F_XZ] + s[stride * K_S_ZY]*r[stride * K_F_YZ] + s[stride * K_S_ZZ]*r[stride * K_F_ZZ];

  	r_s[stride * K_S_XX] = r[stride * K_F_XX]*t[0] + r[stride * K_F_YX]*t[1] + r[stride * K_F_ZX]*t[2];
  	r_s[stride * K_S_YY] = r[stride * K_F_XY]*t[3] + r[stride * K_F_YY]*t[4] + r[stride * K_F_ZY]*t[5];
  	r_s[stride * K_S_ZZ] = r[stride * K_F_XZ]*t[6] + r[stride * K_F_YZ]*t[7] + r[stride * K_F_ZZ]*t[8];

  	r_s[stride * K_S_XY] = r[stride * K_F_XX]*t[3] + r[stride * K_F_YX]*t[4] + r[stride * K_F_ZX]*t[5];
  	r_s[stride * K_S_YZ] = r[stride * K_F_XY]*t[6] + r[stride * K_F_YY]*t[7] + r[stride * K_F_ZY]*t[8];
  	r_s[stride * K_S_ZX] = r[stride * K_F_XZ]*t[0] + r[stride * K_F_YZ]*t[1] + r[stride * K_F_ZZ]*t[2];

}

__global__ void cuda_rotate_tensor_backward(double * rot_new, double * stretch, double * rot_str) {

	int stride =   blockDim.x;
	int r_index =  blockDim.x * blockIdx.x * 9 + threadIdx.x;
	int s_index =  blockDim.x * blockIdx.x * 6 + threadIdx.x;

	double * s = &stretch[s_index];
	double * r =  &rot_new[r_index];
	double * r_s = &rot_str[s_index];

	double t[9];

  	t[0] = s[stride * K_S_XX]*r[stride * K_F_XX] + s[stride * K_S_XY]*r[stride * K_F_XY] + s[stride * K_S_XZ]*r[stride * K_F_XZ];
  	t[1] = s[stride * K_S_YX]*r[stride * K_F_XX] + s[stride * K_S_YY]*r[stride * K_F_XY] + s[stride * K_S_YZ]*r[stride * K_F_XZ];
  	t[2] = s[stride * K_S_ZX]*r[stride * K_F_XX] + s[stride * K_S_ZY]*r[stride * K_F_XY] + s[stride * K_S_ZZ]*r[stride * K_F_XZ];
  	t[3] = s[stride * K_S_XX]*r[stride * K_F_YX] + s[stride * K_S_XY]*r[stride * K_F_YY] + s[stride * K_S_XZ]*r[stride * K_F_YZ];
  	t[4] = s[stride * K_S_YX]*r[stride * K_F_YX] + s[stride * K_S_YY]*r[stride * K_F_YY] + s[stride * K_S_YZ]*r[stride * K_F_YZ];
  	t[5] = s[stride * K_S_ZX]*r[stride * K_F_YX] + s[stride * K_S_ZY]*r[stride * K_F_YY] + s[stride * K_S_ZZ]*r[stride * K_F_YZ];
  	t[6] = s[stride * K_S_XX]*r[stride * K_F_ZX] + s[stride * K_S_XY]*r[stride * K_F_ZY] + s[stride * K_S_XZ]*r[stride * K_F_ZZ];
  	t[7] = s[stride * K_S_YX]*r[stride * K_F_ZX] + s[stride * K_S_YY]*r[stride * K_F_ZY] + s[stride * K_S_YZ]*r[stride * K_F_ZZ];
  	t[8] = s[stride * K_S_ZX]*r[stride * K_F_ZX] + s[stride * K_S_ZY]*r[stride * K_F_ZY] + s[stride * K_S_ZZ]*r[stride * K_F_ZZ];

  	r_s[stride * K_S_XX] = r[stride * K_F_XX]*t[0] + r[stride * K_F_XY]*t[1] + r[stride * K_F_XZ]*t[2];
  	r_s[stride * K_S_YY] = r[stride * K_F_YZ]*t[3] + r[stride * K_F_YY]*t[4] + r[stride * K_F_YZ]*t[5];
  	r_s[stride * K_S_ZZ] = r[stride * K_F_ZX]*t[6] + r[stride * K_F_ZY]*t[7] + r[stride * K_F_ZZ]*t[8];

  	r_s[stride * K_S_XY] = r[stride * K_F_XX]*t[3] + r[stride * K_F_XY]*t[4] + r[stride * K_F_XZ]*t[5];
  	r_s[stride * K_S_YZ] = r[stride * K_F_YX]*t[6] + r[stride * K_F_YY]*t[7] + r[stride * K_F_YZ]*t[8];
  	r_s[stride * K_S_ZX] = r[stride * K_F_ZX]*t[0] + r[stride * K_F_ZY]*t[1] + r[stride * K_F_ZZ]*t[2];

}

__global__ void cuda_comp_hgop(double * m_pos, double * m_vol, double * grad, double * hgop){

	int stride = blockDim.x;
	int disp_index = blockIdx.x * blockDim.x * 24 + threadIdx.x;
	int hgop_index = blockIdx.x * blockDim.x * 32 + threadIdx.x;

	double inv_vol = 1.0 / (12.0 * m_vol[blockIdx.x * blockDim.x + threadIdx.x]);      

  	double hgconst12th[12];
	double x[8], y[8], z[8];

	x[0] = m_pos[disp_index +  0 * stride];
	x[1] = m_pos[disp_index +  1 * stride];
	x[2] = m_pos[disp_index +  2 * stride];
	x[3] = m_pos[disp_index +  3 * stride];
	x[4] = m_pos[disp_index +  4 * stride];
	x[5] = m_pos[disp_index +  5 * stride];
	x[6] = m_pos[disp_index +  6 * stride];
	x[7] = m_pos[disp_index +  7 * stride];

  	double q0 = x[0] - x[1];
  	double q1 = x[2] - x[3];
  	double q2 = x[4] - x[5];
  	double q3 = x[6] - x[7];

  	hgconst12th[0] = ( (x[0]+x[1]) - (x[2]+x[3]) - (x[4]+x[5]) + (x[6]+x[7]) ) * inv_vol;
  	hgconst12th[1] = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[2] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[3] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

	y[0] = m_pos[disp_index +  8 * stride];
	y[1] = m_pos[disp_index +  9 * stride];
	y[2] = m_pos[disp_index + 10 * stride];
	y[3] = m_pos[disp_index + 11 * stride];
	y[4] = m_pos[disp_index + 12 * stride];
	y[5] = m_pos[disp_index + 13 * stride];
	y[6] = m_pos[disp_index + 14 * stride];
	y[7] = m_pos[disp_index + 15 * stride];

  	q0 = (y[0] - y[1]);
  	q1 = (y[2] - y[3]);
  	q2 = (y[4] - y[5]);
  	q3 = (y[6] - y[7]);

  	hgconst12th[4] = ( (y[0]+y[1]) - (y[2]+y[3]) - (y[4]+y[5]) + (y[6]+y[7]) ) * inv_vol;
  	hgconst12th[5] = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[6] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[7] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

	z[0] = m_pos[disp_index + 16 * stride];
	z[1] = m_pos[disp_index + 17 * stride];
	z[2] = m_pos[disp_index + 18 * stride];
	z[3] = m_pos[disp_index + 19 * stride];
	z[4] = m_pos[disp_index + 20 * stride];
	z[5] = m_pos[disp_index + 21 * stride];
	z[6] = m_pos[disp_index + 22 * stride];
	z[7] = m_pos[disp_index + 23 * stride];

  	q0 = (z[0] - z[1]);
  	q1 = (z[2] - z[3]);
  	q2 = (z[4] - z[5]);
  	q3 = (z[6] - z[7]);

  	hgconst12th[8]  = ( (z[0]+z[1]) - (z[2]+z[3]) - (z[4]+z[5]) + (z[6]+z[7]) ) * inv_vol;
  	hgconst12th[9]  = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[10] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[11] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

	/*

	Alternatively:

	const double hgop_arr[32] = {
	 1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0,  1.0,
	 1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0, -1.0,
	 1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0, -1.0,
	-1.0, 1.0,-1.0, 1.0, 1.0,-1.0, 1.0, -1.0};

	for(int i = 0; i < 8; i++){
		hgop(ielem, i   , 1) = hgop_arr[i   ] - (hgconst12th[0] * x[i][1]  + hgconst12th[4] * y[i][1] + hgconst12th[8 ] * z[i][1]);
		hgop(ielem, i+8 , 1) = hgop_arr[i+8 ] - (hgconst12th[1] * x[i][1]  + hgconst12th[5] * y[i][1] + hgconst12th[9 ] * z[i][1]);
		hgop(ielem, i+16, 1) = hgop_arr[i+16] - (hgconst12th[2] * x[i][1]  + hgconst12th[6] * y[i][1] + hgconst12th[10] * z[i][1]);
		hgop(ielem, i+24, 1) = hgop_arr[i+24] - (hgconst12th[3] * x[i][1]  + hgconst12th[7] * y[i][1] + hgconst12th[11] * z[i][1]);
	}

	*/

	double grad_x[8], grad_y[8], grad_z[8];

	#pragma unroll 8
	for(int i = 0; i < 8; i++){
		grad_x[i] = grad[disp_index +  (i +  0) * stride];
		grad_y[i] = grad[disp_index +  (i +  8) * stride];
		grad_z[i] = grad[disp_index +  (i + 16) * stride];
	}

	hgop[hgop_index +  0 * stride] =  1.0 - (hgconst12th[0] * grad_x[0] + hgconst12th[4] * grad_y[0] + hgconst12th[ 8] * grad_z[0]);
	hgop[hgop_index +  1 * stride] =  1.0 - (hgconst12th[0] * grad_x[1] + hgconst12th[4] * grad_y[1] + hgconst12th[ 8] * grad_z[1]);
	hgop[hgop_index +  2 * stride] = -1.0 - (hgconst12th[0] * grad_x[2] + hgconst12th[4] * grad_y[2] + hgconst12th[ 8] * grad_z[2]);
	hgop[hgop_index +  3 * stride] = -1.0 - (hgconst12th[0] * grad_x[3] + hgconst12th[4] * grad_y[3] + hgconst12th[ 8] * grad_z[3]);
	hgop[hgop_index +  4 * stride] = -1.0 - (hgconst12th[0] * grad_x[4] + hgconst12th[4] * grad_y[4] + hgconst12th[ 8] * grad_z[4]);
	hgop[hgop_index +  5 * stride] = -1.0 - (hgconst12th[0] * grad_x[5] + hgconst12th[4] * grad_y[5] + hgconst12th[ 8] * grad_z[5]);
	hgop[hgop_index +  6 * stride] =  1.0 - (hgconst12th[0] * grad_x[6] + hgconst12th[4] * grad_y[6] + hgconst12th[ 8] * grad_z[6]);
	hgop[hgop_index +  7 * stride] =  1.0 - (hgconst12th[0] * grad_x[7] + hgconst12th[4] * grad_y[7] + hgconst12th[ 8] * grad_z[7]);
	hgop[hgop_index +  8 * stride] =  1.0 - (hgconst12th[1] * grad_x[0] + hgconst12th[5] * grad_y[0] + hgconst12th[ 9] * grad_z[0]);
	hgop[hgop_index +  9 * stride] = -1.0 - (hgconst12th[1] * grad_x[1] + hgconst12th[5] * grad_y[1] + hgconst12th[ 9] * grad_z[1]);
	hgop[hgop_index + 10 * stride] = -1.0 - (hgconst12th[1] * grad_x[2] + hgconst12th[5] * grad_y[2] + hgconst12th[ 9] * grad_z[2]);
	hgop[hgop_index + 11 * stride] =  1.0 - (hgconst12th[1] * grad_x[3] + hgconst12th[5] * grad_y[3] + hgconst12th[ 9] * grad_z[3]);
	hgop[hgop_index + 12 * stride] = -1.0 - (hgconst12th[1] * grad_x[4] + hgconst12th[5] * grad_y[4] + hgconst12th[ 9] * grad_z[4]);
	hgop[hgop_index + 13 * stride] =  1.0 - (hgconst12th[1] * grad_x[5] + hgconst12th[5] * grad_y[5] + hgconst12th[ 9] * grad_z[5]);
	hgop[hgop_index + 14 * stride] =  1.0 - (hgconst12th[1] * grad_x[6] + hgconst12th[5] * grad_y[6] + hgconst12th[ 9] * grad_z[6]);
	hgop[hgop_index + 15 * stride] = -1.0 - (hgconst12th[1] * grad_x[7] + hgconst12th[5] * grad_y[7] + hgconst12th[ 9] * grad_z[7]);
	hgop[hgop_index + 16 * stride] =  1.0 - (hgconst12th[2] * grad_x[0] + hgconst12th[6] * grad_y[0] + hgconst12th[10] * grad_z[0]);
	hgop[hgop_index + 17 * stride] = -1.0 - (hgconst12th[2] * grad_x[1] + hgconst12th[6] * grad_y[1] + hgconst12th[10] * grad_z[1]);
	hgop[hgop_index + 18 * stride] =  1.0 - (hgconst12th[2] * grad_x[2] + hgconst12th[6] * grad_y[2] + hgconst12th[10] * grad_z[2]);
	hgop[hgop_index + 19 * stride] = -1.0 - (hgconst12th[2] * grad_x[3] + hgconst12th[6] * grad_y[3] + hgconst12th[10] * grad_z[3]);
	hgop[hgop_index + 20 * stride] =  1.0 - (hgconst12th[2] * grad_x[4] + hgconst12th[6] * grad_y[4] + hgconst12th[10] * grad_z[4]);
	hgop[hgop_index + 21 * stride] = -1.0 - (hgconst12th[2] * grad_x[5] + hgconst12th[6] * grad_y[5] + hgconst12th[10] * grad_z[5]);
	hgop[hgop_index + 22 * stride] =  1.0 - (hgconst12th[2] * grad_x[6] + hgconst12th[6] * grad_y[6] + hgconst12th[10] * grad_z[6]);
	hgop[hgop_index + 23 * stride] = -1.0 - (hgconst12th[2] * grad_x[7] + hgconst12th[6] * grad_y[7] + hgconst12th[10] * grad_z[7]);
	hgop[hgop_index + 24 * stride] = -1.0 - (hgconst12th[3] * grad_x[0] + hgconst12th[7] * grad_y[0] + hgconst12th[11] * grad_z[0]);
	hgop[hgop_index + 25 * stride] =  1.0 - (hgconst12th[3] * grad_x[1] + hgconst12th[7] * grad_y[1] + hgconst12th[11] * grad_z[1]);
	hgop[hgop_index + 26 * stride] = -1.0 - (hgconst12th[3] * grad_x[2] + hgconst12th[7] * grad_y[2] + hgconst12th[11] * grad_z[2]);
	hgop[hgop_index + 27 * stride] =  1.0 - (hgconst12th[3] * grad_x[3] + hgconst12th[7] * grad_y[3] + hgconst12th[11] * grad_z[3]);
	hgop[hgop_index + 28 * stride] =  1.0 - (hgconst12th[3] * grad_x[4] + hgconst12th[7] * grad_y[4] + hgconst12th[11] * grad_z[4]);
	hgop[hgop_index + 29 * stride] = -1.0 - (hgconst12th[3] * grad_x[5] + hgconst12th[7] * grad_y[5] + hgconst12th[11] * grad_z[5]);
	hgop[hgop_index + 30 * stride] =  1.0 - (hgconst12th[3] * grad_x[6] + hgconst12th[7] * grad_y[6] + hgconst12th[11] * grad_z[6]);
	hgop[hgop_index + 31 * stride] = -1.0 - (hgconst12th[3] * grad_x[7] + hgconst12th[7] * grad_y[7] + hgconst12th[11] * grad_z[7]);

}

__global__ void cuda_comp_hgop2(double * m_pos, double * m_vol, double * grad, double * hgop){

	int stride = blockDim.x;
	int disp_index = blockIdx.x * blockDim.x * 24 + threadIdx.x;
	int hgop_index = blockIdx.x * blockDim.x * 32 + threadIdx.x;

	double inv_vol = 1.0 / (12.0 * m_vol[blockIdx.x * blockDim.x + threadIdx.x]);      

  	double hgconst12th[12];

  	double q0 = m_pos[disp_index + stride *0] - m_pos[disp_index + stride *1];
  	double q1 = m_pos[disp_index + stride *2] - m_pos[disp_index + stride *3];
  	double q2 = m_pos[disp_index + stride *4] - m_pos[disp_index + stride *5];
  	double q3 = m_pos[disp_index + stride *6] - m_pos[disp_index + stride *7];

  	hgconst12th[0] = ( 	(m_pos[disp_index + stride *0]+m_pos[disp_index + stride *1]) - 
						(m_pos[disp_index + stride *2]+m_pos[disp_index + stride *3]) - 
						(m_pos[disp_index + stride *4]+m_pos[disp_index + stride *5]) + 
						(m_pos[disp_index + stride *6]+m_pos[disp_index + stride *7]) ) * inv_vol;

  	hgconst12th[1] = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[2] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[3] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

  	q0 = (m_pos[disp_index + stride *0] - m_pos[disp_index + stride *1]);
  	q1 = (m_pos[disp_index + stride *2] - m_pos[disp_index + stride *3]);
  	q2 = (m_pos[disp_index + stride *4] - m_pos[disp_index + stride *5]);
  	q3 = (m_pos[disp_index + stride *6] - m_pos[disp_index + stride *7]);

  	hgconst12th[4] = ( 	(m_pos[disp_index + stride *8]+m_pos[disp_index + stride *9]) - 
						(m_pos[disp_index + stride *10]+m_pos[disp_index + stride *11]) - 
						(m_pos[disp_index + stride *12]+m_pos[disp_index + stride *13]) + 
						(m_pos[disp_index + stride *14]+m_pos[disp_index + stride *15]) ) * inv_vol;

  	hgconst12th[5] = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[6] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[7] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

  	q0 = (m_pos[disp_index + stride *0] - m_pos[disp_index + stride *1]);
  	q1 = (m_pos[disp_index + stride *2] - m_pos[disp_index + stride *3]);
  	q2 = (m_pos[disp_index + stride *4] - m_pos[disp_index + stride *5]);
  	q3 = (m_pos[disp_index + stride *6] - m_pos[disp_index + stride *7]);

  	hgconst12th[8]  = ( (m_pos[disp_index + stride *16]+m_pos[disp_index + stride *17]) - 
						(m_pos[disp_index + stride *18]+m_pos[disp_index + stride *19]) - 
						(m_pos[disp_index + stride *20]+m_pos[disp_index + stride *21]) + 
						(m_pos[disp_index + stride *22]+m_pos[disp_index + stride *23]) ) * inv_vol;

  	hgconst12th[9]  = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[10] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[11] = ( -q0 - q1 + q2 + q3 ) * inv_vol;


	hgop[hgop_index +  0 * stride]= 1.0 - (hgconst12th[0] * grad[disp_index + 0 _M] + hgconst12th[4] * grad[disp_index + 0 _M] + hgconst12th[ 8] * grad[disp_index + 0 _M]);
	hgop[hgop_index +  1 * stride]= 1.0 - (hgconst12th[0] * grad[disp_index + 1 _M] + hgconst12th[4] * grad[disp_index + 1 _M] + hgconst12th[ 8] * grad[disp_index + 1 _M]);
	hgop[hgop_index +  2 * stride]=-1.0 - (hgconst12th[0] * grad[disp_index + 2 _M] + hgconst12th[4] * grad[disp_index + 2 _M] + hgconst12th[ 8] * grad[disp_index + 2 _M]);
	hgop[hgop_index +  3 * stride]=-1.0 - (hgconst12th[0] * grad[disp_index + 3 _M] + hgconst12th[4] * grad[disp_index + 3 _M] + hgconst12th[ 8] * grad[disp_index + 3 _M]);
	hgop[hgop_index +  4 * stride]=-1.0 - (hgconst12th[0] * grad[disp_index + 4 _M] + hgconst12th[4] * grad[disp_index + 4 _M] + hgconst12th[ 8] * grad[disp_index + 4 _M]);
	hgop[hgop_index +  5 * stride]=-1.0 - (hgconst12th[0] * grad[disp_index + 5 _M] + hgconst12th[4] * grad[disp_index + 5 _M] + hgconst12th[ 8] * grad[disp_index + 5 _M]);
	hgop[hgop_index +  6 * stride]= 1.0 - (hgconst12th[0] * grad[disp_index + 6 _M] + hgconst12th[4] * grad[disp_index + 6 _M] + hgconst12th[ 8] * grad[disp_index + 6 _M]);
	hgop[hgop_index +  7 * stride]= 1.0 - (hgconst12th[0] * grad[disp_index + 7 _M] + hgconst12th[4] * grad[disp_index + 7 _M] + hgconst12th[ 8] * grad[disp_index + 7 _M]);
	hgop[hgop_index +  8 * stride]= 1.0 - (hgconst12th[1] * grad[disp_index + 0 _M] + hgconst12th[5] * grad[disp_index + 0 _M] + hgconst12th[ 9] * grad[disp_index + 0 _M]);
	hgop[hgop_index +  9 * stride]=-1.0 - (hgconst12th[1] * grad[disp_index + 1 _M] + hgconst12th[5] * grad[disp_index + 1 _M] + hgconst12th[ 9] * grad[disp_index + 1 _M]);
	hgop[hgop_index + 10 * stride]=-1.0 - (hgconst12th[1] * grad[disp_index + 2 _M] + hgconst12th[5] * grad[disp_index + 2 _M] + hgconst12th[ 9] * grad[disp_index + 2 _M]);
	hgop[hgop_index + 11 * stride]= 1.0 - (hgconst12th[1] * grad[disp_index + 3 _M] + hgconst12th[5] * grad[disp_index + 3 _M] + hgconst12th[ 9] * grad[disp_index + 3 _M]);
	hgop[hgop_index + 12 * stride]=-1.0 - (hgconst12th[1] * grad[disp_index + 4 _M] + hgconst12th[5] * grad[disp_index + 4 _M] + hgconst12th[ 9] * grad[disp_index + 4 _M]);
	hgop[hgop_index + 13 * stride]= 1.0 - (hgconst12th[1] * grad[disp_index + 5 _M] + hgconst12th[5] * grad[disp_index + 5 _M] + hgconst12th[ 9] * grad[disp_index + 5 _M]);
	hgop[hgop_index + 14 * stride]= 1.0 - (hgconst12th[1] * grad[disp_index + 6 _M] + hgconst12th[5] * grad[disp_index + 6 _M] + hgconst12th[ 9] * grad[disp_index + 6 _M]);
	hgop[hgop_index + 15 * stride]=-1.0 - (hgconst12th[1] * grad[disp_index + 7 _M] + hgconst12th[5] * grad[disp_index + 7 _M] + hgconst12th[ 9] * grad[disp_index + 7 _M]);
	hgop[hgop_index + 16 * stride]= 1.0 - (hgconst12th[2] * grad[disp_index + 0 _M] + hgconst12th[6] * grad[disp_index + 0 _M] + hgconst12th[10] * grad[disp_index + 0 _M]);
	hgop[hgop_index + 17 * stride]=-1.0 - (hgconst12th[2] * grad[disp_index + 1 _M] + hgconst12th[6] * grad[disp_index + 1 _M] + hgconst12th[10] * grad[disp_index + 1 _M]);
	hgop[hgop_index + 18 * stride]= 1.0 - (hgconst12th[2] * grad[disp_index + 2 _M] + hgconst12th[6] * grad[disp_index + 2 _M] + hgconst12th[10] * grad[disp_index + 2 _M]);
	hgop[hgop_index + 19 * stride]=-1.0 - (hgconst12th[2] * grad[disp_index + 3 _M] + hgconst12th[6] * grad[disp_index + 3 _M] + hgconst12th[10] * grad[disp_index + 3 _M]);
	hgop[hgop_index + 20 * stride]= 1.0 - (hgconst12th[2] * grad[disp_index + 4 _M] + hgconst12th[6] * grad[disp_index + 4 _M] + hgconst12th[10] * grad[disp_index + 4 _M]);
	hgop[hgop_index + 21 * stride]=-1.0 - (hgconst12th[2] * grad[disp_index + 5 _M] + hgconst12th[6] * grad[disp_index + 5 _M] + hgconst12th[10] * grad[disp_index + 5 _M]);
	hgop[hgop_index + 22 * stride]= 1.0 - (hgconst12th[2] * grad[disp_index + 6 _M] + hgconst12th[6] * grad[disp_index + 6 _M] + hgconst12th[10] * grad[disp_index + 6 _M]);
	hgop[hgop_index + 23 * stride]=-1.0 - (hgconst12th[2] * grad[disp_index + 7 _M] + hgconst12th[6] * grad[disp_index + 7 _M] + hgconst12th[10] * grad[disp_index + 7 _M]);
	hgop[hgop_index + 24 * stride]=-1.0 - (hgconst12th[3] * grad[disp_index + 0 _M] + hgconst12th[7] * grad[disp_index + 0 _M] + hgconst12th[11] * grad[disp_index + 0 _M]);
	hgop[hgop_index + 25 * stride]= 1.0 - (hgconst12th[3] * grad[disp_index + 1 _M] + hgconst12th[7] * grad[disp_index + 1 _M] + hgconst12th[11] * grad[disp_index + 1 _M]);
	hgop[hgop_index + 26 * stride]=-1.0 - (hgconst12th[3] * grad[disp_index + 2 _M] + hgconst12th[7] * grad[disp_index + 2 _M] + hgconst12th[11] * grad[disp_index + 2 _M]);
	hgop[hgop_index + 27 * stride]= 1.0 - (hgconst12th[3] * grad[disp_index + 3 _M] + hgconst12th[7] * grad[disp_index + 3 _M] + hgconst12th[11] * grad[disp_index + 3 _M]);
	hgop[hgop_index + 28 * stride]= 1.0 - (hgconst12th[3] * grad[disp_index + 4 _M] + hgconst12th[7] * grad[disp_index + 4 _M] + hgconst12th[11] * grad[disp_index + 4 _M]);
	hgop[hgop_index + 29 * stride]=-1.0 - (hgconst12th[3] * grad[disp_index + 5 _M] + hgconst12th[7] * grad[disp_index + 5 _M] + hgconst12th[11] * grad[disp_index + 5 _M]);
	hgop[hgop_index + 30 * stride]= 1.0 - (hgconst12th[3] * grad[disp_index + 6 _M] + hgconst12th[7] * grad[disp_index + 6 _M] + hgconst12th[11] * grad[disp_index + 6 _M]);
	hgop[hgop_index + 31 * stride]=-1.0 - (hgconst12th[3] * grad[disp_index + 7 _M] + hgconst12th[7] * grad[disp_index + 7 _M] + hgconst12th[11] * grad[disp_index + 7 _M]);

}

__global__ void cuda_comp_hgop3(double * m_pos, double * m_vol, double * grad, double * hgop){

	double * m_pos_ptr = &m_pos[blockIdx.x * blockDim.x * 24];
	double * grad_ptr = &grad[blockIdx.x * blockDim.x * 24];
	double * hgop_ptr = &hgop[blockIdx.x * blockDim.x * 24];

	double inv_vol = 1.0 / (12.0 * m_vol[blockIdx.x * blockDim.x + threadIdx.x]);      

  	double hgconst12th[12];

  	double q0 = m_pos_ptr[0 _N] - m_pos_ptr[1 _N];
  	double q1 = m_pos_ptr[2 _N] - m_pos_ptr[3 _N];
  	double q2 = m_pos_ptr[4 _N] - m_pos_ptr[5 _N];
  	double q3 = m_pos_ptr[6 _N] - m_pos_ptr[7 _N];

  	hgconst12th[0] = ( 	(m_pos_ptr[0 _N]+m_pos_ptr[1 _N]) - 
						(m_pos_ptr[2 _N]+m_pos_ptr[3 _N]) - 
						(m_pos_ptr[4 _N]+m_pos_ptr[5 _N]) + 
						(m_pos_ptr[6 _N]+m_pos_ptr[7 _N]) ) * inv_vol;

  	hgconst12th[1] = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[2] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[3] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

  	q0 = (m_pos_ptr[0 _N] - m_pos_ptr[1 _N]);
  	q1 = (m_pos_ptr[2 _N] - m_pos_ptr[3 _N]);
  	q2 = (m_pos_ptr[4 _N] - m_pos_ptr[5 _N]);
  	q3 = (m_pos_ptr[6 _N] - m_pos_ptr[7 _N]);

  	hgconst12th[4] = ( 	(m_pos_ptr[8 _N]+m_pos_ptr[9 _N]) - 
						(m_pos_ptr[10 _N]+m_pos_ptr[11 _N]) - 
						(m_pos_ptr[12 _N]+m_pos_ptr[13 _N]) + 
						(m_pos_ptr[14 _N]+m_pos_ptr[15 _N]) ) * inv_vol;

  	hgconst12th[5] = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[6] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[7] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

  	q0 = (m_pos_ptr[0 _N] - m_pos_ptr[1 _N]);
  	q1 = (m_pos_ptr[2 _N] - m_pos_ptr[3 _N]);
  	q2 = (m_pos_ptr[4 _N] - m_pos_ptr[5 _N]);
  	q3 = (m_pos_ptr[6 _N] - m_pos_ptr[7 _N]);

  	hgconst12th[8]  = ( (m_pos_ptr[16 _N]+m_pos_ptr[17 _N]) - 
						(m_pos_ptr[18 _N]+m_pos_ptr[19 _N]) - 
						(m_pos_ptr[20 _N]+m_pos_ptr[21 _N]) + 
						(m_pos_ptr[22 _N]+m_pos_ptr[23 _N]) ) * inv_vol;

  	hgconst12th[9]  = (  q0 - q1 - q2 + q3 ) * inv_vol;
  	hgconst12th[10] = (  q0 + q1 + q2 + q3 ) * inv_vol;
  	hgconst12th[11] = ( -q0 - q1 + q2 + q3 ) * inv_vol;


	hgop_ptr[ 0 _N]= 1.0 - (hgconst12th[0] * grad_ptr[0 _N] + hgconst12th[4] * grad_ptr[0 _N] + hgconst12th[ 8] * grad_ptr[0 _N]);
	hgop_ptr[ 1 _N]= 1.0 - (hgconst12th[0] * grad_ptr[1 _N] + hgconst12th[4] * grad_ptr[1 _N] + hgconst12th[ 8] * grad_ptr[1 _N]);
	hgop_ptr[ 2 _N]=-1.0 - (hgconst12th[0] * grad_ptr[2 _N] + hgconst12th[4] * grad_ptr[2 _N] + hgconst12th[ 8] * grad_ptr[2 _N]);
	hgop_ptr[ 3 _N]=-1.0 - (hgconst12th[0] * grad_ptr[3 _N] + hgconst12th[4] * grad_ptr[3 _N] + hgconst12th[ 8] * grad_ptr[3 _N]);
	hgop_ptr[ 4 _N]=-1.0 - (hgconst12th[0] * grad_ptr[4 _N] + hgconst12th[4] * grad_ptr[4 _N] + hgconst12th[ 8] * grad_ptr[4 _N]);
	hgop_ptr[ 5 _N]=-1.0 - (hgconst12th[0] * grad_ptr[5 _N] + hgconst12th[4] * grad_ptr[5 _N] + hgconst12th[ 8] * grad_ptr[5 _N]);
	hgop_ptr[ 6 _N]= 1.0 - (hgconst12th[0] * grad_ptr[6 _N] + hgconst12th[4] * grad_ptr[6 _N] + hgconst12th[ 8] * grad_ptr[6 _N]);
	hgop_ptr[ 7 _N]= 1.0 - (hgconst12th[0] * grad_ptr[7 _N] + hgconst12th[4] * grad_ptr[7 _N] + hgconst12th[ 8] * grad_ptr[7 _N]);
	hgop_ptr[ 8 _N]= 1.0 - (hgconst12th[1] * grad_ptr[0 _N] + hgconst12th[5] * grad_ptr[0 _N] + hgconst12th[ 9] * grad_ptr[0 _N]);
	hgop_ptr[ 9 _N]=-1.0 - (hgconst12th[1] * grad_ptr[1 _N] + hgconst12th[5] * grad_ptr[1 _N] + hgconst12th[ 9] * grad_ptr[1 _N]);
	hgop_ptr[10 _N]=-1.0 - (hgconst12th[1] * grad_ptr[2 _N] + hgconst12th[5] * grad_ptr[2 _N] + hgconst12th[ 9] * grad_ptr[2 _N]);
	hgop_ptr[11 _N]= 1.0 - (hgconst12th[1] * grad_ptr[3 _N] + hgconst12th[5] * grad_ptr[3 _N] + hgconst12th[ 9] * grad_ptr[3 _N]);
	hgop_ptr[12 _N]=-1.0 - (hgconst12th[1] * grad_ptr[4 _N] + hgconst12th[5] * grad_ptr[4 _N] + hgconst12th[ 9] * grad_ptr[4 _N]);
	hgop_ptr[13 _N]= 1.0 - (hgconst12th[1] * grad_ptr[5 _N] + hgconst12th[5] * grad_ptr[5 _N] + hgconst12th[ 9] * grad_ptr[5 _N]);
	hgop_ptr[14 _N]= 1.0 - (hgconst12th[1] * grad_ptr[6 _N] + hgconst12th[5] * grad_ptr[6 _N] + hgconst12th[ 9] * grad_ptr[6 _N]);
	hgop_ptr[15 _N]=-1.0 - (hgconst12th[1] * grad_ptr[7 _N] + hgconst12th[5] * grad_ptr[7 _N] + hgconst12th[ 9] * grad_ptr[7 _N]);
	hgop_ptr[16 _N]= 1.0 - (hgconst12th[2] * grad_ptr[0 _N] + hgconst12th[6] * grad_ptr[0 _N] + hgconst12th[10] * grad_ptr[0 _N]);
	hgop_ptr[17 _N]=-1.0 - (hgconst12th[2] * grad_ptr[1 _N] + hgconst12th[6] * grad_ptr[1 _N] + hgconst12th[10] * grad_ptr[1 _N]);
	hgop_ptr[18 _N]= 1.0 - (hgconst12th[2] * grad_ptr[2 _N] + hgconst12th[6] * grad_ptr[2 _N] + hgconst12th[10] * grad_ptr[2 _N]);
	hgop_ptr[19 _N]=-1.0 - (hgconst12th[2] * grad_ptr[3 _N] + hgconst12th[6] * grad_ptr[3 _N] + hgconst12th[10] * grad_ptr[3 _N]);
	hgop_ptr[20 _N]= 1.0 - (hgconst12th[2] * grad_ptr[4 _N] + hgconst12th[6] * grad_ptr[4 _N] + hgconst12th[10] * grad_ptr[4 _N]);
	hgop_ptr[21 _N]=-1.0 - (hgconst12th[2] * grad_ptr[5 _N] + hgconst12th[6] * grad_ptr[5 _N] + hgconst12th[10] * grad_ptr[5 _N]);
	hgop_ptr[22 _N]= 1.0 - (hgconst12th[2] * grad_ptr[6 _N] + hgconst12th[6] * grad_ptr[6 _N] + hgconst12th[10] * grad_ptr[6 _N]);
	hgop_ptr[23 _N]=-1.0 - (hgconst12th[2] * grad_ptr[7 _N] + hgconst12th[6] * grad_ptr[7 _N] + hgconst12th[10] * grad_ptr[7 _N]);
	hgop_ptr[24 _N]=-1.0 - (hgconst12th[3] * grad_ptr[0 _N] + hgconst12th[7] * grad_ptr[0 _N] + hgconst12th[11] * grad_ptr[0 _N]);
	hgop_ptr[25 _N]= 1.0 - (hgconst12th[3] * grad_ptr[1 _N] + hgconst12th[7] * grad_ptr[1 _N] + hgconst12th[11] * grad_ptr[1 _N]);
	hgop_ptr[26 _N]=-1.0 - (hgconst12th[3] * grad_ptr[2 _N] + hgconst12th[7] * grad_ptr[2 _N] + hgconst12th[11] * grad_ptr[2 _N]);
	hgop_ptr[27 _N]= 1.0 - (hgconst12th[3] * grad_ptr[3 _N] + hgconst12th[7] * grad_ptr[3 _N] + hgconst12th[11] * grad_ptr[3 _N]);
	hgop_ptr[28 _N]= 1.0 - (hgconst12th[3] * grad_ptr[4 _N] + hgconst12th[7] * grad_ptr[4 _N] + hgconst12th[11] * grad_ptr[4 _N]);
	hgop_ptr[29 _N]=-1.0 - (hgconst12th[3] * grad_ptr[5 _N] + hgconst12th[7] * grad_ptr[5 _N] + hgconst12th[11] * grad_ptr[5 _N]);
	hgop_ptr[30 _N]= 1.0 - (hgconst12th[3] * grad_ptr[6 _N] + hgconst12th[7] * grad_ptr[6 _N] + hgconst12th[11] * grad_ptr[6 _N]);
	hgop_ptr[31 _N]=-1.0 - (hgconst12th[3] * grad_ptr[7 _N] + hgconst12th[7] * grad_ptr[7 _N] + hgconst12th[11] * grad_ptr[7 _N]);

}

__global__ void cuda_comp_moduli(double * dilmod, double * shrmod, double * bulkmod, double * twomu){

	int index = blockIdx.x * blockDim.x + threadIdx.x;

	shrmod[index] = twomu[index];
	dilmod[index] = ONE3RD * (bulkmod[index] + 2.0 * shrmod[index]);

}

__global__ void cuda_comp_facs(	double * gradop12, 
									double * m_vol, 
									double * e_mass, 
									double * e_t_step,
									double * rot_stret, 
									double * shrmod, 
									double * dilmod, 
									double * fac1,
									double * fac2,
									double * bulkq,
									double fac1_pre,
									double linBulkVisc, 
									double quadBulkVisc,
									double hg_visc){

	int stride = blockDim.x;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int r_index = blockIdx.x * blockDim.x * 6 + threadIdx.x;
	//
	//  Compute the aspect ratio.  Aspect ratio is 0.5 * volume / (grad . grad)
	//  With the 12 factors this is actually 6.0 * 12 * volume / (12 * grad . 12 * grad)
	//
	const double aspect = 72.0 * m_vol[index] /
							(gradop12[ 0] * gradop12[ 0] +
							 gradop12[ 1] * gradop12[ 1] +
							 gradop12[ 2] * gradop12[ 2] +
							 gradop12[ 3] * gradop12[ 3] +
							 gradop12[ 4] * gradop12[ 4] +
							 gradop12[ 5] * gradop12[ 5] +
							 gradop12[ 6] * gradop12[ 6] +
							 gradop12[ 7] * gradop12[ 7] +
							 gradop12[ 8] * gradop12[ 8] +
							 gradop12[ 9] * gradop12[ 9] +
							 gradop12[10] * gradop12[10] +
							 gradop12[11] * gradop12[11] +
							 gradop12[12] * gradop12[12] +
							 gradop12[13] * gradop12[13] +
							 gradop12[14] * gradop12[14] +
							 gradop12[15] * gradop12[15] +
							 gradop12[16] * gradop12[16] +
							 gradop12[17] * gradop12[17] +
							 gradop12[18] * gradop12[18] +
							 gradop12[19] * gradop12[19] +
							 gradop12[20] * gradop12[20] +
							 gradop12[21] * gradop12[21] +
							 gradop12[22] * gradop12[22] +
							 gradop12[23] * gradop12[23] );

	const double aspect_inv = 1.0/aspect;
	//
	//  Compute the stable time step and bulk viscosity
	//
	const double dtrial = (sqrt(e_mass[index] * aspect / dilmod[index]));

	const double traced = rot_stret[r_index] + rot_stret[r_index + stride] + rot_stret[r_index];

	double eps;	

	if (traced < 0.0){
		eps = (linBulkVisc - quadBulkVisc * traced * dtrial);
	}
	else{
		eps = linBulkVisc;
  	}

	bulkq[index] = eps * dilmod[index] * dtrial * traced;

	const double cur_time_step = dtrial * ( sqrt( 1.0 + eps * eps) - eps);
	e_t_step[index] = cur_time_step;

	fac1[index] = fac1_pre * shrmod[index] * aspect_inv;
	fac2[index] = hg_visc * sqrt(shrmod[index] * e_mass[index] * aspect_inv);

}

__global__ void cuda_comp_force(	double dt, 
									const double * spin_rate, 
									const double * hgop_for_resist_calc, 
									const double * hgop, 
									const double * vel, 
									const double * fac_1, 
									const double * fac_2, 
									const double * bulk_q, 
									const double * rotated_stress, 
									double * hg_resist_old, 
									double * hg_resist_new, 
									const double * gradop12, 
									double * hg_energy, 
									double * int_energy, 
									double * force_new, 
									bool scaleHGRotation) {

	double bulkq = bulk_q[blockIdx.x * blockDim.x + threadIdx.x];
	double fac1 = fac_1[blockIdx.x * blockDim.x + threadIdx.x];
	double fac2 = fac_2[blockIdx.x * blockDim.x + threadIdx.x];
	const double * spin_rate_ptr = &spin_rate[blockIdx.x * blockDim.x * 3 + threadIdx.x];
	const double * rotated_stress_ptr = &rotated_stress[blockIdx.x * blockDim.x * 6 + threadIdx.x];
	const double * vel_ptr = &vel[blockIdx.x * blockDim.x * 24 + threadIdx.x];
	double * force_new_ptr = &force_new[blockIdx.x * blockDim.x * 24 + threadIdx.x];
	const double * gradop12_ptr = &gradop12[blockIdx.x * blockDim.x * 24 + threadIdx.x];
	double * hg_resist_old_ptr = &hg_resist_old[blockIdx.x * blockDim.x * 12 + threadIdx.x];
	double * hg_resist_new_ptr = &hg_resist_new[blockIdx.x * blockDim.x * 12 + threadIdx.x];
	const double * hgop_ptr = &hgop[blockIdx.x * blockDim.x * 32 + threadIdx.x];
	const double * hgop_for_resist_calc_ptr = &hgop_for_resist_calc[blockIdx.x * blockDim.x * 12 + threadIdx.x];


	double total_stress12th[6];

	total_stress12th[0] = ONE12TH * (rotated_stress_ptr[0 _N] + bulkq);
	total_stress12th[1] = ONE12TH * (rotated_stress_ptr[1 _N] + bulkq);
	total_stress12th[2] = ONE12TH * (rotated_stress_ptr[2 _N] + bulkq);
	total_stress12th[3] = ONE12TH * (rotated_stress_ptr[3 _N]);
	total_stress12th[4] = ONE12TH * (rotated_stress_ptr[4 _N]);
	total_stress12th[5] = ONE12TH * (rotated_stress_ptr[5 _N]);
	//
  	//  NKC, does Presto double need these spin rate terms?  Pronto appears to have dumped them....
  	//
  	const double dwxy = dt * spin_rate_ptr[0 _N];
  	const double dwyz = dt * spin_rate_ptr[1 _N];
  	const double dwzx = dt * spin_rate_ptr[2 _N];
  	//
  	//  Compute new hourglass resitance by the old rotated hourglass resitance plus a hourglass rate term
  	//
  	double hg_resist_total[12];
	double hg_temp[8];

  	if (!scaleHGRotation){
		for(int i = 0; i < 4; ++i) {

			double hg_rate_0 = 0, hg_rate_1 = 0, hg_rate_2 = 0;
	  		double hg_resist_old_0 = hg_resist_old_ptr[((3 * i) + 0) _N];
	  		double hg_resist_old_1 = hg_resist_old_ptr[((3 * i) + 1) _N];
	  		double hg_resist_old_2 = hg_resist_old_ptr[((3 * i) + 2) _N];

			for(int j = 0; j < 8; j++){

				hg_temp[j] = hgop_for_resist_calc_ptr[(i * 8 + j) _N];

			}

			for(int j = 0; j < 8; j++){

				hg_rate_0 += hg_temp[j] * vel_ptr[(j +  0) _N];
				hg_rate_1 += hg_temp[j] * vel_ptr[(j +  8) _N];
				hg_rate_2 += hg_temp[j] * vel_ptr[(j + 16) _N];

			}
			
	  		hg_resist_new_ptr[((i * 3) + 0) _N] = hg_resist_old_0 + dwxy*hg_resist_old_1 - dwzx*hg_resist_old_2 + fac1* hg_rate_0 ;
	  		hg_resist_total[(i * 3) + 0] 	 			= hg_resist_new_ptr[((i * 3) + 0) _N] + fac2* hg_rate_0 ;

		  	hg_resist_new_ptr[((i * 3) + 1) _N]	= hg_resist_old_1 - dwxy*hg_resist_old_0 + dwyz*hg_resist_old_2 + fac1* hg_rate_1 ;
		  	hg_resist_total[(i * 3) + 1]   				= hg_resist_new_ptr[((i * 3) + 1) _N] + fac2* hg_rate_1 ;

		  	hg_resist_new_ptr[((i * 3) + 2) _N]	= hg_resist_old_2 + dwzx*hg_resist_old_0 - dwyz*hg_resist_old_1 + fac1* hg_rate_2 ;
		  	hg_resist_total[(i * 3) + 2]   				= hg_resist_new_ptr[((i * 3) + 2) _N] + fac2* hg_rate_2 ;

		}
  	} 

	else {
		for(int i = 0; i < 4; ++i) {

			double hg_rate_0 = 0, hg_rate_1 = 0, hg_rate_2 = 0;
	  		double hg_resist_old_0 = hg_resist_old_ptr[((3 * i) + 0) _N];
	  		double hg_resist_old_1 = hg_resist_old_ptr[((3 * i) + 1) _N];
	  		double hg_resist_old_2 = hg_resist_old_ptr[((3 * i) + 2) _N];

			for(int j = 0; j < 8; j++){

				hg_temp[j] = hgop_for_resist_calc_ptr[(i * 8 + j) _N];

			}

			for(int j = 0; j < 8; j++){

				hg_rate_0 += hg_temp[j] * vel_ptr[(j +  0) _N];
				hg_rate_1 += hg_temp[j] * vel_ptr[(j +  8) _N];
				hg_rate_2 += hg_temp[j] * vel_ptr[(j + 16) _N];

			}

		  	const double rot_hg_resist_old_0 = hg_resist_old_0 + dwxy*hg_resist_old_1 - dwzx*hg_resist_old_2;
		  	const double rot_hg_resist_old_1 = hg_resist_old_1 - dwxy*hg_resist_old_0 + dwyz*hg_resist_old_2;
		  	const double rot_hg_resist_old_2 = hg_resist_old_2 + dwzx*hg_resist_old_0 - dwyz*hg_resist_old_1;

		  	double fnorm = rot_hg_resist_old_0 *rot_hg_resist_old_0  + rot_hg_resist_old_1 *rot_hg_resist_old_1  + rot_hg_resist_old_2 *rot_hg_resist_old_2 ;

		  	if (fnorm > 1.e-30){
				fnorm = sqrt ( (hg_resist_old_0*hg_resist_old_0 +
			                	hg_resist_old_1*hg_resist_old_1 +
			                	hg_resist_old_2*hg_resist_old_2) / fnorm );
				hg_resist_new_ptr[((i * 3) + 0) _N] = fnorm*rot_hg_resist_old_0 + fac1*hg_rate_0 ;
				hg_resist_new_ptr[((i * 3) + 1) _N] = fnorm*rot_hg_resist_old_1 + fac1*hg_rate_1 ;
				hg_resist_new_ptr[((i * 3) + 2) _N] = fnorm*rot_hg_resist_old_2 + fac1*hg_rate_2 ;
		  	} 

			else {
				hg_resist_new_ptr[((i * 3) + 0) _N] = rot_hg_resist_old_0 + fac1*hg_rate_0 ;
				hg_resist_new_ptr[((i * 3) + 1) _N] = rot_hg_resist_old_1 + fac1*hg_rate_1 ;
				hg_resist_new_ptr[((i * 3) + 2) _N] = rot_hg_resist_old_2 + fac1*hg_rate_2 ;
		  	}

		  	hg_resist_total[(i * 3) + 0] = hg_resist_new_ptr[((i * 3) + 0) _N] + fac2*hg_rate_0 ;
		  	hg_resist_total[(i * 3) + 1] = hg_resist_new_ptr[((i * 3) + 1) _N] + fac2*hg_rate_1 ;
		  	hg_resist_total[(i * 3) + 2] = hg_resist_new_ptr[((i * 3) + 2) _N] + fac2*hg_rate_2 ;

		}
  	}

  	double hg_force_0[8];
  	double hg_force_1[8];
  	double hg_force_2[8];


  	hg_energy[blockIdx.x * blockDim.x + threadIdx.x] = 0.0;
  	int_energy[blockIdx.x * blockDim.x + threadIdx.x] = 0.0;

  	for(int i = 0; i < 8; ++i) {
		hg_force_0[i] =
		   (hg_resist_total[HG_X1] * hgop_ptr[(i +  0) _N] +
		   	hg_resist_total[HG_X2] * hgop_ptr[(i +  8) _N] +
		    hg_resist_total[HG_X3] * hgop_ptr[(i + 16) _N] +
		    hg_resist_total[HG_X4] * hgop_ptr[(i + 24) _N]);

		hg_force_1[i] =
		   (hg_resist_total[HG_Y1] * hgop_ptr[(i +  0) _N] +
		    hg_resist_total[HG_Y2] * hgop_ptr[(i +  8) _N] +
		    hg_resist_total[HG_Y3] * hgop_ptr[(i + 16) _N] +
		    hg_resist_total[HG_Y4] * hgop_ptr[(i + 24) _N]);

		hg_force_2[i] =
		   (hg_resist_total[HG_Z1] * hgop_ptr[(i +  0) _N] +
		    hg_resist_total[HG_Z2] * hgop_ptr[(i +  8) _N] +
		    hg_resist_total[HG_Z3] * hgop_ptr[(i + 16) _N] +
		    hg_resist_total[HG_Z4] * hgop_ptr[(i + 24) _N]);
  	}


  
	for(int i = 0; i < 8; ++i) {
		force_new_ptr[(i +  0) _N] =
		  	total_stress12th[K_S_XX] * gradop12_ptr[(i +  0) _N] +
		  	total_stress12th[K_S_XY] * gradop12_ptr[(i +  8) _N] +
		  	total_stress12th[K_S_XZ] * gradop12_ptr[(i + 16) _N] + hg_force_0[i] ;

		force_new_ptr[(i +  8) _N] =
		  	total_stress12th[K_S_YX] * gradop12_ptr[(i +  0) _N] +
		 	total_stress12th[K_S_YY] * gradop12_ptr[(i +  8) _N] +
		  	total_stress12th[K_S_YZ] * gradop12_ptr[(i + 16) _N] + hg_force_1[i] ;

		force_new_ptr[(i + 16) _N] =
		  	total_stress12th[K_S_ZX] * gradop12_ptr[(i +  0) _N] +
		  	total_stress12th[K_S_ZY] * gradop12_ptr[(i +  8) _N] +
		  	total_stress12th[K_S_ZZ] * gradop12_ptr[(i + 16) _N] + hg_force_2[i] ;

		hg_energy[blockIdx.x * blockDim.x + threadIdx.x]  += 	hg_force_0[i]   *vel_ptr[(i +  0) _N] + \
																hg_force_1[i]   *vel_ptr[(i +  8) _N] + \
																hg_force_2[i]   *vel_ptr[(i + 16) _N];

		int_energy[blockIdx.x * blockDim.x + threadIdx.x] += 	force_new_ptr[(i +  0) _N] * vel_ptr[(i +  0) _N] + \
																force_new_ptr[(i +  8) _N] * vel_ptr[(i +  8) _N] + \
																force_new_ptr[(i + 16) _N] * vel_ptr[(i + 16) _N];

	}

}


inline int cuda_elem_ug3dh8_mi_compute_stretch(	int nelem,
					                              	double dt,
					                              	double * cordel,
													double * m_cordel,
					                              	double * vel,
					                              	double * rotation_old,
					                              	double * mid_vol,
					                              	double * vorticity,
					                              	double * rotation_new,
					                              	double * stretch,
					                              	double * rotated_stretching,
					                              	double * mid_hgop) {
	cudaError_t err;
	int return_value = 0;

/*********************************************************/

	int size_24 = sizeof(double) * nelem * 24;
	int size_9  = sizeof(double) * nelem *  9;
	int size_6  = sizeof(double) * nelem *  6;

	double * gradop12;
		err = cudaMalloc((void**)&gradop12, size_24);
		CheckCudaError(err);
	double * vgrad;
		err = cudaMalloc((void**)&vgrad, size_9);
		CheckCudaError(err);
	double * stretching_tensor;
		err = cudaMalloc((void**)&stretching_tensor, size_6);
		CheckCudaError(err);

/*********************************************************/

	dim3 threads_per_block, blocks;
	
	threads_per_block.x = 256;
	threads_per_block.y = 1;
	threads_per_block.z = 1;

	blocks.x = nelem / threads_per_block.x;
	blocks.y = 1;
	blocks.z = 1;

/*********************************************************/

	cuda_comp_grad_simple<<<blocks, threads_per_block>>>(cordel, gradop12, mid_vol, dt);
	cuda_v_grad<<<blocks, threads_per_block>>>(vel, gradop12, mid_vol, vgrad);
	cuda_additive_decomp<<<blocks, threads_per_block>>>(gradop12, stretching_tensor, vorticity);
	cuda_polar_decomp<<<blocks, threads_per_block>>>(dt, stretching_tensor, stretch, rotation_old, vorticity, rotation_new);
	cuda_rotate_tensor_forward<<<blocks, threads_per_block>>>(rotation_new, stretching_tensor, rotated_stretching);
	cuda_comp_hgop3<<<blocks, threads_per_block>>>(m_cordel, mid_vol, gradop12, mid_hgop);
	cudaThreadSynchronize();

/*********************************************************/

	cudaFree(gradop12);
	cudaFree(vgrad);
	cudaFree(stretching_tensor);


	return return_value;

}

inline int cuda_moduli(	int nelem,	
							double * dilmod, 
							double * shrmod, 
							double * bulkmod, 
							double * twomu){



/*********************************************************/

	dim3 threads_per_block, blocks;
	
	threads_per_block.x = 256;
	threads_per_block.y = 1;
	threads_per_block.z = 1;

	blocks.x = nelem / threads_per_block.x;
	blocks.y = 1;
	blocks.z = 1;

/*********************************************************/

	cuda_comp_moduli<<<blocks, threads_per_block>>>(dilmod, shrmod, bulkmod, twomu);

	return 0;

}


inline int cuda_elem_ug3dh8_mi_compute_divergence_presto(	const int nelem,
						                                 	const double dt,
						                                 	double * cordel,
						                                 	const double * vel,
						                                 	double * rotation,
						                                 	double * stress_new,
						                                 	double * rotated_stretching,
						                                 	const double * spin_rate,
						                                 	const double linBulkVisc,
						                                 	const double quadBulkVisc,
						                                 	double * elem_mass,
						                                 	double * elem_dilmod,
						                                 	double * elem_shrmod,
						                                 	const double hg_stiffness,
						                                 	const double hg_viscosity,
						                                 	double * rotated_stress,
						                                 	double &min_elem_time_step,
						                                 	double *const volume,
						                                 	double *const elem_time_step,
						                                 	double *const hg_resist_old,
						                                 	double *const hg_resist_new,
						                                 	double *const force_new,
						                                 	double *const hg_energy,
						                                 	double *const int_energy,
						                                 	double *const mid_hgop,
						                                 	const bool scaleHGRotation) {

	cudaError_t err;
	int return_value = 0;

  	const double fac1_pre = dt * hg_stiffness * 0.0625;

/*********************************************************/

	int size_32 = sizeof(double) * nelem *  32;
	int size_24 = sizeof(double) * nelem *  24;
	int size_6  = sizeof(double) * nelem *   6;
	int size_1  = sizeof(double) * nelem *   1;

	double * gradop12;
		err = cudaMalloc((void**)&gradop12, size_24);
		CheckCudaError(err);
	double * hgop;
		err = cudaMalloc((void**)&hgop, size_32);
		CheckCudaError(err);
	double * fac1;
		err = cudaMalloc((void**)&fac1, size_1);
		CheckCudaError(err);
	double * fac2;
		err = cudaMalloc((void**)&fac2, size_1);
		CheckCudaError(err);
	double * bulkq;
		err = cudaMalloc((void**)&bulkq, size_1);
		CheckCudaError(err);
	double * total_stress12th;
		err = cudaMalloc((void**)&total_stress12th, size_6);
		CheckCudaError(err);

/*********************************************************/

	dim3 threads_per_block, blocks;
	
	threads_per_block.x = 256;
	threads_per_block.y = 1;
	threads_per_block.z = 1;

	blocks.x = nelem / threads_per_block.x;
	blocks.y = 1;
	blocks.z = 1;


/*********************************************************/

	cuda_comp_grad_simple<<<blocks, threads_per_block>>> (	cordel, 
															gradop12, 
															volume, 
															dt);


	cuda_comp_facs<<<blocks, threads_per_block>>>(	gradop12, 
													volume, 
													elem_mass, 
													elem_time_step, 
													rotated_stretching, 
													elem_shrmod, 
													elem_dilmod, 
													fac1, 
													fac2, 
													bulkq,
													fac1_pre, 
													linBulkVisc, 
													quadBulkVisc, 
													hg_viscosity );


	cuda_rotate_tensor_backward<<<blocks, threads_per_block>>>(	rotation,
																stress_new, 
																rotated_stress);

	cuda_comp_hgop3<<<blocks, threads_per_block>>>(	cordel, 
													volume, 
													gradop12, 
													hgop);


	cuda_comp_force<<<blocks, threads_per_block>>>(	dt, 
													spin_rate, 
													mid_hgop, 
													hgop, 
													vel, 
													fac1, 
													fac2, 
													bulkq, 
													rotated_stress, 
													hg_resist_old, 
													hg_resist_new, 
													gradop12, 
													hg_energy, 
													int_energy, 
													force_new, 
													scaleHGRotation);

	cudaThreadSynchronize();

/*********************************************************/

	cudaFree(gradop12);
	cudaFree(hgop);
	cudaFree(fac1);
	cudaFree(fac2);
	cudaFree(bulkq);
	cudaFree(total_stress12th);

	return return_value;

}


double cuda_internalForce(		const int num_elements, 
                               	const double dt,
                               	double current_stable_time_step,
                               	double * const element_time_step,
                               	double * const coordinates,
								double * mid_coordinates,
                               	double * const velocity,
                               	double * rotation_old, 
								double * rotation_new, 
                               	double * const midstep_volume, 
                               	double * const vorticity_tensor, 
                               	double * const stretch,
                               	double * const strain_rate, 
                               	double * const mid_hgop, 
                               	double * stress_new,
                               	double * rotated_stress,
                               	double * const material_eff_bulk_mod,
                               	double * const material_eff_twomu,
                               	double * shrmod,
                               	double * dilmod,
                               	double * element_mass,
                               	double * const force_new,
                               	double * const hourglass_energy,
                               	double * const internal_energy,
                               	double * const hg_resistance_old, 
								double * const hg_resistance_new){

/*********************************************************/

	timeval start,stop,result;


	gettimeofday(&start, NULL);

	cuda_fill<<<num_elements / 256, 256>>>(coordinates, velocity, rotation_old);

	gettimeofday(&stop, NULL);
	timersub(&stop, &start, &result);
	double time = (result.tv_sec + result.tv_usec/1000000.0);

/*********************************************************/

	int err = cuda_elem_ug3dh8_mi_compute_stretch(	num_elements, 
													dt, 
													coordinates, 
													mid_coordinates,
													velocity, 
													rotation_old, 
													midstep_volume, 
													vorticity_tensor, 
													rotation_new, 
													stretch, 
													strain_rate, 
													mid_hgop);

/*********************************************************/

	err = cuda_moduli(	num_elements,
						dilmod, 
						shrmod, 
						material_eff_bulk_mod, 
						material_eff_twomu);

/*********************************************************/

    double dilatationalHGParam = 0.05;
    double deviatoricHGParam   = 0.0;
    const bool scaleHGRotation            = false; // default?
    const double linear_bulk_viscosity    = 0.0; // default?
    const double quadratic_bulk_viscosity = 0.0; // default?

/*********************************************************/

	err = cuda_elem_ug3dh8_mi_compute_divergence_presto(	num_elements, 
															dt, 
															coordinates, 
															velocity, 
															rotation_new, 
															stress_new, 
															strain_rate, 
															vorticity_tensor, 
															linear_bulk_viscosity, 
															quadratic_bulk_viscosity, 
															element_mass, 
															dilmod, 
															shrmod, 
															dilatationalHGParam, 
															deviatoricHGParam, 
															rotated_stress, 
															current_stable_time_step, 
															midstep_volume, 
															element_time_step, 
															hg_resistance_old, 
															hg_resistance_new, 
															force_new, 
															hourglass_energy, 
															internal_energy, 
															mid_hgop, 
															scaleHGRotation);


/*********************************************************/


	return time;

}


double run_cuda_kernel(int n){

	if(n > 2400000){ //capacity of Tesla memory
		std::cout << "Nope.avi\n";
		return 1337;
	}

	timeval start,stop,result;
	cudaError_t err = cudaSuccess;

	int size_32 = sizeof(double) * n * 32;
	int size_24 = sizeof(double) * n * 24;
	int size_12 = sizeof(double) * n * 12;
	int size_9 = sizeof(double)  * n *  9;
	int size_6 = sizeof(double)  * n *  6;
	int size_3 = sizeof(double)  * n *  3;
	int size_1 = sizeof(double)  * n *  1;

	int nelem = n;

/*********************************************************/

  	double * position;	
		err = cudaMalloc((void**)&position, size_24);
		CheckCudaError(err);
  	double * m_position;	
		err = cudaMalloc((void**)&m_position, size_24);
		CheckCudaError(err);	
	double * velocity;
		err = cudaMalloc((void**)&velocity, size_24);
		CheckCudaError(err);
	double * mid_pos;
		err = cudaMalloc((void**)&mid_pos, size_24);
		CheckCudaError(err);
  	double * rot_old;
		err = cudaMalloc((void**)&rot_old, size_9);
		CheckCudaError(err);
  	double * rot_new;
		err = cudaMalloc((void**)&rot_new, size_9);
		CheckCudaError(err);			
	double * gradop12;
		err = cudaMalloc((void**)&gradop12, size_24);
		CheckCudaError(err);	
	double * vel_grad;
		err = cudaMalloc((void**)&vel_grad, size_9);
		CheckCudaError(err);	
	double * stretch;
		err = cudaMalloc((void**)&stretch, size_6);
		CheckCudaError(err);	
	double * vorticity;
		err = cudaMalloc((void**)&vorticity, size_3);
		CheckCudaError(err);	
	double * rot_stret;
		err = cudaMalloc((void**)&rot_stret, size_6);
		CheckCudaError(err);	
	double * mid_vol;
		err = cudaMalloc((void**)&mid_vol, size_1);
		CheckCudaError(err);	
	double * hgop;		
		err = cudaMalloc((void**)&hgop, size_32);
		CheckCudaError(err);	
	double * strain_rate;		
		err = cudaMalloc((void**)&strain_rate, size_6);
		CheckCudaError(err);
	double * shrmod;		
		err = cudaMalloc((void**)&shrmod, size_1);
		CheckCudaError(err);
	double * dilmod;
		err = cudaMalloc((void**)&dilmod, size_1);
		CheckCudaError(err);
	double * elem_mass;
		err = cudaMalloc((void**)&elem_mass, size_1);
		CheckCudaError(err);
	double * elem_t_step;
		err = cudaMalloc((void**)&elem_t_step, size_1);
		CheckCudaError(err);
	double * force_new;
		err = cudaMalloc((void**)&force_new, size_24);
		CheckCudaError(err);
	double * hg_energy;
		err = cudaMalloc((void**)&hg_energy, size_1);
		CheckCudaError(err);
	double * hg_resist_n;
		err = cudaMalloc((void**)&hg_resist_n, size_12);
		CheckCudaError(err);
	double * hg_resist_o;
		err = cudaMalloc((void**)&hg_resist_o, size_12);
		CheckCudaError(err);
	double * intern_energy;
		err = cudaMalloc((void**)&intern_energy, size_1);
		CheckCudaError(err);
	double * stress_new;	
		err = cudaMalloc((void**)&stress_new, size_6);
		CheckCudaError(err);		
	double * rot_stress;
		err = cudaMalloc((void**)&rot_stress, size_6);
		CheckCudaError(err);
	double * bulk_modulus;
		err = cudaMalloc((void**)&bulk_modulus, size_1);
		CheckCudaError(err);
	double * two_mu;
		err = cudaMalloc((void**)&two_mu, size_1);
		CheckCudaError(err);
	double * s_temp;
		err = cudaMalloc((void**)&s_temp, size_6);
		CheckCudaError(err);

/*********************************************************/

	double dt = 0.25;
	double	stable_time_step = 0;

/*********************************************************/

	gettimeofday(&start, NULL);

	double fill_time = cuda_internalForce(	nelem, 
											dt, 
											stable_time_step, 
											elem_t_step, 
											position,
											m_position, 
											velocity, 
											rot_old, 
											rot_new, 
											mid_vol, 
											vorticity, 
											stretch, 
											strain_rate,
											hgop, 
											stress_new,
											rot_stret,
											bulk_modulus,
											two_mu,
											shrmod,
											dilmod,
											elem_mass,
											force_new,
											hg_energy,
											intern_energy,
											hg_resist_o,
											hg_resist_n);

	gettimeofday(&stop, NULL);
	timersub(&stop, &start, &result);
	double time = (result.tv_sec + result.tv_usec/1000000.0);

/***********************************************************/

	cudaFree(position);
	cudaFree(velocity);
	cudaFree(mid_pos);
	cudaFree(rot_old);
	cudaFree(rot_new);
	cudaFree(gradop12);
	cudaFree(vel_grad);
	cudaFree(stretch);
	cudaFree(vorticity);
	cudaFree(rot_stret);
	cudaFree(mid_vol);
	cudaFree(hgop);
	cudaFree(strain_rate);
	cudaFree(shrmod);
	cudaFree(dilmod);
	cudaFree(elem_mass);
	cudaFree(elem_t_step);
	cudaFree(force_new);
	cudaFree(hg_energy);
	cudaFree(hg_resist_o);
	cudaFree(hg_resist_n);
	cudaFree(intern_energy);
	cudaFree(stress_new);
	cudaFree(rot_stress);
	cudaFree(bulk_modulus);
	cudaFree(two_mu);
	cudaFree(s_temp);

	return time - fill_time;

}
