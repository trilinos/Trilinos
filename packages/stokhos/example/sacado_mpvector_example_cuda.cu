// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "stdio.h"
#include "Stokhos_Sacado_Kokkos.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "sacado_mpvector_example.hpp"

// The function to compute the polynomial chaos expansion of,
// written as a template function
template <class ScalarType>
__host__ __device__
void simple_function(const ScalarType& x, ScalarType& y) {
  ScalarType u = x*x;
  ScalarType v = std::pow(std::log(u),2.0);
  y = v/(x + 1.0);
}

template <typename vector_type>
__global__ 
void mpkernel(int n, int sz, double *dev_x, double *dev_y, 
	      bool reset, bool print) 
{
  typedef typename vector_type::storage_type storage_type;

  int threads_per_block = blockDim.x*blockDim.y*blockDim.z;

  int block = blockIdx.x + (blockIdx.y + blockIdx.z*gridDim.y)*gridDim.x;
  int local_thread = 
    threadIdx.x + (threadIdx.y + threadIdx.z*blockDim.y)*blockDim.x;

  // multi-point expansions
  vector_type x(sz, 0.0), y(sz, 0.0);

  // Loop over elements, propagating sz samples simultaneously
  for (int e=0; e<n; e++) {

    // Initialize x
    int offset = (block*n+e)*sz*threads_per_block + local_thread;
    int stride = threads_per_block;

    if (reset && storage_type::supports_reset) {
      storage_type& x_s = x.storage();
      storage_type& y_s = y.storage();
      x_s.shallowReset(dev_x+offset, sz, stride, false);
      y_s.shallowReset(dev_y+offset, sz, stride, false);
    }
    else {
      for (int i=0; i<sz; i++)
	x.fastAccessCoeff(i) = dev_x[offset+i*stride];
    }

    simple_function(x,y);

    // Print x and y
    if (print) {
      for (int tidx = 0; tidx<threads_per_block; tidx++) {
	if (local_thread == tidx) {
	  printf("x(%i) = [ ",tidx);
	  for (int i=0; i<sz; i++)
	    printf("%g ", x.coeff(i));
	  printf("]\n\n");
	}
	__syncthreads();
      }
      
      for (int tidx = 0; tidx<threads_per_block; tidx++) {
	if (local_thread == tidx) {
	  printf("y(%i) = [ ",tidx);
	  for (int i=0; i<sz; i++)
	    printf("%g ", y.coeff(i));
	  printf("]\n\n");
	}
	__syncthreads();
      }
    }

    // Return result
    if (!(reset && vector_type::storage_type::supports_reset)) {
      for (int i=0; i<sz; i++)
    	dev_y[offset+i*stride] = y.fastAccessCoeff(i);
    }

  }
}

__global__ 
void kernel(int n, int sz, double *dev_x, double *dev_y) 
{
  int threads_per_block = blockDim.x*blockDim.y*blockDim.z;
  int block = blockIdx.x + (blockIdx.y + blockIdx.z*gridDim.y)*gridDim.x;
  int local_thread = 
    threadIdx.x + (threadIdx.y + threadIdx.z*blockDim.y)*blockDim.x;

  // Loop over elements
  double x, y;
  for (int e=0; e<n; e++) {
    int offset = (block*n+e)*sz*threads_per_block + local_thread;
    int stride = threads_per_block;

    for (int i=0; i<sz; i++) {

      // Initialize x
      x = dev_x[offset+i*stride];

      // Compute function
      simple_function(x,y);

      // Return result
      dev_y[offset+i*stride] = y;

    }

  }
}

// Partial specialization of vector example runner for CUDA
template <int MaxSize>
struct MPVectorExample<MaxSize,KokkosArray::Cuda> {
   typedef KokkosArray::Cuda node_type;

  static bool 
  run(Storage_Method storage_method, int n, int sz, int nblocks, int nthreads, 
      bool reset, bool print) {
    typedef MPVectorTypes<MaxSize, node_type> MPT;

    bool status;
    if (storage_method == STATIC)
      status = run_impl<typename MPT::static_vector>(
    	n, sz, nblocks, nthreads, reset, print);
    else if (storage_method == STATIC_FIXED)
      status = run_impl<typename MPT::static_fixed_vector>(
	n, sz, nblocks, nthreads, reset, print);
    else if (storage_method == LOCAL)
      status = run_impl<typename MPT::local_vector>(
    	n, sz, nblocks, nthreads, reset, print);
    else if (storage_method == DYNAMIC)
      status = run_impl<typename MPT::dynamic_vector>(
    	n, sz, nblocks, nthreads, reset, print);
    else if (storage_method == DYNAMIC_STRIDED)
      status = run_impl<typename MPT::dynamic_strided_vector>(
    	n, sz, nblocks, nthreads, reset, print);
    else if (storage_method == DYNAMIC_THREADED)
      status = run_impl<typename MPT::dynamic_threaded_vector>(
    	n, sz, nblocks, nthreads, reset, print);
    
    return status;
  }

private:

  template <typename vector_type>
  static bool 
  run_impl(int n, int sz, int nblocks, int nthreads, bool reset, bool print) {

    // Setup CUDA thread blocks
    dim3 grid(nblocks, 1, 1);
    dim3 block(nthreads, 1, 1);
    
    // Allocate memory inputs and outpus
    int len = nblocks*nthreads*sz*n;
    std::cout << "total size = " << len << std::endl;
    
    double *x = new double[len];
    double *y = new double[len];
    double *y_mp = new double[len];
    double *dev_x, *dev_y;
    cudaMalloc( (void**)&dev_x, len*sizeof(double) );
    cudaMalloc( (void**)&dev_y, len*sizeof(double) );
    
    // Initialize x
    for (int i=0; i<len; i++)
      x[i] = static_cast<double>(i+1)/static_cast<double>(len);
    
    // Transfer input data to device
    {
      TEUCHOS_FUNC_TIME_MONITOR("Host -> Device Transfer");
      cudaMemcpy(dev_x, x, len*sizeof(double), cudaMemcpyHostToDevice);
    }
    
    // Invoke kernel
    {
      TEUCHOS_FUNC_TIME_MONITOR("Device calculation");
      kernel<<<grid,block>>>(n, sz, dev_x, dev_y);
      cudaDeviceSynchronize();
    }
    
    // Transfer results from device
    {
      TEUCHOS_FUNC_TIME_MONITOR("Device -> Host Transfer");
      cudaMemcpy(y, dev_y, len*sizeof(double), cudaMemcpyDeviceToHost);
    }
    
    // Invoke kernel
    {
      TEUCHOS_FUNC_TIME_MONITOR("Device calculation (MP)");
      
      mpkernel<vector_type><<<grid,block>>>(n, sz, dev_x, dev_y, reset, print);
      cudaDeviceSynchronize();
    }
    
    // Transfer results from device
    {
      TEUCHOS_FUNC_TIME_MONITOR("Device -> Host Transfer (MP)");
      cudaMemcpy(y_mp, dev_y, len*sizeof(double), cudaMemcpyDeviceToHost);
    }
    
    // Check results agree
    double rtol = 1e-15;
    double atol = 1e-15;
    bool agree = true;
    for (int i=0; i<len; i++) {
      if (std::abs(y[i]-y_mp[i]) > std::abs(y[i])*rtol + atol) {
	agree = false;
	break;
      }
    }
    
    if (print) {
      std::cout << "x    = [ ";
      for (int i=0; i<len; i++)
	std::cout << x[i] << " ";
      std::cout << "]" << std::endl;
      
      std::cout << "y      [ ";
      for (int i=0; i<len; i++)
	std::cout << y[i] << " ";
      std::cout << "]" << std::endl;
      
      std::cout << "y_mp = [ ";
      for (int i=0; i<len; i++)
	std::cout << y_mp[i] << " ";
      std::cout << "]" << std::endl;
    }
    
    // Clean up memory
    delete [] x;
    delete [] y;
    delete [] y_mp;
    cudaFree(dev_x);
    cudaFree(dev_y);
  
    return agree;
  }

};

template struct MPVectorExample<MaxSize, KokkosArray::Cuda>;
