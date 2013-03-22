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

#include "Stokhos_Sacado_Kokkos.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "sacado_mpvector_example.hpp"

// The function to compute the polynomial chaos expansion of,
// written as a template function
template <class ScalarType>
void simple_function(const ScalarType& x, ScalarType& y) {
  ScalarType u = x*x;
  ScalarType v = std::pow(std::log(u),2.0);
  y = v/(x + 1.0);
}

template <typename vector_type>
void mpkernel(int offset, int stride, int n, int sz, 
	      double *dev_x, double *dev_y, 
	      bool reset, bool print) 
{
  typedef typename vector_type::storage_type storage_type;

  // multi-point expansions
  vector_type x(sz, 0.0), y(sz, 0.0);

  // Loop over elements, propagating sz samples simultaneously
  for (int e=0; e<n; e++) {

    // Initialize x
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
      std::cout << "x(0) = [ ";
      for (int i=0; i<sz; i++)
	std::cout << x.coeff(i) << " ";
      std::cout << "]" << std::endl << std::endl;
      
      std::cout << "y(0) = [ ";
      for (int i=0; i<sz; i++)
	std::cout << y.coeff(i) << " ";
      std::cout << "]" << std::endl << std::endl;
    }

    // Return result
    if (!(reset && vector_type::storage_type::supports_reset)) {
      for (int i=0; i<sz; i++)
    	dev_y[offset+i*stride] = y.fastAccessCoeff(i);
    }

    offset += stride*sz;

  }
}

void kernel(int offset, int stride, int n, int sz, double *dev_x, double *dev_y)
{
  // Loop over elements
  double x, y;
  for (int e=0; e<n; e++) {
    
    for (int i=0; i<sz; i++) {

      // Initialize x
      x = dev_x[offset+i*stride];

      // Compute function
      simple_function(x,y);

      // Return result
      dev_y[offset+i*stride] = y;

    }

    offset += stride*sz;

  }
}

// Partial specialization of vector example runner for CUDA
template <int MaxSize>
struct MPVectorExample<MaxSize,KokkosArray::Host> {
  typedef KokkosArray::Host node_type;

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
    else if (storage_method == DYNAMIC_THREADED) {
      std::cout << "Host node doesn't support dynamic-threaded storage!"
    		<< std::endl;
      status = false;
    }
    
    return status;
  }

private:

  template <typename vector_type>
  static bool 
  run_impl(int n, int sz, int nblocks, int nthreads, bool reset, bool print) {

    // Allocate memory inputs and outpus
    int len = nblocks*nthreads*sz*n;
    std::cout << "total size = " << len << std::endl;

    int stride = 1;
    
    double *x = new double[len];
    double *y = new double[len];
    double *y_mp = new double[len];
    
    // Initialize x
    for (int i=0; i<len; i++)
      x[i] = static_cast<double>(i+1)/static_cast<double>(len);
    
    // Invoke kernel
    {
      TEUCHOS_FUNC_TIME_MONITOR("Host calculation");
      for (int offset=0; offset<len; offset += sz*n)
	kernel(offset, stride, n, sz, x, y);
    }
    
    // Invoke kernel
    {
      TEUCHOS_FUNC_TIME_MONITOR("Host calculation (MP)");
      for (int offset=0; offset<len; offset += sz*n)
	mpkernel<vector_type>(offset, stride, n, sz, x, y_mp, reset, print);
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
  
    return agree;
  }

};

template struct MPVectorExample<MaxSize, KokkosArray::Host>;
