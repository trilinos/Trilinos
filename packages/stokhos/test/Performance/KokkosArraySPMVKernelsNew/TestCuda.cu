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

#include "TestStochastic.hpp"

#include "KokkosArray_Cuda.hpp"

#include "Stokhos_Cuda_CrsMatrix.hpp"
#include "Stokhos_Cuda_BlockCrsMatrix.hpp"
#include "Stokhos_Cuda_StochasticProductTensor.hpp"
#include "Stokhos_Cuda_CrsProductTensor.hpp"
#include "Stokhos_Cuda_FlatSparse3Tensor.hpp"
#include "Stokhos_Cuda_FlatSparse3Tensor_kji.hpp"

namespace unit_test {

template<typename Scalar>
struct performance_test_driver<Scalar,KokkosArray::Cuda> {
  static void run(bool test_flat, bool test_orig, bool test_block, 
		  bool symmetric) {
    typedef KokkosArray::Cuda Device;
    
    int nGrid;
    int nIter; 

    // All methods compared against flat-original
    // if (test_flat) {
    //   nGrid = 5 ;
    //   nIter = 1 ; 
    //   performance_test_driver_all<Scalar,Device>( 
    // 	3 , 1 ,  9 , nGrid , nIter , test_block );
    //   performance_test_driver_all<Scalar,Device>( 
    // 	5 , 1 ,  5 , nGrid , nIter , test_block );
    // }
    
    // Just polynomial methods compared against original
    if (test_orig) {
      nGrid = 32 ;
      nIter = 1 ; 
      performance_test_driver_poly<Scalar,Device,Stokhos::DefaultSparseMatOps>( 
	3 , 1 , 12 , nGrid , nIter , test_block , symmetric );
      performance_test_driver_poly<Scalar,Device,Stokhos::DefaultSparseMatOps>( 
	5 , 1 ,  6 , nGrid , nIter , test_block , symmetric );
    }
    
  }

};

}

template <typename Scalar>
int mainCuda(bool test_flat, bool test_orig, bool test_block, bool symmetric, 
	     int device_id)
{
  typedef unsigned long long int IntType ;

  KokkosArray::Cuda::initialize( KokkosArray::Cuda::SelectDevice(0) );

  cudaSetDevice(device_id);
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, device_id);
  std::cout << std::endl 
	    << "Device " << device_id << ": " << deviceProp.name 
	    << std::endl;

  std::cout << std::endl << "\"Cuda Performance\"" << std::endl ;
  unit_test::performance_test_driver<Scalar,KokkosArray::Cuda>::run(
    test_flat, test_orig, test_block, symmetric);

  KokkosArray::Cuda::finalize();

  return 0 ;
}

template int mainCuda<float>(bool, bool, bool, bool, int);
template int mainCuda<double>(bool, bool, bool, bool, int);
