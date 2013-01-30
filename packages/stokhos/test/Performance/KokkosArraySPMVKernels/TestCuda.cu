/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <TestStochastic.hpp>

#include <KokkosArray_Cuda.hpp>
#include <Host/KokkosArray_Host_ProductTensor.hpp>
#include <Cuda/KokkosArray_Cuda_SymmetricDiagonalSpec.hpp>
#include <Cuda/KokkosArray_Cuda_ProductTensor.hpp>
#include <Cuda/KokkosArray_Cuda_CrsProductTensorLegendre.hpp>
#include <Cuda/KokkosArray_Cuda_StochasticProductTensor.hpp>
#include <Cuda/KokkosArray_Cuda_BlockCrsMatrix.hpp>
#include <Cuda/KokkosArray_Cuda_CrsMatrix.hpp>

namespace unit_test {

template<typename Scalar>
struct performance_test_driver<Scalar,KokkosArray::Cuda> {
  static void run(bool test_flat, bool test_orig, bool test_block) {
    typedef KokkosArray::Cuda Device;
    
    int nGrid;
    int nIter; 

    // All methods compared against flat-original
    if (test_flat) {
      nGrid = 5 ;
      nIter = 1 ; 
      performance_test_driver_all<Scalar,Device>( 
	3 , 1 ,  9 , nGrid , nIter , test_block );
      performance_test_driver_all<Scalar,Device>( 
	5 , 1 ,  5 , nGrid , nIter , test_block );
    }
    
    // Just polynomial methods compared against original
    if (test_orig) {
      nGrid = 32 ;
      nIter = 1 ; 
      performance_test_driver_poly<Scalar,Device>( 
	3 , 1 , 12 , nGrid , nIter , test_block );
      performance_test_driver_poly<Scalar,Device>( 
	5 , 1 ,  6 , nGrid , nIter , test_block );
    }
    
  }

};

}

template <typename Scalar>
int mainCuda(bool test_flat, bool test_orig, bool test_block, int device_id)
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
    test_flat, test_orig, test_block);

  KokkosArray::Cuda::finalize();

  return 0 ;
}

template int mainCuda<float>(bool, bool, bool, int);
template int mainCuda<double>(bool, bool, bool, int);
