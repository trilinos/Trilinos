// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TestStochastic.hpp"

#include "Kokkos_Core.hpp"

#include "Stokhos_Cuda_CrsMatrix.hpp"
#include "Stokhos_Cuda_BlockCrsMatrix.hpp"
#include "Stokhos_Cuda_StochasticProductTensor.hpp"
#include "Stokhos_Cuda_SymmetricDiagonalSpec.hpp"
#include "Stokhos_Cuda_CrsProductTensor.hpp"
#include "Stokhos_Cuda_TiledCrsProductTensor.hpp"
#include "Stokhos_Cuda_SimpleTiledCrsProductTensor.hpp"
#include "Stokhos_Cuda_CooProductTensor.hpp"
#include "Stokhos_Cuda_FlatSparse3Tensor.hpp"
#include "Stokhos_Cuda_FlatSparse3Tensor_kji.hpp"
#include "Stokhos_Cuda_LexicographicBlockSparse3Tensor.hpp"
#include "Stokhos_Cuda_LinearSparse3Tensor.hpp"

namespace unit_test {

template<typename Scalar>
struct performance_test_driver<Scalar,Kokkos::Cuda> {
  static void run(bool test_flat, bool test_orig, bool test_lin,
                  bool test_block, bool symmetric) {
    typedef Kokkos::Cuda Device;

    int nGrid;
    int nIter;

    // All methods compared against flat-original
    if (test_flat) {
      nGrid = 5 ;
      nIter = 1 ;
      performance_test_driver_all<Scalar,Device>(
        3 , 1 ,  9 , nGrid , nIter , test_block , symmetric );
      performance_test_driver_all<Scalar,Device>(
        5 , 1 ,  5 , nGrid , nIter , test_block , symmetric );
    }

    // Just polynomial methods compared against original
    if (test_orig) {
      nGrid = 32 ;
      nIter = 1 ;
      performance_test_driver_poly<Scalar,Device,Stokhos::DefaultMultiply>(
        3 , 1 , 12 , nGrid , nIter , test_block , symmetric );
      performance_test_driver_poly<Scalar,Device,Stokhos::DefaultMultiply>(
        5 , 1 ,  6 , nGrid , nIter , test_block , symmetric );
    }

     // Just polynomial methods compared against original
    if (test_lin) {
      nGrid = 32 ;
      nIter = 10 ;
      performance_test_driver_linear<Scalar,Device,Stokhos::DefaultMultiply>(
        31 ,  255 , 32 , nGrid , nIter , test_block , symmetric );
    }

  }

};

}

template <typename Scalar>
int mainCuda(bool test_flat, bool test_orig, bool test_lin, bool test_block,
             bool symmetric, int device_id)
{
  typedef unsigned long long int IntType ;

  Kokkos::InitializationSettings init_args;
  init_args.set_device_id(device_id);
  Kokkos::initialize( init_args );
  Kokkos::print_configuration( std::cout );

  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, device_id);
  std::cout << std::endl
            << "Device " << device_id << ": " << deviceProp.name
            << std::endl;

  cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

  std::cout << std::endl << "\"Cuda Performance\"" << std::endl ;
  unit_test::performance_test_driver<Scalar,Kokkos::Cuda>::run(
    test_flat, test_orig, test_lin, test_block, symmetric);

  Kokkos::finalize();

  cudaDeviceReset();

  return 0 ;
}

template int mainCuda<float>(bool, bool, bool, bool, bool, int);
template int mainCuda<double>(bool, bool, bool, bool, bool, int);
