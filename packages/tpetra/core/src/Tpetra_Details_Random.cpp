// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_Random.hpp"
#include "Teuchos_TestForException.hpp"

namespace Tpetra {
namespace Details {

namespace { // (anonymous)

unsigned int getSeedFromRank(int mpi_rank) {
  // Seed the pseudorandom number generator using the calling
  // process' rank.  This helps decorrelate different process'
  // pseudorandom streams.  It's not perfect but it's effective and
  // doesn't require MPI communication.  The seed also includes bits
  // from the standard library's rand().  
  uint64_t myRank =static_cast<uint64_t>(mpi_rank);
  uint64_t seed64 = static_cast<uint64_t> (std::rand ()) + myRank + 17311uLL;
  unsigned int seed = static_cast<unsigned int> (seed64&0xffffffff);
  return seed;
}

#ifdef KOKKOS_ENABLE_CUDA
Kokkos::Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space> * cuda_pool_=nullptr;

void finalize_cuda_pool() {
  if(cuda_pool_ != nullptr) {
    delete cuda_pool_;
    cuda_pool_ = nullptr;
  }    
}
#endif // KOKKOS_ENABLE_CUDA


#ifdef KOKKOS_ENABLE_HIP
Kokkos::Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space> * hip_pool_=nullptr;

void finalize_hip_pool() {
  if(hip_pool_ != nullptr) {
    delete hip_pool_;
    hip_pool_ = nullptr;
  }    
}
#endif // KOKKOS_ENABLE_HIP

#ifdef KOKKOS_ENABLE_SYCL
Kokkos::Random_XorShift64_Pool<typename Kokkos::Experimental::SYCLDeviceUSMSpace::execution_space> * sycl_pool_=nullptr;

void finalize_sycl_pool() {
  if(sycl_pool_ != nullptr) {
    delete sycl_pool_;
    sycl_pool_ = nullptr;
  }    
}
#endif // KOKKOS_ENABLE_SYCL


#ifdef KOKKOS_ENABLE_OPENMP
Kokkos::Random_XorShift64_Pool<Kokkos::OpenMP> * openmp_pool_=nullptr;

void finalize_openmp_pool() {
  if(openmp_pool_ != nullptr) {
    delete openmp_pool_;
    openmp_pool_ = nullptr;
  }   
}
#endif // KOKKOS_ENABLE_OPENMP


#ifdef KOKKOS_ENABLE_SERIAL
Kokkos::Random_XorShift64_Pool<Kokkos::Serial> * serial_pool_=nullptr;

void finalize_serial_pool() {
  if(serial_pool_ != nullptr) {
    delete serial_pool_;
    serial_pool_ = nullptr;
  }   
}
#endif // KOKKOS_ENABLE_SERIAL

} // namespace (anonymous)


/********************************************************************************/
#ifdef KOKKOS_ENABLE_CUDA
void 
Static_Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space>::
resetPool(int mpi_rank) {
  using pool_type = Kokkos::Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space>;

  if(isSet())
    delete cuda_pool_;
  else
    Kokkos::push_finalize_hook(finalize_cuda_pool);

  cuda_pool_ = new pool_type(getSeedFromRank(mpi_rank));
} 

bool 
Static_Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space>::
isSet() {
  return cuda_pool_!=nullptr;
}

Kokkos::Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space> &
Static_Random_XorShift64_Pool<typename Kokkos::CudaSpace::execution_space>::
getPool() {
  TEUCHOS_TEST_FOR_EXCEPTION(!isSet(),std::runtime_error,"Tpetra::Details::Static_Random_XorShift64_Pool: resetPool() must be called before getPool");
  return *cuda_pool_;
}
#endif // KOKKOS_ENABLE_CUDA


/********************************************************************************/
#ifdef KOKKOS_ENABLE_HIP
void 
Static_Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space>::
resetPool(int mpi_rank) {
  using pool_type = Kokkos::Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space>;

  if(isSet())
    delete hip_pool_;
  else
    Kokkos::push_finalize_hook(finalize_hip_pool);

  hip_pool_ = new pool_type(getSeedFromRank(mpi_rank));
} 

bool 
Static_Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space>::
isSet() {
  return hip_pool_!=nullptr;
}

Kokkos::Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space> &
Static_Random_XorShift64_Pool<typename Kokkos::HIPSpace::execution_space>::
getPool() {
  TEUCHOS_TEST_FOR_EXCEPTION(!isSet(),std::runtime_error,"Tpetra::Details::Static_Random_XorShift64_Pool: resetPool() must be called before getPool");
  return *hip_pool_;
}
#endif // KOKKOS_ENABLE_HIP


/********************************************************************************/
#ifdef KOKKOS_ENABLE_SYCL
void 
Static_Random_XorShift64_Pool<typename Kokkos::Experimental::SYCLDeviceUSMSpace::execution_space>::
resetPool(int mpi_rank) {
  using pool_type = Kokkos::Random_XorShift64_Pool<typename Kokkos::Experimental::SYCLDeviceUSMSpace::execution_space>;

  if(isSet())
    delete sycl_pool_;
  else
    Kokkos::push_finalize_hook(finalize_sycl_pool);

  sycl_pool_ = new pool_type(getSeedFromRank(mpi_rank));
} 

bool 
Static_Random_XorShift64_Pool<typename Kokkos::Experimental::SYCLDeviceUSMSpace::execution_space>::
isSet() {
  return sycl_pool_!=nullptr;
}

Kokkos::Random_XorShift64_Pool<typename Kokkos::Experimental::SYCLDeviceUSMSpace::execution_space> &
Static_Random_XorShift64_Pool<typename Kokkos::Experimental::SYCLDeviceUSMSpace::execution_space>::
getPool() {
  TEUCHOS_TEST_FOR_EXCEPTION(!isSet(),std::runtime_error,"Tpetra::Details::Static_Random_XorShift64_Pool: resetPool() must be called before getPool");
  return *sycl_pool_;
}
#endif // KOKKOS_ENABLE_SYCL


/********************************************************************************/
#ifdef KOKKOS_ENABLE_OPENMP
void 
Static_Random_XorShift64_Pool<Kokkos::OpenMP>::
resetPool(int mpi_rank) {
  using pool_type = Kokkos::Random_XorShift64_Pool<Kokkos::OpenMP>;

  if(isSet())
    delete openmp_pool_;
  else
    Kokkos::push_finalize_hook(finalize_openmp_pool);

  openmp_pool_ = new pool_type(getSeedFromRank(mpi_rank));
} 

bool 
Static_Random_XorShift64_Pool<Kokkos::OpenMP>::
isSet() {
  return openmp_pool_!=nullptr;
}

Kokkos::Random_XorShift64_Pool<Kokkos::OpenMP> &
Static_Random_XorShift64_Pool<Kokkos::OpenMP>::
getPool() {
  TEUCHOS_TEST_FOR_EXCEPTION(!isSet(),std::runtime_error,"Tpetra::Details::Static_Random_XorShift64_Pool: resetPool() must be called before getPool");
  return *openmp_pool_;
}
#endif // KOKKOS_ENABLE_OPENMP


/********************************************************************************/
#ifdef KOKKOS_ENABLE_SERIAL
void 
Static_Random_XorShift64_Pool<Kokkos::Serial>::
resetPool(int mpi_rank) {
  using pool_type = Kokkos::Random_XorShift64_Pool<Kokkos::Serial>;

  if(isSet())
    delete serial_pool_;
  else
    Kokkos::push_finalize_hook(finalize_serial_pool);

  serial_pool_ = new pool_type(getSeedFromRank(mpi_rank));
} 

bool 
Static_Random_XorShift64_Pool<Kokkos::Serial>::
isSet() {
  return serial_pool_!=nullptr;
}

Kokkos::Random_XorShift64_Pool<Kokkos::Serial> &
Static_Random_XorShift64_Pool<Kokkos::Serial>::
getPool() {
  TEUCHOS_TEST_FOR_EXCEPTION(!isSet(),std::runtime_error,"Tpetra::Details::Static_Random_XorShift64_Pool: resetPool() must be called before getPool");
  return *serial_pool_;
}
#endif // KOKKOS_ENABLE_SERIAL


} // namespace Details
} // namespace Tpetra
