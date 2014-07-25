/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/



#include <Kokkos_Core.hpp>
#include <cstdio>
#include <cstdlib>

#ifndef KOKKOS_RANDOM_HPP
#define KOKKOS_RANDOM_HPP

// These generators are based on Vigna, Sebastiano (2014). "An experimental exploration of Marsaglia's xorshift generators, scrambled"
// See: http://arxiv.org/abs/1402.6246

namespace Kokkos {
  template<class DeviceType>
  class Random_XorShift64_Pool;

  template<class DeviceType>
  class Random_XorShift64 {
  private:
    uint64_t state_;
    const int chunk_num_;
    friend class Random_XorShift64_Pool<DeviceType>;
  public:

    enum {MAX_URAND = (1<<31)-1};
    enum {MAX_URAND64 = (1<<63)-1};

    KOKKOS_INLINE_FUNCTION
    Random_XorShift64 (uint64_t state, int chunk_num)
     : state_(state),chunk_num_(chunk_num){}

    KOKKOS_INLINE_FUNCTION
    uint32_t urand() {
      state_ ^= state_ >> 12;
      state_ ^= state_ << 25;
      state_ ^= state_ >> 27;

      uint64_t tmp = state_ * 2685821657736338717LL;
      return static_cast<uint32_t>(tmp>>32);
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64() {
      state_ ^= state_ >> 12;
      state_ ^= state_ << 25;
      state_ ^= state_ >> 27;
      return (state_ * 2685821657736338717LL) - 1;
    }
  };

  template<class DeviceType = Kokkos::Impl::DefaultDeviceType>
  class Random_XorShift64_Pool {
  private:
    typedef View<int*,DeviceType> lock_type;
    typedef View<int*,DeviceType> pool_chunk_count_type;
    typedef View<uint64_t*,DeviceType> state_data_type;
    lock_type locks_;
    pool_chunk_count_type chunk_count_;
    state_data_type state_;
    int num_states_;
    int chunk_size_;
    int num_chunks_;

  public:
    typedef Random_XorShift64<DeviceType> generator_type;

    Random_XorShift64_Pool(unsigned int seed) {
      num_states_ = 0;
      chunk_size_ = 0;
      num_chunks_ = 0;
      init(seed,DeviceType::max_hardware_threads(),1);
    }

    void init(unsigned int seed, int num_states, int chunk_size) {
      num_states_ = num_states;
      chunk_size_ = chunk_size;
      num_chunks_ = num_states/chunk_size;

      locks_ = lock_type("Kokkos::Random_XorShift64::locks",num_chunks_);
      chunk_count_ = pool_chunk_count_type("Kokkos::Random_XorShift64::chunk_count",num_chunks_);
      state_ = state_data_type("Kokkos::Random_XorShift64::state",num_states_);

      typename state_data_type::HostMirror h_state = create_mirror_view(state_);
      typename lock_type::HostMirror h_lock = create_mirror_view(locks_);
      typename pool_chunk_count_type::HostMirror h_chunk_count = create_mirror_view(chunk_count_);
      srand(seed);
      for(int i = 0; i < num_states_; i++) {
        int n1 = rand();
        int n2 = rand();
        int n3 = rand();
        int n4 = rand();
        h_state(i) = (((static_cast<uint64_t>(n1)) & 0xffff)<<00) |
                     (((static_cast<uint64_t>(n2)) & 0xffff)<<16) |
                     (((static_cast<uint64_t>(n3)) & 0xffff)<<32) |
                     (((static_cast<uint64_t>(n4)) & 0xffff)<<48);
      }
      deep_copy(state_,h_state);

      for(int i = 0; i < num_chunks_; i++) {
        h_lock(i) = 0;
        h_chunk_count(i) = 0;
      }
      deep_copy(locks_,h_lock);
      deep_copy(chunk_count_,h_chunk_count);
    }

    KOKKOS_INLINE_FUNCTION
    Random_XorShift64<DeviceType> get_state() const {
      const int i = DeviceType::hardware_thread_id();;
      return Random_XorShift64<DeviceType>(state_(i),i);
    }

    KOKKOS_INLINE_FUNCTION
    void free_state(const Random_XorShift64<DeviceType>& state) const {
      state_(state.chunk_num_) = state.state_;
    }
  };


  template<class DeviceType>
  class Random_XorShift1024_Pool;

  template<class DeviceType>
  class Random_XorShift1024 {
  private:
    int p_;
    const int chunk_num_;
    uint64_t state_[16];
    friend class Random_XorShift1024_Pool<DeviceType>;
  public:

    enum {MAX_URAND = (1<<31)-1};
    enum {MAX_URAND64 = (1<<63)-1};

    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024 (uint64_t* state, int p, int chunk_num):
      p_(p),chunk_num_(chunk_num){
      for(int i=0 ; i<16; i++)
        state_[i] = state[i];
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand() {
      uint64_t state_0 = state_[ p_ ];
      uint64_t state_1 = state_[ p_ = ( p_ + 1 ) & 15 ];
      state_1 ^= state_1 << 31;
      state_1 ^= state_1 >> 11;
      state_0 ^= state_0 >> 30;
      const uint64_t tmp = ( state_[ p_ ] = state_0 ^ state_1 ) * 1181783497276652981LL;
      return static_cast<uint32_t>(tmp>>32) - 1;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64() {
      uint64_t state_0 = state_[ p_ ];
      uint64_t state_1 = state_[ p_ = ( p_ + 1 ) & 15 ];
      state_1 ^= state_1 << 31;
      state_1 ^= state_1 >> 11;
      state_0 ^= state_0 >> 30;
      return (( state_[ p_ ] = state_0 ^ state_1 ) * 1181783497276652981LL) - 1;
    }
  };




  template<class DeviceType = Kokkos::Impl::DefaultDeviceType>
  class Random_XorShift1024_Pool {
  private:
    typedef View<int*,DeviceType> int_view_type;
    typedef View<uint64_t*[16],DeviceType> state_data_type;

    int_view_type locks_;
    int_view_type chunk_count_;
    state_data_type state_;
    int_view_type p_;
    int num_states_;
    int chunk_size_;
    int num_chunks_;

  public:
    typedef Random_XorShift1024<DeviceType> generator_type;

    Random_XorShift1024_Pool(unsigned int seed){
      num_states_ = 0;
      chunk_size_ = 0;
      num_chunks_ = 0;
      init(seed,DeviceType::max_hardware_threads(),1);
    }

    void init(unsigned int seed, int num_states, int chunk_size) {
      num_states_ = num_states;
      chunk_size_ = chunk_size;
      num_chunks_ = num_states/chunk_size;

      locks_ = int_view_type("Kokkos::Random_XorShift1024::locks",num_chunks_);
      chunk_count_ = int_view_type("Kokkos::Random_XorShift1024::chunk_count",num_chunks_);
      state_ = state_data_type("Kokkos::Random_XorShift1024::state",num_states_);
      p_ = int_view_type("Kokkos::Random_XorShift1024::p",num_states_);

      typename state_data_type::HostMirror h_state = create_mirror_view(state_);
      typename int_view_type::HostMirror h_lock = create_mirror_view(locks_);
      typename int_view_type::HostMirror h_chunk_count = create_mirror_view(chunk_count_);
      typename int_view_type::HostMirror h_p = create_mirror_view(p_);
      srand(seed);
      for(int i = 0; i < num_states_; i++) {
        for(int j = 0; j < 16 ; j++) {
          int n1 = rand();
          int n2 = rand();
          int n3 = rand();
          int n4 = rand();
          h_state(i,j) = (((static_cast<uint64_t>(n1)) & 0xffff)<<00) |
                         (((static_cast<uint64_t>(n2)) & 0xffff)<<16) |
                         (((static_cast<uint64_t>(n3)) & 0xffff)<<32) |
                         (((static_cast<uint64_t>(n4)) & 0xffff)<<48);
        }
        h_p(i) = 0;
      }
      deep_copy(state_,h_state);

      for(int i = 0; i < num_chunks_; i++) {
        h_lock(i) = 0;
        h_chunk_count(i) = 0;
      }
      deep_copy(locks_,h_lock);
      deep_copy(chunk_count_,h_chunk_count);
    }

    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024<DeviceType> get_state() const {
      const int i = DeviceType::hardware_thread_id();
      return Random_XorShift1024<DeviceType>(&state_(i,0),p_(i),i);
    };

    KOKKOS_INLINE_FUNCTION
    void free_state(const Random_XorShift1024<DeviceType>& state) const {
      for(int i = 0; i<16; i++)
        state_(state.chunk_num_,i) = state.state_[i];
      p_(state.chunk_num_) = state.p_;
    }
  };

#ifdef KOKKOS_HAVE_CUDA
  template<>
  class Random_XorShift1024<Kokkos::Cuda> {
  private:
    uint64_t* state_;
    int p_;
    const int stride_ ;
    const int chunk_num_;
    friend class Random_XorShift1024_Pool<Kokkos::Cuda>;
  public:
    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024 (uint64_t* state, int p, const int stride, const int chunk_num);

    KOKKOS_INLINE_FUNCTION
    uint32_t urand();

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64();
  };


template<>
Random_XorShift64_Pool<Kokkos::Cuda>::Random_XorShift64_Pool(unsigned int seed) {
  num_states_ = 0;
  chunk_size_ = 0;
  num_chunks_ = 0;
  init(seed,4*32768,1024);
}

template<>
KOKKOS_INLINE_FUNCTION
Random_XorShift64<Kokkos::Cuda> Random_XorShift64_Pool<Kokkos::Cuda>::get_state() const {
#ifdef __CUDA_ARCH__
  const int chunk_offset = (threadIdx.x*blockDim.y + threadIdx.y)*blockDim.z+threadIdx.z;
  __shared__ int block_i;
  if(chunk_offset == 0) {
    int i = ((blockIdx.x*gridDim.y+blockIdx.y)*blockDim.z+blockIdx.z)%num_chunks_;

    while(Kokkos::atomic_compare_exchange(&locks_(i),0,1)) {
      i++;
      if(i==num_chunks_) i = 0;
    }

    chunk_count_(i) += blockDim.x*blockDim.y*blockDim.z;
    block_i = i;
  }
  __syncthreads();

  return Random_XorShift64<Kokkos::Cuda>(state_(block_i*chunk_size_ + chunk_offset),block_i);
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
void Random_XorShift64_Pool<Kokkos::Cuda>::free_state(const Random_XorShift64<Kokkos::Cuda> &state) const {
#ifdef __CUDA_ARCH__
  const int chunk_offset = (threadIdx.x*blockDim.y + threadIdx.y)*blockDim.z+threadIdx.z;
  state_(state.chunk_num_*chunk_size_ + chunk_offset) = state.state_;
  int count = Kokkos::atomic_fetch_add(&chunk_count_(state.chunk_num_),-1);
  if(count == 1)
    locks_(state.chunk_num_) = 0;
  return;
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
Random_XorShift1024<Kokkos::Cuda>::Random_XorShift1024(uint64_t* state, int p, const int stride, const int chunk_num):
  state_ (state), p_(p), stride_ (stride), chunk_num_ (chunk_num) {
}

template<>
KOKKOS_INLINE_FUNCTION
uint32_t Random_XorShift1024<Kokkos::Cuda>::urand() {
  uint64_t state_0 = state_[ p_*stride_ ];
  uint64_t state_1 = state_[ (p_ = ( p_ + 1 ) & 15)*stride_ ];
  state_1 ^= state_1 << 31;
  state_1 ^= state_1 >> 11;
  state_0 ^= state_0 >> 30;
  const uint64_t tmp = ( state_[ p_*stride_ ] = state_0 ^ state_1 ) * 1181783497276652981LL;
  return static_cast<uint32_t>(tmp>>32) - 1;
}

template<>
KOKKOS_INLINE_FUNCTION
uint64_t Random_XorShift1024<Kokkos::Cuda>::urand64() {
  uint64_t state_0 = state_[ p_*stride_ ];
  uint64_t state_1 = state_[ (p_ = ( p_ + 1 ) & 15)*stride_ ];
  state_1 ^= state_1 << 31;
  state_1 ^= state_1 >> 11;
  state_0 ^= state_0 >> 30;
  return ( state_[ p_*stride_ ] = state_0 ^ state_1 ) * 1181783497276652981LL;
}


template<>
Random_XorShift1024_Pool<Kokkos::Cuda>::Random_XorShift1024_Pool(unsigned int seed) {
  num_states_ = 0;
  chunk_size_ = 0;
  num_chunks_ = 0;
  init(seed,4*32768,1024);
}

template<>
KOKKOS_INLINE_FUNCTION
Random_XorShift1024<Kokkos::Cuda> Random_XorShift1024_Pool<Kokkos::Cuda>::get_state() const {
#ifdef __CUDA_ARCH__
  const int chunk_offset = (threadIdx.x*blockDim.y + threadIdx.y)*blockDim.z+threadIdx.z;
  __shared__ int block_i;
  if(chunk_offset == 0) {
    int i = ((blockIdx.x*gridDim.y+blockIdx.y)*blockDim.z+blockIdx.z)%num_chunks_;

    while(Kokkos::atomic_compare_exchange(&locks_(i),0,1)) {
      i++;
      if(i==num_chunks_) i = 0;
    }

    chunk_count_(i) += blockDim.x*blockDim.y*blockDim.z;
    block_i = i;
  }
  __syncthreads();

  return Random_XorShift1024<Kokkos::Cuda>(&state_(block_i*chunk_size_ + chunk_offset,0),
                                            p_(block_i*chunk_size_ + chunk_offset),
                                            int(&state_(block_i*chunk_size_ + chunk_offset,1)-
                                                &state_(block_i*chunk_size_ + chunk_offset,0)),
                                            block_i);
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
void Random_XorShift1024_Pool<Kokkos::Cuda>::free_state(const Random_XorShift1024<Kokkos::Cuda> &state) const {
#ifdef __CUDA_ARCH__
  const int chunk_offset = (threadIdx.x*blockDim.y + threadIdx.y)*blockDim.z+threadIdx.z;
  p_(state.chunk_num_*chunk_size_ + chunk_offset) = state.p_;
  int count = Kokkos::atomic_fetch_add(&chunk_count_(state.chunk_num_),-1);
  if(count == 1)
    locks_(state.chunk_num_) = 0;
  return;
#endif
}

#endif
}


#endif
