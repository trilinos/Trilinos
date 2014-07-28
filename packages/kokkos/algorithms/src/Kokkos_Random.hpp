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
#ifdef KOKKOS_HAVE_CUDA
#include <Kokkos_Cuda.hpp>
#endif
#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifndef KOKKOS_RANDOM_HPP
#define KOKKOS_RANDOM_HPP

// These generators are based on Vigna, Sebastiano (2014). "An experimental exploration of Marsaglia's xorshift generators, scrambled"
// See: http://arxiv.org/abs/1402.6246

namespace Kokkos {

  //Template functions to get equidistributed random numbers from a generator

  //Maximum value for the rand() call depending on Scalar type

  template<class Generator, class Scalar>
  struct rand;


  template<class Generator>
  struct rand<Generator,int> {
    KOKKOS_INLINE_FUNCTION
    static int max(){return Generator::MAX_RAND;}
    KOKKOS_INLINE_FUNCTION
    static int draw(Generator& gen)
                          {return gen.rand();}
    KOKKOS_INLINE_FUNCTION
    static int draw(Generator& gen, const int& range)
                          {return gen.rand(range);}
    KOKKOS_INLINE_FUNCTION
    static int draw(Generator& gen, const int& start, const int& end)
                          {return gen.rand(start,end);}

  };

  template<class Generator>
  struct rand<Generator,unsigned int> {
    KOKKOS_INLINE_FUNCTION
    static unsigned int max(){return Generator::MAX_URAND;}
    KOKKOS_INLINE_FUNCTION
    static unsigned int draw(Generator& gen)
                          {return gen.urand();}
    KOKKOS_INLINE_FUNCTION
    static unsigned int draw(Generator& gen, const unsigned int& range)
                          {return gen.urand(range);}
    KOKKOS_INLINE_FUNCTION
    static unsigned int draw(Generator& gen, const unsigned int& start, const unsigned int& end)
                          {return gen.urand(start,end);}

  };

  template<class Generator>
  struct rand<Generator,int64_t> {
    KOKKOS_INLINE_FUNCTION
    static int64_t max(){return Generator::MAX_RAND64;}
    KOKKOS_INLINE_FUNCTION
    static int64_t draw(Generator& gen)
                          {return gen.rand64();}
    KOKKOS_INLINE_FUNCTION
    static int64_t draw(Generator& gen, const int64_t& range)
                          {return gen.rand64(range);}
    KOKKOS_INLINE_FUNCTION
    static int64_t draw(Generator& gen, const int64_t& start, const int64_t& end)
                          {return gen.rand64(start,end);}

  };

  template<class Generator>
  struct rand<Generator,uint64_t> {
    KOKKOS_INLINE_FUNCTION
    static uint64_t max(){return Generator::MAX_URAND64;}
    KOKKOS_INLINE_FUNCTION
    static uint64_t draw(Generator& gen)
                          {return gen.urand64();}
    KOKKOS_INLINE_FUNCTION
    static uint64_t draw(Generator& gen, const uint64_t& range)
                          {return gen.urand64(range);}
    KOKKOS_INLINE_FUNCTION
    static uint64_t draw(Generator& gen, const uint64_t& start, const uint64_t& end)
                          {return gen.urand64(start,end);}

  };

  template<class Generator>
  struct rand<Generator,float> {
    KOKKOS_INLINE_FUNCTION
    static float max(){return 1.0f;}
    KOKKOS_INLINE_FUNCTION
    static float draw(Generator& gen)
                          {return gen.frand();}
    KOKKOS_INLINE_FUNCTION
    static float draw(Generator& gen, const float& range)
                          {return gen.frand(range);}
    KOKKOS_INLINE_FUNCTION
    static float draw(Generator& gen, const float& start, const float& end)
                          {return gen.frand(start,end);}

  };

  template<class Generator>
  struct rand<Generator,double> {
    KOKKOS_INLINE_FUNCTION
    static double max(){return 1.0;}
    KOKKOS_INLINE_FUNCTION
    static double draw(Generator& gen)
                          {return gen.drand();}
    KOKKOS_INLINE_FUNCTION
    static double draw(Generator& gen, const double& range)
                          {return gen.drand(range);}
    KOKKOS_INLINE_FUNCTION
    static double draw(Generator& gen, const double& start, const double& end)
                          {return gen.drand(start,end);}

  };

  template<class DeviceType>
  class Random_XorShift64_Pool;

  template<class DeviceType>
  class Random_XorShift64 {
  private:
    uint64_t state_;
    const int state_idx_;
    friend class Random_XorShift64_Pool<DeviceType>;
  public:

    typedef DeviceType device_type;

    enum {MAX_URAND = 0xffffffffU};
    enum {MAX_URAND64 = 0xffffffffffffffffULL-1};
    enum {MAX_RAND = static_cast<int>(0xffffffff/2)};
    enum {MAX_RAND64 = static_cast<int64_t>(0xffffffffffffffffLL/2-1)};

    KOKKOS_INLINE_FUNCTION
    Random_XorShift64 (uint64_t state, int state_idx)
     : state_(state),state_idx_(state_idx){}

    KOKKOS_INLINE_FUNCTION
    uint32_t urand() {
      state_ ^= state_ >> 12;
      state_ ^= state_ << 25;
      state_ ^= state_ >> 27;

      uint64_t tmp = state_ * 2685821657736338717ULL;
      tmp = tmp>>16;
      return static_cast<uint32_t>(tmp&MAX_URAND);
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64() {
      state_ ^= state_ >> 12;
      state_ ^= state_ << 25;
      state_ ^= state_ >> 27;
      return (state_ * 2685821657736338717ULL) - 1;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& range) {
      const uint32_t max_val = (MAX_URAND/range)*range;
      uint32_t tmp = urand();
      while(tmp>=max_val)
        urand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& start, const uint32_t& end ) {
      return urand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& range) {
      const uint64_t max_val = (MAX_URAND64/range)*range;
      uint64_t tmp = urand64();
      while(tmp>=max_val)
        urand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& start, const uint64_t& end ) {
      return urand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int rand() {
      return static_cast<int>(urand()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& range) {
      const int max_val = (MAX_RAND/range)*range;
      int tmp = rand();
      while(tmp>=max_val)
        rand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& start, const int& end ) {
      return rand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64() {
      return static_cast<int64_t>(urand64()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& range) {
      const int64_t max_val = (MAX_RAND64/range)*range;
      int64_t tmp = rand64();
      while(tmp>=max_val)
        rand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& start, const int64_t& end ) {
      return rand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    float frand() {
      return 1.0f * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& start, const float& end ) {
      return frand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    double drand() {
      return 1.0 * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& start, const double& end ) {
      return frand(end-start)+start;
    }

    //Marsaglia polar method for drawing a standard normal distributed random number
    KOKKOS_INLINE_FUNCTION
    double normal() {
      double S = 2.0;
      double U;
      while(S>=1.0) {
        U = drand();
        const double V = drand();
        S = U*U+V*V;
      }
      return U*sqrt(-2.0*log(S)/S);
    }

    KOKKOS_INLINE_FUNCTION
    double normal(const double& mean, const double& std_dev=1.0) {
      return mean + normal()*std_dev;
    }

  };

  template<class DeviceType = Kokkos::Impl::DefaultDeviceType>
  class Random_XorShift64_Pool {
  private:
    typedef View<int*,DeviceType> lock_type;
    typedef View<uint64_t*,DeviceType> state_data_type;
    lock_type locks_;
    state_data_type state_;
    int num_states_;

  public:
    typedef Random_XorShift64<DeviceType> generator_type;
    typedef DeviceType device_type;

    Random_XorShift64_Pool(unsigned int seed) {
      num_states_ = 0;
      init(seed,DeviceType::max_hardware_threads());
    }

    void init(unsigned int seed, int num_states) {
      num_states_ = num_states;

      locks_ = lock_type("Kokkos::Random_XorShift64::locks",num_states_);
      state_ = state_data_type("Kokkos::Random_XorShift64::state",num_states_);

      typename state_data_type::HostMirror h_state = create_mirror_view(state_);
      typename lock_type::HostMirror h_lock = create_mirror_view(locks_);
      srand(seed);
      for(int i = 0; i < num_states_; i++) {
        int n1 = ::rand();
        int n2 = ::rand();
        int n3 = ::rand();
        int n4 = ::rand();
        h_state(i) = (((static_cast<uint64_t>(n1)) & 0xffff)<<00) |
                     (((static_cast<uint64_t>(n2)) & 0xffff)<<16) |
                     (((static_cast<uint64_t>(n3)) & 0xffff)<<32) |
                     (((static_cast<uint64_t>(n4)) & 0xffff)<<48);
        h_lock(i) = 0;
      }
      deep_copy(state_,h_state);
      deep_copy(locks_,h_lock);
    }

    KOKKOS_INLINE_FUNCTION
    Random_XorShift64<DeviceType> get_state() const {
      const int i = DeviceType::hardware_thread_id();;
      return Random_XorShift64<DeviceType>(state_(i),i);
    }

    KOKKOS_INLINE_FUNCTION
    void free_state(const Random_XorShift64<DeviceType>& state) const {
      state_(state.state_idx_) = state.state_;
    }
  };


  template<class DeviceType>
  class Random_XorShift1024_Pool;

  template<class DeviceType>
  class Random_XorShift1024 {
  private:
    int p_;
    const int state_idx_;
    uint64_t state_[16];
    friend class Random_XorShift1024_Pool<DeviceType>;
  public:

    typedef DeviceType device_type;

    enum {MAX_URAND = 0xffffffffU};
    enum {MAX_URAND64 = 0xffffffffffffffffULL-1};
    enum {MAX_RAND = static_cast<int>(0xffffffffU/2)};
    enum {MAX_RAND64 = static_cast<int64_t>(0xffffffffffffffffULL/2-1)};

    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024 (uint64_t* state, int p, int state_idx):
      p_(p),state_idx_(state_idx){
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
      uint64_t tmp = ( state_[ p_ ] = state_0 ^ state_1 ) * 1181783497276652981ULL;
      tmp = tmp>>16;
      return static_cast<uint32_t>(tmp&MAX_URAND);
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

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& range) {
      const uint32_t max_val = (MAX_URAND/range)*range;
      uint32_t tmp = urand();
      while(tmp>=max_val)
        urand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& start, const uint32_t& end ) {
      return urand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& range) {
      const uint64_t max_val = (MAX_URAND64/range)*range;
      uint64_t tmp = urand64();
      while(tmp>=max_val)
        urand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& start, const uint64_t& end ) {
      return urand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int rand() {
      return static_cast<int>(urand()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& range) {
      const int max_val = (MAX_RAND/range)*range;
      int tmp = rand();
      while(tmp>=max_val)
        rand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& start, const int& end ) {
      return rand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64() {
      return static_cast<int64_t>(urand64()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& range) {
      const int64_t max_val = (MAX_RAND64/range)*range;
      int64_t tmp = rand64();
      while(tmp>=max_val)
        rand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& start, const int64_t& end ) {
      return rand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    float frand() {
      return 1.0f * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& start, const float& end ) {
      return frand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    double drand() {
      return 1.0 * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& start, const double& end ) {
      return frand(end-start)+start;
    }

    //Marsaglia polar method for drawing a standard normal distributed random number
    KOKKOS_INLINE_FUNCTION
    double normal() {
      double S = 2.0;
      double U;
      while(S>=1.0) {
        U = drand();
        const double V = drand();
        S = U*U+V*V;
      }
      return U*sqrt(-2.0*log(S)/S);
    }

    KOKKOS_INLINE_FUNCTION
    double normal(const double& mean, const double& std_dev=1.0) {
      return mean + normal()*std_dev;
    }
  };

#ifdef KOKKOS_HAVE_CUDA
  template<>
  class Random_XorShift1024<Kokkos::Cuda> {
  private:
    int p_;
    const int state_idx_;
    uint64_t* state_;
    friend class Random_XorShift1024_Pool<Kokkos::Cuda>;
  public:

    typedef Kokkos::Cuda device_type;

    enum {MAX_URAND = 0xffffffffU};
    enum {MAX_URAND64 = 0xffffffffffffffffULL-1};
    enum {MAX_RAND = static_cast<int>(0xffffffffU/2)};
    enum {MAX_RAND64 = static_cast<int64_t>(0xffffffffffffffffULL/2-1)};

    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024 (uint64_t* state, int p, int state_idx):
      p_(p),state_idx_(state_idx),state_(state){
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand() {
      uint64_t state_0 = state_[ p_ ];
      uint64_t state_1 = state_[ p_ = ( p_ + 1 ) & 15 ];
      state_1 ^= state_1 << 31;
      state_1 ^= state_1 >> 11;
      state_0 ^= state_0 >> 30;
      uint64_t tmp = ( state_[ p_ ] = state_0 ^ state_1 ) * 1181783497276652981ULL;
      tmp = tmp>>16;
      return static_cast<uint32_t>(tmp&MAX_URAND);
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

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& range) {
      const uint32_t max_val = (MAX_URAND/range)*range;
      uint32_t tmp = urand();
      while(tmp>=max_val)
        urand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint32_t urand(const uint32_t& start, const uint32_t& end ) {
      return urand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& range) {
      const uint64_t max_val = (MAX_URAND64/range)*range;
      uint64_t tmp = urand64();
      while(tmp>=max_val)
        urand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    uint64_t urand64(const uint64_t& start, const uint64_t& end ) {
      return urand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int rand() {
      return static_cast<int>(urand()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& range) {
      const int max_val = (MAX_RAND/range)*range;
      int tmp = rand();
      while(tmp>=max_val)
        rand();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int rand(const int& start, const int& end ) {
      return rand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64() {
      return static_cast<int64_t>(urand64()/2);
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& range) {
      const int64_t max_val = (MAX_RAND64/range)*range;
      int64_t tmp = rand64();
      while(tmp>=max_val)
        rand64();
      return tmp%range;
    }

    KOKKOS_INLINE_FUNCTION
    int64_t rand64(const int64_t& start, const int64_t& end ) {
      return rand64(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    float frand() {
      return 1.0f * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    float frand(const float& start, const float& end ) {
      return frand(end-start)+start;
    }

    KOKKOS_INLINE_FUNCTION
    double drand() {
      return 1.0 * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& range) {
      return range * urand64()/MAX_URAND64;
    }

    KOKKOS_INLINE_FUNCTION
    double drand(const double& start, const double& end ) {
      return frand(end-start)+start;
    }

    //Marsaglia polar method for drawing a standard normal distributed random number
    KOKKOS_INLINE_FUNCTION
    double normal() {
      double S = 2.0;
      double U;
      while(S>=1.0) {
        U = drand();
        const double V = drand();
        S = U*U+V*V;
      }
      return U*sqrt(-2.0*log(S)/S);
    }

    KOKKOS_INLINE_FUNCTION
    double normal(const double& mean, const double& std_dev=1.0) {
      return mean + normal()*std_dev;
    }
  };

#endif

  template<class DeviceType = Kokkos::Impl::DefaultDeviceType>
  class Random_XorShift1024_Pool {
  private:
    typedef View<int*,DeviceType> int_view_type;
    typedef View<uint64_t*[16],DeviceType> state_data_type;

    int_view_type locks_;
    state_data_type state_;
    int_view_type p_;
    int num_states_;

  public:
    typedef Random_XorShift1024<DeviceType> generator_type;

    typedef DeviceType device_type;

    Random_XorShift1024_Pool(unsigned int seed){
      num_states_ = 0;
      init(seed,DeviceType::max_hardware_threads());
    }

    void init(unsigned int seed, int num_states) {
      num_states_ = num_states;

      locks_ = int_view_type("Kokkos::Random_XorShift1024::locks",num_states_);
      state_ = state_data_type("Kokkos::Random_XorShift1024::state",num_states_);
      p_ = int_view_type("Kokkos::Random_XorShift1024::p",num_states_);

      typename state_data_type::HostMirror h_state = create_mirror_view(state_);
      typename int_view_type::HostMirror h_lock = create_mirror_view(locks_);
      typename int_view_type::HostMirror h_p = create_mirror_view(p_);
      srand(seed);
      for(int i = 0; i < num_states_; i++) {
        for(int j = 0; j < 16 ; j++) {
          int n1 = ::rand();
          int n2 = ::rand();
          int n3 = ::rand();
          int n4 = ::rand();
          h_state(i,j) = (((static_cast<uint64_t>(n1)) & 0xffff)<<00) |
                         (((static_cast<uint64_t>(n2)) & 0xffff)<<16) |
                         (((static_cast<uint64_t>(n3)) & 0xffff)<<32) |
                         (((static_cast<uint64_t>(n4)) & 0xffff)<<48);
        }
        h_p(i) = 0;
        h_lock(i) = 0;
      }
      deep_copy(state_,h_state);
      deep_copy(locks_,h_lock);
    }

    KOKKOS_INLINE_FUNCTION
    Random_XorShift1024<DeviceType> get_state() const {
      const int i = DeviceType::hardware_thread_id();
      return Random_XorShift1024<DeviceType>(&state_(i,0),p_(i),i);
    };

    KOKKOS_INLINE_FUNCTION
    void free_state(const Random_XorShift1024<DeviceType>& state) const {
      for(int i = 0; i<16; i++)
        state_(state.state_idx_,i) = state.state_[i];
      p_(state.state_idx_) = state.p_;
    }
  };

#ifdef KOKKOS_HAVE_CUDA
#ifdef __CUDACC__
template<>
Random_XorShift64_Pool<Kokkos::Cuda>::Random_XorShift64_Pool(unsigned int seed) {
  num_states_ = 0;
  init(seed,4*32768);
}

template<>
KOKKOS_INLINE_FUNCTION
Random_XorShift64<Kokkos::Cuda> Random_XorShift64_Pool<Kokkos::Cuda>::get_state() const {
#ifdef __CUDA_ARCH__
  const int i_offset = (threadIdx.x*blockDim.y + threadIdx.y)*blockDim.z+threadIdx.z;
  int i = ((blockIdx.x*gridDim.y+blockIdx.y)*gridDim.z + blockIdx.z) *
           blockDim.x*blockDim.y*blockDim.z + i_offset;
  while(Kokkos::atomic_compare_exchange(&locks_(i),0,1)) {
      i+=blockDim.x*blockDim.y*blockDim.z;
      if(i>=num_states_) {i = i_offset;}
  }

  return Random_XorShift64<Kokkos::Cuda>(state_(i),i);
#else
  return Random_XorShift64<Kokkos::Cuda>(state_(0),0);
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
void Random_XorShift64_Pool<Kokkos::Cuda>::free_state(const Random_XorShift64<Kokkos::Cuda> &state) const {
#ifdef __CUDA_ARCH__
  state_(state.state_idx_) = state.state_;
  locks_(state.state_idx_) = 0;
  return;
#endif
}


template<>
Random_XorShift1024_Pool<Kokkos::Cuda>::Random_XorShift1024_Pool(unsigned int seed) {
  num_states_ = 0;
  init(seed,4*32768);
}

template<>
KOKKOS_INLINE_FUNCTION
Random_XorShift1024<Kokkos::Cuda> Random_XorShift1024_Pool<Kokkos::Cuda>::get_state() const {
#ifdef __CUDA_ARCH__
  const int i_offset = (threadIdx.x*blockDim.y + threadIdx.y)*blockDim.z+threadIdx.z;
  int i = ((blockIdx.x*gridDim.y+blockIdx.y)*gridDim.z + blockIdx.z) *
           blockDim.x*blockDim.y*blockDim.z + i_offset;
  while(Kokkos::atomic_compare_exchange(&locks_(i),0,1)) {
      i+=blockDim.x*blockDim.y*blockDim.z;
      if(i>=num_states_) {i = i_offset;}
  }

  return Random_XorShift1024<Kokkos::Cuda>(&state_(i,0), p_(i), i);
#else
  return Random_XorShift1024<Kokkos::Cuda>(&state_(0,0), p_(0), 0);
#endif
}

template<>
KOKKOS_INLINE_FUNCTION
void Random_XorShift1024_Pool<Kokkos::Cuda>::free_state(const Random_XorShift1024<Kokkos::Cuda> &state) const {
#ifdef __CUDA_ARCH__
  for(int i=0; i<16; i++)
    state_(state.state_idx_,i) = state.state_[i];
  locks_(state.state_idx_) = 0;
  return;
#endif
}
#endif
#endif
}


#endif
