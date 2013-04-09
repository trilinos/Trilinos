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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_ATOMIC_FETCH_ADD_H_
#define KOKKOSARRAY_ATOMIC_FETCH_ADD_H_


template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,int>::value, int  >::type
atomic_fetch_add(T* dest, T val) {
  #if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL)
	return __sync_fetch_and_add(dest,val);
  #endif

  #ifdef OMP31_ATOMICS
	int retval;
    #pragma omp atomic capture
	{
	  retval = dest[0];
	  dest[0] += val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    return atomicAdd(dest,val);
  #endif
}


template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,long long int>::value, long long int  >::type
atomic_fetch_add(T* dest, T val) {
  #if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL)
	return __sync_fetch_and_add(dest,val);
  #endif

  #ifdef OMP31_ATOMICS
	long long int retval;
    #pragma omp atomic capture
	{
	  retval = dest[0];
	  dest[0] += val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    unsigned long long int return_value =
    atomicAdd(reinterpret_cast<unsigned long long int*>(dest),*reinterpret_cast<unsigned long long int*>(&val));
    return *reinterpret_cast<unsigned long long int*>(&return_value);
  #endif
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,long int>::value && (sizeof(T) == sizeof(int)), long int >::type
atomic_fetch_add(T* dest, T val) {
  return (long int) atomic_fetch_add((int*) dest,*reinterpret_cast<int*>(&val));
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,long int>::value && (sizeof(T) == sizeof(long long int)), long int >::type
atomic_fetch_add(T* dest, T val) {
  return (long int) atomic_fetch_add((long long int*) dest,*reinterpret_cast<long long int*>(&val));
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,unsigned int>::value, unsigned int  >::type
atomic_fetch_add(T* dest, T val) {
  #ifdef GNU_ATOMICS_GCC
	return __sync_fetch_and_add(dest,val);
  #endif
  #ifdef GNU_ATOMICS_INTEL
	  volatile int* address_as_int = reinterpret_cast<volatile int*>(dest);
	  int old = *address_as_int;
	  int assumed;
	  unsigned int newval;
	  do {
	    assumed = old;
	    newval = val + *(reinterpret_cast<unsigned int*>(&assumed));
	    old = __sync_val_compare_and_swap(address_as_int, assumed, *reinterpret_cast<int*>(&newval));
	  } while (assumed != old);
	  return *reinterpret_cast<unsigned int*>(&old);
  #endif

  #ifdef OMP31_ATOMICS
	unsigned int retval;
    #pragma omp atomic capture
	{
	  retval = dest[0];
	  dest[0] += val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    return atomicAdd(dest,val);
  #endif
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T, unsigned long long int >::value, unsigned long long int  >::type
atomic_fetch_add(T* dest, T val) {
  #ifdef GNU_ATOMICS_GCC
	return __sync_fetch_and_add(dest,val);
  #endif
  #ifdef GNU_ATOMICS_INTEL
	  volatile long long int* address_as_int = reinterpret_cast<volatile long long int*>(dest);
	  long long int old = *address_as_int;
	  long long int assumed;
	  unsigned long long int newval;
	  do {
	    assumed = old;
	    newval = val + *(reinterpret_cast<unsigned long long int*>(&assumed));
	    old = __sync_val_compare_and_swap(address_as_int, assumed, *reinterpret_cast<long long int*>(&newval));
	  } while (assumed != old);
	  return *reinterpret_cast<unsigned long long int*>(&old);
  #endif

  #ifdef OMP31_ATOMICS
	unsigned long long int retval;
    #pragma omp atomic capture
	{
	  retval = dest[0];
	  dest[0] += val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    return atomicAdd(dest,val);
  #endif
}


template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,unsigned long int>::value && (sizeof(T) == sizeof(int)), unsigned long int >::type
atomic_fetch_add(T* dest, T val) {
  return (unsigned long int) atomic_fetch_add((unsigned int*) dest,*reinterpret_cast<unsigned int*>(&val));
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,unsigned long int>::value && (sizeof(T) == sizeof(long long int)), unsigned long int >::type
atomic_fetch_add(volatile T* dest, T val) {
  return (unsigned long int) atomic_fetch_add((unsigned long long int*) dest,*reinterpret_cast<unsigned long long int*>(&val));
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T, float >::value, float  >::type
atomic_fetch_add(T* dest, T val) {
  #if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL)
	  int32_t* address_as_int = reinterpret_cast<int32_t*>(dest);
	  int32_t old = *address_as_int;
	  int32_t assumed;
	  float newval;
	  do {
	    assumed = old;
	    newval = val + *(reinterpret_cast<float*>(&assumed));
	    old = __sync_val_compare_and_swap(address_as_int, assumed, *reinterpret_cast<int32_t*>(&newval));
	  } while (assumed != old);
	  return *reinterpret_cast<float*>(&old);
  #endif

  #ifdef OMP31_ATOMICS
	float retval;
    #pragma omp atomic capture
	{
	  retval = dest[0];
	  dest[0] += val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    return atomicAdd(dest,val);
  #endif
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T, double >::value, double  >::type
atomic_fetch_add(T* dest, T val) {
  #if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL)
	  int64_t* address_as_int = reinterpret_cast<int64_t*>(dest);
	  int64_t old = *address_as_int;
	  int64_t assumed;
	  double newval;
	  do {
	    assumed = old;
	    newval = val + *(reinterpret_cast<double*>(&assumed));
	    old = __sync_val_compare_and_swap(address_as_int, assumed, *reinterpret_cast<int64_t*>(&newval));
	  } while (assumed != old);
	  return *reinterpret_cast<double*>(&old);
  #endif

  #ifdef OMP31_ATOMICS
	double retval;
    #pragma omp atomic capture
	{
	  retval = dest[0];
	  dest[0] += val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
      unsigned long long int* address_as_int = (unsigned long long int*) dest;
      unsigned long long int old = *address_as_int;
	  unsigned long long int assumed;
	  double newval;
	  do {
	    assumed = old;
	    newval = val + *(reinterpret_cast<double*>(&assumed));
	    old = atomicCAS(address_as_int, assumed, *reinterpret_cast<unsigned long long int*>(&newval));
	  } while (assumed != old);
	  return *reinterpret_cast<double*>(&old);
  #endif
}

#endif /* KOKKOSARRAY_ATOMIC_FETCH_ADD_H_ */
