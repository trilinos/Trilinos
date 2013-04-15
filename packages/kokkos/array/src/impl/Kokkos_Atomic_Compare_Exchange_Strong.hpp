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


template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,int>::value, int  >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  #if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL)
	return __sync_val_compare_and_swap(dest,compare,val);
  #endif

  #ifdef OMP31_ATOMICS
	T retval;
    #pragma omp critical
	{
	  retval = dest[0];
	  if(retval == compare)
		dest[0] = val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    return atomicCAS((T*) dest,compare,val);
  #endif
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,long long int>::value, long long int  >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  #if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL)
	return __sync_val_compare_and_swap(dest,compare,val);
  #endif

  #ifdef OMP31_ATOMICS
	T retval;
    #pragma omp critical
	{
	  retval = dest[0];
	  if(retval == compare)
		dest[0] = val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    return atomicCAS((unsigned long long int*) reinterpret_cast<volatile unsigned long long int*>(dest),
    		*reinterpret_cast<unsigned long long int*>(&compare),*reinterpret_cast<unsigned long long int*>(&val));
  #endif
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,long int>::value && (sizeof(T) == sizeof(int)), long int >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  return (long int) atomic_compare_exchange((int*) dest,*reinterpret_cast<int*>(&compare),*reinterpret_cast<int*>(&val));
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,long int>::value && (sizeof(T) == sizeof(long long int)), long int >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  return (long int) atomic_compare_exchange((long long int*) dest,*reinterpret_cast<long long int*>(&compare),*reinterpret_cast<long long int*>(&val));
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,unsigned int>::value, unsigned int  >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  #ifdef GNU_ATOMICS_GCC
	return __sync_val_compare_and_swap(dest,compare,val);
  #endif

  #ifdef GNU_ATOMICS_INTEL
	volatile int* dest_int = reinterpret_cast<volatile int*>(dest);
	int compare_int = *reinterpret_cast<int*> (&compare);
	int val_int = *reinterpret_cast<int*> (&val);
	int return_value = __sync_val_compare_and_swap(dest_int,compare_int,val_int);
	return *reinterpret_cast<T*>(&return_value);
  #endif

  #ifdef OMP31_ATOMICS
	T retval;
    #pragma omp critical
	{
	  retval = dest[0];
	  if(retval == compare)
		dest[0] = val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    return atomicCAS((T*) dest,compare,val);
  #endif
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,unsigned long long int>::value, unsigned long long int  >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  #ifdef GNU_ATOMICS_GCC
	return __sync_val_compare_and_swap(dest,compare,val);
  #endif

  #ifdef GNU_ATOMICS_INTEL
	volatile long long int* dest_int = reinterpret_cast<volatile long long int*>(dest);
	long long int compare_int = *reinterpret_cast<long long int*> (&compare);
	long long int val_int = *reinterpret_cast<long long int*> (&val);
	long long int return_value = __sync_val_compare_and_swap(dest_int,compare_int,val_int);
	return *reinterpret_cast<T*>(&return_value);
  #endif

  #ifdef OMP31_ATOMICS
	T retval;
    #pragma omp critical
	{
	  retval = dest[0];
	  if(retval == compare)
		dest[0] = val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    return atomicCAS((T*) dest,compare,val);
  #endif
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,unsigned long int>::value && (sizeof(T) == sizeof(unsigned int)), unsigned long int >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  return (unsigned long int) atomic_compare_exchange((unsigned int*) dest,*reinterpret_cast<unsigned int*>(&compare),*reinterpret_cast<unsigned int*>(&val));
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,unsigned long int>::value && (sizeof(T) == sizeof(unsigned long long int)), unsigned long int >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  return (unsigned long int) atomic_compare_exchange((unsigned long long int*) dest,*reinterpret_cast<unsigned long long int*>(&compare),*reinterpret_cast<unsigned long long int*>(&val));
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,float>::value, float  >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  #if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL)
	volatile int32_t* dest_int = reinterpret_cast<volatile int32_t*>(dest);
	int32_t compare_int = *reinterpret_cast<int32_t*> (&compare);
	int32_t val_int = *reinterpret_cast<int32_t*> (&val);
	int32_t return_value = __sync_val_compare_and_swap(dest_int,compare_int,val_int);
	return *reinterpret_cast<float*> (&return_value);
  #endif

  #ifdef OMP31_ATOMICS
	float retval;
    #pragma omp critical
	{
	  retval = dest[0];
	  if(retval == compare)
		dest[0] = val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    volatile unsigned int* dest_int = reinterpret_cast<volatile unsigned int*>(dest);
    unsigned int compare_int = *reinterpret_cast<unsigned int*> (&compare);
    unsigned int val_int = *reinterpret_cast<unsigned int*> (&val);
    unsigned int return_value = atomicCAS((unsigned int*) dest_int,compare_int,val_int);
	return *reinterpret_cast<float*> (&return_value);
  #endif
}

template<typename T>
KOKKOSARRAY_INLINE_FUNCTION
typename KokkosArray::Impl::enable_if<KokkosArray::Impl::is_same<T,double>::value, double  >::type
atomic_compare_exchange(volatile T* dest, T compare, T val) {
  #if defined(GNU_ATOMICS_GCC) || defined(GNU_ATOMICS_INTEL)
    volatile int64_t* dest_int = reinterpret_cast<volatile int64_t*>(dest);
    int64_t compare_int = *reinterpret_cast<int64_t*> (&compare);
    int64_t val_int = *reinterpret_cast<int64_t*> (&val);
	int64_t return_value = __sync_val_compare_and_swap(dest_int,compare_int,val_int);
	return *reinterpret_cast<double*> (&return_value);
  #endif

  #ifdef OMP31_ATOMICS
	double retval;
    #pragma omp critical
	{
	  retval = dest[0];
	  if(retval == compare)
		dest[0] = val;
	}
    return retval;
  #endif

  #ifdef CUDA_ATOMICS
    volatile unsigned long long int* dest_int = reinterpret_cast<volatile unsigned long long int*>(dest);
    unsigned long long int compare_int = *reinterpret_cast<unsigned long long int*> (&compare);
    unsigned long long int val_int = *reinterpret_cast<unsigned long long int*> (&val);
    unsigned long long int return_value = atomicCAS((unsigned long long int*) dest_int,compare_int,val_int);
	return *reinterpret_cast<double*> (&return_value);
  #endif
}
