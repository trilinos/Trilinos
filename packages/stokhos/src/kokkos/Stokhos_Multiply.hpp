// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_MULTIPLY_HPP
#define STOKHOS_MULTIPLY_HPP

//#include "Kokkos_Macros.hpp"
//#include "Kokkos_Pair.hpp"
//#include "impl/Kokkos_Traits.hpp"

#include "Kokkos_Core.hpp"

#include <vector> // for std::vector (needed below)

namespace Stokhos {

class DefaultMultiply {};

template <unsigned> class IntegralRank {};

template <typename T> struct ViewRank {
  typedef IntegralRank< T::Rank > type;
};

template <typename T> struct ViewRank< std::vector<T> > {
  typedef IntegralRank< T::Rank > type;
};

template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          typename ColumnIndicesType = void,
          typename VectorRank = typename ViewRank<InputVectorType>::type,
          typename ImplTag = DefaultMultiply
          > class Multiply;

template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType>
void multiply(const MatrixType& A,
              const InputVectorType& x,
              OutputVectorType& y) {
  typedef Multiply<MatrixType,InputVectorType,OutputVectorType> multiply_type;
  multiply_type::apply( A, x, y );
}

namespace { // (anonymous)

// Work-around for CWG 1558.  See
// https://en.cppreference.com/w/cpp/types/void_t
template<class... Ts> struct make_void { typedef void type; };
template<class... Ts>
using replace_me_with_void_t_in_cxx17 =
  typename make_void<Ts...>::type;

template<class T, class = replace_me_with_void_t_in_cxx17<> >
struct const_type_impl {
  using type = T;
};

template<class T>
struct const_type_impl<T,
  replace_me_with_void_t_in_cxx17<typename T::const_type> > {
  using type = typename T::const_type;
};

template<class T>
using const_type_t = typename const_type_impl<T>::type;

} // namespace (anonymous)

template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType>
void multiply(const MatrixType& A,
              const InputVectorType& x,
              OutputVectorType& y,
              DefaultMultiply tag) {
  // mfh 29 Jul 2019: Not sure why, but std::vector claims to be a
  // Kokkos::View using Kokkos::is_view.  This is why I check instead
  // whether the class has a const_type typedef.
  using input_vector_type = const_type_t<InputVectorType>;
  using multiply_type =
    Multiply<MatrixType, input_vector_type, OutputVectorType>;
  multiply_type::apply( A, x, y );
}

template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          typename ColumnIndicesType>
void multiply(const MatrixType& A,
              const InputVectorType& x,
              OutputVectorType& y,
              const ColumnIndicesType& col) {
  typedef Multiply<MatrixType,InputVectorType,OutputVectorType,ColumnIndicesType> multiply_type;
  multiply_type::apply( A, x, y, col );
}

template <typename MatrixType,
          typename InputVectorType,
          typename OutputVectorType,
          typename ColumnIndicesType>
void multiply(const MatrixType& A,
              const InputVectorType& x,
              OutputVectorType& y,
              const ColumnIndicesType& col,
              DefaultMultiply tag) {
  typedef Multiply<MatrixType,InputVectorType,OutputVectorType,ColumnIndicesType> multiply_type;
  multiply_type::apply( A, x, y, col );
}

template <typename BlockSpec> class BlockMultiply;

namespace details {

/*
 * Compute work range = (begin, end) such that adjacent threads/blocks write to
 * separate cache lines
 */
template <typename scalar_type, typename execution_space, typename size_type>
KOKKOS_INLINE_FUNCTION
Kokkos::pair<size_type, size_type>
compute_work_range( const execution_space device,
                    const size_type work_count,
                    const size_type thread_count,
                    const size_type thread_rank)
{
#if defined( KOKKOS_ENABLE_CUDA )
  enum { cache_line =
         std::is_same<execution_space,Kokkos::Cuda>::value ? 128 : 64 };
#else
  enum { cache_line = 64 };
#endif

  enum { work_align = cache_line / sizeof(scalar_type) };
  enum { work_shift = Kokkos::Impl::power_of_two< work_align >::value };
  enum { work_mask  = work_align - 1 };

  const size_type work_per_thread =
    ( ( ( ( work_count + work_mask ) >> work_shift ) + thread_count - 1 ) /
      thread_count ) << work_shift ;

  size_type work_begin = thread_rank * work_per_thread;
  size_type work_end = work_begin + work_per_thread;
  if (work_begin > work_count)
    work_begin = work_count;
  if (work_end > work_count)
    work_end = work_count;

  return Kokkos::make_pair(work_begin, work_end);
}

// Functor implementing assignment update for multiply kernels
struct MultiplyAssign {
  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION
  void operator()(Scalar& y, const Scalar& x) const { y = x; }
};

// Functor implementing += update for multiply kernels
struct MultiplyUpdate {
  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION
  void operator()(Scalar& y, const Scalar& x) const { y += x; }
};

// Functor implementing scaled assignment update for multiply kernels
template <typename Value>
struct MultiplyScaledAssign {
  const Value a;
  MultiplyScaledAssign(const Value& a_) : a(a_) {}
  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION
  void operator()(Scalar& y, const Scalar& x) const { y = a*x; }
};

// Functor implementing += update for multiply kernels
template <typename Value>
struct MultiplyScaledUpdate {
  const Value a;
  MultiplyScaledUpdate(const Value& a_) : a(a_) {}
  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION
  void operator()(Scalar& y, const Scalar& x) const { y += a*x; }
};

// Functor implementing saxpby update for multiply kernels
template <typename Value>
struct MultiplyScaledUpdate2 {
  const Value a;
  const Value b;
  MultiplyScaledUpdate2(const Value& a_, const Value& b_) : a(a_), b(b_) {}
  template <typename Scalar>
  KOKKOS_INLINE_FUNCTION
  void operator()(Scalar& y, const Scalar& x) const { y = a*x + b*y; }
};

} // namespace details

} // namespace Stokhos

#endif
