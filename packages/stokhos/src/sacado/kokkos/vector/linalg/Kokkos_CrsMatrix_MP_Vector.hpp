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

#ifndef KOKKOS_CRSMATRIX_MP_VECTOR_HPP
#define KOKKOS_CRSMATRIX_MP_VECTOR_HPP

#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "Kokkos_ArithTraits_MP_Vector.hpp"
#include "KokkosSparse_spmv.hpp"
#include "Kokkos_Blas1_MP_Vector.hpp" // for some utilities

#include "Kokkos_Core.hpp"
#include "Stokhos_Multiply.hpp"

#include "Teuchos_TestForException.hpp"

//----------------------------------------------------------------------------
// Specializations of KokkosSparse::CrsMatrix for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

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

namespace Stokhos {

namespace details {

template <typename Matrix, typename InputVector, typename OutputVector,
          typename Update = MultiplyAssign,
          typename Enabled = void>
class MPMultiply {};

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< const Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP,
          typename Update>
class MPMultiply< KokkosSparse::CrsMatrix< const Sacado::MP::Vector<MatrixStorage>,
                                           MatrixOrdinal,
                                           MatrixDevice,
                                           MatrixMemory,
                                           MatrixSize>,
                  Kokkos::View< const Sacado::MP::Vector<InputStorage>*,
                                InputP... >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                OutputP... >,
                  Update
#ifdef KOKKOS_ENABLE_CUDA
                  , typename std::enable_if<
                      !std::is_same<typename MatrixDevice::execution_space,Kokkos::Cuda>::value >::type
#endif
                  >
{

public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;
  typedef OutputVectorValue scalar_type;

  typedef typename MatrixDevice::execution_space execution_space;
  typedef typename execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix< const MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize > matrix_type;
  typedef Kokkos::View< const InputVectorValue*,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputP... > output_vector_type;
  typedef Update update_type;

  const matrix_type  m_A;
  const input_vector_type  m_x;
  const output_vector_type  m_y;
  const update_type m_update;

  MPMultiply( const matrix_type & A,
              const input_vector_type & x,
              const output_vector_type & y,
              const update_type& update )
    : m_A( A )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    // Compute mat-vec for this row
    const size_type iEntryBegin = m_A.graph.row_map[iRow];
    const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];
    scalar_type sum = 0.0;
    for (size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry) {
      size_type iCol = m_A.graph.entries(iEntry);
      sum += m_A.values(iEntry) * m_x(iCol);
    }
    m_update( m_y(iRow), sum );
  } // operator()

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    const size_type row_count = A.graph.row_map.extent(0)-1;
    Kokkos::parallel_for( row_count, MPMultiply(A,x,y,update) );
  }
};

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>**,...>,
//   x and y are rank 2, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP,
          typename Update>
class MPMultiply< KokkosSparse::CrsMatrix< const Sacado::MP::Vector<MatrixStorage>,
                                           MatrixOrdinal,
                                           MatrixDevice,
                                           MatrixMemory,
                                           MatrixSize >,
                  Kokkos::View< const Sacado::MP::Vector<InputStorage>**,
                                InputP... >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                                OutputP... >,
                  Update
#ifdef KOKKOS_ENABLE_CUDA
                  , typename std::enable_if<
                    !std::is_same<typename MatrixDevice::execution_space,Kokkos::Cuda>::value >::type
#endif
                  >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;
  typedef OutputVectorValue scalar_type;

  typedef typename MatrixDevice::execution_space execution_space;
  typedef typename execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix< const MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize > matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;
  typedef Update update_type;

  const matrix_type  m_A;
  const input_vector_type  m_x;
  const output_vector_type  m_y;
  const update_type m_update;

  MPMultiply( const matrix_type & A,
              const input_vector_type & x,
              const output_vector_type & y,
              const update_type& update )
    : m_A( A )
    , m_x( x )
    , m_y( y )
    , m_update( update )
    {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type iRow ) const
  {
    scalar_type sum;
    // Loop over columns of x, y
    const size_type num_col = m_y.extent(1);
    for (size_type col=0; col<num_col; ++col) {
      // Compute mat-vec for this row
      const size_type iEntryBegin = m_A.graph.row_map[iRow];
      const size_type iEntryEnd   = m_A.graph.row_map[iRow+1];
      sum = 0.0;
      for (size_type iEntry = iEntryBegin; iEntry < iEntryEnd; ++iEntry) {
        size_type iCol = m_A.graph.entries(iEntry);
        sum += m_A.values(iEntry) * m_x(iCol,col);
      }
      m_update( m_y(iRow,col), sum );
    } // x, y column loop
  } // operator()

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    const size_type row_count = A.graph.row_map.extent(0)-1;
    Kokkos::parallel_for( row_count, MPMultiply(A,x,y,update) );
  }
};

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP,
          typename Update>
class MPMultiply< KokkosSparse::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                           MatrixOrdinal,
                                           MatrixDevice,
                                           MatrixMemory,
                                           MatrixSize>,
                  Kokkos::View< const Sacado::MP::Vector<InputStorage>*,
                                InputP... >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                OutputP... >,
                  Update
                  >
{

public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;
  typedef OutputVectorValue scalar_type;

  typedef typename MatrixDevice::execution_space execution_space;
  typedef typename execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix< MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize > matrix_type;
  typedef typename matrix_type::const_type const_matrix_type;

  typedef Kokkos::View< const InputVectorValue*,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputP... > output_vector_type;
  typedef Update update_type;

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    const_matrix_type cA = A;
    MPMultiply< const_matrix_type, input_vector_type, output_vector_type, update_type >::apply(
      cA, x, y, update);
  }
};

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>**,...>,
//   x and y are rank 2, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP,
          typename Update>
class MPMultiply< KokkosSparse::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                           MatrixOrdinal,
                                           MatrixDevice,
                                           MatrixMemory,
                                           MatrixSize>,
                  Kokkos::View< const Sacado::MP::Vector<InputStorage>**,
                                InputP... >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                                OutputP... >,
                  Update
                  >
{

public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;
  typedef OutputVectorValue scalar_type;

  typedef typename MatrixDevice::execution_space execution_space;
  typedef typename execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix< MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize > matrix_type;
  typedef typename matrix_type::const_type const_matrix_type;

  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;
  typedef Update update_type;

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y,
                     const update_type & update )
  {
    const_matrix_type cA = A;
    MPMultiply< const_matrix_type, input_vector_type, output_vector_type, update_type >::apply(
      cA, x, y, update);
  }
};

} // namespace details

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class Multiply< KokkosSparse::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
                                         MatrixMemory,
                                         MatrixSize >,
                Kokkos::View< const Sacado::MP::Vector<InputStorage>*,
                              InputP... >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                              OutputP... >
                >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef typename MatrixDevice::execution_space execution_space;
  typedef typename execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix< MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize > matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< const InputVectorValue*,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputP... > output_vector_type;

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y )
  {
    typedef details::MPMultiply<matrix_type,input_vector_type,output_vector_type> multiply_type;
    multiply_type::apply(A,x,y,details::MultiplyAssign());
  }
};

// Kernel implementing y = A * x where
//   A == KokkosSparse::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>**,...>,
//   x and y are rank 2, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename MatrixDevice,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename ... InputP,
          typename OutputStorage,
          typename ... OutputP>
class Multiply< KokkosSparse::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                         MatrixOrdinal,
                                         MatrixDevice,
                                         MatrixMemory,
                                         MatrixSize >,
                Kokkos::View< const Sacado::MP::Vector<InputStorage>**,
                              InputP... >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                              OutputP... >
                >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef typename MatrixDevice::execution_space execution_space;
  typedef typename execution_space::size_type size_type;

  typedef KokkosSparse::CrsMatrix< MatrixValue,
                                   MatrixOrdinal,
                                   MatrixDevice,
                                   MatrixMemory,
                                   MatrixSize > matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< const InputVectorValue**,
                        InputP... > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputP... > output_vector_type;

public:

  static void apply( const matrix_type & A,
                     const input_vector_type & x,
                     const output_vector_type & y )
  {
    typedef details::MPMultiply<matrix_type,input_vector_type,output_vector_type> multiply_type;
    multiply_type::apply(A,x,y,details::MultiplyAssign());
  }
};

} // namespace Stokhos

namespace KokkosSparse {

template <typename AlphaType,
          typename BetaType,
          typename MatrixType,
          typename InputType,
          typename ... InputP,
          typename OutputType,
          typename ... OutputP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View< InputType, InputP... > >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View< OutputType, OutputP... > >::value
  >::type
spmv(
  const char mode[],
  const AlphaType& a,
  const MatrixType& A,
  const Kokkos::View< InputType, InputP... >& x,
  const BetaType& b,
  const Kokkos::View< OutputType, OutputP... >& y,
  const RANK_ONE)
{
  typedef Kokkos::View< OutputType, OutputP... > OutputVectorType;
  typedef Kokkos::View< InputType, InputP... > InputVectorType;
  using input_vector_type = const_type_t<InputVectorType>;
  typedef typename InputVectorType::array_type::non_const_value_type value_type;

  if(mode[0]!='N') {
    Kokkos::Impl::raise_error(
      "Stokhos spmv not implemented for transposed or conjugated matrix-vector multiplies");
  }

  if (!Sacado::is_constant(a) || !Sacado::is_constant(b)) {
    Kokkos::Impl::raise_error(
      "MV_Multiply not implemented for non-constant a or b");
  }

  value_type aa = Sacado::Value<AlphaType>::eval(a);
  value_type bb = Sacado::Value<BetaType>::eval(b);
  if (bb == value_type(0)) {
    if (aa == value_type(1)) {
      // y = A*x
      typedef Stokhos::details::MultiplyAssign UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        input_vector_type, OutputVectorType,
        UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType() );
    }
    else {
      // y = a*A*x
      typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        input_vector_type, OutputVectorType,
        UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa) );
    }
  }
  else if (bb == value_type(1)) {
    if (aa == value_type(1)) {
      // y += A*x
      typedef Stokhos::details::MultiplyUpdate UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        input_vector_type, OutputVectorType,
        UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType() );
    }
    else {
      // y += a*A*x
      typedef Stokhos::details::MultiplyScaledUpdate<value_type> UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        input_vector_type, OutputVectorType,
        UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa) );
    }
  }
  else {
    // y = a*A*x + b*y
    typedef Stokhos::details::MultiplyScaledUpdate2<value_type> UpdateType;
    typedef Stokhos::details::MPMultiply<MatrixType,
      input_vector_type, OutputVectorType,
      UpdateType> multiply_type;
    multiply_type::apply( A, x, y, UpdateType(aa,bb) );
  }
}

template <typename AlphaType,
          typename BetaType,
          typename MatrixType,
          typename InputType,
          typename ... InputP,
          typename OutputType,
          typename ... OutputP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View< InputType, InputP... > >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View< OutputType, OutputP... > >::value
  >::type
spmv(
  KokkosKernels::Experimental::Controls,
  const char mode[],
  const AlphaType& a,
  const MatrixType& A,
  const Kokkos::View< InputType, InputP... >& x,
  const BetaType& b,
  const Kokkos::View< OutputType, OutputP... >& y,
  const RANK_ONE)
{
  spmv(mode, a, A, x, b, y, RANK_ONE());
}

template <typename AlphaType,
          typename BetaType,
          typename MatrixType,
          typename InputType,
          typename ... InputP,
          typename OutputType,
          typename ... OutputP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View< InputType, InputP... > >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View< OutputType, OutputP... > >::value
  >::type
spmv(
  const char mode[],
  const AlphaType& a,
  const MatrixType& A,
  const Kokkos::View< InputType, InputP... >& x,
  const BetaType& b,
  const Kokkos::View< OutputType, OutputP... >& y,
  const RANK_TWO)
{
  if(mode[0]!='N') {
    Kokkos::Impl::raise_error(
      "Stokhos spmv not implemented for transposed or conjugated matrix-vector multiplies");
  }
  if (y.extent(1) == 1) {
    auto y_1D = subview(y, Kokkos::ALL(), 0);
    auto x_1D = subview(x, Kokkos::ALL(), 0);
    spmv(mode, a, A, x_1D, b, y_1D, RANK_ONE());
  }
  else {
    typedef Kokkos::View< OutputType, OutputP... > OutputVectorType;
    typedef Kokkos::View< InputType, InputP... > InputVectorType;
    using input_vector_type = const_type_t<InputVectorType>;
    typedef typename InputVectorType::array_type::non_const_value_type value_type;

    if (!Sacado::is_constant(a) || !Sacado::is_constant(b)) {
      Kokkos::Impl::raise_error(
        "Stokhos spmv not implemented for non-constant a or b");
    }

    value_type aa = Sacado::Value<AlphaType>::eval(a);
    value_type bb = Sacado::Value<BetaType>::eval(b);
    if (bb == value_type(0)) {
      if (aa == value_type(1)) {
        // y = A*x
        typedef Stokhos::details::MultiplyAssign UpdateType;
        typedef Stokhos::details::MPMultiply<MatrixType,
          input_vector_type, OutputVectorType,
          UpdateType> multiply_type;
        multiply_type::apply( A, x, y, UpdateType() );
      }
      else {
        // y = a*A*x
        typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
        typedef Stokhos::details::MPMultiply<MatrixType,
          input_vector_type, OutputVectorType,
          UpdateType> multiply_type;
        multiply_type::apply( A, x, y, UpdateType(aa) );
      }
    }
    else if (bb == value_type(1)) {
      if (aa == value_type(1)) {
        // y += A*x
        typedef Stokhos::details::MultiplyUpdate UpdateType;
        typedef Stokhos::details::MPMultiply<MatrixType,
          input_vector_type, OutputVectorType,
          UpdateType> multiply_type;
        multiply_type::apply( A, x, y, UpdateType() );
      }
      else {
        // y += a*A*x
        typedef Stokhos::details::MultiplyScaledUpdate<value_type> UpdateType;
        typedef Stokhos::details::MPMultiply<MatrixType,
          input_vector_type, OutputVectorType,
          UpdateType> multiply_type;
        multiply_type::apply( A, x, y, UpdateType(aa) );
      }
    }
    else {
      // y = a*A*x + b*y
      typedef Stokhos::details::MultiplyScaledUpdate2<value_type> UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        input_vector_type, OutputVectorType,
        UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa,bb) );
    }
  }
}

template <typename AlphaType,
          typename BetaType,
          typename MatrixType,
          typename InputType,
          typename ... InputP,
          typename OutputType,
          typename ... OutputP>
typename std::enable_if<
  Kokkos::is_view_mp_vector< Kokkos::View< InputType, InputP... > >::value &&
  Kokkos::is_view_mp_vector< Kokkos::View< OutputType, OutputP... > >::value
  >::type
spmv(
  KokkosKernels::Experimental::Controls,
  const char mode[],
  const AlphaType& a,
  const MatrixType& A,
  const Kokkos::View< InputType, InputP... >& x,
  const BetaType& b,
  const Kokkos::View< OutputType, OutputP... >& y,
  const RANK_TWO)
{
  spmv(mode, a, A, x, b, y, RANK_TWO());
}

}

#endif /* #ifndef KOKKOS_CRSMATRIX_MP_VECTOR_HPP */
