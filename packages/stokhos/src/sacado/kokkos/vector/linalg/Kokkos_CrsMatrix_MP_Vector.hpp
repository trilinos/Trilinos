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
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_MV_MP_Vector.hpp" // for some utilities

#include "Kokkos_Core.hpp"
#include "Stokhos_Multiply.hpp"

#include "Teuchos_TestForException.hpp"

//----------------------------------------------------------------------------
// Specializations of Kokkos::CrsMatrix for Sacado::MP::Vector scalar type
//----------------------------------------------------------------------------

namespace Stokhos {

namespace details {

template <typename Matrix, typename InputVector, typename OutputVector,
          typename Update = MultiplyAssign>
class MPMultiply {};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory,
          typename Update>
class MPMultiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                     MatrixOrdinal,
                                     Device,
                                     MatrixMemory,
                                     MatrixSize>,
                  Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                                InputLayout,
                                Device,
                                InputMemory >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                                OutputLayout,
                                Device,
                                OutputMemory >,
                  Update
                >
{

public:

  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;
  typedef OutputVectorValue scalar_type;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef Kokkos::View< InputVectorValue*,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;
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
    const size_type row_count = A.graph.row_map.dimension_0()-1;
    Kokkos::parallel_for( row_count, MPMultiply(A,x,y,update) );
  }
};

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>**,...>,
//   x and y are rank 2, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory,
          typename Update>
class MPMultiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                     MatrixOrdinal,
                                     Device,
                                     MatrixMemory,
                                     MatrixSize>,
                  Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                                InputLayout,
                                Device,
                                InputMemory >,
                  Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                                OutputLayout,
                                Device,
                                OutputMemory >,
                  Update
                  >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;
  typedef OutputVectorValue scalar_type;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue**,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;
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
    const size_type num_col = m_y.dimension_1();
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
    const size_type row_count = A.graph.row_map.dimension_0()-1;
    Kokkos::parallel_for( row_count, MPMultiply(A,x,y,update) );
  }
};

} // namespace details

// Kernel implementing y = A * x where
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>*,...>,
//   x and y are rank 1, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                              InputLayout,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>*,
                              OutputLayout,
                              Device,
                              OutputMemory >
                >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue*,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue*,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;

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
//   A == Kokkos::CrsMatrix< Sacado::MP::Vector<...>,...>,
//   x, y == Kokkos::View< Sacado::MP::Vector<...>**,...>,
//   x and y are rank 2, any layout
// We spell everything out here to make sure the ranks and devices match.
//
// This implementation uses overloaded operators for MP::Vector.
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
class Multiply< Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                                   MatrixOrdinal,
                                   Device,
                                   MatrixMemory,
                                   MatrixSize>,
                Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                              InputLayout,
                              Device,
                              InputMemory >,
                Kokkos::View< Sacado::MP::Vector<OutputStorage>**,
                              OutputLayout,
                              Device,
                              OutputMemory >
                >
{
public:
  typedef Sacado::MP::Vector<MatrixStorage> MatrixValue;
  typedef Sacado::MP::Vector<InputStorage> InputVectorValue;
  typedef Sacado::MP::Vector<OutputStorage> OutputVectorValue;

  typedef Device device_type;
  typedef typename device_type::size_type size_type;

  typedef Kokkos::CrsMatrix< MatrixValue,
                             MatrixOrdinal,
                             Device,
                             MatrixMemory,
                             MatrixSize> matrix_type;
  typedef typename matrix_type::values_type matrix_values_type;
  typedef Kokkos::View< InputVectorValue**,
                        InputLayout,
                        Device,
                        InputMemory > input_vector_type;
  typedef Kokkos::View< OutputVectorValue**,
                        OutputLayout,
                        Device,
                        OutputMemory > output_vector_type;

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

namespace Kokkos {

// Overload of Kokkos::MV_Multiply for Sacado::MP::Vector scalar types
template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
                                          OutputLayout,
                                          Device,
                                          OutputMemory >& y,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
    OutputLayout, Device, OutputMemory > OutputVectorType;
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
    MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
  typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*,
    InputLayout, Device, InputMemory > InputVectorType;
  typedef Stokhos::details::MPMultiply<MatrixType,InputVectorType,
    OutputVectorType> multiply_type;

  multiply_type::apply( A, x, y, Stokhos::details::MultiplyAssign() );
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
                     OutputLayout,
                     Device,
                     OutputMemory >& y,
  const Sacado::MP::Vector<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  typedef typename InputStorage::value_type value_type;
  typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
    OutputLayout, Device, OutputMemory > OutputVectorType;
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
    MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
  typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*,
    InputLayout, Device, InputMemory > InputVectorType;

  if (!Sacado::is_constant(a)) {
    Impl::raise_error(
      "MV_Multiply not implemented for non-constant a");
  }

  value_type aa = a.fastAccessCoeff(0);
  if (aa == value_type(1)) {
    // y = A*x
    typedef Stokhos::details::MultiplyAssign UpdateType;
    typedef Stokhos::details::MPMultiply<MatrixType,InputVectorType,
      OutputVectorType,UpdateType> multiply_type;
    multiply_type::apply( A, x, y, UpdateType() );
  }
  else {
    // y = a*A*x
    typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
    typedef Stokhos::details::MPMultiply<MatrixType,InputVectorType,
      OutputVectorType,UpdateType> multiply_type;
    multiply_type::apply( A, x, y, UpdateType(aa) );
  }
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Sacado::MP::Vector<InputStorage>& b,
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
                     OutputLayout,
                     Device,
                     OutputMemory >& y,
  const Sacado::MP::Vector<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>*,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  typedef typename InputStorage::value_type value_type;
  typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*,
    OutputLayout, Device, OutputMemory > OutputVectorType;
  typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
    MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
  typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*,
    InputLayout, Device, InputMemory > InputVectorType;

  if (!Sacado::is_constant(a) || !Sacado::is_constant(b)) {
    Impl::raise_error(
      "MV_Multiply not implemented for non-constant a or b");
  }

  value_type aa = a.fastAccessCoeff(0);
  value_type bb = b.fastAccessCoeff(0);
  if (bb == value_type(0)) {
    if (aa == value_type(1)) {
      // y = A*x
      typedef Stokhos::details::MultiplyAssign UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        InputVectorType,OutputVectorType,UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType() );
    }
    else {
      // y = a*A*x
      typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        InputVectorType,OutputVectorType,UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa) );
    }
  }
  else if (bb == value_type(1)) {
    if (aa == value_type(1)) {
      // y += A*x
      typedef Stokhos::details::MultiplyUpdate UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        InputVectorType,OutputVectorType,UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType() );
    }
    else {
      // y += a*A*x
      typedef Stokhos::details::MultiplyScaledUpdate<value_type> UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        InputVectorType,OutputVectorType,UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa) );
    }
  }
  else {
    // y = a*A*x + b*y
    typedef Stokhos::details::MultiplyScaledUpdate2<value_type> UpdateType;
    typedef Stokhos::details::MPMultiply<MatrixType,
        InputVectorType,OutputVectorType,UpdateType> multiply_type;
    multiply_type::apply( A, x, y, UpdateType(aa,bb) );
  }
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
                      OutputLayout,
                      Device,
                      OutputMemory >& y,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  if (y.dimension_1() == 1) {
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*, OutputLayout,
      Device,OutputMemory > OutputView1D;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*, InputLayout,
      Device, InputMemory > InputView1D;
    OutputView1D y_1D = subview<OutputView1D>(y, ALL(), 0);
    InputView1D x_1D = subview<InputView1D>(x, ALL(), 0);
    MV_Multiply(y_1D, A, x_1D);
  }
  else {
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
      OutputLayout, Device, OutputMemory > OutputVectorType;
    typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
      MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>**,
      InputLayout, Device, InputMemory > InputVectorType;

    // y = A*x
    typedef Stokhos::details::MultiplyAssign UpdateType;
    typedef Stokhos::details::MPMultiply<MatrixType,
      InputVectorType,OutputVectorType,UpdateType> multiply_type;
    multiply_type::apply( A, x, y, UpdateType() );
  }
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
                      OutputLayout,
                      Device,
                      OutputMemory >& y,
  const Sacado::MP::Vector<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  if (y.dimension_1() == 1) {
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*, OutputLayout,
      Device,OutputMemory > OutputView1D;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*, InputLayout,
      Device, InputMemory > InputView1D;
    OutputView1D y_1D = subview<OutputView1D>(y, ALL(), 0);
    InputView1D x_1D = subview<InputView1D>(x, ALL(), 0);
    MV_Multiply(y_1D, a, A, x_1D);
  }
  else {
    typedef typename InputStorage::value_type value_type;
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
      OutputLayout, Device, OutputMemory > OutputVectorType;
    typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
      MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>**,
      InputLayout, Device, InputMemory > InputVectorType;

    if (!Sacado::is_constant(a)) {
      Impl::raise_error(
        "MV_Multiply not implemented for non-constant a");
    }

    value_type aa = a.fastAccessCoeff(0);
    if (aa == value_type(1)) {
      // y = A*x
      typedef Stokhos::details::MultiplyAssign UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        InputVectorType,OutputVectorType,UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType() );
    }
    else {
      // y = a*A*x
      typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        InputVectorType,OutputVectorType,UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa) );
    }
  }
}

template <typename Device,
          typename MatrixStorage,
          typename MatrixOrdinal,
          typename MatrixMemory,
          typename MatrixSize,
          typename InputStorage,
          typename InputLayout,
          typename InputMemory,
          typename OutputStorage,
          typename OutputLayout,
          typename OutputMemory>
void
MV_Multiply(
  const Sacado::MP::Vector<InputStorage>& b,
  const Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
                      OutputLayout,
                      Device,
                      OutputMemory >& y,
  const Sacado::MP::Vector<InputStorage>& a,
  const Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
                           MatrixOrdinal,
                           Device,
                           MatrixMemory,
                           MatrixSize>& A,
  const Kokkos::View< Sacado::MP::Vector<InputStorage>**,
                      InputLayout,
                      Device,
                      InputMemory >& x)
{
  if (y.dimension_1() == 1) {
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>*, OutputLayout,
      Device,OutputMemory > OutputView1D;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>*, InputLayout,
      Device, InputMemory > InputView1D;
    OutputView1D y_1D = subview<OutputView1D>(y, ALL(), 0);
    InputView1D x_1D = subview<InputView1D>(x, ALL(), 0);
    MV_Multiply(b, y_1D, a, A, x_1D);
  }
  else {
    typedef typename InputStorage::value_type value_type;
    typedef Kokkos::View< Sacado::MP::Vector< OutputStorage>**,
      OutputLayout, Device, OutputMemory > OutputVectorType;
    typedef Kokkos::CrsMatrix< Sacado::MP::Vector<MatrixStorage>,
      MatrixOrdinal, Device, MatrixMemory, MatrixSize> MatrixType;
    typedef Kokkos::View< Sacado::MP::Vector<InputStorage>**,
      InputLayout, Device, InputMemory > InputVectorType;

    if (!Sacado::is_constant(a) || !Sacado::is_constant(b)) {
      Impl::raise_error(
        "MV_Multiply not implemented for non-constant a or b");
    }

    value_type aa = a.fastAccessCoeff(0);
    value_type bb = b.fastAccessCoeff(0);
    if (bb == value_type(0)) {
      if (aa == value_type(1)) {
        // y = A*x
        typedef Stokhos::details::MultiplyAssign UpdateType;
        typedef Stokhos::details::MPMultiply<MatrixType,
          InputVectorType,OutputVectorType,UpdateType> multiply_type;
        multiply_type::apply( A, x, y, UpdateType() );
      }
      else {
        // y = a*A*x
        typedef Stokhos::details::MultiplyScaledAssign<value_type> UpdateType;
        typedef Stokhos::details::MPMultiply<MatrixType,
          InputVectorType,OutputVectorType,UpdateType> multiply_type;
        multiply_type::apply( A, x, y, UpdateType(aa) );
      }
    }
    else if (bb == value_type(1)) {
      if (aa == value_type(1)) {
        // y += A*x
        typedef Stokhos::details::MultiplyUpdate UpdateType;
        typedef Stokhos::details::MPMultiply<MatrixType,
          InputVectorType,OutputVectorType,UpdateType> multiply_type;
        multiply_type::apply( A, x, y, UpdateType() );
      }
      else {
        // y += a*A*x
        typedef Stokhos::details::MultiplyScaledUpdate<value_type> UpdateType;
        typedef Stokhos::details::MPMultiply<MatrixType,
          InputVectorType,OutputVectorType,UpdateType> multiply_type;
        multiply_type::apply( A, x, y, UpdateType(aa) );
      }
    }
    else {
      // y = a*A*x + b*y
      typedef Stokhos::details::MultiplyScaledUpdate2<value_type> UpdateType;
      typedef Stokhos::details::MPMultiply<MatrixType,
        InputVectorType,OutputVectorType,UpdateType> multiply_type;
      multiply_type::apply( A, x, y, UpdateType(aa,bb) );
    }
  }
}

}

#endif /* #ifndef KOKKOS_CRSMATRIX_MP_VECTOR_HPP */
