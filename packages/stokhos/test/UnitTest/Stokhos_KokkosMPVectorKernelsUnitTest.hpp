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

#ifndef STOKHOS_KOKKOS_CORE_KERNELS_UNIT_TEST_HPP
#define STOKHOS_KOKKOS_CORE_KERNELS_UNIT_TEST_HPP

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado_Kokkos.hpp"
#include "Kokkos_CrsMatrix_MP_Vector.hpp"

namespace KokkosMPVectorKernelsUnitTest {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::Array;
using Teuchos::ParameterList;

template< typename IntType >
inline
IntType map_fem_graph_coord( const IntType& N,
                             const IntType& i,
                             const IntType& j,
                             const IntType& k )
{
  return k + N * ( j + N * i );
}

template < typename ordinal >
inline
ordinal generate_fem_graph( ordinal N,
                            std::vector< std::vector<ordinal> >& graph )
{
  graph.resize( N * N * N, std::vector<ordinal>() );

  ordinal total = 0;

  for ( int i = 0; i < (int) N; ++i ) {
    for ( int j = 0; j < (int) N; ++j ) {
      for ( int k = 0; k < (int) N; ++k ) {

        const ordinal row = map_fem_graph_coord((int)N,i,j,k);

        graph[row].reserve(27);

        for ( int ii = -1; ii < 2; ++ii ) {
          for ( int jj = -1; jj < 2; ++jj ) {
            for ( int kk = -1; kk < 2; ++kk ) {
              if ( 0 <= i + ii && i + ii < (int) N &&
                   0 <= j + jj && j + jj < (int) N &&
                   0 <= k + kk && k + kk < (int) N ) {
                ordinal col = map_fem_graph_coord((int)N,i+ii,j+jj,k+kk);

                graph[row].push_back(col);
              }
            }}}
        total += graph[row].size();
      }}}

  return total;
}

template <typename scalar, typename ordinal>
inline
scalar generate_matrix_coefficient( const ordinal nFEM,
                                    const ordinal nStoch,
                                    const ordinal iRowFEM,
                                    const ordinal iColFEM,
                                    const ordinal iStoch )
{
  const scalar A_fem = ( 10.0 + scalar(iRowFEM) / scalar(nFEM) ) +
    (  5.0 + scalar(iColFEM) / scalar(nFEM) );

  const scalar A_stoch = ( 1.0 + scalar(iStoch) / scalar(nStoch) );

  return A_fem + A_stoch;
  //return 1.0;
}

template <typename scalar, typename ordinal>
inline
scalar generate_vector_coefficient( const ordinal nFEM,
                                    const ordinal nStoch,
                                    const ordinal iColFEM,
                                    const ordinal iStoch )
{
  const scalar X_fem = 100.0 + scalar(iColFEM) / scalar(nFEM);
  const scalar X_stoch =  1.0 + scalar(iStoch) / scalar(nStoch);
  return X_fem + X_stoch;
  //return 1.0;
}

template <typename Scalar> struct ScalarTol {};
template <> struct ScalarTol<float> {  static float  tol() { return 1e-5;  } };
template <> struct ScalarTol<double> { static double tol() { return 1e-12; } };

template <typename array_type, typename scalar_type>
bool compare_rank_2_views(const array_type& y,
                          const array_type& y_exp,
                          const scalar_type rel_tol,
                          const scalar_type abs_tol,
                          Teuchos::FancyOStream& out)
{
  typedef typename array_type::size_type size_type;
  typename array_type::HostMirror hy = Kokkos::create_mirror_view(y);
  typename array_type::HostMirror hy_exp = Kokkos::create_mirror_view(y_exp);
  Kokkos::deep_copy(hy, y);
  Kokkos::deep_copy(hy_exp, y_exp);

  size_type num_rows = y.dimension_0();
  size_type num_cols = y.dimension_1();
  bool success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      scalar_type diff = std::abs( hy(i,j) - hy_exp(i,j) );
      scalar_type tol = rel_tol*std::abs(hy_exp(i,j)) + abs_tol;
      bool s = diff < tol;
      out << "y_expected(" << i << "," << j << ") - "
          << "y(" << i << "," << j << ") = " << hy_exp(i,j)
          << " - " << hy(i,j) << " == "
          << diff << " < " << tol << " : ";
      if (s)
        out << "passed";
      else
        out << "failed";
      out << std::endl;
      success = success && s;
    }
  }

  return success;
}

template <typename VectorType, typename MultiplyTag>
bool test_embedded_vector(const typename VectorType::ordinal_type nGrid,
                          const typename VectorType::ordinal_type stoch_length,
                          Kokkos::DeviceConfig dev_config,
                          MultiplyTag tag,
                          Teuchos::FancyOStream& out)
{
  typedef typename VectorType::ordinal_type ordinal_type;
  typedef typename VectorType::value_type scalar_type;
  typedef typename VectorType::storage_type storage_type;
  typedef typename storage_type::device_type device_type;
  typedef Kokkos::LayoutRight Layout;
  typedef Kokkos::View< VectorType*, Layout, device_type > block_vector_type;
  typedef Kokkos::CrsMatrix< VectorType, ordinal_type, device_type > block_matrix_type;
  typedef typename block_matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename block_matrix_type::values_type matrix_values_type;

  // Generate FEM graph:
  ordinal_type fem_length = nGrid * nGrid * nGrid;
  std::vector< std::vector<ordinal_type> > fem_graph;
  ordinal_type fem_graph_length = generate_fem_graph( nGrid, fem_graph );

  //------------------------------
  // Generate input multivector:

  block_vector_type x =
    block_vector_type(Kokkos::allocate_without_initializing,
                      "x", fem_length, stoch_length);
  block_vector_type y =
    block_vector_type(Kokkos::allocate_without_initializing,
                      "y", fem_length, stoch_length);

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror_view( x );
  typename block_vector_type::HostMirror hy = Kokkos::create_mirror_view( y );

  // View the block vector as an array of the embedded intrinsic type.
  typename block_vector_type::HostMirror::array_type hax = hx ;
  typename block_vector_type::HostMirror::array_type hay = hy ;

  for (ordinal_type iRowFEM=0; iRowFEM<fem_length; ++iRowFEM) {
    for (ordinal_type iRowStoch=0; iRowStoch<stoch_length; ++iRowStoch) {
      hax(iRowFEM,iRowStoch) =
        generate_vector_coefficient<scalar_type>(
          fem_length, stoch_length, iRowFEM, iRowStoch );
      hay(iRowFEM,iRowStoch) = 0.0;
    }
  }

  Kokkos::deep_copy( x, hx );
  Kokkos::deep_copy( y, hy );

  //------------------------------
  // Generate block matrix

  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("test crs graph"), fem_graph);
  matrix_values_type matrix_values =
    matrix_values_type(
      Kokkos::allocate_without_initializing,
      "matrix", fem_graph_length, stoch_length);
  block_matrix_type matrix(
    "block_matrix", fem_length, matrix_values, matrix_graph);
  matrix.dev_config = dev_config;

  typename matrix_values_type::HostMirror hM =
    Kokkos::create_mirror_view( matrix.values );

  typename matrix_values_type::HostMirror::array_type haM = hM ;

  for (ordinal_type iRowFEM=0, iEntryFEM=0; iRowFEM<fem_length; ++iRowFEM) {
    const ordinal_type row_size = fem_graph[iRowFEM].size();
    for (ordinal_type iRowEntryFEM=0; iRowEntryFEM<row_size;
         ++iRowEntryFEM, ++iEntryFEM) {
      const ordinal_type iColFEM = fem_graph[iRowFEM][iRowEntryFEM];

      for (ordinal_type k=0; k<stoch_length; ++k) {
        haM(iEntryFEM,k) =
          generate_matrix_coefficient<scalar_type>(
            fem_length, stoch_length, iRowFEM, iColFEM, k);
      }
    }
  }

  Kokkos::deep_copy( matrix.values, hM );

  //------------------------------
  // multiply

  Stokhos::multiply( matrix, x, y, tag );

  //------------------------------
  // generate correct answer

  typedef typename block_vector_type::array_type array_type;
  array_type ay_expected =
    array_type("ay_expected", fem_length, stoch_length);
  typename array_type::HostMirror hay_expected =
    Kokkos::create_mirror_view(ay_expected);
  for (ordinal_type iRowFEM=0, iEntryFEM=0; iRowFEM<fem_length; ++iRowFEM) {
    const ordinal_type row_size = fem_graph[iRowFEM].size();
    for (ordinal_type iRowEntryFEM=0; iRowEntryFEM<row_size;
         ++iRowEntryFEM, ++iEntryFEM) {
      const ordinal_type iColFEM = fem_graph[iRowFEM][iRowEntryFEM];
      for (ordinal_type k=0; k<stoch_length; ++k) {
        hay_expected(iRowFEM, k) +=
          generate_matrix_coefficient<scalar_type>(
            fem_length, stoch_length, iRowFEM, iColFEM, k) *
          generate_vector_coefficient<scalar_type>(
            fem_length, stoch_length, iColFEM, k );
      }
    }
  }
  Kokkos::deep_copy( ay_expected, hay_expected );

  //------------------------------
  // check

  typename block_vector_type::array_type ay = y;
  scalar_type rel_tol = ScalarTol<scalar_type>::tol();
  scalar_type abs_tol = ScalarTol<scalar_type>::tol();
  bool success = compare_rank_2_views(ay, ay_expected, rel_tol, abs_tol, out);

  return success;
}



}

#endif
