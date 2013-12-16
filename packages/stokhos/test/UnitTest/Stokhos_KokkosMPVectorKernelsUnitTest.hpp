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

template <typename OrdinalType, typename ScalarType, typename Device>
struct UnitTestSetup {
  typedef OrdinalType ordinal_type;
  typedef ScalarType value_type;
  typedef typename Device::host_mirror_device_type host_device;
  typedef Kokkos::CrsMatrix<value_type,ordinal_type,host_device> matrix_type;
  typedef Kokkos::View<value_type[],host_device> vector_type;

  // The *local* vector size per thread.  Total size is the local
  // vector size * number of threads per team
  static const ordinal_type local_vec_size = 3;

  ordinal_type nGrid, fem_length, fem_graph_length, stoch_length;
  value_type rel_tol, abs_tol;
  std::vector< std::vector<ordinal_type> > fem_graph;
  std::vector< matrix_type > A;
  std::vector< vector_type > x, y;

  void setup(const ordinal_type threads_per_vector) {

    nGrid = 10;
    rel_tol = ScalarTol<ScalarType>::tol();
    abs_tol = ScalarTol<ScalarType>::tol();

    // Generate FEM graph:
    fem_length = nGrid * nGrid * nGrid;
    fem_graph_length = generate_fem_graph( nGrid, fem_graph );

    // Generate FEM matrices and vectors
    stoch_length = local_vec_size * threads_per_vector;
    A.resize(stoch_length);
    x.resize(stoch_length);
    y.resize(stoch_length);
    typedef typename matrix_type::StaticCrsGraphType matrix_graph_type;
    typedef typename matrix_type::values_type matrix_values_type;
    for (ordinal_type block=0; block<stoch_length; block++) {
      matrix_graph_type matrix_graph =
        Kokkos::create_staticcrsgraph<matrix_graph_type>(
          std::string("testing"), fem_graph);
      matrix_values_type matrix_values =
        matrix_values_type(
          Kokkos::allocate_without_initializing, "matrix", fem_graph_length);
      A[block] = matrix_type("matrix", fem_length, matrix_values, matrix_graph);
      x[block] = vector_type(
        Kokkos::allocate_without_initializing, "x", fem_length);
      y[block] = vector_type(
        Kokkos::allocate_without_initializing, "y", fem_length);

      for (ordinal_type iRowFEM=0, iEntryFEM=0; iRowFEM<fem_length; ++iRowFEM) {
        const ordinal_type row_size = fem_graph[iRowFEM].size();
        for (ordinal_type iRowEntryFEM=0; iRowEntryFEM<row_size;
             ++iRowEntryFEM, ++iEntryFEM) {
          const ordinal_type iColFEM = fem_graph[iRowFEM][iRowEntryFEM];

          value_type v =
            generate_matrix_coefficient<value_type>(
              fem_length, stoch_length, iRowFEM, iColFEM, block);
          A[block].values(iEntryFEM) = v;
        }

        x[block](iRowFEM) =
          generate_vector_coefficient<value_type>(
            fem_length, stoch_length, iRowFEM, block );
        y[block](iRowFEM) = 0.0;
      }
    }

    // Multiply
    for (ordinal_type block=0; block<stoch_length; block++)
      Kokkos::MV_Multiply(y[block], A[block], x[block]);
  }

  void cleanup() {
    x.resize(0);
    y.resize(0);
    A.resize(0);
  }

  template <typename vec_type>
  bool test_original(const vec_type& yy,
                     Teuchos::FancyOStream& out) const
  {
    typename vec_type::array_type ayy = yy ;

    bool success = true;
    for (ordinal_type block=0; block<stoch_length; ++block) {
      for (ordinal_type i=0; i<fem_length; ++i) {
        value_type diff = std::abs( y[block][i] - ayy(block,i) );
        value_type tol = rel_tol*std::abs(y[block][i]) + abs_tol;
        bool s = diff < tol;
        if (!s) {
          out << "y_expected[" << block << "][" << i << "] - "
              << "y(" << i << "," << block << ") = " << y[block][i]
              << " - " << ayy(block,i) << " == "
              << diff << " < " << tol << " : failed"
              << std::endl;
        }
        success = success && s;
      }
    }

    return success;
  }

  template <typename vec_type>
  bool test_commuted(const vec_type& yy,
                     Teuchos::FancyOStream& out) const
  {
    typename vec_type::array_type ayy = yy ;

    bool success = true;
    for (ordinal_type block=0; block<stoch_length; ++block) {
      for (ordinal_type i=0; i<fem_length; ++i) {
        value_type diff = std::abs( y[block][i] - ayy(i,block) );
        value_type tol = rel_tol*std::abs(y[block][i]) + abs_tol;
        bool s = diff < tol;
        if (!s) {
          out << "y_expected[" << block << "][" << i << "] - "
              << "y(" << i << "," << block << ") = " << y[block][i]
              << " - " << ayy(i,block) << " == "
              << diff << " < " << tol << " : failed"
              << std::endl;
        }
        success = success && s;
      }
    }

    return success;
  }

};

template <typename VectorType,
          typename OrdinalType, typename ScalarType, typename Device,
          typename MultiplyTag>
bool test_embedded_vector(
  const UnitTestSetup<OrdinalType,ScalarType,Device>& setup,
  Kokkos::DeviceConfig dev_config,
  MultiplyTag tag,
  Teuchos::FancyOStream& out)
{
  typedef OrdinalType ordinal_type;
  typedef ScalarType value_type;
  typedef typename VectorType::storage_type storage_type;
  typedef typename storage_type::device_type device_type;
  typedef Kokkos::LayoutRight Layout;
  typedef Kokkos::View< VectorType*, Layout, device_type > block_vector_type;
  typedef Kokkos::CrsMatrix< VectorType, ordinal_type, device_type > block_matrix_type;
  typedef typename block_matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename block_matrix_type::values_type matrix_values_type;

  //------------------------------
  // Generate input multivector:

  block_vector_type x =
    block_vector_type(Kokkos::allocate_without_initializing,
                      "x", setup.fem_length, setup.stoch_length);
  block_vector_type y =
    block_vector_type(Kokkos::allocate_without_initializing,
                      "y", setup.fem_length, setup.stoch_length);

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror_view( x );
  typename block_vector_type::HostMirror hy = Kokkos::create_mirror_view( y );

  // View the block vector as an array of the embedded intrinsic type.
  typename block_vector_type::HostMirror::array_type hax = hx ;
  typename block_vector_type::HostMirror::array_type hay = hy ;

  for (ordinal_type iRowFEM=0; iRowFEM<setup.fem_length; ++iRowFEM) {
    for (ordinal_type iRowStoch=0; iRowStoch<setup.stoch_length; ++iRowStoch) {
      hax(iRowFEM,iRowStoch) =
        generate_vector_coefficient<ScalarType>(
          setup.fem_length, setup.stoch_length, iRowFEM, iRowStoch );
      hay(iRowFEM,iRowStoch) = 0.0;
    }
  }

  Kokkos::deep_copy( x, hx );
  Kokkos::deep_copy( y, hy );

  //------------------------------

  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("test crs graph"), setup.fem_graph);
  matrix_values_type matrix_values =
    matrix_values_type(
      Kokkos::allocate_without_initializing,
      "matrix", setup.fem_graph_length, setup.stoch_length);
  block_matrix_type matrix(
    "block_matrix", setup.fem_length, matrix_values, matrix_graph);
  matrix.dev_config = dev_config;

  typename matrix_values_type::HostMirror hM =
    Kokkos::create_mirror_view( matrix.values );

  typename matrix_values_type::HostMirror::array_type haM = hM ;

  for (ordinal_type iRowFEM=0, iEntryFEM=0; iRowFEM<setup.fem_length; ++iRowFEM) {
    const ordinal_type row_size = setup.fem_graph[iRowFEM].size();
    for (ordinal_type iRowEntryFEM=0; iRowEntryFEM<row_size;
         ++iRowEntryFEM, ++iEntryFEM) {
      const ordinal_type iColFEM = setup.fem_graph[iRowFEM][iRowEntryFEM];

      for (ordinal_type k=0; k<setup.stoch_length; ++k) {
        haM(iEntryFEM,k) =
          generate_matrix_coefficient<ScalarType>(
            setup.fem_length, setup.stoch_length, iRowFEM, iColFEM, k);
      }
    }
  }

  Kokkos::deep_copy( matrix.values, hM );

  //------------------------------

  Stokhos::multiply( matrix, x, y, tag );

  Kokkos::deep_copy( hy, y );

  bool success = setup.test_commuted(hy, out);
  return success;
}



}

#endif
