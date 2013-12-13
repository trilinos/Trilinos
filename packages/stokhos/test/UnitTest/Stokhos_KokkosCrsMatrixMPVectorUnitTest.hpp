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

#include "Teuchos_UnitTestHelpers.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado_Kokkos.hpp"
#include "Kokkos_CrsMatrix_MP_Vector.hpp"

// For computing DeviceConfig
#include "Kokkos_hwloc.hpp"
#include "Kokkos_Cuda.hpp"

// Helper functions
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

// Helper function to build a diagonal matrix
template <typename MatrixType>
MatrixType
buildDiagonalMatrix(typename MatrixType::ordinal_type nrow,
                    typename MatrixType::ordinal_type mp_vector_size) {
  typedef typename MatrixType::ordinal_type ordinal_type;
  typedef typename MatrixType::StaticCrsGraphType matrix_graph_type;
  typedef typename MatrixType::values_type matrix_values_type;

  std::vector< std::vector<ordinal_type> > graph(nrow);
  for (ordinal_type i=0; i<nrow; ++i)
    graph[i] = std::vector<ordinal_type>(1, i);
  ordinal_type graph_length = nrow;

  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>("graph", graph);
  matrix_values_type matrix_values =
    matrix_values_type("values", graph_length, mp_vector_size);

  MatrixType matrix("matrix", nrow, matrix_values, matrix_graph);
  return matrix;
}

//
// Tests
//

// Kernel to set diagonal of a matrix to prescribed values
template <typename MatrixType>
struct ReplaceDiagonalValuesKernel {
  typedef typename MatrixType::device_type device_type;
  typedef typename MatrixType::size_type size_type;
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

  const MatrixType m_matrix;
  ReplaceDiagonalValuesKernel(const MatrixType matrix) : m_matrix(matrix) {};

  // Replace diagonal entry for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    const ordinal_type row = i;
    const ordinal_type col = i;
    scalar_type val = scalar_type(row);
    m_matrix.replaceValues(row, &col, 1, &val, true);
  }

  // Kernel launch
  static void apply(const MatrixType matrix) {
    const size_type nrow = matrix.numRows();
    Kokkos::parallel_for( nrow, ReplaceDiagonalValuesKernel(matrix) );
  }

  // Check the result is as expected
  static bool check(const MatrixType matrix,
                    Teuchos::FancyOStream& out) {
    typedef typename MatrixType::values_type matrix_values_type;
    typename matrix_values_type::HostMirror host_matrix_values =
      Kokkos::create_mirror_view(matrix.values);
    Kokkos::deep_copy(host_matrix_values, matrix.values);
    const ordinal_type nrow = matrix.numRows();
    bool success = true;
    for (ordinal_type row=0; row<nrow; ++row) {
      bool s = compareVecs(host_matrix_values(row),
                           "matrix_values(row)",
                           scalar_type(row),
                           "scalar_type(row)",
                           0.0, 0.0, out);
      success = success && s;
    }
    return success;
  }
};

// Kernel to add values to the diagonal of a matrix
template <typename MatrixType>
struct AddDiagonalValuesKernel {
  typedef typename MatrixType::device_type device_type;
  typedef typename MatrixType::size_type size_type;
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

  const MatrixType m_matrix;
  AddDiagonalValuesKernel(const MatrixType matrix) : m_matrix(matrix) {};

  // Replace diagonal entry for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    const ordinal_type row = i;
    const ordinal_type col = i;
    scalar_type val = scalar_type(row);
    m_matrix.sumIntoValues(row, &col, 1, &val, true);
  }

  // Kernel launch
  static void apply(const MatrixType matrix) {
    const size_type nrow = matrix.numRows();
    Kokkos::parallel_for( nrow, AddDiagonalValuesKernel(matrix) );
  }

  // Check the result is as expected
  static bool check(const MatrixType matrix,
                    Teuchos::FancyOStream& out) {
    typedef typename MatrixType::values_type matrix_values_type;
    typename matrix_values_type::HostMirror host_matrix_values =
      Kokkos::create_mirror_view(matrix.values);
    Kokkos::deep_copy(host_matrix_values, matrix.values);
    const ordinal_type nrow = matrix.numRows();
    bool success = true;
    for (ordinal_type row=0; row<nrow; ++row) {
      bool s = compareVecs(host_matrix_values(row),
                           "matrix_values(row)",
                           scalar_type(row),
                           "scalar_type(row)",
                           0.0, 0.0, out);
      success = success && s;
    }
    return success;
  }
};

// Kernel to add values to the diagonal of a matrix where each thread
// adds to the same row (checks atomic really works)
template <typename MatrixType>
struct AddDiagonalValuesAtomicKernel {
  typedef typename MatrixType::device_type device_type;
  typedef typename MatrixType::size_type size_type;
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

  const MatrixType m_matrix;
  AddDiagonalValuesAtomicKernel(const MatrixType matrix) : m_matrix(matrix) {};

  // Replace diagonal entry for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    const ordinal_type row = 0;
    const ordinal_type col = 0;
    scalar_type val = scalar_type(i);
    m_matrix.sumIntoValues(row, &col, 1, &val, true);
  }

  // Kernel launch
  static void apply(const MatrixType matrix) {
    const size_type nrow = matrix.numRows();
    Kokkos::parallel_for( nrow, AddDiagonalValuesAtomicKernel(matrix) );
  }

  // Check the result is as expected
  static bool check(const MatrixType matrix,
                    Teuchos::FancyOStream& out) {
    typedef typename MatrixType::values_type matrix_values_type;
    typename matrix_values_type::HostMirror host_matrix_values =
      Kokkos::create_mirror_view(matrix.values);
    Kokkos::deep_copy(host_matrix_values, matrix.values);
    const ordinal_type nrow = matrix.numRows();
    bool success = true;
    for (ordinal_type row=0; row<nrow; ++row) {
      scalar_type val;
      if (row == 0)
        val = scalar_type( nrow*(nrow-1)/2 );
      else
        val = scalar_type(0.0);
      bool s = compareVecs(host_matrix_values(row),
                           "matrix_values(row)",
                           val,
                           "val",
                           0.0, 0.0, out);
      success = success && s;
    }
    return success;
  }
};

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, ReplaceValues, Scalar, Ordinal, Device )
{
  const Ordinal VectorSize = 3;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> MatrixScalar;
  typedef Kokkos::CrsMatrix<MatrixScalar,Ordinal,Device> Matrix;

  // Build diagonal matrix
  Ordinal nrow = 10;
  Matrix matrix = buildDiagonalMatrix<Matrix>(nrow, VectorSize);

  // Launch our kernel
  typedef ReplaceDiagonalValuesKernel<Matrix> kernel;
  kernel::apply(matrix);

  // Check the result
  success = kernel::check(matrix, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, SumIntoValues, Scalar, Ordinal, Device )
{
  const Ordinal VectorSize = 3;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> MatrixScalar;
  typedef Kokkos::CrsMatrix<MatrixScalar,Ordinal,Device> Matrix;

  // Build diagonal matrix
  Ordinal nrow = 10;
  Matrix matrix = buildDiagonalMatrix<Matrix>(nrow, VectorSize);

  // Launch our kernel
  typedef AddDiagonalValuesKernel<Matrix> kernel;
  kernel::apply(matrix);

  // Check the result
  success = kernel::check(matrix, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, SumIntoValuesAtomic, Scalar, Ordinal, Device )
{
  const Ordinal VectorSize = 3;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> MatrixScalar;
  typedef Kokkos::CrsMatrix<MatrixScalar,Ordinal,Device> Matrix;

  // Build diagonal matrix
  Ordinal nrow = 10;
  Matrix matrix = buildDiagonalMatrix<Matrix>(nrow, VectorSize);

  // Launch our kernel
  typedef AddDiagonalValuesAtomicKernel<Matrix> kernel;
  kernel::apply(matrix);

  // Check the result
  success = kernel::check(matrix, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_CrsMatrix_MP, SpMv, Scalar, Ordinal, Device )
{
  const Ordinal VectorSize = 3;
  typedef Stokhos::StaticFixedStorage<Ordinal,Scalar,VectorSize,Device> Storage;
  typedef Sacado::MP::Vector<Storage> MatrixScalar;
  typedef Kokkos::CrsMatrix<MatrixScalar,Ordinal,Device> Matrix;
  typedef Kokkos::View<MatrixScalar*, Kokkos::LayoutRight, Device> Vector;

  typedef typename Matrix::StaticCrsGraphType matrix_graph_type;
  typedef typename Matrix::values_type matrix_values_type;
  typedef typename matrix_values_type::HostMirror host_matrix_values_type;
  typedef typename Vector::HostMirror HostVector;

  // Generate FEM graph:
  Ordinal nGrid = 5;
  Ordinal fem_length = nGrid * nGrid * nGrid;
  std::vector< std::vector<Ordinal> > fem_graph;
  Ordinal fem_graph_length = generate_fem_graph( nGrid, fem_graph );
  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>("graph", fem_graph);

  // Create matrix
  matrix_values_type matrix_values =
    matrix_values_type("values", fem_graph_length, VectorSize);
  host_matrix_values_type host_matrix_values =
    Kokkos::create_mirror_view(matrix_values);
  for (Ordinal rowFEM=0, entryFEM=0; rowFEM<fem_length; ++rowFEM) {
    const Ordinal row_size = fem_graph[rowFEM].size();
    for (Ordinal rowEntryFEM=0; rowEntryFEM<row_size; ++rowEntryFEM,++entryFEM){
      const Ordinal colFEM = fem_graph[rowFEM][rowEntryFEM];
      for (Ordinal k=0; k<VectorSize; ++k) {
        host_matrix_values(entryFEM,k) =
          generate_matrix_coefficient<Scalar>(
            fem_length, VectorSize, rowFEM, colFEM, k);
      }
    }
  }
  Kokkos::deep_copy(matrix_values, host_matrix_values);
  Matrix A("matrix", fem_length, matrix_values, matrix_graph);

  // Create x and y vectors
  Vector x("x", fem_length, VectorSize);
  Vector y("y", fem_length, VectorSize);
  HostVector host_x = Kokkos::create_mirror_view(x);
  for (Ordinal i=0; i<fem_length; ++i)
    for (Ordinal j=0; j<VectorSize; ++j)
      host_x(i,j) = generate_vector_coefficient<Scalar>(
        fem_length, VectorSize, i, j);
  Kokkos::deep_copy(x, host_x);

  // Setup DeviceConfig
  if (Kokkos::Impl::is_same<Device,Kokkos::Cuda>::value) {
    // For GPU, 1 thread per vector and 16 rows per block
    A.dev_config = Kokkos::DeviceConfig(0, 1, 16);
  }
  else {
    // Use hwloc to get number of cores, hypertheads for CPU/MIC.
    size_t num_cores =
      Kokkos::hwloc::get_available_numa_count() *
      Kokkos::hwloc::get_available_cores_per_numa();
    size_t num_hyper_threads =
      Kokkos::hwloc::get_available_threads_per_core();
    size_t vector_threads = 1;
    size_t row_threads = num_hyper_threads / vector_threads;
    A.dev_config =
      Kokkos::DeviceConfig(num_cores, vector_threads, row_threads);
  }

  // Multiply -- Using our kernel.  Besides not being optimal, the default one
  // doesn't work because the implementation thinks its a 2-D view
  Kokkos::MV_Multiply(y, A, x);

  // Check
  HostVector host_y = Kokkos::create_mirror_view(y);
  Kokkos::deep_copy(host_y, y);
  success = true;
  for (Ordinal rowFEM=0, entryFEM=0; rowFEM<fem_length; ++rowFEM) {
    MatrixScalar y_expected = MatrixScalar(0.0);
    const Ordinal row_size = fem_graph[rowFEM].size();
    for (Ordinal rowEntryFEM=0; rowEntryFEM<row_size; ++rowEntryFEM,++entryFEM){
      const Ordinal colFEM = fem_graph[rowFEM][rowEntryFEM];
      for (Ordinal k=0; k<VectorSize; ++k) {
        y_expected.fastAccessCoeff(k) +=
          generate_matrix_coefficient<Scalar>(
            fem_length, VectorSize, rowFEM, colFEM, k) *
          generate_vector_coefficient<Scalar>(
            fem_length, VectorSize, colFEM, k);
      }
    }
    bool s = compareVecs(host_y(rowFEM), "host_y(rowFEM)", y_expected, "y_expected",
                         0.0, 0.0, out);
    success = success && s;
  }
}

#define CRSMATRIX_MP_VECTOR_TESTS_SCALAR_ORDINAL_DEVICE(SCALAR, ORDINAL, DEVICE)\
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
  Kokkos_CrsMatrix_MP, ReplaceValues, SCALAR, ORDINAL, DEVICE )       \
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
  Kokkos_CrsMatrix_MP, SumIntoValues, SCALAR, ORDINAL, DEVICE )       \
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
  Kokkos_CrsMatrix_MP, SumIntoValuesAtomic, SCALAR, ORDINAL, DEVICE ) \
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(                                 \
  Kokkos_CrsMatrix_MP, SpMv, SCALAR, ORDINAL, DEVICE )
