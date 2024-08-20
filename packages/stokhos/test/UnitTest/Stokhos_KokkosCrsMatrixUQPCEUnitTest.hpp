// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHelpers.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE.hpp"
#include "Kokkos_CrsMatrix_UQ_PCE_Cuda.hpp"
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

// For computing DeviceConfig
#include "Kokkos_Core.hpp"

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

template <typename kokkos_cijk_type, typename ordinal_type>
kokkos_cijk_type build_cijk(ordinal_type stoch_dim,
                            ordinal_type poly_ord)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;

  typedef typename kokkos_cijk_type::value_type value_type;
  typedef typename kokkos_cijk_type::execution_space execution_space;
  typedef Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> one_d_basis;
  typedef Stokhos::LegendreBasis<ordinal_type,value_type> legendre_basis;
  typedef Stokhos::CompletePolynomialBasis<ordinal_type,value_type> product_basis;
  typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk;

  // Create product basis
  Array< RCP<const one_d_basis> > bases(stoch_dim);
  for (ordinal_type i=0; i<stoch_dim; i++)
    bases[i] = rcp(new legendre_basis(poly_ord, true));
  RCP<const product_basis> basis = rcp(new product_basis(bases));

  // Triple product tensor
  RCP<Cijk> cijk = basis->computeTripleProductTensor();

  // Kokkos triple product tensor
  kokkos_cijk_type kokkos_cijk =
    Stokhos::create_product_tensor<execution_space>(*basis, *cijk);

  return kokkos_cijk;
}

// Reasonable tolerances for common precisions
template <typename Scalar> struct ScalarTol {};
template <> struct ScalarTol<float> {  static float  tol() { return 1e-4;  } };
template <> struct ScalarTol<double> { static double tol() { return 1e-10; } };

// Compare two rank-2 views for equality, to given precision
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

  size_type num_rows = y.extent(0);
  size_type num_cols = y.extent(1);
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

template <typename vector_type, typename scalar_type>
bool compareRank1(const vector_type& y,
                          const vector_type& y_exp,
                          const scalar_type rel_tol,
                          const scalar_type abs_tol,
                          Teuchos::FancyOStream& out)
{
  typedef typename vector_type::size_type size_type;
  typename vector_type::HostMirror hy = Kokkos::create_mirror_view(y);
  typename vector_type::HostMirror hy_exp = Kokkos::create_mirror_view(y_exp);
  Kokkos::deep_copy(hy, y);
  Kokkos::deep_copy(hy_exp, y_exp);

  size_type num_rows = y.extent(0);
  bool success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<Kokkos::dimension_scalar(y); ++j) {
      scalar_type diff = std::abs( hy(i).fastAccessCoeff(j) - hy_exp(i).fastAccessCoeff(j) );
      scalar_type tol = rel_tol*std::abs(hy_exp(i).fastAccessCoeff(j)) + abs_tol;
      bool s = diff < tol;
      out << "y_expected(" << i << ").coeff(" << j << ") - "
          << "y(" << i << ").coeff(" << j << ") = " << hy_exp(i).fastAccessCoeff(j)
          << " - " << hy(i).fastAccessCoeff(j) << " == "
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

template <typename vector_type, typename scalar_type>
bool compareRank2(const vector_type& y,
                          const vector_type& y_exp,
                          const scalar_type rel_tol,
                          const scalar_type abs_tol,
                          Teuchos::FancyOStream& out)
{
  typedef typename vector_type::size_type size_type;
  typename vector_type::HostMirror hy = Kokkos::create_mirror_view(y);
  typename vector_type::HostMirror hy_exp = Kokkos::create_mirror_view(y_exp);
  Kokkos::deep_copy(hy, y);
  Kokkos::deep_copy(hy_exp, y_exp);

  size_type num_rows = y.extent(0);
  size_type num_cols = y.extent(1);
  bool success = true;

 for (size_type col = 0; col < num_cols; ++col){
 for (size_type i=0; i<num_rows; ++i) {
   for (size_type j=0; j<Kokkos::dimension_scalar(y); ++j) {
      scalar_type diff = std::abs( hy(i,col).fastAccessCoeff(j) - hy_exp(i,col).fastAccessCoeff(j) );
      scalar_type tol = rel_tol*std::abs(hy_exp(i,col).fastAccessCoeff(j)) + abs_tol;
      bool s = diff < tol;
      out << "y_expected(" << i << ").coeff(" << j << ") - "
          << "y(" << i << ").coeff(" << j << ") = " << hy_exp(i,col).fastAccessCoeff(j)
          << " - " << hy(i,col).fastAccessCoeff(j) << " == "
          << diff << " < " << tol << " : ";
      if (s)
        out << "passed";
      else
        out << "failed";
      out << std::endl;
      success = success && s;
    }
  }
  }


  return success;
}


// Helper function to build a diagonal matrix
template <typename MatrixType, typename CijkType>
MatrixType
buildDiagonalMatrix(typename MatrixType::ordinal_type nrow,
                    typename MatrixType::ordinal_type pce_size,
                    const CijkType& cijk) {
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
    Kokkos::make_view<matrix_values_type>("values", cijk, graph_length, pce_size);

  MatrixType matrix("matrix", nrow, matrix_values, matrix_graph);
  return matrix;
}

//
// Tests
//

// Kernel to set diagonal of a matrix to prescribed values
template <typename MatrixType>
struct ReplaceDiagonalValuesKernel {
  typedef typename MatrixType::execution_space execution_space;
  typedef typename MatrixType::size_type size_type;
  typedef typename MatrixType::value_type value_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

  const MatrixType m_matrix;
  ReplaceDiagonalValuesKernel(const MatrixType matrix) : m_matrix(matrix) {};

  // Replace diagonal entry for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    const ordinal_type row = i;
    const ordinal_type col = i;
    value_type val = value_type(row);
    m_matrix.replaceValues(row, &col, 1, &val, false, true);
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
    value_type val_expected(Kokkos::cijk(matrix.values));
    for (ordinal_type row=0; row<nrow; ++row) {
      val_expected = row;
      bool s = compareVecs(host_matrix_values(row),
                           "matrix_values(row)",
                           val_expected,
                           "val_expected",
                           0.0, 0.0, out);
      success = success && s;
    }
    return success;
  }
};

// Kernel to add values to the diagonal of a matrix
template <typename MatrixType>
struct AddDiagonalValuesKernel {
  typedef typename MatrixType::execution_space execution_space;
  typedef typename MatrixType::size_type size_type;
  typedef typename MatrixType::value_type value_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

  const MatrixType m_matrix;
  AddDiagonalValuesKernel(const MatrixType matrix) : m_matrix(matrix) {};

  // Replace diagonal entry for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    const ordinal_type row = i;
    const ordinal_type col = i;
    value_type val = value_type(row);
    m_matrix.sumIntoValues(row, &col, 1, &val, false, true);
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
    value_type val_expected(Kokkos::cijk(matrix.values));
    for (ordinal_type row=0; row<nrow; ++row) {
      val_expected = row;
      bool s = compareVecs(host_matrix_values(row),
                           "matrix_values(row)",
                           val_expected,
                           "val_expected",
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
  typedef typename MatrixType::execution_space execution_space;
  typedef typename MatrixType::size_type size_type;
  typedef typename MatrixType::value_type value_type;
  typedef typename MatrixType::ordinal_type ordinal_type;

  const MatrixType m_matrix;
  AddDiagonalValuesAtomicKernel(const MatrixType matrix) : m_matrix(matrix) {};

  // Replace diagonal entry for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    const ordinal_type row = 0;
    const ordinal_type col = 0;
    value_type val = value_type(i);
    m_matrix.sumIntoValues(row, &col, 1, &val, false, true);
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
    value_type val_expected(Kokkos::cijk(matrix.values));
    for (ordinal_type row=0; row<nrow; ++row) {
      value_type val;
      if (row == 0)
        val_expected = nrow*(nrow-1)/2;
      else
        val_expected = 0.0;
      bool s = compareVecs(host_matrix_values(row),
                           "matrix_values(row)",
                           val_expected,
                           "val_expected",
                           0.0, 0.0, out);
      success = success && s;
    }
    return success;
  }
};

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(
  Kokkos_CrsMatrix_PCE, ReplaceValues, MatrixScalar )
{
  typedef typename MatrixScalar::ordinal_type Ordinal;
  typedef typename MatrixScalar::execution_space Device;
  typedef typename MatrixScalar::cijk_type Cijk;
  typedef KokkosSparse::CrsMatrix<MatrixScalar,Ordinal,Device> Matrix;

  // Build Cijk tensor
  const Ordinal stoch_dim = 2;
  const Ordinal poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  // Build diagonal matrix
  const Ordinal nrow = 10;
  const Ordinal pce_size = cijk.dimension();
  Matrix matrix = buildDiagonalMatrix<Matrix>(nrow, pce_size, cijk);

  // Launch our kernel
  typedef ReplaceDiagonalValuesKernel<Matrix> kernel;
  kernel::apply(matrix);

  // Check the result
  success = kernel::check(matrix, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(
  Kokkos_CrsMatrix_PCE, SumIntoValues, MatrixScalar )
{
  typedef typename MatrixScalar::ordinal_type Ordinal;
  typedef typename MatrixScalar::execution_space Device;
  typedef typename MatrixScalar::cijk_type Cijk;
  typedef KokkosSparse::CrsMatrix<MatrixScalar,Ordinal,Device> Matrix;

  // Build Cijk tensor
  const Ordinal stoch_dim = 2;
  const Ordinal poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  // Build diagonal matrix
  const Ordinal nrow = 10;
  const Ordinal pce_size = cijk.dimension();
  Matrix matrix = buildDiagonalMatrix<Matrix>(nrow, pce_size, cijk);

  // Launch our kernel
  typedef AddDiagonalValuesKernel<Matrix> kernel;
  kernel::apply(matrix);

  // Check the result
  success = kernel::check(matrix, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(
  Kokkos_CrsMatrix_PCE, SumIntoValuesAtomic, MatrixScalar )
{
  typedef typename MatrixScalar::ordinal_type Ordinal;
  typedef typename MatrixScalar::execution_space Device;
  typedef typename MatrixScalar::cijk_type Cijk;
  typedef KokkosSparse::CrsMatrix<MatrixScalar,Ordinal,Device> Matrix;

  // Build Cijk tensor
  const Ordinal stoch_dim = 2;
  const Ordinal poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  // Build diagonal matrix
  const Ordinal nrow = 10;
  const Ordinal pce_size = cijk.dimension();
  Matrix matrix = buildDiagonalMatrix<Matrix>(nrow, pce_size, cijk);

  // Launch our kernel
  typedef AddDiagonalValuesAtomicKernel<Matrix> kernel;
  kernel::apply(matrix);

  // Check the result
  success = kernel::check(matrix, out);
}

template <typename PCEType, typename Multiply>
bool test_embedded_pce(const typename PCEType::ordinal_type nGrid,
                       const typename PCEType::ordinal_type stoch_dim,
                       const typename PCEType::ordinal_type poly_ord,
                       KokkosSparse::DeviceConfig dev_config,
                       Multiply multiply_op,
                       Teuchos::FancyOStream& out)
{
  typedef typename PCEType::ordinal_type ordinal_type;
  typedef typename PCEType::value_type scalar_type;
  typedef typename PCEType::storage_type storage_type;
  typedef typename PCEType::cijk_type cijk_type;
  typedef typename storage_type::execution_space execution_space;
  typedef Kokkos::LayoutLeft Layout;
  typedef Kokkos::View< PCEType*, Layout, execution_space > block_vector_type;
  typedef KokkosSparse::CrsMatrix< PCEType, ordinal_type, execution_space > block_matrix_type;
  typedef typename block_matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename block_matrix_type::values_type matrix_values_type;

  // Build Cijk tensor
  cijk_type cijk = build_cijk<cijk_type>(stoch_dim, poly_ord);
  const ordinal_type stoch_length = cijk.dimension();
  // const ordinal_type align = 8;
  // const ordinal_type stoch_length_aligned = (stoch_length+align-1) & ~(align-1);
  const ordinal_type stoch_length_aligned = stoch_length;

  // Check pce_length == storage_type::static_size for static storage
  TEUCHOS_TEST_FOR_EXCEPTION(
    storage_type::is_static && storage_type::static_size != stoch_length,
    std::logic_error,
    "Static storage size must equal pce size");

  // Generate FEM graph:
  const ordinal_type fem_length = nGrid * nGrid * nGrid;
  std::vector< std::vector<ordinal_type> > fem_graph;
  const ordinal_type fem_graph_length = generate_fem_graph( nGrid, fem_graph );

  //------------------------------
  // Generate input/output multivectors -- Sacado dimension is always last,
  // regardless of LayoutLeft/Right

  block_vector_type x =
    Kokkos::make_view<block_vector_type>("x", cijk, fem_length, stoch_length_aligned);
  block_vector_type y =
    Kokkos::make_view<block_vector_type>("y", cijk, fem_length, stoch_length_aligned);

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror_view( x );
  typename block_vector_type::HostMirror hy = Kokkos::create_mirror_view( y );

  // View the block vector as an array of the embedded intrinsic type.
  typename block_vector_type::HostMirror::array_type hax = hx ;
  typename block_vector_type::HostMirror::array_type hay = hy ;

  for (ordinal_type iRowFEM=0; iRowFEM<fem_length; ++iRowFEM) {
    for (ordinal_type iRowStoch=0; iRowStoch<stoch_length; ++iRowStoch) {
      hax(iRowStoch, iRowFEM) =
        generate_vector_coefficient<scalar_type>(
          fem_length, stoch_length, iRowFEM, iRowStoch );
      hay(iRowStoch, iRowFEM) = 0.0;
    }
  }

  Kokkos::deep_copy( x, hx );
  Kokkos::deep_copy( y, hy );

  //------------------------------
  // Generate block matrix -- it is always LayoutRight (currently)

  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("test crs graph"), fem_graph);
  matrix_values_type matrix_values =
    Kokkos::make_view<matrix_values_type>(
      Kokkos::ViewAllocateWithoutInitializing("matrix"), cijk, fem_graph_length, stoch_length_aligned);
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
        haM(iEntryFEM, k) =
          generate_matrix_coefficient<scalar_type>(
            fem_length, stoch_length, iRowFEM, iColFEM, k);
      }
    }
  }

  Kokkos::deep_copy( matrix.values, hM );

  //------------------------------
  // multiply

  multiply_op( matrix, x, y );

  //------------------------------
  // generate correct answer

  typedef typename block_vector_type::array_type array_type;
  array_type ay_expected =
    array_type("ay_expected", stoch_length_aligned, fem_length);
  typename array_type::HostMirror hay_expected =
    Kokkos::create_mirror_view(ay_expected);
  typename cijk_type::HostMirror host_cijk =
    Kokkos::create_mirror_view(cijk);
  Kokkos::deep_copy(host_cijk, cijk);
  for (ordinal_type iRowFEM=0, iEntryFEM=0; iRowFEM<fem_length; ++iRowFEM) {
    const ordinal_type row_size = fem_graph[iRowFEM].size();
    for (ordinal_type iRowEntryFEM=0; iRowEntryFEM<row_size;
         ++iRowEntryFEM, ++iEntryFEM) {
      const ordinal_type iColFEM = fem_graph[iRowFEM][iRowEntryFEM];
      for (ordinal_type i=0; i<stoch_length; ++i) {
        const ordinal_type num_entry = host_cijk.num_entry(i);
        const ordinal_type entry_beg = host_cijk.entry_begin(i);
        const ordinal_type entry_end = entry_beg + num_entry;
        scalar_type tmp = 0;
        for (ordinal_type entry = entry_beg; entry < entry_end; ++entry) {
          const ordinal_type j = host_cijk.coord(entry,0);
          const ordinal_type k = host_cijk.coord(entry,1);
          const scalar_type a_j =
            generate_matrix_coefficient<scalar_type>(
              fem_length, stoch_length, iRowFEM, iColFEM, j);
          const scalar_type a_k =
            generate_matrix_coefficient<scalar_type>(
              fem_length, stoch_length, iRowFEM, iColFEM, k);
          const scalar_type x_j =
            generate_vector_coefficient<scalar_type>(
              fem_length, stoch_length, iColFEM, j);
          const scalar_type x_k =
            generate_vector_coefficient<scalar_type>(
              fem_length, stoch_length, iColFEM, k);
          tmp += host_cijk.value(entry) * ( a_j * x_k + a_k * x_j );
        }
        hay_expected(i, iRowFEM) += tmp;
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

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(
  Kokkos_CrsMatrix_PCE, Multiply, Scalar, MultiplyOp )
{
  typedef typename Scalar::ordinal_type Ordinal;

  const Ordinal nGrid = 5;
  const Ordinal stoch_dim = 2;
  const Ordinal poly_ord = 3;
  KokkosSparse::DeviceConfig dev_config;

  success = test_embedded_pce<Scalar>(
    nGrid, stoch_dim, poly_ord, dev_config, MultiplyOp(), out);
}

struct Kokkos_MV_Multiply_Op {
  template <typename Matrix, typename InputVector, typename OutputVector>
  void operator() (const Matrix& A,
                   const InputVector& x,
                   OutputVector& y) const {
    KokkosSparse::spmv("N", typename Matrix::value_type(1.0) , A, x, typename Matrix::value_type(0.0), y);
  }
};

template <typename Tag>
struct Stokhos_MV_Multiply_Op {
  Tag tag;
  Stokhos_MV_Multiply_Op(const Tag& tg = Tag()) : tag(tg) {}

  template <typename Matrix, typename InputVector, typename OutputVector>
  void operator() (const Matrix& A,
                   const InputVector& x,
                   OutputVector& y) const {
    Stokhos::multiply(A, x, y, tag);
  }
};

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(
  Kokkos_CrsMatrix_PCE, MeanMultiplyRank1, Scalar )
{
  typedef typename Scalar::ordinal_type Ordinal;

  const Ordinal nGrid = 5;
  const Ordinal stoch_dim = 5;
  const Ordinal poly_ord = 3;
  KokkosSparse::DeviceConfig dev_config;

  typedef typename Scalar::ordinal_type ordinal_type;
  typedef typename Scalar::value_type scalar_type;
  typedef typename Scalar::storage_type storage_type;
  typedef typename Scalar::cijk_type cijk_type;
  typedef typename storage_type::execution_space execution_space;
  typedef Kokkos::LayoutLeft Layout;
  typedef Kokkos::View< Scalar*, Layout, execution_space > block_vector_type;
  typedef KokkosSparse::CrsMatrix< Scalar, ordinal_type, execution_space > block_matrix_type;
  typedef typename block_matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename block_matrix_type::values_type matrix_values_type;

  // Build Cijk tensor
  cijk_type cijk = build_cijk<cijk_type>(stoch_dim, poly_ord);
  cijk_type mean_cijk =
    Stokhos::create_mean_based_product_tensor<execution_space, ordinal_type, scalar_type>();
  const ordinal_type stoch_length = cijk.dimension();
  const ordinal_type align = 8;
  const ordinal_type stoch_length_aligned = (stoch_length+align-1) & ~(align-1);

  // Generate FEM graph:
  const ordinal_type fem_length = nGrid * nGrid * nGrid;
  std::vector< std::vector<ordinal_type> > fem_graph;
  const ordinal_type fem_graph_length = generate_fem_graph( nGrid, fem_graph );

  block_vector_type x =
    Kokkos::make_view<block_vector_type>("x", cijk, fem_length, stoch_length_aligned);
  block_vector_type y =
    Kokkos::make_view<block_vector_type>("y", cijk, fem_length, stoch_length_aligned);

  block_vector_type y_expected =
    Kokkos::make_view<block_vector_type>("y", cijk, fem_length, stoch_length_aligned);

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror_view( x );
  typename block_vector_type::HostMirror hy = Kokkos::create_mirror_view( y );
  typename block_vector_type::HostMirror hy_expected =
    Kokkos::create_mirror_view( y_expected );

  // View the block vector as an array of the embedded intrinsic type.
  typename block_vector_type::HostMirror::array_type hax = hx ;
  typename block_vector_type::HostMirror::array_type hay = hy ;
  typename block_vector_type::HostMirror::array_type hay_expected =
    hy_expected ;

  for (ordinal_type iRowFEM=0; iRowFEM<fem_length; ++iRowFEM) {
    for (ordinal_type iRowStoch=0; iRowStoch<stoch_length; ++iRowStoch) {
      hax(iRowStoch,iRowFEM) =
        generate_vector_coefficient<scalar_type>(
          fem_length, stoch_length, iRowFEM, iRowStoch );
      hay(iRowStoch,iRowFEM) = 0.0;
      hay_expected(iRowStoch,iRowFEM) = 0.0;
    }
  }
  Kokkos::deep_copy( x, hx );
  Kokkos::deep_copy( y, hy );
  Kokkos::deep_copy( y_expected, hy_expected );

  //------------------------------
  // Generate block matrix -- it is always LayoutRight (currently)
  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("test crs graph"), fem_graph);
  matrix_values_type matrix_values =
    Kokkos::make_view<matrix_values_type>(
      Kokkos::ViewAllocateWithoutInitializing("matrix"), mean_cijk, fem_graph_length, ordinal_type(1)); //instead of stoch_length
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

      haM(iEntryFEM, 0) =
        generate_matrix_coefficient<scalar_type>(
          fem_length, 1, iRowFEM, iColFEM, 0);
      for (ordinal_type iRowStoch=0; iRowStoch<stoch_length; ++iRowStoch) {
        hay_expected(iRowStoch,iRowFEM) +=
          haM(iEntryFEM, 0) * hax(iRowStoch,iColFEM);
      }
    }
  }
  Kokkos::deep_copy( matrix.values, hM );
  Kokkos::deep_copy( y_expected, hy_expected );

  /*
  //Generate same matrix with stochastic dim = Kokkos::dimension_scalar(x) (i.e. not = 1)
  matrix_values_type full_matrix_values =
    Kokkos::make_view<matrix_values_type>(
      Kokkos::ViewAllocateWithoutInitializing("matrix"), cijk, fem_graph_length, stoch_length_aligned);
  block_matrix_type full_matrix(
    "block_matrix", fem_length, full_matrix_values, matrix_graph);
  matrix.dev_config = dev_config;

  typename matrix_values_type::HostMirror full_hM =
    Kokkos::create_mirror_view( full_matrix.values );
  typename matrix_values_type::HostMirror::array_type full_haM = full_hM ;

  for (ordinal_type iRowFEM=0, iEntryFEM=0; iRowFEM<fem_length; ++iRowFEM) {
    const ordinal_type row_size = fem_graph[iRowFEM].size();
    for (ordinal_type iRowEntryFEM=0; iRowEntryFEM<row_size;
         ++iRowEntryFEM, ++iEntryFEM) {
      const ordinal_type iColFEM = fem_graph[iRowFEM][iRowEntryFEM];

      for (ordinal_type k=0; k<stoch_length; ++k) {
        if (k == 0)
          full_haM(iEntryFEM, k) =
            generate_matrix_coefficient<scalar_type>(
              fem_length, 1, iRowFEM, iColFEM, k);
        else
          full_haM(iEntryFEM, k) = 0.0;
      }
    }
  }

  Kokkos::deep_copy( full_matrix.values, full_hM );
  */

  //------------------------------
  // multiply

  KokkosSparse::spmv("N", Scalar(1.0) , matrix, x, Scalar(0.0), y);

  //------------------------------
  // multiply with same matrix but with sacado_size = x.sacado_size

  //Kokkos::MV_Multiply( y_expected, full_matrix, x );

  //------------------------------
  // check
  scalar_type rel_tol = ScalarTol<scalar_type>::tol();
  scalar_type abs_tol = ScalarTol<scalar_type>::tol();
  success = compareRank1(y, y_expected, rel_tol, abs_tol, out);
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(
  Kokkos_CrsMatrix_PCE, MeanMultiplyRank2, Scalar )
{
  typedef typename Scalar::ordinal_type ordinal_type;
  typedef typename Scalar::value_type scalar_type;
  typedef typename Scalar::storage_type storage_type;
  typedef typename Scalar::cijk_type cijk_type;
  typedef typename storage_type::execution_space execution_space;
  typedef Kokkos::LayoutLeft Layout;
  typedef Kokkos::View< Scalar**, Layout, execution_space > block_vector_type;
  typedef KokkosSparse::CrsMatrix< Scalar, ordinal_type, execution_space > block_matrix_type;
  typedef typename block_matrix_type::StaticCrsGraphType matrix_graph_type;
  typedef typename block_matrix_type::values_type matrix_values_type;


  const ordinal_type nGrid = 5;
  const ordinal_type stoch_dim = 2;
  const ordinal_type poly_ord = 3;
  KokkosSparse::DeviceConfig dev_config;

  // Build Cijk tensor
  cijk_type cijk = build_cijk<cijk_type>(stoch_dim, poly_ord);
  cijk_type mean_cijk =
    Stokhos::create_mean_based_product_tensor<execution_space, ordinal_type, scalar_type>();
  const ordinal_type stoch_length = cijk.dimension();
  const ordinal_type align = 8;
  const ordinal_type stoch_length_aligned = (stoch_length+align-1) & ~(align-1);
  const ordinal_type num_cols = 2;
  // Generate FEM graph:
  const ordinal_type fem_length = nGrid * nGrid * nGrid;
  std::vector< std::vector<ordinal_type> > fem_graph;
  const ordinal_type fem_graph_length = generate_fem_graph( nGrid, fem_graph );

  block_vector_type x =
    Kokkos::make_view<block_vector_type>("x", cijk, fem_length, num_cols, stoch_length_aligned);
  block_vector_type y =
    Kokkos::make_view<block_vector_type>("y", cijk, fem_length, num_cols, stoch_length_aligned);

  block_vector_type y_expected =
    Kokkos::make_view<block_vector_type>("y_expected", cijk, fem_length, num_cols, stoch_length_aligned);

  typename block_vector_type::HostMirror hx = Kokkos::create_mirror_view( x );
  typename block_vector_type::HostMirror hy = Kokkos::create_mirror_view( y );
  typename block_vector_type::HostMirror hy_expected =
    Kokkos::create_mirror_view( y_expected );

  for (ordinal_type i=0; i<num_cols; ++i){
    for (ordinal_type iRowFEM=0; iRowFEM<fem_length; ++iRowFEM) {
      for (ordinal_type iRowStoch=0; iRowStoch<stoch_length; ++iRowStoch) {
        hx(iRowFEM,i).fastAccessCoeff(iRowStoch) =
          generate_vector_coefficient<scalar_type>(
            fem_length, stoch_length, iRowFEM, iRowStoch );
        hy(iRowFEM,i).fastAccessCoeff(iRowStoch) = 0.0;
        hy_expected(iRowFEM,i).fastAccessCoeff(iRowStoch) = 0.0;
       }
    }
  }
  Kokkos::deep_copy( x, hx );
  Kokkos::deep_copy( y, hy );
  Kokkos::deep_copy( y_expected, hy_expected );

 //------------------------------
  // Generate matrix with stochastic dimension 1
  matrix_graph_type matrix_graph =
    Kokkos::create_staticcrsgraph<matrix_graph_type>(
      std::string("test crs graph"), fem_graph);
  matrix_values_type matrix_values =
    Kokkos::make_view<matrix_values_type>(
      Kokkos::ViewAllocateWithoutInitializing("matrix"), mean_cijk, fem_graph_length, ordinal_type(1));
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

      haM(iEntryFEM, 0) =
        generate_matrix_coefficient<scalar_type>(
          fem_length, 1, iRowFEM, iColFEM, 0);
      for (ordinal_type i=0; i<num_cols; ++i){
        for (ordinal_type iRowStoch=0; iRowStoch<stoch_length; ++iRowStoch) {
          hy_expected(iRowFEM,i).fastAccessCoeff(iRowStoch) +=
            haM(iEntryFEM, 0) * hx(iColFEM,i).fastAccessCoeff(iRowStoch);
        }
      }
    }
  }

  Kokkos::deep_copy( matrix.values, hM );
  Kokkos::deep_copy( y_expected, hy_expected );

  /*
  //Generate same matrix with stochastic dim = Kokkos::dimension_scalar(x) (i.e. not = 1)
  matrix_values_type full_matrix_values =
    Kokkos::make_view<matrix_values_type>(
      Kokkos::ViewAllocateWithoutInitializing("matrix"), cijk, fem_graph_length, stoch_length_aligned);
  block_matrix_type full_matrix(
    "block_matrix", fem_length, full_matrix_values, matrix_graph);
  matrix.dev_config = dev_config;

  typename matrix_values_type::HostMirror full_hM =
    Kokkos::create_mirror_view( full_matrix.values );

  typename matrix_values_type::HostMirror::array_type full_haM = full_hM ;

  for (ordinal_type iRowFEM=0, iEntryFEM=0; iRowFEM<fem_length; ++iRowFEM) {
    const ordinal_type row_size = fem_graph[iRowFEM].size();
    for (ordinal_type iRowEntryFEM=0; iRowEntryFEM<row_size;
         ++iRowEntryFEM, ++iEntryFEM) {
      const ordinal_type iColFEM = fem_graph[iRowFEM][iRowEntryFEM];

      for (ordinal_type k=0; k<stoch_length; ++k) {
        if (k == 0)
          full_haM(iEntryFEM, k) =
          generate_matrix_coefficient<scalar_type>(
            fem_length, 1, iRowFEM, iColFEM, k);
        else
          full_haM(iEntryFEM, k) = 0.0;
      }
    }
  }

  Kokkos::deep_copy( full_matrix.values, full_hM );
  */

  //------------------------------
  // multiply

  KokkosSparse::spmv("N", Scalar(1.0) , matrix, x, Scalar(0.0), y);

  //------------------------------
  // multiply with full matrix

  //Kokkos::MV_Multiply( y_expected, full_matrix, x );

  //------------------------------
  // check

  scalar_type rel_tol = ScalarTol<scalar_type>::tol();
  scalar_type abs_tol = ScalarTol<scalar_type>::tol();
  success = compareRank2(y, y_expected, rel_tol, abs_tol, out);
}


typedef Kokkos_MV_Multiply_Op KokkosMultiply;

#define CRSMATRIX_UQ_PCE_TESTS_MATRIXSCALAR( SCALAR )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                 \
    Kokkos_CrsMatrix_PCE, ReplaceValues, SCALAR )                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                 \
    Kokkos_CrsMatrix_PCE, SumIntoValues, SCALAR )                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                 \
    Kokkos_CrsMatrix_PCE, SumIntoValuesAtomic, SCALAR )                 
#define CRSMATRIX_UQ_PCE_MEAN_MULTIPLY_TESTS( SCALAR )                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                 \
    Kokkos_CrsMatrix_PCE, MeanMultiplyRank1, SCALAR )                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                 \
    Kokkos_CrsMatrix_PCE, MeanMultiplyRank2, SCALAR )                        
#define CRS_MATRIX_UQ_PCE_MULTIPLY_TESTS_SCALAR_OP( SCALAR, OP )   \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                              \
    Kokkos_CrsMatrix_PCE, Multiply, SCALAR, OP )

#define CRS_MATRIX_UQ_PCE_MULTIPLY_TESTS_SCALAR( SCALAR ) \
  CRS_MATRIX_UQ_PCE_MULTIPLY_TESTS_SCALAR_OP( SCALAR, KokkosMultiply )
  
#define CRSMATRIX_UQ_PCE_TESTS_STORAGE( STORAGE )                       \
  typedef Sacado::UQ::PCE<STORAGE> UQ_PCE_ ## STORAGE;                  \
  CRSMATRIX_UQ_PCE_TESTS_MATRIXSCALAR( UQ_PCE_ ## STORAGE )             \
  CRS_MATRIX_UQ_PCE_MULTIPLY_TESTS_SCALAR( UQ_PCE_ ## STORAGE )         \
  CRSMATRIX_UQ_PCE_MEAN_MULTIPLY_TESTS( UQ_PCE_ ## STORAGE )            
   
#define CRSMATRIX_UQ_PCE_TESTS_ORDINAL_SCALAR_DEVICE( ORDINAL, SCALAR, DEVICE ) \
  typedef Stokhos::DynamicStorage<ORDINAL,SCALAR,DEVICE> DS;            \
  CRSMATRIX_UQ_PCE_TESTS_STORAGE( DS )

#define CRSMATRIX_UQ_PCE_TESTS_DEVICE( DEVICE )                         \
  CRSMATRIX_UQ_PCE_TESTS_ORDINAL_SCALAR_DEVICE( int, double, DEVICE )
