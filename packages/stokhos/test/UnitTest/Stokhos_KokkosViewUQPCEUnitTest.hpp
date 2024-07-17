// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestHelpers.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

//
// Tests various View< Sacado::UQ::PCE<...>,...> operations work
// as expected
//

// Helper functions

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

template <typename scalar, typename ordinal>
inline
scalar generate_pce_coefficient( const ordinal nFEM,
                                    const ordinal nStoch,
                                    const ordinal iColFEM,
                                    const ordinal iStoch )
{
  const scalar X_fem = 100.0 + scalar(iColFEM) / scalar(nFEM);
  const scalar X_stoch =  1.0 + scalar(iStoch) / scalar(nStoch);
  return X_fem + X_stoch;
  //return 1.0;
}

template <typename ViewType>
bool
checkPCEView(const ViewType& v, Teuchos::FancyOStream& out) {
  typedef ViewType view_type;
  typedef typename view_type::size_type size_type;
  typedef typename view_type::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename host_array_type::value_type scalar_type;

  // Copy to host
  host_view_type h_v = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(h_v, v);
  host_array_type h_a = h_v;

  size_type num_rows, num_cols;

  // For layout left, sacado dimension becomes first dimension
  // instead of last
  bool is_right = std::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right) {
    num_rows = h_a.extent(0);
    num_cols = h_a.extent(1);
  }
  else {
    num_rows = h_a.extent(1);
    num_cols = h_a.extent(0);
  }
  bool success = true;
  if (is_right) {
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        scalar_type val = h_a(i,j);
        scalar_type val_expected =
          generate_pce_coefficient<scalar_type>(
            num_rows, num_cols, i, j);
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }
  else {
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        scalar_type val = h_a(j,i);
        scalar_type val_expected =
          generate_pce_coefficient<scalar_type>(
            num_rows, num_cols, i, j);
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }

  return success;
}

template <typename ViewType>
bool
checkConstantPCEView(const ViewType& v,
                     const typename ViewType::value_type& v_expected,
                     Teuchos::FancyOStream& out) {
  typedef ViewType view_type;
  typedef typename view_type::size_type size_type;
  typedef typename view_type::HostMirror host_view_type;
  typedef typename Kokkos::IntrinsicScalarType<host_view_type>::type scalar_type;

  // Copy to host
  host_view_type h_v = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(h_v, v);

  const size_type num_rows = h_v.extent(0);
  const size_type num_cols = Kokkos::dimension_scalar(h_v);

  bool success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      scalar_type val = h_v(i).coeff(j);
      scalar_type val_expected = v_expected.coeff(j);
      TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
    }
  }

  return success;
}

template <typename DataType, typename LayoutType, typename ExecutionSpace>
struct ApplyView {
  typedef Kokkos::View<DataType,LayoutType,ExecutionSpace> type;
};

struct NoLayout {};
template <typename DataType, typename ExecutionSpace>
struct ApplyView<DataType,NoLayout,ExecutionSpace> {
  typedef Kokkos::View<DataType,ExecutionSpace> type;
};

//
// Tests
//

const int global_num_rows = 11;
const int global_num_cols = 9;

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, Size, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  //typedef size_t size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  TEUCHOS_TEST_EQUALITY(v.size(), num_rows, out, success);
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;

  bool is_right = std::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_v);

  success = checkPCEView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_NonContiguous, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef Kokkos::HostSpace HostDevice;
  typedef Kokkos::View<PCE*,typename ViewType::array_layout,HostDevice,Kokkos::MemoryUnmanaged> HostViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);

  Teuchos::Array<PCE> a(num_rows);
  for (size_type i=0; i<num_rows; ++i) {
    a[i].reset(cijk);
    for (size_type j=0; j<num_cols; ++j)
      a[i].fastAccessCoeff(j) = generate_pce_coefficient<Scalar>(
        num_rows, num_cols, i, j);
  }
  // HostViewType ha(a.getRawPtr(), cijk, num_rows, num_cols);
  HostViewType ha =
    Kokkos::make_view<HostViewType>(a.getRawPtr(), cijk, num_rows, num_cols);
  Kokkos::deep_copy(v, ha);

  success = checkPCEView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_ConstantScalar, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  Scalar val = 1.2345;

  Kokkos::deep_copy( v, val );

  PCE pce_val(cijk); pce_val.fastAccessCoeff(0) = val;
  success = checkConstantPCEView(v, pce_val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_ConstantPCE, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  Scalar val = 1.2345;

  Kokkos::deep_copy( v, PCE(val) );

  PCE pce_val(cijk); pce_val.fastAccessCoeff(0) = val;
  success = checkConstantPCEView(v, pce_val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_ConstantPCE2, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  PCE val(cijk);
  for (size_type j=0; j<num_cols; ++j)
    val.fastAccessCoeff(j) =
      generate_pce_coefficient<Scalar>(num_rows, num_cols, size_type(0), j);

  Kokkos::deep_copy( v, val );

  success = checkConstantPCEView(v, val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Kokkos_View_PCE, DeepCopy_Subview_Range, Storage )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE**,Kokkos::LayoutLeft,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows1 = global_num_rows;
  const size_type num_rows2 = global_num_rows*2;
  const size_type num_cols = 5;
  const size_type num_pce = cijk.dimension();
  ViewType v1 = Kokkos::make_view<ViewType>("view1", cijk, num_rows1, num_cols);
  ViewType v2 = Kokkos::make_view<ViewType>("view2", cijk, num_rows2, num_cols);

  for (size_type j=0; j<num_cols; ++j) {
    std::pair<size_type,size_type> rows( 0, num_rows1 );
    ViewType v1s = Kokkos::subview( v1, rows, std::pair<size_t,size_t> (j,j+1) );
    ViewType v2s = Kokkos::subview( v2, rows, std::pair<size_t,size_t> (j,j+1) );
    Kokkos::deep_copy( v1s, Scalar(j+1) );
    Kokkos::deep_copy( v2s, v1s );
  }

  // Check
  success = true;
  host_view_type hv2 = Kokkos::create_mirror_view( v2 );
  Kokkos::deep_copy( hv2, v2 );
  for (size_type j=0; j<num_cols; ++j) {
    for (size_type i=0; i<num_rows1; ++i) {
      for (size_type k=0; k<num_pce; ++k) {
        Scalar val = hv2(i,j).fastAccessCoeff(k);
        Scalar val_expected = k == 0 ? Scalar(j+1) : Scalar(0);
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
    for (size_type i=num_rows1; i<num_rows2; ++i) {
      for (size_type k=0; k<num_pce; ++k) {
        Scalar val = hv2(i,j).fastAccessCoeff(k);
        Scalar val_expected = 0;
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_HostArray, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  host_array_type h_a = Kokkos::create_mirror_view(v);

  bool is_right = std::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_a);

  success = checkPCEView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_DeviceArray, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename ViewType::array_type array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;
  array_type a = v;

  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_a(i,j) = generate_pce_coefficient<Scalar>(
        num_rows, num_cols, i, j);
  Kokkos::deep_copy(a, h_v);

  success = checkPCEView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, Unmanaged, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;

  bool is_right = std::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_v);

  // Create unmanaged view
  ViewType v2 =
    Kokkos::make_view<ViewType>( v.data(), cijk, num_rows, num_cols );

  success = checkPCEView(v2, out);
}


namespace Test {

template< class ViewType >
struct PCEAtomicFunctor {
  typedef typename ViewType::execution_space execution_space ;

  typedef typename ViewType::value_type pce_type ;
  typedef typename pce_type::storage_type::value_type scalar_type ;

  scalar_type m_s ;
  ViewType m_v ;

  PCEAtomicFunctor( const ViewType & v , const scalar_type & s ) :
    m_v( v ), m_s( s )
  {
    Kokkos::parallel_for( m_v.extent(0) , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( int i ) const
  {
    pce_type v( m_s );
    atomic_assign( & m_v(i) , v );
  }
};

}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeviceAtomic, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v = Kokkos::make_view<ViewType>("view", cijk, num_rows, num_cols);
  Scalar val = 1.2345;

  (void) Test::PCEAtomicFunctor<ViewType>( v , val );

  success = checkConstantPCEView(v, val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, AssignData, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  //typedef size_t size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v1 = Kokkos::make_view<ViewType>("view1", cijk, num_rows, num_cols);
  ViewType v2 = Kokkos::make_view<ViewType>("view2", cijk, num_rows, num_cols);
  Scalar val1 = 1.234;
  Scalar val2 = 5.678;
  Kokkos::deep_copy(v1, val1);
  Kokkos::deep_copy(v2, val2);

  auto s1 = Kokkos::subview(v1, std::make_pair(0, 1));
  auto s2 = Kokkos::subview(v2, std::make_pair(0, 1));

  s1.assign_data(s2.data());

  success = checkConstantPCEView(s1, val2, out);
}

/*
*/

#define VIEW_UQ_PCE_TESTS_STORAGE_LAYOUT( STORAGE, LAYOUT )             \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, Size, STORAGE, LAYOUT )                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy, STORAGE, LAYOUT )                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_NonContiguous, STORAGE, LAYOUT )          \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_ConstantScalar, STORAGE, LAYOUT )         \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_ConstantPCE, STORAGE, LAYOUT )            \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_ConstantPCE2, STORAGE, LAYOUT )           \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, Unmanaged, STORAGE, LAYOUT )                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, AssignData, STORAGE, LAYOUT )

// Some tests the fail, or fail to compile

  /*
    // These don't compile as there are no deep_copy overloads between
    // static view spec and its array_type.  That could be done, but
    // would require some additional work to ensure the shapes match.
    // It is simple enough to create an array_type view, so deep copying
    // between matching views isn't much more trouble.
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_HostArray, STORAGE, LAYOUT )               \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                       \
    Kokkos_View_PCE, DeepCopy_DeviceArray, STORAGE, LAYOUT )
  */

#define VIEW_UQ_PCE_TESTS_STORAGE( STORAGE )                            \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_UQ_PCE_TESTS_STORAGE_LAYOUT(STORAGE, NoLayout)                   \
  VIEW_UQ_PCE_TESTS_STORAGE_LAYOUT(STORAGE, LayoutLeft)                 \
  VIEW_UQ_PCE_TESTS_STORAGE_LAYOUT(STORAGE, LayoutRight)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_Subview_Range, STORAGE )

#define VIEW_UQ_PCE_TESTS_ORDINAL_SCALAR_DEVICE( ORDINAL, SCALAR, DEVICE ) \
  typedef Stokhos::DynamicStorage<ORDINAL,SCALAR,DEVICE> DS;            \
  VIEW_UQ_PCE_TESTS_STORAGE( DS )

#define VIEW_UQ_PCE_TESTS_DEVICE( DEVICE )                              \
  VIEW_UQ_PCE_TESTS_ORDINAL_SCALAR_DEVICE( int, double, DEVICE )
