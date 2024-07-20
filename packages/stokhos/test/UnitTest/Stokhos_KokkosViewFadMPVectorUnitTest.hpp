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

#include "Kokkos_View_Fad_Fwd.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Stokhos_Ensemble_Sizes.hpp"
#include "Sacado.hpp"
#include "Sacado_Fad_MP_Vector.hpp"

// For computing DeviceConfig
#include "Kokkos_Core.hpp"

//
// Tests various View< Fad< Sacado::MP::Vector<...> >,...> operations work
// as expected
//

// Helper functions

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

template <typename ViewType>
bool
checkVectorView(const ViewType& v,
                Teuchos::FancyOStream& out) {
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

  // For static, layout left, sacado dimension becomes first dimension
  // instead of last
  bool is_right = std::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right || !view_type::is_contiguous) {
    num_rows = h_a.extent(0);
    num_cols = h_a.extent(1);
  }
  else {
    num_rows = h_a.extent(1);
    num_cols = h_a.extent(0);
  }
  bool success = true;
  if (is_right || !view_type::is_contiguous) {
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        scalar_type val = h_a(i,j);
        scalar_type val_expected =
          generate_vector_coefficient<scalar_type>(
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
          generate_vector_coefficient<scalar_type>(
            num_rows, num_cols, i, j);
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }

  return success;
}

template <typename ViewType>
bool
checkConstantFadVectorView(const ViewType& view,
                           const typename ViewType::value_type& v,
                           Teuchos::FancyOStream& out) {
  typedef ViewType view_type;
  typedef typename view_type::value_type fad_vector_type;
  typedef typename fad_vector_type::value_type vector_type;
  typedef typename vector_type::storage_type storage_type;
  typedef typename view_type::size_type size_type;
  typedef typename view_type::HostMirror host_view_type;

  // Copy to host
  host_view_type h_view = Kokkos::create_mirror_view(view);
  Kokkos::deep_copy(h_view, view);

  const size_type num_rows = h_view.extent(0);
  const size_type num_fad = Kokkos::dimension_scalar(h_view)-1;
  const size_type num_ensemble = storage_type::static_size;
  bool success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type k=0; k<num_ensemble; ++k) {
      TEUCHOS_TEST_EQUALITY(h_view(i).val().coeff(k), v.val().coeff(k), out, success);
      for (size_type j=0; j<num_fad; ++j) {
        TEUCHOS_TEST_EQUALITY(h_view(i).dx(j).coeff(k), v.dx(j).coeff(k), out, success);
      }
    }
  }

  return success;
}

template <typename ViewType>
bool
checkConstantFadVectorView2(const ViewType& view,
                            const typename ViewType::value_type& v,
                            Teuchos::FancyOStream& out) {
  typedef ViewType view_type;
  typedef typename view_type::value_type fad_vector_type;
  typedef typename fad_vector_type::value_type vector_type;
  typedef typename vector_type::storage_type storage_type;
  typedef typename view_type::size_type size_type;
  typedef typename view_type::HostMirror host_view_type;

  // Copy to host
  host_view_type h_view = Kokkos::create_mirror_view(view);
  Kokkos::deep_copy(h_view, view);

  bool success = true;
  const size_type num_fad = Kokkos::dimension_scalar(h_view)-1;
  const size_type num_ensemble = storage_type::static_size;
  for (size_type i0=0; i0<h_view.extent(0); ++i0) {
  for (size_type i1=0; i1<h_view.extent(1); ++i1) {
  for (size_type i2=0; i2<h_view.extent(2); ++i2) {
  for (size_type i3=0; i3<h_view.extent(3); ++i3) {
  for (size_type i4=0; i4<h_view.extent(4); ++i4) {
  for (size_type i5=0; i5<h_view.extent(5); ++i5) {
  for (size_type i6=0; i6<h_view.extent(6); ++i6) {
    for (size_type k=0; k<num_ensemble; ++k)
      TEUCHOS_TEST_EQUALITY(h_view.access(i0,i1,i2,i3,i4,i5,i6,0).val().coeff(k),
                            v.val().coeff(k), out, success);
    for (size_type j=0; j<num_fad; ++j) {
      for (size_type k=0; k<num_ensemble; ++k)
        TEUCHOS_TEST_EQUALITY(h_view.access(i0,i1,i2,i3,i4,i5,i6,0).dx(j).coeff(k),
                              v.dx(j).coeff(k), out, success);
    }
  }}}}}}}

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
const int global_ensemble_size = STOKHOS_DEFAULT_ENSEMBLE_SIZE;
const int global_fad_size = 5;

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_Fad_MP, Size, Scalar, Layout )
{
  typedef typename Scalar::value_type Vector;
  typedef typename Vector::execution_space Device;
  typedef typename ApplyView<Scalar*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_fad_size+1;
  ViewType v("view", num_rows, num_cols);
  TEUCHOS_TEST_EQUALITY(v.size(), num_rows, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_Fad_MP, DeepCopy_ConstantScalar, Scalar, Layout )
{
  typedef typename Scalar::value_type Vector;
  typedef typename Vector::value_type BaseScalar;
  typedef typename Vector::execution_space Device;
  typedef typename ApplyView<Scalar*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_fad_size+1;
  ViewType v("view", num_rows, num_cols);
  BaseScalar val = 1.2345;

  Kokkos::deep_copy( v, val );

  success = checkConstantFadVectorView(v, Scalar(val), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_Fad_MP, Rank6, Scalar, Layout )
{
  // Try and create a rank-6 view
  typedef typename Scalar::value_type Vector;
  typedef typename Vector::value_type BaseScalar;
  typedef typename Vector::execution_space Device;
  typedef typename ApplyView<Scalar******,Layout,Device>::type ViewType;

  ViewType v("view", 1, 2, 3, 4, 4, 3, global_fad_size+1);
  BaseScalar val = 1.2345;

  Kokkos::deep_copy( v, val );

  success = checkConstantFadVectorView2(v, Scalar(val), out);
}

#define VIEW_FAD_MP_VECTOR_TESTS_SCALAR_LAYOUT( SCALAR, LAYOUT )        \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_Fad_MP, Size, SCALAR, LAYOUT )                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_Fad_MP, DeepCopy_ConstantScalar, SCALAR, LAYOUT )       \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_Fad_MP, Rank6, SCALAR, LAYOUT )

#define VIEW_FAD_MP_VECTOR_TESTS_SCALAR( SCALAR )                       \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_FAD_MP_VECTOR_TESTS_SCALAR_LAYOUT(SCALAR, NoLayout)              \
  VIEW_FAD_MP_VECTOR_TESTS_SCALAR_LAYOUT(SCALAR, LayoutLeft)            \
  VIEW_FAD_MP_VECTOR_TESTS_SCALAR_LAYOUT(SCALAR, LayoutRight)

#define VIEW_FAD_MP_VECTOR_TESTS_DEVICE( DEVICE )                       \
  typedef Stokhos::StaticFixedStorage<int,double,global_ensemble_size,DEVICE> SFS; \
  typedef Sacado::MP::Vector< SFS > MP_SFS;                             \
  typedef Sacado::Fad::DFad< MP_SFS > Fad_MP_SFS;                       \
  VIEW_FAD_MP_VECTOR_TESTS_SCALAR( Fad_MP_SFS )
