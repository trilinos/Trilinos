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

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestHelpers.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"

// For computing DeviceConfig
#include "Kokkos_Core.hpp"

//
// Tests various View< Sacado::MP::Vector<...>,...> operations work
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
  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
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
checkConstantVectorView(const ViewType& v,
                        const typename ViewType::value_type& v_expected,
                        Teuchos::FancyOStream& out) {
  typedef ViewType view_type;
  typedef typename view_type::size_type size_type;
  typedef typename view_type::HostMirror host_view_type;
  typedef typename host_view_type::array_type::value_type scalar_type;

  // Copy to host
  host_view_type h_v = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(h_v, v);

  const size_type num_rows = h_v.extent(0);
  const size_type num_cols = Kokkos::dimension_scalar(h_v);
  bool success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      scalar_type val = h_v(i).fastAccessCoeff(j);
      scalar_type val_expected = v_expected.fastAccessCoeff(j);
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
const int global_num_cols = 8;  // Currently must be a multiple of 8 based on
                                // alignment assumptions for SFS

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, Size, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  TEUCHOS_TEST_EQUALITY(v.size(), num_rows, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;

  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_vector_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_vector_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_v);

  success = checkVectorView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_ConstantScalar, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  Scalar val = 1.2345;

  Kokkos::deep_copy( v, val );

  success = checkConstantVectorView(v, Vector(num_cols, val), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_ConstantVector, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  Scalar val = 1.2345;

  Kokkos::deep_copy( v, Vector(val) );

  success = checkConstantVectorView(v, Vector(num_cols, val), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_ConstantVector2, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  Vector val(num_cols, 0.0);
  for (size_type j=0; j<num_cols; ++j)
    val.fastAccessCoeff(j) =
      generate_vector_coefficient<Scalar>(num_rows, num_cols, size_type(0), j);

  Kokkos::deep_copy( v, val );

  success = checkConstantVectorView(v, val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Kokkos_View_MP, DeepCopy_Subview_Range, Storage )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector**,Kokkos::LayoutLeft,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows1 = global_num_rows;
  const size_type num_rows2 = global_num_rows*2;
  const size_type num_cols = 5;
  const size_type num_vec =
    Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v1("view1", num_rows1, num_cols, num_vec);
  ViewType v2("view2", num_rows2, num_cols, num_vec);

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
      for (size_type k=0; k<num_vec; ++k) {
        Scalar val = hv2(i,j).fastAccessCoeff(k);
        Scalar val_expected = j+1;
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
    for (size_type i=num_rows1; i<num_rows2; ++i) {
      for (size_type k=0; k<num_vec; ++k) {
        Scalar val = hv2(i,j).fastAccessCoeff(k);
        Scalar val_expected = 0;
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_HostArray, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  host_array_type h_a = Kokkos::create_mirror_view(v);

  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_vector_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_vector_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_a);

  success = checkVectorView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_DeviceArray, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename ViewType::array_type array_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;
  array_type a = v;

  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_a(i,j) = generate_vector_coefficient<Scalar>(
        num_rows, num_cols, i, j);
  Kokkos::deep_copy(a, h_v);

  success = checkVectorView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, Unmanaged, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;

  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_vector_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_vector_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_v);

  // Create unmanaged view
  ViewType v2(v.data(), num_rows, num_cols);

  success = checkVectorView(v2, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, PartitionHost, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::static_size;
  ViewType v("view", num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);

  const size_type num_cols_part = num_cols/2;
  auto h_v1 = Kokkos::partition<num_cols_part>(h_v, 0);
  auto h_v2 = Kokkos::partition<num_cols_part>(h_v, num_cols_part);

  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols_part; ++j) {
      h_v1(i).fastAccessCoeff(j) = generate_vector_coefficient<Scalar>(
        num_rows, num_cols, i, j);
      h_v2(i).fastAccessCoeff(j) = generate_vector_coefficient<Scalar>(
        num_rows, num_cols, i, j+num_cols_part);
    }
  }
  Kokkos::deep_copy(v, h_v);

  success = checkVectorView(v, out);
}

/*
// This test does not work because we can't call deep_copy on partitioned views
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, PartitionDevice, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  // This test requires static storage, so it always passes if it isn't
  if (!Storage::is_static) {
    success = true;
    return;
  }

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::static_size;
  ViewType v("view", num_rows, num_cols);

  const size_type num_cols_part = num_cols/2;
  auto v1 = Kokkos::partition<num_cols_part>(v, 0);
  auto v2 = Kokkos::partition<num_cols_part>(v, num_cols_part);

  typename decltype(v1)::HostMirror h_v1 = Kokkos::create_mirror_view(v1);
  typename decltype(v2)::HostMirror h_v2 = Kokkos::create_mirror_view(v2);

  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols_part; ++j) {
      h_v1(i).fastAccessCoeff(j) = generate_vector_coefficient<Scalar>(
        num_rows, num_cols, i, j);
      h_v2(i).fastAccessCoeff(j) = generate_vector_coefficient<Scalar>(
        num_rows, num_cols, i, j+num_cols_part);
    }
  }
  Kokkos::deep_copy(v1, h_v1);
  Kokkos::deep_copy(v2, h_v2);

  success = checkVectorView(v, out);
}
*/

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, Flatten, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename Kokkos::FlatArrayType<ViewType>::type flat_view_type;
  typedef typename flat_view_type::HostMirror host_flat_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);

  // Create flattened view
  flat_view_type flat_v = v;
  host_flat_view_type h_flat_v = Kokkos::create_mirror_view(flat_v);
  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_flat_v(i*num_cols+j) = generate_vector_coefficient<Scalar>(
        num_rows, num_cols, i, j);
  Kokkos::deep_copy(flat_v, h_flat_v);

  success = checkVectorView(v, out);
}

namespace Test {

template< class ViewType >
struct MPVectorAtomicFunctor {
  typedef typename ViewType::execution_space execution_space ;

  typedef typename ViewType::value_type vector_type ;
  typedef typename vector_type::storage_type::value_type scalar_type ;

  scalar_type m_s ;
  ViewType m_v ;

  MPVectorAtomicFunctor( const ViewType & v , const scalar_type & s ) : m_v( v ), m_s( s )
  {
    Kokkos::parallel_for( m_v.extent(0) , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( int i ) const
  {
    vector_type v( m_s );
    atomic_assign( & m_v(i) , v );
  }
};

}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeviceAtomic, Storage, Layout )
{
  typedef typename Storage::execution_space Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  Scalar val = 1.2345;

  (void) Test::MPVectorAtomicFunctor<ViewType>( v , val );

  success = checkConstantVectorView(v, val, out);
}


#define VIEW_MP_VECTOR_TESTS_STORAGE_LAYOUT( STORAGE, LAYOUT )          \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, Size, STORAGE, LAYOUT )                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy, STORAGE, LAYOUT )                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy_ConstantScalar, STORAGE, LAYOUT )          \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy_ConstantVector, STORAGE, LAYOUT )          \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy_ConstantVector2, STORAGE, LAYOUT )         \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, Unmanaged, STORAGE, LAYOUT )                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, Flatten, STORAGE, LAYOUT )

// Some tests the fail, or fail to compile

  /*
    // These don't compile as there are no deep_copy overloads between
    // static view spec and its array_type.  That could be done, but
    // would require some additional work to ensure the shapes match.
    // It is simple enough to create an array_type view, so deep copying
    // between matching views isn't much more trouble.
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy_HostArray, STORAGE, LAYOUT )               \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                       \
    Kokkos_View_MP, DeepCopy_DeviceArray, STORAGE, LAYOUT )
  */

#define VIEW_MP_VECTOR_TESTS_STORAGE( STORAGE )                         \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_MP_VECTOR_TESTS_STORAGE_LAYOUT(STORAGE, NoLayout)                \
  VIEW_MP_VECTOR_TESTS_STORAGE_LAYOUT(STORAGE, LayoutLeft)              \
  VIEW_MP_VECTOR_TESTS_STORAGE_LAYOUT(STORAGE, LayoutRight)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy_Subview_Range, STORAGE )

#define VIEW_MP_VECTOR_TESTS_ORDINAL_SCALAR_DEVICE( ORDINAL, SCALAR, DEVICE ) \
  typedef Stokhos::StaticFixedStorage<ORDINAL,SCALAR,global_num_cols,DEVICE> SFS;     \
  typedef Stokhos::DynamicStorage<ORDINAL,SCALAR,DEVICE> DS;            \
  VIEW_MP_VECTOR_TESTS_STORAGE( SFS )                                   \
  VIEW_MP_VECTOR_TESTS_STORAGE( DS )                                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, PartitionHost, SFS, NoLayout )                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, PartitionHost, SFS, LayoutLeft )                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, PartitionHost, SFS, LayoutRight )

#define VIEW_MP_VECTOR_TESTS_DEVICE( DEVICE )                           \
  VIEW_MP_VECTOR_TESTS_ORDINAL_SCALAR_DEVICE( int, double, DEVICE )
