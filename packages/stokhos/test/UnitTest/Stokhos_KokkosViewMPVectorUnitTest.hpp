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

#include "Stokhos_Sacado_Kokkos.hpp"

// For computing DeviceConfig
#include "Kokkos_hwloc.hpp"
#include "Kokkos_Cuda.hpp"

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

  size_type num_rows = v.dimension_0();
  size_type num_cols = v.dimension_1();
  bool success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      scalar_type val = h_a(i,j);
      scalar_type val_expected =
        generate_vector_coefficient<scalar_type>(
          num_rows, num_cols, i, j);
      TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
    }
  }

  return success;
}

template <typename ViewType>
bool
checkConstantVectorView(const ViewType& v,
                        const typename ViewType::array_type::value_type& val_expected,
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

  size_type num_rows = v.dimension_0();
  size_type num_cols = v.dimension_1();
  bool success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      scalar_type val = h_a(i,j);
      TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
    }
  }

  return success;
}

template <typename DataType, typename LayoutType, typename DeviceType>
struct ApplyView {
  typedef Kokkos::View<DataType,LayoutType,DeviceType> type;
};

struct NoLayout {};
template <typename DataType, typename DeviceType>
struct ApplyView<DataType,NoLayout,DeviceType> {
  typedef Kokkos::View<DataType,DeviceType> type;
};

//
// Tests
//

const int global_num_rows = 11;
const int global_num_cols = 9;

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy, Storage, Layout )
{
  typedef typename Storage::device_type Device;
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

  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_a(i,j) = generate_vector_coefficient<Scalar>(
        num_rows, num_cols, i, j);
  Kokkos::deep_copy(v, h_v);

  success = checkVectorView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_ConstantScalar, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  Scalar val = 1.2345;

  Kokkos::deep_copy( v, val );

  success = checkConstantVectorView(v, val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_ConstantVector, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::MP::Vector<Storage> Vector;
  typedef typename ApplyView<Vector*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = Storage::is_static ? Storage::static_size : global_num_cols;
  ViewType v("view", num_rows, num_cols);
  Scalar val = 1.2345;

  Kokkos::deep_copy( v, Vector(val) );

  success = checkConstantVectorView(v, val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_HostArray, Storage, Layout )
{
  typedef typename Storage::device_type Device;
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

  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_a(i,j) = generate_vector_coefficient<Scalar>(
        num_rows, num_cols, i, j);
  Kokkos::deep_copy(v, h_a);

  success = checkVectorView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_MP, DeepCopy_DeviceArray, Storage, Layout )
{
  typedef typename Storage::device_type Device;
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

namespace Test {

template< class ViewType >
struct MPVectorAtomicFunctor {
  typedef typename ViewType::device_type device_type ;

  typedef typename ViewType::value_type vector_type ;
  typedef typename vector_type::storage_type::value_type scalar_type ;

  scalar_type m_s ;
  ViewType m_v ;

  MPVectorAtomicFunctor( const ViewType & v , const scalar_type & s ) : m_v( v ), m_s( s )
  {
    Kokkos::parallel_for( m_v.dimension_0() , *this );
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
  typedef typename Storage::device_type Device;
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
  Scalar val = 1.2345;

  (void) Test::MPVectorAtomicFunctor<ViewType>( v , val );

  success = checkConstantVectorView(v, val, out);
}


#define VIEW_MP_VECTOR_TESTS_STORAGE_LAYOUT( STORAGE, LAYOUT )          \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy, STORAGE, LAYOUT )                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy_ConstantScalar, STORAGE, LAYOUT )

// Some tests the fail, or fail to compile

  /*
   // This doesn't compile because deep_copy() only has an overload for
   // scalar_type, and not value_type
   TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                \
    Kokkos_View_MP, DeepCopy_ConstantVector, STORAGE, LAYOUT )
  */

  /*
    // This compiles but gives the wrong answer for LayoutLeft or LayoutRight
    // (but is correct for default layout)
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_MP, DeepCopy_HostArray, STORAGE, LAYOUT )               \
  */

  /*
    // This doesn't compile, I believe for the same reason the above generates
    // the wrong answer
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                       \
    Kokkos_View_MP, DeepCopy_DeviceArray, STORAGE, LAYOUT )
  */

#define VIEW_MP_VECTOR_TESTS_STORAGE( STORAGE )                         \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_MP_VECTOR_TESTS_STORAGE_LAYOUT(STORAGE, NoLayout)                \
  VIEW_MP_VECTOR_TESTS_STORAGE_LAYOUT(STORAGE, LayoutLeft)              \
  VIEW_MP_VECTOR_TESTS_STORAGE_LAYOUT(STORAGE, LayoutRight)

#define VIEW_MP_VECTOR_TESTS_SCALAR_ORDINAL_DEVICE( SCALAR, ORDINAL, DEVICE ) \
  typedef Stokhos::StaticFixedStorage<SCALAR,ORDINAL,global_num_cols,DEVICE> SFS;     \
  VIEW_MP_VECTOR_TESTS_STORAGE( SFS )

#define VIEW_MP_VECTOR_TESTS_DEVICE( DEVICE ) \
  VIEW_MP_VECTOR_TESTS_SCALAR_ORDINAL_DEVICE( int, double, DEVICE )
