// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER
#include "Teuchos_TestingHelpers.hpp"

#include "Sacado_Kokkos.hpp"

template <typename FadType1, typename FadType2>
bool checkFads(const FadType1& x, const FadType2& x2,
               Teuchos::FancyOStream& out, double tol = 1.0e-15)
{
  bool success = true;

  // Check sizes match
  TEUCHOS_TEST_EQUALITY(x.size(), x2.size(), out, success);

  // Check values match
  TEUCHOS_TEST_EQUALITY(x.val(), x2.val(), out, success);

  // Check derivatives match
  for (int i=0; i<x.size(); ++i)
    TEUCHOS_TEST_FLOATING_EQUALITY(x.dx(i), x2.dx(i), tol, out, success);

  return success;
}

template <typename fadtype, typename ordinal>
inline
fadtype generate_fad( const ordinal num_rows,
                      const ordinal num_cols,
                      const ordinal fad_size,
                      const ordinal row,
                      const ordinal col )
{
  typedef typename fadtype::value_type scalar;
  fadtype x(fad_size, scalar(0.0));

  const scalar x_row = 100.0 + scalar(num_rows) / scalar(row+1);
  const scalar x_col =  10.0 + scalar(num_cols) / scalar(col+1);
  x.val() = x_row + x_col;
  for (ordinal i=0; i<fad_size; ++i) {
    const scalar x_fad = 1.0 + scalar(fad_size) / scalar(i+1);
    x.fastAccessDx(i) = x_row + x_col + x_fad;
  }
  return x;
}

template <typename DataType, typename LayoutType, typename DeviceType,
          typename Memory = void>
struct ApplyView {
  typedef Kokkos::View<DataType,LayoutType,DeviceType,Memory> type;
};

struct NoLayout {};
template <typename DataType, typename DeviceType, typename Memory>
struct ApplyView<DataType,NoLayout,DeviceType,Memory> {
  typedef Kokkos::View<DataType,DeviceType,Memory> type;
};

const int global_num_rows = 11;
const int global_num_cols = 7;
const int global_fad_size = 5;

// Kernel to multiply two views
template <typename InputViewType, typename OutputViewType = InputViewType>
struct MultiplyKernel {
  typedef typename InputViewType::device_type device_type;
  typedef typename InputViewType::size_type size_type;

  const InputViewType  m_v1, m_v2;
  const OutputViewType m_v3;

  MultiplyKernel(const InputViewType  v1,
                 const InputViewType  v2,
                 const OutputViewType v3) :
    m_v1(v1), m_v2(v2), m_v3(v3) {};

  // Multiply entries for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    m_v3(i) = m_v1(i)*m_v2(i);
  }

  // Kernel launch
  static void apply(const InputViewType  v1,
                    const InputViewType  v2,
                    const OutputViewType v3) {
    const size_type nrow = v1.dimension_0();
    Kokkos::parallel_for( nrow, MultiplyKernel(v1,v2,v3) );
  }
};

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DeepCopy, FadType, Layout, Device )
{
  typedef typename ApplyView<FadType**,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_v(i,j) = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
  Kokkos::deep_copy(v, h_v);

  // Copy back
  host_view_type h_v2 = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(h_v2, v);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      success = success && checkFads(f, h_v2(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DeepCopy_ConstantScalar, FadType, Layout, Device )
{
  typedef typename ApplyView<FadType**,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadType::value_type value_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant scalar
  value_type a = 2.3456;
  Kokkos::deep_copy( v, a );

  // Copy to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, v);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
      FadType f = FadType(fad_size, a);
#else
      FadType f = a;
#endif
      success = success && checkFads(f, hv(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DeepCopy_ConstantZero, FadType, Layout, Device )
{
  typedef typename ApplyView<FadType**,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadType::value_type value_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant scalar
  value_type a = 0.0;
  Kokkos::deep_copy( v, a );

  // Copy to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, v);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
      FadType f = FadType(fad_size, a);
#else
      FadType f = a;
#endif
      success = success && checkFads(f, hv(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DeepCopy_ConstantFad, FadType, Layout, Device )
{
  typedef typename ApplyView<FadType**,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant scalar
  FadType a = 2.3456;
  Kokkos::deep_copy( v, a );

  // Copy to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, v);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
      FadType f = FadType(fad_size, a.val());
#else
      FadType f = a;
#endif
      success = success && checkFads(f, hv(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DeepCopy_ConstantFadFull, FadType, Layout, Device )
{
  typedef typename ApplyView<FadType**,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant Fad
  FadType a(fad_size, 2.3456);
  for (size_type i=0; i<fad_size; ++i)
    a.fastAccessDx(i) = 7.89 + (i+1);
  Kokkos::deep_copy( v, a );

  // Copy to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, v);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      success = success && checkFads(a, hv(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Multiply, FadType, Layout, Device )
{
  typedef typename ApplyView<FadType*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  ViewType v1("view1", num_rows, fad_size+1), v2("view2", num_rows, fad_size+1);
  host_view_type h_v1 = Kokkos::create_mirror_view(v1);
  host_view_type h_v2 = Kokkos::create_mirror_view(v2);
  for (size_type i=0; i<num_rows; ++i) {
    h_v1(i) = generate_fad<FadType>(
      num_rows, size_type(2), fad_size, i, size_type(0));
    h_v2(i) = generate_fad<FadType>(
      num_rows, size_type(2), fad_size, i, size_type(1));
  }
  Kokkos::deep_copy(v1, h_v1);
  Kokkos::deep_copy(v2, h_v2);

  // Launch kernel
  ViewType v3("view3", num_rows, fad_size+1);
  MultiplyKernel<ViewType>::apply(v1,v2,v3);

  // Copy back
  host_view_type h_v3 = Kokkos::create_mirror_view(v3);
  Kokkos::deep_copy(h_v3, v3);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    FadType f1 =
      generate_fad<FadType>(num_rows, size_type(2), fad_size, i, size_type(0));
    FadType f2 =
      generate_fad<FadType>(num_rows, size_type(2), fad_size, i, size_type(1));
    FadType f3 = f1*f2;
    success = success && checkFads(f3, h_v3(i), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, MultiplyConst, FadType, Layout, Device )
{
  typedef typename ApplyView<const FadType*,Layout,Device,Kokkos::MemoryUnmanaged>::type ConstViewType;
  typedef typename ApplyView<FadType*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  ViewType v1("view1", num_rows, fad_size+1);
  ViewType v2("view2", num_rows, fad_size+1);
  host_view_type h_v1 = Kokkos::create_mirror_view(v1);
  host_view_type h_v2 = Kokkos::create_mirror_view(v2);
  for (size_type i=0; i<num_rows; ++i) {
    h_v1(i) = generate_fad<FadType>(
      num_rows, size_type(2), fad_size, i, size_type(0));
    h_v2(i) = generate_fad<FadType>(
      num_rows, size_type(2), fad_size, i, size_type(1));
  }
  Kokkos::deep_copy(v1, h_v1);
  Kokkos::deep_copy(v2, h_v2);

  ConstViewType cv1 = v1;
  ConstViewType cv2 = v2;

  // Launch kernel
  ViewType v3("view3", num_rows, fad_size+1);
  MultiplyKernel<ConstViewType,ViewType>::apply(cv1,cv2,v3);

  // Copy back
  host_view_type h_v3 = Kokkos::create_mirror_view(v3);
  Kokkos::deep_copy(h_v3, v3);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    FadType f1 =
      generate_fad<FadType>(num_rows, size_type(2), fad_size, i, size_type(0));
    FadType f2 =
      generate_fad<FadType>(num_rows, size_type(2), fad_size, i, size_type(1));
    FadType f3 = f1*f2;
    success = success && checkFads(f3, h_v3(i), out);
  }
}

// Tests that require view spec

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ShmemSize, FadType, Layout, Device )
{
  typedef typename ApplyView<FadType**,Layout,Device>::type ViewType;
  typedef typename FadType::value_type value_type;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Compute shared memory size for View
  const size_type shmem_size =
    ViewType::shmem_size(num_rows, num_cols, fad_size+1);

  // Check
  static const size_type align = 8;
  static const size_type mask  = align - 1;
  const size_type shmem_size_expected =
    ( sizeof(value_type) * global_num_rows * global_num_cols * (fad_size+1) +
      mask ) & ~mask;
  TEUCHOS_TEST_EQUALITY(shmem_size, shmem_size_expected, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Unmanaged, FadType, Layout, Device )
{
  typedef typename FadType::value_type scalar_type;
  typedef typename ApplyView<scalar_type***,Layout,Device>::type ViewType;
  typedef typename ApplyView<FadType**,Layout,Device,Kokkos::MemoryUnmanaged>::type FadViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadViewType::HostMirror fad_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j,0) = f.val();
      for (size_type k=0; k<fad_size; k++)
        h_v(i,j,k+1) = f.dx(k);
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Create unmanaged view
  FadViewType v_fad(v.ptr_on_device(), num_rows, num_cols, fad_size+1);

  // Copy back -- can't use create_mirror_view() because v_fad is unmanaged
  fad_host_view_type h_v_fad("host_view_fad", num_rows, num_cols, fad_size+1);
  Kokkos::deep_copy(h_v_fad, v_fad);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      success = success && checkFads(f, h_v_fad(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Unmanaged2, FadType, Layout, Device )
{
  typedef typename FadType::value_type scalar_type;
  typedef typename ApplyView<scalar_type***,Layout,Device>::type ViewType;
  typedef typename ApplyView<FadType**,Layout,Device>::type FadViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadViewType::HostMirror fad_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j,0) = f.val();
      for (size_type k=0; k<fad_size; k++)
        h_v(i,j,k+1) = f.dx(k);
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Create unmanaged view
  FadViewType v_fad( v.ptr_on_device(), num_rows, num_cols, fad_size+1);

  // Copy back -- can't use create_mirror_view() because v_fad is unmanaged
  fad_host_view_type h_v_fad("host_view_fad", num_rows, num_cols, fad_size+1);
  Kokkos::deep_copy(h_v_fad, v_fad);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      success = success && checkFads(f, h_v_fad(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, UnmanagedConst, FadType, Layout, Device )
{
  typedef typename FadType::value_type scalar_type;
  typedef typename ApplyView<scalar_type***,Layout,Device>::type ViewType;
  typedef typename ApplyView<const scalar_type***,Layout,Device>::type ConstViewType;
  typedef typename ApplyView<FadType**,Layout,Device,Kokkos::MemoryUnmanaged>::type FadViewType;
  typedef typename ApplyView<const FadType**,Layout,Device,Kokkos::MemoryUnmanaged>::type ConstFadViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadViewType::HostMirror fad_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j,0) = f.val();
      for (size_type k=0; k<fad_size; k++)
        h_v(i,j,k+1) = f.dx(k);
    }
  }
  Kokkos::deep_copy(v, h_v);
  ConstViewType v_const = v;

  // Create unmanaged view
  ConstFadViewType v_fad(
    v_const.ptr_on_device(), num_rows, num_cols, fad_size+1);

  // Copy back -- can't use create_mirror_view() because v_fad is unmanaged
  fad_host_view_type h_v_fad("host_view_fad", num_rows, num_cols, fad_size+1);
  Kokkos::deep_copy(h_v_fad, v_fad);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      success = success && checkFads(f, h_v_fad(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, UnmanagedConst2, FadType, Layout, Device )
{
  typedef typename FadType::value_type scalar_type;
  typedef typename ApplyView<scalar_type***,Layout,Device>::type ViewType;
  typedef typename ApplyView<const scalar_type***,Layout,Device>::type ConstViewType;
  typedef typename ApplyView<FadType**,Layout,Device>::type FadViewType;
  typedef typename ApplyView<const FadType**,Layout,Device>::type ConstFadViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadViewType::HostMirror fad_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j,0) = f.val();
      for (size_type k=0; k<fad_size; k++)
        h_v(i,j,k+1) = f.dx(k);
    }
  }
  Kokkos::deep_copy(v, h_v);
  ConstViewType v_const = v;

  // Create unmanaged view
  ConstFadViewType v_fad(v_const.ptr_on_device(), num_rows, num_cols, fad_size+1);

  // Copy back -- can't use create_mirror_view() because v_fad is unmanaged
  fad_host_view_type h_v_fad("host_view_fad", num_rows, num_cols, fad_size+1);
  Kokkos::deep_copy(h_v_fad, v_fad);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      success = success && checkFads(f, h_v_fad(i,j), out);
    }
  }
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ShmemSize, FadType, Layout, Device )
{
  typedef typename ApplyView<FadType**,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;

  // Compute shared memory size for View
  const size_type shmem_size =
    ViewType::shmem_size(num_rows, num_cols);

  // Check
  static const size_type align = 8;
  static const size_type mask  = align - 1;
  const size_type shmem_size_expected =
    ( sizeof(FadType) * global_num_rows * global_num_cols + mask ) & ~mask;
  TEUCHOS_TEST_EQUALITY(shmem_size, shmem_size_expected, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Unmanaged, FadType, Layout, Device ) {}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Unmanaged2, FadType, Layout, Device ) {}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, UnmanagedConst, FadType, Layout, Device ) {}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, UnmanagedConst2, FadType, Layout, Device ) {}

#endif

#define VIEW_FAD_TESTS_FLD( F, L, D )                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy_ConstantScalar, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy_ConstantZero, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy_ConstantFad, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy_ConstantFadFull, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Unmanaged, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Unmanaged2, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, UnmanagedConst, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, UnmanagedConst2, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Multiply, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, MultiplyConst, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, ShmemSize, F, L, D )

#define VIEW_FAD_TESTS_FD( F, D )                                       \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_FAD_TESTS_FLD( F, NoLayout, D)                                   \
  VIEW_FAD_TESTS_FLD( F, LayoutLeft, D)                                 \
  VIEW_FAD_TESTS_FLD( F, LayoutRight, D)

typedef Sacado::Fad::DFad<double> DFadType;
typedef Sacado::Fad::SLFad<double,2*global_fad_size> SLFadType;
typedef Sacado::Fad::SFad<double,global_fad_size> SFadType;

typedef Sacado::ELRFad::DFad<double> ELRDFadType;
typedef Sacado::ELRFad::SLFad<double,2*global_fad_size> ELRSLFadType;
typedef Sacado::ELRFad::SFad<double,global_fad_size> ELRSFadType;

typedef Sacado::CacheFad::DFad<double> CacheDFadType;
typedef Sacado::CacheFad::SLFad<double,2*global_fad_size> CacheSLFadType;
typedef Sacado::CacheFad::SFad<double,global_fad_size> CacheSFadType;

typedef Sacado::ELRCacheFad::DFad<double> ELRCacheDFadType;
typedef Sacado::ELRCacheFad::SLFad<double,2*global_fad_size> ELRCacheSLFadType;
typedef Sacado::ELRCacheFad::SFad<double,global_fad_size> ELRCacheSFadType;

// We can't use DFad unless we use the View specialization
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
#define VIEW_FAD_TESTS_D( D )                           \
  VIEW_FAD_TESTS_FD( SFadType, D )                      \
  VIEW_FAD_TESTS_FD( SLFadType, D )                     \
  VIEW_FAD_TESTS_FD( DFadType, D )                      \
  VIEW_FAD_TESTS_FD( ELRSFadType, D )                   \
  VIEW_FAD_TESTS_FD( ELRSLFadType, D )                  \
  VIEW_FAD_TESTS_FD( ELRDFadType, D )                   \
  VIEW_FAD_TESTS_FD( CacheSFadType, D )                 \
  VIEW_FAD_TESTS_FD( CacheSLFadType, D )                \
  VIEW_FAD_TESTS_FD( CacheDFadType, D )                 \
  VIEW_FAD_TESTS_FD( ELRCacheSFadType, D )              \
  VIEW_FAD_TESTS_FD( ELRCacheSLFadType, D )             \
  VIEW_FAD_TESTS_FD( ELRCacheDFadType, D )
#else
#define VIEW_FAD_TESTS_D( D )                        \
  VIEW_FAD_TESTS_FD( SFadType, D )                   \
  VIEW_FAD_TESTS_FD( SLFadType, D )                  \
  VIEW_FAD_TESTS_FD( ELRSFadType, D )                \
  VIEW_FAD_TESTS_FD( ELRSLFadType, D )               \
  VIEW_FAD_TESTS_FD( CacheSFadType, D )              \
  VIEW_FAD_TESTS_FD( CacheSLFadType, D )             \
  VIEW_FAD_TESTS_FD( ELRCacheSFadType, D )           \
  VIEW_FAD_TESTS_FD( ELRCacheSLFadType, D )
#endif
