// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_TestingHelpers.hpp"

#include "Sacado.hpp"
#include "Kokkos_DynRankView_Fad.hpp"

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

template <typename FadType1, typename FadType2>
bool checkNestedFads(const FadType1& x, const FadType2& x2,
                     Teuchos::FancyOStream& out, double tol = 1.0e-15)
{
  bool success = true;

  // Check sizes match
  TEUCHOS_TEST_EQUALITY(x.size(), x2.size(), out, success);

  // Check values match
  success = success && checkFads(x.val(), x2.val(), out, tol);

  // Check derivatives match
  for (int i=0; i<x.size(); ++i)
    success = success && checkFads(x.dx(i), x2.dx(i), out, tol);

  return success;
}

template <typename fadfadtype, typename ordinal>
inline
fadfadtype generate_nested_fad( const ordinal num_rows,
                                const ordinal num_cols,
                                const ordinal outer_fad_size,
                                const ordinal inner_fad_size,
                                const ordinal row,
                                const ordinal col )
{
  typedef typename fadfadtype::value_type fadtype;
  typedef typename fadtype::value_type scalar;
  fadfadtype x(outer_fad_size, scalar(0.0));
  fadtype y(inner_fad_size, scalar(0.0));

  const scalar x_row = 1000.0 + scalar(num_rows) / scalar(row+1);
  const scalar x_col =  100.0 + scalar(num_cols) / scalar(col+1);
  y.val() = x_row + x_col;
  for (ordinal j=0; j<inner_fad_size; ++j) {
    const scalar y_fad = 1.0 + scalar(inner_fad_size) / scalar(j+1);
      y.fastAccessDx(j) = x_row + x_col + y_fad;
  }
  x.val() = y;
  for (ordinal i=0; i<outer_fad_size; ++i) {
    const scalar x_fad = 10.0 + scalar(outer_fad_size) / scalar(i+1);
    y.val() = x_fad;
    for (ordinal j=0; j<inner_fad_size; ++j) {
      const scalar y_fad = 1.0 + scalar(inner_fad_size) / scalar(j+1);
      y.fastAccessDx(j) = x_row + x_col + x_fad + y_fad;
    }
    x.fastAccessDx(i) = y;
  }
  return x;
}

const int global_num_rows = 11;
const int global_num_cols = 7;
const int global_outer_fad_size = 5;
const int global_inner_fad_size = 3;

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_FadFad, DeepCopy, FadFadType, Layout, Device )
{
  typedef Kokkos::View<FadFadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type outer_fad_size = global_outer_fad_size;
  const size_type inner_fad_size = global_inner_fad_size;

  // Create and fill view
  ViewType v1("view1", num_rows, num_cols, outer_fad_size+1);
  host_view_type h_v1 = Kokkos::create_mirror_view(v1);
  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_v1(i,j) = generate_nested_fad<FadFadType>(num_rows,
                                                  num_cols,
                                                  outer_fad_size,
                                                  inner_fad_size,
                                                  i, j);
  Kokkos::deep_copy(v1, h_v1);

  // Deep copy
  ViewType v2("view2", num_rows, num_cols, outer_fad_size+1);
  Kokkos::deep_copy(v2, v1);

  // Copy back
  host_view_type h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadFadType f = generate_nested_fad<FadFadType>(num_rows,
                                                     num_cols,
                                                     outer_fad_size,
                                                     inner_fad_size,
                                                     i, j);
      success = success && checkNestedFads(f, h_v2(i,j), out);
    }
  }
}

#ifdef HAVE_SACADO_KOKKOS

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_DynRankView_FadFad, DeepCopy, FadFadType, Layout, Device )
{
  typedef Kokkos::DynRankView<FadFadType,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type outer_fad_size = global_outer_fad_size;
  const size_type inner_fad_size = global_inner_fad_size;

  // Create and fill view
  ViewType v1("view1", num_rows, num_cols, outer_fad_size+1);
  host_view_type h_v1 = Kokkos::create_mirror_view(v1);
  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_v1(i,j) = generate_nested_fad<FadFadType>(num_rows,
                                                  num_cols,
                                                  outer_fad_size,
                                                  inner_fad_size,
                                                  i, j);
  Kokkos::deep_copy(v1, h_v1);

  // Deep copy
  ViewType v2("view2", num_rows, num_cols, outer_fad_size+1);
  Kokkos::deep_copy(v2, v1);

  // Copy back
  host_view_type h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadFadType f = generate_nested_fad<FadFadType>(num_rows,
                                                     num_cols,
                                                     outer_fad_size,
                                                     inner_fad_size,
                                                     i, j);
      success = success && checkNestedFads(f, h_v2(i,j), out);
    }
  }
}

// To test DynRankView - View interoperabitlity
// Deep copy of DynRankView to View
// Assignment operator
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_DynRankView_FadFad, Interop, FadFadType, Layout, Device )
{
  typedef Kokkos::DynRankView<FadFadType,Layout,Device> DRViewType;
  typedef typename DRViewType::size_type size_type;
  typedef typename DRViewType::HostMirror host_view_type;

  typedef Kokkos::View<FadFadType**,Layout,Device> NoDynViewType;
  typedef typename NoDynViewType::HostMirror host_nondynrankview_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type outer_fad_size = global_outer_fad_size;
  const size_type inner_fad_size = global_inner_fad_size;

  // Create and fill view
  DRViewType v1("drview1", num_rows, num_cols, outer_fad_size+1);
  host_view_type h_v1 = Kokkos::create_mirror_view(v1);

  NoDynViewType ndv2("nodview2", num_rows, num_cols, outer_fad_size+1);
  host_nondynrankview_type h_ndv2 = Kokkos::create_mirror_view(ndv2);

  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_v1(i,j) = generate_nested_fad<FadFadType>(num_rows,
                                                  num_cols,
                                                  outer_fad_size,
                                                  inner_fad_size,
                                                  i, j);
  Kokkos::deep_copy(v1, h_v1); //v1 unused here

  // Deep copy DynRankView to View on device
  Kokkos::deep_copy(ndv2, h_v1);
  // Assign View to DynRankView
  DRViewType v2("drview2", num_rows, num_cols, outer_fad_size+1);
  v2 = ndv2 ;

  // Copy back
  host_view_type h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadFadType f = generate_nested_fad<FadFadType>(num_rows,
                                                     num_cols,
                                                     outer_fad_size,
                                                     inner_fad_size,
                                                     i, j);
      success = success && checkNestedFads(f, h_v2(i,j), out);
    }
  }
}


// Deep copy of DynRankView to View
// Copy ctor 
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_DynRankView_FadFad, Interop2, FadFadType, Layout, Device )
{
  typedef Kokkos::DynRankView<FadFadType,Layout,Device> DRViewType;
  typedef typename DRViewType::size_type size_type;
  typedef typename DRViewType::HostMirror host_view_type;

  typedef Kokkos::View<FadFadType**,Layout,Device> NoDynViewType;
  typedef typename NoDynViewType::HostMirror host_nondynrankview_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type outer_fad_size = global_outer_fad_size;
  const size_type inner_fad_size = global_inner_fad_size;

  // Create and fill view
  DRViewType v1("drview1", num_rows, num_cols, outer_fad_size+1);
  host_view_type h_v1 = Kokkos::create_mirror_view(v1);

  NoDynViewType ndv2("nodview2", num_rows, num_cols, outer_fad_size+1);
  host_nondynrankview_type h_ndv2 = Kokkos::create_mirror_view(ndv2);

  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_v1(i,j) = generate_nested_fad<FadFadType>(num_rows,
                                                  num_cols,
                                                  outer_fad_size,
                                                  inner_fad_size,
                                                  i, j);
  Kokkos::deep_copy(v1, h_v1); //v1 unused here

  // Deep copy DynRankView to View on device
  Kokkos::deep_copy(ndv2, h_v1);
  // Copy construct DynRankView from View
  DRViewType v2(ndv2) ;

  // Copy back
  host_view_type h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadFadType f = generate_nested_fad<FadFadType>(num_rows,
                                                     num_cols,
                                                     outer_fad_size,
                                                     inner_fad_size,
                                                     i, j);
      success = success && checkNestedFads(f, h_v2(i,j), out);
    }
  }
}



#else

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_DynRankView_FadFad, DeepCopy, FadFadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_DynRankView_FadFad, Interop, FadFadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_DynRankView_FadFad, Interop2, FadFadType, Layout, Device ) {}

#endif

#define VIEW_FAD_TESTS_FLD( F, L, D )                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_FadFad, DeepCopy, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_DynRankView_FadFad, DeepCopy, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_DynRankView_FadFad, Interop, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_DynRankView_FadFad, Interop2, F, L, D )

#define VIEW_FAD_TESTS_FD( F, D )                                       \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_FAD_TESTS_FLD( F, LayoutLeft, D)                                 \
  VIEW_FAD_TESTS_FLD( F, LayoutRight, D)

// We've unified the implementation for the different Fad variants, so
// there is no reason to test ELRFad, CacheFad, and ELRCacheFad.
typedef Sacado::Fad::SFad<double,global_inner_fad_size> InnerFadType;
typedef Sacado::Fad::DFad<InnerFadType> DFadType;
typedef Sacado::Fad::SLFad<InnerFadType,2*global_outer_fad_size> SLFadType;
typedef Sacado::Fad::SFad<InnerFadType,global_outer_fad_size> SFadType;

// These tests are only relevant when we have the experimental view
// specialization
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
#if SACADO_TEST_DFAD
#define VIEW_FAD_TESTS_D( D )                            \
  VIEW_FAD_TESTS_FD( SFadType, D )                       \
  VIEW_FAD_TESTS_FD( SLFadType, D )                      \
  VIEW_FAD_TESTS_FD( DFadType, D )
#else
#define VIEW_FAD_TESTS_D( D )                            \
  VIEW_FAD_TESTS_FD( SFadType, D )                       \
  VIEW_FAD_TESTS_FD( SLFadType, D )
#endif
#else
#define VIEW_FAD_TESTS_D( D ) /* */
#endif
