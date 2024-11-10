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

template <typename T>
struct is_dfad {
  static const bool value = false;
};

template <typename T>
struct is_dfad< Sacado::Fad::Exp::DFad<T> > {
  static const bool value = true;
};

template <typename FadType1, typename FadType2>
bool checkFads(const FadType1& x, const FadType2& x2,
               Teuchos::FancyOStream& out, double tol = 1.0e-15)
{
  bool success = true;

  // Check sizes match
  TEUCHOS_TEST_EQUALITY(x.size(), x2.size(), out, success);

  // Check values match
  TEUCHOS_TEST_FLOATING_EQUALITY(x.val(), x2.val(), tol, out, success);

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

#ifndef GLOBAL_FAD_SIZE
#define GLOBAL_FAD_SIZE 5
#endif
const int global_num_rows = 11;
const int global_num_cols = 7;
const int global_fad_size = GLOBAL_FAD_SIZE;

struct AddTag {
  static double init() { return 0.0; }
  template <typename T1, typename T2>
  static auto apply(const T1& a, const T2& b) -> decltype(a+b)
  {
    return a+b;
  }
};
struct SubTag {
  static double init() { return 0.0; }
  template <typename T1, typename T2>
  static auto apply(const T1& a, const T2& b) -> decltype(a-b)
  {
    return a-b;
  }
};
struct MulTag {
  static double init() { return 1.0; }
  template <typename T1, typename T2>
  static auto apply(const T1& a, const T2& b) -> decltype(a*b)
  {
    return a*b;
  }
};
struct DivTag {
  static double init() { return 1.0; }
  template <typename T1, typename T2>
  static auto apply(const T1& a, const T2& b) -> decltype(a/b)
  {
    return a/b;
  }
};
struct MaxTag {
  static double init() { return 1.0; }
  template <typename T1, typename T2>
  static auto apply(const T1& a, const T2& b) -> decltype(max(a,b))
  {
    return max(a,b);
  }
};
struct MinTag {
  static double init() { return 1.0; }
  template <typename T1, typename T2>
  static auto apply(const T1& a, const T2& b) -> decltype(min(a,b))
  {
    return min(a,b);
  }
};

// Kernel to test atomic_add
template <typename ViewType, typename ScalarViewType, bool OperFetch>
struct AtomicKernel {
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::size_type size_type;
  typedef typename Kokkos::TeamPolicy< execution_space>::member_type team_handle;
  typedef typename Kokkos::ThreadLocalScalarType<ViewType>::type local_scalar_type;
  static const size_type stride = Kokkos::ViewScalarStride<ViewType>::stride;

  const ViewType m_v;
  const ScalarViewType m_s;

  AtomicKernel(const ViewType& v, const ScalarViewType& s) :
    m_v(v), m_s(s) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (AddTag tag, const size_type i) const {
    local_scalar_type x = m_v(i);
    if (OperFetch)
      Kokkos::atomic_add_fetch(&(m_s()), x);
    else
      Kokkos::atomic_fetch_add(&(m_s()), x);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (SubTag tag, const size_type i) const {
    local_scalar_type x = m_v(i);
    if (OperFetch)
      Kokkos::atomic_sub_fetch(&(m_s()), x);
    else
      Kokkos::atomic_fetch_sub(&(m_s()), x);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (MulTag tag, const size_type i) const {
    local_scalar_type x = m_v(i);
    if (OperFetch)
      Kokkos::atomic_mul_fetch(&(m_s()), x);
    else
      Kokkos::atomic_fetch_mul(&(m_s()), x);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (DivTag tag, const size_type i) const {
    local_scalar_type x = m_v(i);
    if (OperFetch)
      Kokkos::atomic_div_fetch(&(m_s()), x);
    else
      Kokkos::atomic_fetch_div(&(m_s()), x);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (MaxTag tag, const size_type i) const {
    local_scalar_type x = m_v(i);
    if (OperFetch)
      Kokkos::atomic_max_fetch(&(m_s()), x);
    else
      Kokkos::atomic_fetch_max(&(m_s()), x);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (MinTag tag, const size_type i) const {
    local_scalar_type x = m_v(i);
    if (OperFetch)
      Kokkos::atomic_min_fetch(&(m_s()), x);
    else
      Kokkos::atomic_fetch_min(&(m_s()), x);
  }

  template <typename Tag>
  KOKKOS_INLINE_FUNCTION
  void operator()( Tag tag, const team_handle& team ) const
  {
    const size_type i = team.league_rank()*team.team_size() + team.team_rank();
    if (i < m_v.extent(0))
      (*this)(tag, i);
  }

  // Kernel launch
  template <typename Tag>
  static void apply(Tag tag, const ViewType& v, const ScalarViewType& s) {
    const size_type nrow = v.extent(0);

#if defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      Kokkos::is_view_fad_contiguous<ViewType>::value &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      Kokkos::is_view_fad_contiguous<ViewType>::value &&
      is_dfad<typename ViewType::non_const_value_type>::value;
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      Kokkos::is_view_fad_contiguous<ViewType>::value &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      Kokkos::is_view_fad_contiguous<ViewType>::value &&
      is_dfad<typename ViewType::non_const_value_type>::value;
#else
    const bool use_team = false;
#endif

    if (use_team) {
      const size_type team_size = 256 / stride;
      Kokkos::TeamPolicy<execution_space, Tag> policy(
        (nrow+team_size-1)/team_size, team_size, stride );
      Kokkos::parallel_for( policy, AtomicKernel(v,s) );
    }
    else {
      Kokkos::RangePolicy<execution_space, Tag> policy( 0, nrow );
      Kokkos::parallel_for( policy, AtomicKernel(v,s) );
    }
  }
};

template <typename FadType, typename Layout, typename Device, bool OperFetch,
          typename TagType>
bool testAtomic(const TagType& tag, Teuchos::FancyOStream& out)
{
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef Kokkos::View<FadType,Layout,Device> ScalarViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename ScalarViewType::HostMirror host_scalar_view_type;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
  ScalarViewType s0;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v  = ViewType ("view", num_rows);
  s0 = ScalarViewType ("");
#else
  v  = ViewType ("view", num_rows, fad_size+1);
  s0 = ScalarViewType ("", fad_size+1);
#endif
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i)
    h_v(i) =
      generate_fad<FadType>(num_rows, size_type(1), fad_size, i, size_type(0));
  Kokkos::deep_copy(v, h_v);

  Kokkos::deep_copy(s0, tag.init());

  // Create scalar view
  ScalarViewType s;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  s = ScalarViewType ("scalar view");
#else
  s = ScalarViewType ("scalar view", fad_size+1);
#endif
  Kokkos::deep_copy( s, tag.init() );

  // Call atomic_add kernel, which adds up entries in v
  AtomicKernel<ViewType,ScalarViewType,OperFetch>::apply( tag, v, s );

  // Copy to host
  host_scalar_view_type hs = Kokkos::create_mirror_view(s);
  Kokkos::deep_copy(hs, s);

  // Compute correct result
  auto b = Kokkos::create_mirror_view(s0);
  Kokkos::deep_copy(b, s0);

  for (size_type i=0; i<num_rows; ++i)
    b() = tag.apply(b(), h_v(i));

  // Check
  bool success = checkFads(b(), hs(), out);

  return success;
}

// Test atomic_oper_fetch form

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicAddFetch, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, true>(AddTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicSubFetch, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, true>(SubTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicMulFetch, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, true>(MulTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicDivFetch, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, true>(DivTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicMaxFetch, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, true>(MaxTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicMinFetch, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, true>(MinTag(), out);
}

// Test atomic_fetch_oper form

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicFetchAdd, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, false>(AddTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicFetchSub, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, false>(SubTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicFetchMul, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, false>(MulTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicFetchDiv, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, false>(DivTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicFetchMax, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, false>(MaxTag(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicFetchMin, FadType, Layout, Device )
{
  success = testAtomic<FadType, Layout, Device, false>(MinTag(), out);
}

#define VIEW_FAD_TESTS_FLD( F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicAddFetch, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicSubFetch, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicMulFetch, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicDivFetch, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicMaxFetch, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicMinFetch, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicFetchAdd, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicFetchSub, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicFetchMul, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicFetchDiv, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicFetchMax, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicFetchMin, F, L, D )

using Kokkos::LayoutLeft;
using Kokkos::LayoutRight;
typedef Kokkos::LayoutContiguous<Kokkos::LayoutLeft> LeftContiguous;
typedef Kokkos::LayoutContiguous<Kokkos::LayoutRight> RightContiguous;

#define VIEW_FAD_TESTS_FD( F, D )                                       \
  VIEW_FAD_TESTS_FLD( F, LayoutLeft, D )                                \
  VIEW_FAD_TESTS_FLD( F, LayoutRight, D )                               \
  VIEW_FAD_TESTS_FLD( F, LeftContiguous, D )                            \
  VIEW_FAD_TESTS_FLD( F, RightContiguous, D )

// Full set of atomics only implemented for new design
#if SACADO_ENABLE_NEW_DESIGN
typedef Sacado::Fad::Exp::DFad<double> DFadType;
typedef Sacado::Fad::Exp::SLFad<double,2*global_fad_size> SLFadType;
typedef Sacado::Fad::Exp::SFad<double,global_fad_size> SFadType;

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
