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

template <typename T>
struct is_sfad {
  static const bool value = false;
};

template <typename T, int N>
struct is_sfad< Sacado::Fad::SFad<T,N> > {
  static const bool value = true;
};

template <typename T>
struct is_dfad {
  static const bool value = false;
};

template <typename T>
struct is_dfad< Sacado::Fad::DFad<T> > {
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

// Kernel to multiply two views
template <typename InputViewType1,
          typename InputViewType2 = InputViewType1,
          typename OutputViewType = InputViewType1>
struct MultiplyKernel {
  typedef typename InputViewType1::execution_space execution_space;
  typedef typename InputViewType1::size_type size_type;
  typedef Kokkos::RangePolicy< execution_space> range_policy_type;
  typedef Kokkos::TeamPolicy< execution_space> team_policy_type;
  typedef typename team_policy_type::member_type team_handle;

  const InputViewType1 m_v1;
  const InputViewType2 m_v2;
  const OutputViewType m_v3;
  const bool m_update;

  MultiplyKernel(const InputViewType1 v1,
                 const InputViewType2 v2,
                 const OutputViewType v3,
                 const bool update) :
    m_v1(v1), m_v2(v2), m_v3(v3), m_update(update) {};

  // Multiply entries for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    if (m_update)
      m_v3(i) += m_v1(i)*m_v2(i);
    else
      m_v3(i) = m_v1(i)*m_v2(i);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_handle& team ) const
  {
    const size_type i = team.league_rank()*team.team_size() + team.team_rank();
    if (i < m_v1.extent(0))
      (*this)(i);
  }

  // Kernel launch
  static void apply(const InputViewType1 v1,
                    const InputViewType2 v2,
                    const OutputViewType v3,
                    const bool update = false) {
    const size_type nrow = v1.extent(0);

#if defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const size_type stride = Kokkos::ViewScalarStride<InputViewType1>::stride;
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<InputViewType1>::value ||
        Kokkos::is_dynrankview_fad_contiguous<InputViewType1>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const size_type stride = team_policy_type::vector_length_max(); // 32
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<InputViewType1>::value ||
        Kokkos::is_dynrankview_fad_contiguous<InputViewType1>::value ) &&
      is_dfad<typename InputViewType1::non_const_value_type>::value;
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const size_type stride = Kokkos::ViewScalarStride<InputViewType1>::stride;
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<InputViewType1>::value ||
        Kokkos::is_dynrankview_fad_contiguous<InputViewType1>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const size_type stride = team_policy_type::vector_length_max(); // 64
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<InputViewType1>::value ||
        Kokkos::is_dynrankview_fad_contiguous<InputViewType1>::value ) &&
      is_dfad<typename InputViewType1::non_const_value_type>::value;
#else
    const size_type stride = 1;
    const bool use_team = false;
#endif

    if (use_team) {
      const size_type team_size = 256 / stride;
      team_policy_type policy( (nrow+team_size-1)/team_size, team_size, stride );
      Kokkos::parallel_for( policy, MultiplyKernel(v1,v2,v3,update) );
    }
    else {
      range_policy_type policy( 0, nrow );
      Kokkos::parallel_for( policy, MultiplyKernel(v1,v2,v3,update) );
    }
  }
};

// Kernel to assign a constant to a view
template <typename ViewType>
struct ScalarAssignKernel {
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::value_type::value_type ScalarType;
  typedef Kokkos::TeamPolicy< execution_space> team_policy_type;
  typedef Kokkos::RangePolicy< execution_space> range_policy_type;
  typedef typename team_policy_type::member_type team_handle;
  static const size_type stride = Kokkos::ViewScalarStride<ViewType>::stride;

  const ViewType   m_v;
  const ScalarType m_s;

  ScalarAssignKernel(const ViewType& v, const ScalarType& s) :
    m_v(v), m_s(s) {};

  // Multiply entries for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    m_v(i) = m_s;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_handle& team ) const
  {
    const size_type i = team.league_rank()*team.team_size() + team.team_rank();
    if (i < m_v.extent(0))
      (*this)(i);
  }

  // Kernel launch
  static void apply(const ViewType& v, const ScalarType& s) {
    const size_type nrow = v.extent(0);

#if defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      is_dfad<typename ViewType::non_const_value_type>::value;
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      is_dfad<typename ViewType::non_const_value_type>::value;
#else
    const bool use_team = false;
#endif

    if (use_team) {
      const size_type team_size = 256 / stride;
      team_policy_type policy( (nrow+team_size-1)/team_size, team_size, stride );
      Kokkos::parallel_for( policy, ScalarAssignKernel(v,s) );
    }
    else {
      range_policy_type policy( 0, nrow );
      Kokkos::parallel_for( policy, ScalarAssignKernel(v,s) );
    }
  }
};

// Kernel to assign a constant to a view
template <typename ViewType, typename ScalarViewType>
struct ValueAssignKernel {
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::value_type ValueType;
  typedef Kokkos::TeamPolicy< execution_space> team_policy_type;
  typedef Kokkos::RangePolicy< execution_space> range_policy_type;
  typedef typename team_policy_type::member_type team_handle;
  typedef typename Kokkos::ThreadLocalScalarType<ViewType>::type local_scalar_type;
  static const size_type stride = Kokkos::ViewScalarStride<ViewType>::stride;

  const ViewType m_v;
  const ScalarViewType m_s;

  ValueAssignKernel(const ViewType& v, const ScalarViewType& s) :
    m_v(v), m_s(s) {};

  // Multiply entries for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    local_scalar_type s = Sacado::partition_scalar<stride>(m_s());    
    m_v(i) = s;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_handle& team ) const
  {
    const size_type i = team.league_rank()*team.team_size() + team.team_rank();
    if (i < m_v.extent(0))
      (*this)(i);
  }

  // Kernel launch
  static void apply(const ViewType& v, const ScalarViewType& s) {
    const size_type nrow = v.extent(0);

#if defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      is_dfad<typename ViewType::non_const_value_type>::value;
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      is_dfad<typename ViewType::non_const_value_type>::value;
#else
    const bool use_team = false;
#endif

    if (use_team) {
      const size_type team_size = 256 / stride;
      team_policy_type policy( (nrow+team_size-1)/team_size, team_size, stride );
      Kokkos::parallel_for( policy, ValueAssignKernel(v,s) );
    }
    else {
      range_policy_type policy( 0, nrow );
      Kokkos::parallel_for( policy, ValueAssignKernel(v,s) );
    }
  }
};

// Kernel to assign a column of a rank-2 to a rank-1
template <typename InputViewType,
          typename OutputViewType,
          typename Enabled = void>
struct AssignRank2Rank1Kernel {
  typedef typename InputViewType::execution_space execution_space;
  typedef typename InputViewType::size_type size_type;
  typedef Kokkos::TeamPolicy< execution_space> team_policy_type;
  typedef Kokkos::RangePolicy< execution_space> range_policy_type;
  typedef typename team_policy_type::member_type team_handle;
  static const size_type stride = Kokkos::ViewScalarStride<InputViewType>::stride;

  const InputViewType m_v1;
  const OutputViewType m_v2;
  const size_type m_col;

  AssignRank2Rank1Kernel(const InputViewType v1,
                         const OutputViewType v2,
                         const size_type col) :
    m_v1(v1), m_v2(v2), m_col(col) {
    static_assert( unsigned(InputViewType::rank) == 2 ,
                   "Require rank-2 input view" );
    static_assert( unsigned(OutputViewType::rank) == 1 ,
                   "Require rank-1 output view" );
  };

  // Multiply entries for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    m_v2(i) = m_v1(i,m_col);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_handle& team ) const
  {
    const size_type i = team.league_rank()*team.team_size() + team.team_rank();
    if (i < m_v1.extent(0))
      (*this)(i);
  }

  // Kernel launch
  static void apply(const InputViewType v1,
                    const OutputViewType v2,
                    const size_type col) {
    const size_type nrow = v1.extent(0);

#if defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<InputViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<InputViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<InputViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<InputViewType>::value ) &&
      is_dfad<typename InputViewType::non_const_value_type>::value;
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<InputViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<InputViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<InputViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<InputViewType>::value ) &&
      is_dfad<typename InputViewType::non_const_value_type>::value;
#else
    const bool use_team = false;
#endif

    if (use_team) {
      const size_type team_size = 256 / stride;
      team_policy_type policy( (nrow+team_size-1)/team_size, team_size, stride );
      Kokkos::parallel_for( policy, AssignRank2Rank1Kernel(v1,v2,col) );
    }
    else {
      range_policy_type policy( 0, nrow );
      Kokkos::parallel_for( policy, AssignRank2Rank1Kernel(v1,v2,col) );
    }
  }
};

// Kernel to test atomic_add
template <typename ViewType, typename ScalarViewType>
struct AtomicAddKernel {
  typedef typename ViewType::execution_space execution_space;
  typedef typename ViewType::size_type size_type;
  typedef Kokkos::TeamPolicy< execution_space> team_policy_type;
  typedef Kokkos::RangePolicy< execution_space> range_policy_type;
  typedef typename team_policy_type::member_type team_handle;
  typedef typename Kokkos::ThreadLocalScalarType<ViewType>::type local_scalar_type;
  static const size_type stride = Kokkos::ViewScalarStride<ViewType>::stride;

  const ViewType m_v;
  const ScalarViewType m_s;

  AtomicAddKernel(const ViewType& v, const ScalarViewType& s) :
    m_v(v), m_s(s) {};

  // Multiply entries for row 'i' with a value
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type i) const {
    local_scalar_type x = m_v(i);
    Kokkos::atomic_add(&(m_s()), x);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_handle& team ) const
  {
    const size_type i = team.league_rank()*team.team_size() + team.team_rank();
    if (i < m_v.extent(0))
      (*this)(i);
  }

  // Kernel launch
  static void apply(const ViewType& v, const ScalarViewType& s) {
    const size_type nrow = v.extent(0);

#if defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::Cuda>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      is_dfad<typename ViewType::non_const_value_type>::value;
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
    const bool use_team =
      std::is_same<execution_space, Kokkos::HIP>::value &&
      ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
        Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      is_dfad<typename ViewType::non_const_value_type>::value;
#else
    const bool use_team = false;
#endif

    if (use_team) {
      const size_type team_size = 256 / stride;
      team_policy_type policy( (nrow+team_size-1)/team_size, team_size, stride );
      Kokkos::parallel_for( policy, AtomicAddKernel(v,s) );
    }
    else {
      range_policy_type policy( 0, nrow );
      Kokkos::parallel_for( policy, AtomicAddKernel(v,s) );
    }
  }
};

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Size, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
    v = ViewType("view", num_rows);
#else
    const size_type fad_size = global_fad_size;
    v = ViewType("view", num_rows, fad_size+1);
#endif
  TEUCHOS_TEST_EQUALITY(v.size(), num_rows, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DeepCopy, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
    v = ViewType ("view", num_rows, num_cols);
#else
    v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
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
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadType::value_type value_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols);
#else
  const size_type fad_size = global_fad_size;
  v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
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
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadType::value_type value_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols);
#else
  const size_type fad_size = global_fad_size;
  v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
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
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols);
#else
  const size_type fad_size = global_fad_size;
  v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
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
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols);
#else
  v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant Fad
  FadType a(fad_size, 2.3456);
  for (size_type i=0; i<fad_size; ++i)
    a.fastAccessDx(i) = 7.89 + (i+1);

  // Copy to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, a);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      success = success && checkFads(a, hv(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, LocalDeepCopy, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType***,Layout,Device> ViewType;
  typedef Kokkos::View<FadType,Layout,Device> ScalarViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename ScalarViewType::HostMirror host_scalar_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type num_slices = 10;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols, num_slices);
#else
  v = ViewType ("view", num_rows, num_cols, num_slices, fad_size+1);
#endif
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant Fad to the device
  // Can't deep_copy directly because that doesn't work with DFad
  FadType a(fad_size, 2.3456);
  for (size_type i=0; i<fad_size; ++i)
    a.fastAccessDx(i) = 7.89 + (i+1);
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  ScalarViewType a_view("a");
#else
  ScalarViewType a_view("a", fad_size+1);
#endif
  host_scalar_view_type ha_view = Kokkos::create_mirror_view(a_view);
  ha_view() = a;
  Kokkos::deep_copy( a_view, ha_view );

  // Excersize local_deep_copy by setting each row of s to a
  Kokkos::parallel_for(Kokkos::RangePolicy<Device>(0,num_rows),
                       KOKKOS_LAMBDA(const int i)
  {
    auto s = Kokkos::subview(v,i,Kokkos::ALL,Kokkos::ALL);
    Kokkos::Experimental::local_deep_copy(s,a_view());
  });

  // Copy back to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, a);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      for (size_type k=0; k<num_slices; ++k) {
        success = success && checkFads(a, hv(i,j,k), out);
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, LocalDeepCopyTeam, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType***,Layout,Device> ViewType;
  typedef Kokkos::View<FadType,Layout,Device> ScalarViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename ScalarViewType::HostMirror host_scalar_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type num_slices = 10;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols, num_slices);
#else
  v = ViewType ("view", num_rows, num_cols, num_slices, fad_size+1);
#endif
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant Fad to the device
  // Can't deep_copy directly because that doesn't work with DFad
  FadType a(fad_size, 2.3456);
  for (size_type i=0; i<fad_size; ++i)
    a.fastAccessDx(i) = 7.89 + (i+1);
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  ScalarViewType a_view("a");
#else
  ScalarViewType a_view("a", fad_size+1);
#endif
  host_scalar_view_type ha_view = Kokkos::create_mirror_view(a_view);
  ha_view() = a;
  Kokkos::deep_copy( a_view, ha_view );

  // Excersize local_deep_copy by setting each row of s to a
  typedef Kokkos::TeamPolicy<Device> Policy;
  static const size_type stride = Kokkos::ViewScalarStride<ViewType>::stride;
  Kokkos::parallel_for(Policy(num_rows,Kokkos::AUTO,stride),
                       KOKKOS_LAMBDA(const typename Policy::member_type& team)
  {
    int i = team.league_rank();
    auto s = Kokkos::subview(v,i,Kokkos::ALL,Kokkos::ALL);
    Kokkos::Experimental::local_deep_copy(team,s,a_view());
  });

  // Copy back to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, a);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      for (size_type k=0; k<num_slices; ++k) {
        success = success && checkFads(a, hv(i,j,k), out);
      }
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ScalarAssign, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadType::value_type value_type;

  const size_type num_rows = global_num_rows;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows);
#else
  const size_type fad_size = global_fad_size;
  v = ViewType ("view", num_rows, fad_size+1);
#endif
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant scalar
  value_type a = 2.3456;
  ScalarAssignKernel<ViewType>::apply( v, a );

  // Copy to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, v);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
    FadType f = FadType(fad_size, a);
#else
    FadType f = a;
#endif
    success = success && checkFads(f, hv(i), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ValueAssign, FadType, Layout, Device )
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
  ScalarViewType a;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows);
  a = ScalarViewType ("fad");
#else
  v = ViewType ("view", num_rows, fad_size+1);
  a = ScalarViewType ("fad", fad_size+1);
#endif
  typename ViewType::array_type va = v;
  Kokkos::deep_copy( va, 1.0 );

  // Deep copy a constant scalar
  Kokkos::deep_copy(a, 2.3456);

  Kokkos::parallel_for(Kokkos::RangePolicy< Device>(0, fad_size), KOKKOS_LAMBDA(const int i) {
    a().fastAccessDx(i) = 7.89 + i;
  });
  Kokkos::fence();

  ValueAssignKernel<ViewType, ScalarViewType>::apply( v, a );
  Kokkos::fence();

  // Copy to host
  host_view_type hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, v);

  host_scalar_view_type ha = Kokkos::create_mirror_view(a);
  Kokkos::deep_copy(ha, a);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    success = success && checkFads(ha(), hv(i), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Resize, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
    v = ViewType ("view", num_rows, num_cols);
#else
    v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_v(i,j) = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
  Kokkos::deep_copy(v, h_v);

  // Resize
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  Kokkos::resize(v, num_rows, num_cols+1);
#else
  Kokkos::resize(v, num_rows, num_cols+1, fad_size+1);
#endif

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
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
    FadType f = 0.0;
#else
    FadType f(fad_size, 0.0);
#endif
    success = success && checkFads(f, h_v2(i,num_cols), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Multiply, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  ViewType v1, v2;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v1 = ViewType ("view1", num_rows);
  v2 = ViewType ("view2", num_rows);
#else
  v1 = ViewType ("view1", num_rows, fad_size+1);
  v2 = ViewType ("view2", num_rows, fad_size+1);
#endif
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
  ViewType v3;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v3 = ViewType ("view3", num_rows);
#else
  v3 = ViewType ("view3", num_rows, fad_size+1);
#endif
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
  Kokkos_View_Fad, MultiplyUpdate, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  ViewType v1, v2;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v1 = ViewType ("view1", num_rows);
  v2 = ViewType ("view2", num_rows);
#else
  v1 = ViewType ("view1", num_rows, fad_size+1);
  v2 = ViewType ("view2", num_rows, fad_size+1);
#endif
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
  ViewType v3;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v3 = ViewType ("view3", num_rows);
#else
  v3 = ViewType ("view3", num_rows, fad_size+1);
#endif
  Kokkos::deep_copy(v3, 1.0);
  MultiplyKernel<ViewType>::apply(v1,v2,v3,true);

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
    FadType f3 = 1.0 + f1*f2;
    success = success && checkFads(f3, h_v3(i), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, MultiplyConst, FadType, Layout, Device )
{
  typedef Kokkos::View<const FadType*,Layout,Device,Kokkos::MemoryUnmanaged> ConstViewType;
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  ViewType v1, v2;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v1 = ViewType ("view1", num_rows);
  v2 = ViewType ("view2", num_rows);
#else
  v1 = ViewType ("view1", num_rows, fad_size+1);
  v2 = ViewType ("view2", num_rows, fad_size+1);
#endif
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

  // Launch kernel
  ViewType v3;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v3 = ViewType ("view3", num_rows);
#else
  v3 = ViewType ("view3", num_rows, fad_size+1);
#endif
  MultiplyKernel<ConstViewType,ViewType,ViewType>::apply(cv1,v2,v3);

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
  Kokkos_View_Fad, MultiplyMixed, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = 2;
  const size_type fad_size = global_fad_size;

  // Create and fill views -- do everything on the host for this test
  FadType f0 = generate_fad<FadType>(
    num_rows, size_type(2), fad_size, size_type(0), size_type(0));
  FadType f1 = generate_fad<FadType>(
    num_rows, size_type(2), fad_size, size_type(1), size_type(0));
  host_view_type h_v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  h_v = host_view_type ("view1", num_rows);
#else
  h_v = host_view_type ("view1", num_rows, fad_size+1);
#endif
  h_v(0) = f0;
  h_v(1) = f1;

  FadType f2 = f0 * h_v(1);

  // Check
  FadType f3 = f0 * f1;
  success = checkFads(f3, f2, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AtomicAdd, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef Kokkos::View<FadType,Layout,Device> ScalarViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ScalarViewType::HostMirror host_scalar_view_type;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows);
#else
  v = ViewType ("view", num_rows, fad_size+1);
#endif
  Kokkos::deep_copy(v, 2.3456);

  Kokkos::parallel_for(Kokkos::RangePolicy<Device>(0, num_rows), KOKKOS_LAMBDA(const size_type i) {
    for (size_type j = 0; j < fad_size; ++j)
      v(i).fastAccessDx(j) = 7.89 + j;
  });

  // Create scalar view
  ScalarViewType s;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  s = ScalarViewType ("scalar view");
#else
  s = ScalarViewType ("scalar view", fad_size+1);
#endif

  //  Call atomic_add kernel, which adds up entries in v
  AtomicAddKernel<ViewType,ScalarViewType>::apply( v, s );

  // Copy to host
  host_scalar_view_type hs = Kokkos::create_mirror_view(s);
  Kokkos::deep_copy(hs, s);

  // Check
  auto hv = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hv, v);

  FadType b = num_rows*hv(0);
  success = checkFads(b, hs(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Rank8, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*******,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", 100, 1, 2, 3, 4, 5, 6);
#else
  v = ViewType ("view", 100, 1, 2, 3, 4, 5, 6, fad_size+1);
#endif
  host_view_type h_v = Kokkos::create_mirror_view(v);
  typename host_view_type::array_type h_a = h_v;
  Kokkos::deep_copy(h_a, 1.0);

  FadType f1 = FadType(fad_size, 2.0);
  h_v(99,0,1,2,3,4,5) = f1;
  FadType f2 = h_v(99,0,1,2,3,4,5);

  // Check
  success = checkFads(f1, f2, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Roger, FadType, Layout, Device )
{
  Kokkos::View<FadType*,Layout,Device> a;
  Kokkos::View<FadType**,Layout,Device> b;
  Kokkos::View<FadType***,Layout,Device> c;
  Kokkos::View<FadType****,Layout,Device> d;
  Kokkos::View<FadType*****,Layout,Device> e;
  Kokkos::View<FadType******,Layout,Device> f;
  Kokkos::View<FadType*******,Layout,Device> g;

#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  a = Kokkos::View<FadType*,Layout,Device>("a",4);
  b = Kokkos::View<FadType**,Layout,Device> ("b",4,4);
  c = Kokkos::View<FadType***,Layout,Device> ("c",4,4,4);
  d = Kokkos::View<FadType****,Layout,Device> ("d",4,4,4,4);
  e = Kokkos::View<FadType*****,Layout,Device> ("e",4,4,4,4,4);
  f = Kokkos::View<FadType******,Layout,Device> ("f",4,4,4,4,4,4);
  g = Kokkos::View<FadType*******,Layout,Device> ("g",4,4,4,4,4,4,4);
#else
  const unsigned fad_size = global_fad_size;
  a = Kokkos::View<FadType*,Layout,Device>("a",4,fad_size+1);
  b = Kokkos::View<FadType**,Layout,Device> ("b",4,4,fad_size+1);
  c = Kokkos::View<FadType***,Layout,Device> ("c",4,4,4,fad_size+1);
  d = Kokkos::View<FadType****,Layout,Device> ("d",4,4,4,4,fad_size+1);
  e = Kokkos::View<FadType*****,Layout,Device> ("e",4,4,4,4,4,fad_size+1);
  f = Kokkos::View<FadType******,Layout,Device> ("f",4,4,4,4,4,4,fad_size+1);
  g = Kokkos::View<FadType*******,Layout,Device> ("g",4,4,4,4,4,4,4,fad_size+1);
#endif

  typedef typename Device::memory_space memory_space;
  const bool is_accessible =
    Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                    memory_space>::accessible;
  if (is_accessible) {
    a(0) = FadType(1.0);
    f(0,0,0,0,0,0) = FadType(1.0);
    g(0,0,0,0,0,0,0) = FadType(1.0);
  }

  // Check
  success = true;
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AssignDifferentStrides, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Layout,Device> ViewType1;
  typedef Kokkos::View<FadType*,Layout,Device> ViewType2;
  typedef typename ViewType1::size_type size_type;
  typedef typename ViewType1::HostMirror host_view_type1;
  typedef typename ViewType2::HostMirror host_view_type2;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  ViewType1 v1;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v1 = ViewType1 ("view1", num_rows, num_cols);
#else
  v1 = ViewType1 ("view1", num_rows, num_cols, fad_size+1);
#endif
  host_view_type1 h_v1 = Kokkos::create_mirror_view(v1);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      h_v1(i,j) = generate_fad<FadType>(
        num_rows, num_cols, fad_size, i, j);
    }
  }
  Kokkos::deep_copy(v1, h_v1);

  // Launch kernel
  ViewType2 v2;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v2 = ViewType2 ("view2", num_rows);
#else
  v2 = ViewType2 ("view2", num_rows, fad_size+1);
#endif
  AssignRank2Rank1Kernel<ViewType1,ViewType2>::apply(v1,v2,1);

  // Copy back
  host_view_type2 h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    FadType f =
      generate_fad<FadType>(num_rows, num_cols, fad_size, i, size_type(1));
    success = success && checkFads(f, h_v2(i), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ScalarValue, FadType, Layout, Device )
{
  typedef typename Sacado::ScalarType<FadType>::type ScalarType;
  typedef Kokkos::View<FadType,Layout,Device> ViewType1;
  typedef Kokkos::View<ScalarType,Layout,Device> ViewType2;
  typedef typename ViewType1::size_type size_type;
  typedef typename ViewType1::HostMirror host_view_type1;
  typedef typename ViewType2::HostMirror host_view_type2;

  const int fad_size = global_fad_size;

  // Create and fill views
  ViewType1 v1;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v1 = ViewType1 ("view1");
#else
  v1 = ViewType1 ("view1", fad_size+1);
#endif
  host_view_type1 h_v1 = Kokkos::create_mirror_view(v1);
  h_v1() = generate_fad<FadType>(1, 1, fad_size, 0, 0);
  Kokkos::deep_copy(v1, h_v1);

  // Launch kernel
  ViewType2 v2 = ViewType2 ("view2");
  Kokkos::parallel_for(Kokkos::RangePolicy<Device>(0,1),
                       KOKKOS_LAMBDA(const size_type i)
  {
    v2() = Sacado::scalarValue(v1());
  });

  // Copy back
  host_view_type2 h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check
  success = true;
  TEUCHOS_TEST_EQUALITY(h_v1().val(), h_v2(), out, success);
}

#if defined(HAVE_SACADO_KOKKOS) && defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankDimensionScalar, FadType, Layout, Device )
{
  typedef Kokkos::DynRankView<double,Layout,Device> DoubleViewType;
  typedef Kokkos::DynRankView<FadType,Layout,Device> FadViewType;
  typedef typename FadViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create views
  DoubleViewType v1("view1", num_rows);
  FadViewType v2 ("view2", num_rows, fad_size+1);

  // Check dimension scalar works
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(v1), 0, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(v2), fad_size+1, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankAssignStatic0, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType,Layout,Device> StaticViewType;
  typedef Kokkos::DynRankView<FadType,Layout,Device> DynamicViewType;
  typedef typename StaticViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  StaticViewType v1("view", fad_size+1);
  auto h_v1 = Kokkos::create_mirror_view(v1);
  h_v1() = generate_fad<FadType>(num_rows, num_cols, fad_size, size_type(0), size_type(0));
  Kokkos::deep_copy(v1, h_v1);

  // Assign static to dynamic
  DynamicViewType v2 = v1;

  // Copy back
  auto h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check dimensions are correct
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(v2), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(v2.stride_0(), v1.stride_0(), out, success);

  // Check values
  FadType f =
    generate_fad<FadType>(num_rows, num_cols, fad_size, size_type(0), size_type(0));
  success = success && checkFads(f, h_v2(), out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankAssignStatic1, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*,Layout,Device> StaticViewType;
  typedef Kokkos::DynRankView<FadType,Layout,Device> DynamicViewType;
  typedef typename StaticViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  StaticViewType v1("view", num_rows, fad_size+1);
  auto h_v1 = Kokkos::create_mirror_view(v1);
  for (size_type i=0; i<num_rows; ++i)
    h_v1(i) =
      generate_fad<FadType>(num_rows, num_cols, fad_size, i, size_type(0));
  Kokkos::deep_copy(v1, h_v1);

  // Assign static to dynamic
  DynamicViewType v2 = v1;

  // Copy back
  auto h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check dimensions are correct
  TEUCHOS_TEST_EQUALITY(v2.extent(0), num_rows, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(v2), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(v2.stride_0(), v1.stride_0(), out, success);
  TEUCHOS_TEST_EQUALITY(v2.stride_1(), v1.stride_1(), out, success);

  // Check values
  for (size_type i=0; i<num_rows; ++i) {
    FadType f =
      generate_fad<FadType>(num_rows, num_cols, fad_size, i, size_type(0));
    success = success && checkFads(f, h_v2(i), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankAssignStatic2, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Layout,Device> StaticViewType;
  typedef Kokkos::DynRankView<FadType,Layout,Device> DynamicViewType;
  typedef typename StaticViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill views
  StaticViewType v1("view", num_rows, num_cols, fad_size+1);
  auto h_v1 = Kokkos::create_mirror_view(v1);
  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_v1(i,j) = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
  Kokkos::deep_copy(v1, h_v1);

  // Assign static to dynamic
  DynamicViewType v2 = v1;

  // Copy back
  auto h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check dimensions are correct
  TEUCHOS_TEST_EQUALITY(v2.extent(0), num_rows, out, success);
  TEUCHOS_TEST_EQUALITY(v2.extent(1), num_cols, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(v2), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(v2.stride_0(), v1.stride_0(), out, success);
  TEUCHOS_TEST_EQUALITY(v2.stride_1(), v1.stride_1(), out, success);
  TEUCHOS_TEST_EQUALITY(v2.stride_2(), v1.stride_2(), out, success);

  // Check values
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      success = success && checkFads(f, h_v2(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankMultiply, FadType, Layout, Device )
{
  typedef Kokkos::DynRankView<FadType,Layout,Device> ViewType;
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

  // Launch kernel
  ViewType v3("view3", num_rows, fad_size+1);
  MultiplyKernel<ViewType>::apply(v1,v2,v3);

  // Copy back
  host_view_type h_v3 = Kokkos::create_mirror_view(v3);
  Kokkos::deep_copy(h_v3, v3);

  // Check
  success = true;
  TEUCHOS_TEST_EQUALITY(v3.rank(), 1, out, success);
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
  Kokkos_View_Fad, SubdynrankviewCol, FadType, Layout, Device )
{
  typedef Kokkos::DynRankView<FadType,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j) = f;
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Create subview of first column
  size_type col = 1;
  auto s = Kokkos::subdynrankview(v, Kokkos::ALL(), col);

  // Copy back
  typedef decltype(s) SubviewType;
  typedef typename SubviewType::HostMirror HostSubviewType;

  // Note:  don't create h_s through create_mirror_view and deep_copy
  // since Kokkos doesn't support deep_copy of non-contiguous views
  //HostSubviewType h_s = Kokkos::create_mirror_view(s);
  //Kokkos::deep_copy(h_s, s);
  HostSubviewType h_s = Kokkos::subdynrankview(h_v, Kokkos::ALL(), col);

  // Check
  success = true;
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(s), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(h_s), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(h_s.extent(0), num_rows, out, success);
  TEUCHOS_TEST_EQUALITY(h_s.extent(1), 1, out, success);
  TEUCHOS_TEST_EQUALITY(h_s.extent(7), 1, out, success);

  for (size_type i=0; i<num_rows; ++i) {
    FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, col);
    success = success && checkFads(f, h_s(i), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, SubdynrankviewRow, FadType, Layout, Device )
{
  typedef Kokkos::DynRankView<FadType,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type num_planes = 9;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, num_planes, fad_size+1);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      for (size_type k=0; k<num_planes; ++k) {
        h_v(i,j,k) = (k+1)*f;
      }
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Create subview of first column
  size_type row = 2;
  auto s = Kokkos::subdynrankview(v, row, Kokkos::ALL(),  Kokkos::ALL());

  // Copy back
  typedef decltype(s) SubviewType;
  typedef typename SubviewType::HostMirror HostSubviewType;

  // Note:  don't create h_s through create_mirror_view and deep_copy
  // since Kokkos doesn't support deep_copy of non-contiguous views
  //HostSubviewType h_s = Kokkos::create_mirror_view(s);
  //Kokkos::deep_copy(h_s, s);
  HostSubviewType h_s =
    Kokkos::subdynrankview(h_v, row, Kokkos::ALL(),  Kokkos::ALL());

  // Check
  success = true;
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(s), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(h_s), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(h_s.extent(0), num_cols, out, success);
  TEUCHOS_TEST_EQUALITY(h_s.extent(1), num_planes, out, success);
  TEUCHOS_TEST_EQUALITY(h_s.extent(2), 1, out, success);
  TEUCHOS_TEST_EQUALITY(h_s.extent(7), 1, out, success);

  for (size_type j=0; j<num_cols; ++j) {
    FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, row, j);
    for (size_type k=0; k<num_planes; ++k) {
      FadType g = (k+1)*f;
      success = success && checkFads(g, h_s(j,k), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, SubdynrankviewScalar, FadType, Layout, Device )
{
  typedef Kokkos::DynRankView<FadType,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols, fad_size+1);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j) = f;
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Create subview of first column
  size_type row = 3;
  size_type col = 1;
  auto s = Kokkos::subdynrankview(v, row, col);

  // Copy back
  typedef decltype(s) SubviewType;
  typedef typename SubviewType::HostMirror HostSubviewType;

  // Note:  don't create h_s through create_mirror_view and deep_copy
  // since Kokkos doesn't support deep_copy of non-contiguous views
  //HostSubviewType h_s = Kokkos::create_mirror_view(s);
  //Kokkos::deep_copy(h_s, s);
  HostSubviewType h_s = Kokkos::subdynrankview(h_v, row, col);

  // Check
  success = true;
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(s), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(h_s), fad_size+1, out, success);
  FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, row, col);
  success = success && checkFads(f, h_s(), out);
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankDimensionScalar, FadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankAssignStatic0, FadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankAssignStatic1, FadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankAssignStatic2, FadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, DynRankMultiply, FadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, SubdynrankviewCol, FadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, SubdynrankviewRow, FadType, Layout, Device ) {}
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, SubdynrankviewScalar, FadType, Layout, Device ) {}

#endif

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Subview, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols);
#else
  v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j) = f;
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Create subview of first column
  size_type col = 1;
  auto s = Kokkos::subview(v, Kokkos::ALL(), col);

  // Copy back
  typedef decltype(s) SubviewType;
  typedef typename SubviewType::HostMirror HostSubviewType;

  // Note:  don't create h_s through create_mirror_view and deep_copy
  // since Kokkos doesn't support deep_copy of non-contiguous views
  //HostSubviewType h_s = Kokkos::create_mirror_view(s);
  //Kokkos::deep_copy(h_s, s);
  HostSubviewType h_s =  Kokkos::subview(h_v, Kokkos::ALL(), col);

  // Check
  success = true;
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(s), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(h_s), fad_size+1, out, success);
#endif
  for (size_type i=0; i<num_rows; ++i) {
    FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, col);
    success = success && checkFads(f, h_s(i), out);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Subview2, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType***,Layout,Device> ViewType;
  typedef typename ViewType::HostMirror host_view_type;

  // Test various subview operations to check the resulting indexing is correct.
  // We only need to run these tests on the host because the indexing does
  // not depend on the device.

  const int num_cell = 5;
  const int num_qp = 4;
  const int num_dim = 3;
  const int num_deriv = 2;

  // Create and fill view
  host_view_type v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = host_view_type ("view", num_cell, num_qp, num_dim);
#else
  v = host_view_type ("view", num_cell, num_qp, num_dim, num_deriv+1);
#endif
  for (int cell=0; cell < num_cell; ++cell) {
    for (int qp=0; qp < num_qp; ++qp) {
      for (int dim = 0; dim < num_dim; ++dim) {
        v(cell,qp,dim).val() = 100.*cell + 10.*qp + 1.*dim;
        for (int deriv = 0; deriv < num_deriv; ++deriv) {
          v(cell,qp,dim).fastAccessDx(deriv) = v(cell,qp,dim).val() + (1.0*deriv)/10.;
        }
      }
    }
  }

  success = true;

  out << "checking subview(v,ALL,*,*)..." << std::endl;
  for (int qp=0; qp < num_qp; ++qp) {
    for (int dim=0; dim < num_dim; ++dim) {
      auto v_tmp = subview(v,Kokkos::ALL(),qp,dim);
      for (int cell=0; cell < num_cell; ++cell) {
        out << "\tChecking (" << cell << "," << qp << "," << dim << ")" << std::endl;
        success = success && checkFads(v(cell,qp,dim), v_tmp(cell), out);
      }
    }
  }

  out << "checking subview(v,*,ALL,*)..." << std::endl;
  for (int cell=0; cell < num_cell; ++cell) {
    for (int dim=0; dim < num_dim; ++dim) {
      auto v_tmp = subview(v,cell,Kokkos::ALL(),dim);
      for (int qp=0; qp < num_qp; ++qp) {
        out << "\tChecking (" << cell << "," << qp << "," << dim << ")" << std::endl;
        success = success && checkFads(v(cell,qp,dim), v_tmp(qp), out);
      }
    }
  }

  out << "checking subview(v,*,*,ALL)..." << std::endl;
  for (int cell=0; cell < num_cell; ++cell) {
    for (int qp=0; qp < num_qp; ++qp) {
      auto v_tmp = subview(v,cell,qp,Kokkos::ALL());
      for (int dim=0; dim < num_dim; ++dim) {
        out << "\tChecking (" << cell << "," << qp << "," << dim << ")" << std::endl;
        success = success && checkFads(v(cell,qp,dim), v_tmp(dim), out);
      }
    }
  }

  out << "checking subview(v,ALL,ALL,*)..." << std::endl;
  for (int dim=0; dim < num_dim; ++dim) {
    auto v_tmp = subview(v,Kokkos::ALL(),Kokkos::ALL(),dim);
    for (int cell=0; cell < num_cell; ++cell) {
      for (int qp=0; qp < num_qp; ++qp) {
        out << "\tChecking (" << cell << "," << qp << "," << dim << ")" << std::endl;
        success = success && checkFads(v(cell,qp,dim), v_tmp(cell,qp), out);
      }
    }
  }

  out << "checking subview(v,*,ALL,ALL)..." << std::endl;
  for (int cell=0; cell < num_cell; ++cell) {
    auto v_tmp = subview(v,cell,Kokkos::ALL(),Kokkos::ALL());
    for (int qp=0; qp < num_qp; ++qp) {
      for (int dim=0; dim < num_dim; ++dim) {
        out << "\tChecking (" << cell << "," << qp << "," << dim << ")" << std::endl;
        success = success && checkFads(v(cell,qp,dim), v_tmp(qp,dim), out);
      }
    }
  }

  out << "checking subview(v,ALL,*,ALL)..." << std::endl;
  for (int qp=0; qp < num_qp; ++qp) {
    auto v_tmp = subview(v,Kokkos::ALL(),qp,Kokkos::ALL());
    for (int cell=0; cell < num_cell; ++cell) {
      for (int dim=0; dim < num_dim; ++dim) {
        out << "\tChecking (" << cell << "," << qp << "," << dim << ")" << std::endl;
        success = success && checkFads(v(cell,qp,dim), v_tmp(cell,dim), out);
      }
    }
  }

  out << "checking subview(v,range,range,range)..." << std::endl;
  auto v_tmp = subview(v,std::make_pair(1,5),std::make_pair(1,4),std::make_pair(1,3));
  for (int cell=1; cell < num_cell; ++cell) {
    for (int qp=1; qp < num_qp; ++qp) {
      for (int dim=1; dim < num_dim; ++dim) {
        out << "\tChecking (" << cell << "," << qp << "," << dim << ")" << std::endl;
        success = success && checkFads(v(cell,qp,dim), v_tmp(cell-1,qp-1,dim-1), out);
      }
    }
  }
}

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ConstViewAssign, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType*,Layout,Device> ViewType;
  typedef Kokkos::View<const FadType,Layout,Device> ConstViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename ViewType::execution_space exec_space;

  const size_type num_rows = global_num_rows;
  const size_type fad_size = global_fad_size;

  // Create and fill view
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  ViewType v1("view1", num_rows);
#else
  ViewType v1("view1", num_rows, fad_size+1);
#endif
  host_view_type h_v1 = Kokkos::create_mirror_view(v1);
  for (size_type i=0; i<num_rows; ++i) {
    FadType f = generate_fad<FadType>(num_rows, size_type(1), fad_size, i,
                                      size_type(0));
    h_v1(i) = f;
  }
  Kokkos::deep_copy(v1, h_v1);

#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  ViewType v2("view2", num_rows);
#else
  ViewType v2("view2", num_rows, fad_size+1);
#endif

  static const size_type stride = Kokkos::ViewScalarStride<ViewType>::stride;
#if defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
  const bool use_team =
    std::is_same<exec_space, Kokkos::Cuda>::value &&
    ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
      Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_CUDA) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
  const bool use_team =
    std::is_same<exec_space, Kokkos::Cuda>::value &&
    ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
      Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
    is_dfad<typename ViewType::non_const_value_type>::value;
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL)
  const bool use_team =
    std::is_same<exec_space, Kokkos::HIP>::value &&
    ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
      Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
      ( stride > 1 );
#elif defined (KOKKOS_ENABLE_HIP) && defined (SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
  const bool use_team =
    std::is_same<exec_space, Kokkos::HIP>::value &&
    ( Kokkos::is_view_fad_contiguous<ViewType>::value ||
      Kokkos::is_dynrankview_fad_contiguous<ViewType>::value ) &&
    is_dfad<typename ViewType::non_const_value_type>::value;
#else
  const bool use_team = false;
#endif

  if (use_team) {
    typedef Kokkos::TeamPolicy<exec_space> team_policy;
    Kokkos::parallel_for(team_policy(num_rows, 1, stride),
                         KOKKOS_LAMBDA(typename team_policy::member_type team)
    {
      const int i = team.league_rank();
      typename ConstViewType::reference_type x = v1(i);
      v2(i) = x;
    });
  }
  else {
    Kokkos::parallel_for(Kokkos::RangePolicy<exec_space>(0,num_rows),
                         KOKKOS_LAMBDA(const int i)
    {
      typename ConstViewType::reference_type x = v1(i);
      v2(i) = x;
    });
  }

  // Copy back
  host_view_type h_v2 = Kokkos::create_mirror_view(v2);
  Kokkos::deep_copy(h_v2, v2);

  // Check
  success = true;
  for (size_type i=0; i<num_rows; ++i) {
    FadType f = generate_fad<FadType>(num_rows, size_type(1), fad_size, i,
                                      size_type(0));
    success = success && checkFads(f, h_v2(i), out);
  }
}
#else
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ConstViewAssign, FadType, Layout, Device ) {}
#endif

// Tests that require view spec

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ShmemSize, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename FadType::value_type value_type;
  typedef typename ViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Compute shared memory size for View
  const size_type shmem_size =
    ViewType::shmem_size(num_rows, num_cols, fad_size+1);

  // Check
  const size_type align = 8;
  const size_type mask  = align - 1;
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols);
#else
  v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
  const size_type shmem_size_expected =
    (( sizeof(value_type) * global_num_rows * global_num_cols * (fad_size+1) + mask ) & ~mask) + sizeof(typename ViewType::traits::value_type);
  TEUCHOS_TEST_EQUALITY(shmem_size, shmem_size_expected, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Unmanaged, FadType, Layout, Device )
{
  // For LayoutContiguous or LayoutNatural, strip out the layout they are templated on
  typedef typename Kokkos::inner_layout<Layout>::type TestLayout;

  typedef typename FadType::value_type scalar_type;
  typedef Kokkos::View<scalar_type***,TestLayout,Device> ViewType;
  typedef Kokkos::View<FadType**,TestLayout,Device,Kokkos::MemoryUnmanaged> FadViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadViewType::HostMirror fad_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
  host_view_type h_v;
  if (Kokkos::is_view_fad_contiguous<FadViewType>::value &&
      std::is_same<TestLayout, Kokkos::LayoutLeft >::value) {
    v = ViewType ("view", fad_size+1, num_rows, num_cols);
    h_v = Kokkos::create_mirror_view(v);
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
        for (size_type k=0; k<fad_size; k++)
          h_v(k,i,j) = f.dx(k);
        h_v(fad_size,i,j) = f.val();
      }
    }
  }
  else {
    v = ViewType ("view", num_rows, num_cols, fad_size+1);
    h_v = Kokkos::create_mirror_view(v);
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
        for (size_type k=0; k<fad_size; k++)
          h_v(i,j,k) = f.dx(k);
        h_v(i,j,fad_size) = f.val();
      }
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Create unmanaged view
  FadViewType v_fad;//( v.data(), num_rows, num_cols, fad_size+1);
  fad_host_view_type h_v_fad;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v_fad = FadViewType ( v.data(), num_rows, num_cols);
  h_v_fad = fad_host_view_type ("host_view_fad", num_rows, num_cols);
#else
  v_fad = FadViewType ( v.data(), num_rows, num_cols, fad_size+1);
  h_v_fad = fad_host_view_type ("host_view_fad", num_rows, num_cols, fad_size+1);
#endif

  // Copy back -- can't use create_mirror_view() because v_fad is unmanaged
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
  // For LayoutContiguous or LayoutNatural, strip out the layout they are templated on
  typedef typename Kokkos::inner_layout<Layout>::type TestLayout;

  typedef typename FadType::value_type scalar_type;
  typedef Kokkos::View<scalar_type***,TestLayout,Device> ViewType;
  typedef Kokkos::View<FadType**,TestLayout,Device> FadViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadViewType::HostMirror fad_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
  host_view_type h_v;
  if (Kokkos::is_view_fad_contiguous<FadViewType>::value &&
      std::is_same<TestLayout, Kokkos::LayoutLeft >::value) {
    v = ViewType ("view", fad_size+1, num_rows, num_cols);
    h_v = Kokkos::create_mirror_view(v);
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
        for (size_type k=0; k<fad_size; k++)
          h_v(k,i,j) = f.dx(k);
        h_v(fad_size,i,j) = f.val();
      }
    }
  }
  else {
    v = ViewType ("view", num_rows, num_cols, fad_size+1);
    h_v = Kokkos::create_mirror_view(v);
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
        for (size_type k=0; k<fad_size; k++)
          h_v(i,j,k) = f.dx(k);
        h_v(i,j,fad_size) = f.val();
      }
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Create unmanaged view
  FadViewType v_fad;//( v.data(), num_rows, num_cols, fad_size+1);
  fad_host_view_type h_v_fad;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v_fad = FadViewType ( v.data(), num_rows, num_cols);
  h_v_fad = fad_host_view_type ("host_view_fad", num_rows, num_cols);
#else
  v_fad = FadViewType ( v.data(), num_rows, num_cols, fad_size+1);
  h_v_fad = fad_host_view_type ("host_view_fad", num_rows, num_cols, fad_size+1);
#endif

  // Copy back -- can't use create_mirror_view() because v_fad is unmanaged
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
  // For LayoutContiguous or LayoutNatural, strip out the layout they are templated on
  typedef typename Kokkos::inner_layout<Layout>::type TestLayout;

  typedef typename FadType::value_type scalar_type;
  typedef Kokkos::View<scalar_type***,TestLayout,Device> ViewType;
  typedef Kokkos::View<const scalar_type***,TestLayout,Device> ConstViewType;
  typedef Kokkos::View<FadType**,TestLayout,Device,Kokkos::MemoryUnmanaged> FadViewType;
  typedef Kokkos::View<const FadType**,TestLayout,Device,Kokkos::MemoryUnmanaged> ConstFadViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadViewType::HostMirror fad_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
  host_view_type h_v;
  if (Kokkos::is_view_fad_contiguous<FadViewType>::value &&
      std::is_same<TestLayout, Kokkos::LayoutLeft >::value) {
    v = ViewType ("view", fad_size+1, num_rows, num_cols);
    h_v = Kokkos::create_mirror_view(v);
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
        for (size_type k=0; k<fad_size; k++)
          h_v(k,i,j) = f.dx(k);
        h_v(fad_size,i,j) = f.val();
      }
    }
  }
  else {
    v = ViewType ("view", num_rows, num_cols, fad_size+1);
    h_v = Kokkos::create_mirror_view(v);
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
        for (size_type k=0; k<fad_size; k++)
          h_v(i,j,k) = f.dx(k);
        h_v(i,j,fad_size) = f.val();
      }
    }
  }
  Kokkos::deep_copy(v, h_v);
  ConstViewType v_const = v;

  // Create unmanaged view

  ConstFadViewType v_fad;//( v.data(), num_rows, num_cols, fad_size+1);
  fad_host_view_type h_v_fad;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v_fad = ConstFadViewType ( v_const.data(), num_rows, num_cols);
  h_v_fad = fad_host_view_type ("host_view_fad", num_rows, num_cols);
#else
  v_fad = ConstFadViewType ( v_const.data(), num_rows, num_cols, fad_size+1);
  h_v_fad = fad_host_view_type ("host_view_fad", num_rows, num_cols, fad_size+1);
#endif

  // Copy back -- can't use create_mirror_view() because v_fad is unmanaged
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
  // For LayoutContiguous or LayoutNatural, strip out the layout they are templated on
  typedef typename Kokkos::inner_layout<Layout>::type TestLayout;
  typedef typename FadType::value_type scalar_type;
  typedef Kokkos::View<scalar_type***,TestLayout,Device> ViewType;
  typedef Kokkos::View<const scalar_type***,TestLayout,Device> ConstViewType;
  typedef Kokkos::View<FadType**,TestLayout,Device> FadViewType;
  typedef Kokkos::View<const FadType**,TestLayout,Device> ConstFadViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename FadViewType::HostMirror fad_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
  host_view_type h_v;
  if (Kokkos::is_view_fad_contiguous<FadViewType>::value &&
      std::is_same<TestLayout, Kokkos::LayoutLeft >::value) {
    v = ViewType ("view", fad_size+1, num_rows, num_cols);
    h_v = Kokkos::create_mirror_view(v);
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
        for (size_type k=0; k<fad_size; k++)
          h_v(k,i,j) = f.dx(k);
        h_v(fad_size,i,j) = f.val();
      }
    }
  }
  else {
    v = ViewType ("view", num_rows, num_cols, fad_size+1);
    h_v = Kokkos::create_mirror_view(v);
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
        for (size_type k=0; k<fad_size; k++)
          h_v(i,j,k) = f.dx(k);
        h_v(i,j,fad_size) = f.val();
      }
    }
  }
  Kokkos::deep_copy(v, h_v);
  ConstViewType v_const = v;

  // Create unmanaged view
  ConstFadViewType v_fad;
  fad_host_view_type h_v_fad;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v_fad = ConstFadViewType (v_const.data(), num_rows, num_cols);
  h_v_fad = fad_host_view_type ("host_view_fad", num_rows, num_cols);
#else
  v_fad = ConstFadViewType (v_const.data(), num_rows, num_cols, fad_size+1);
  h_v_fad = fad_host_view_type ("host_view_fad", num_rows, num_cols, fad_size+1);
#endif

  // Copy back -- can't use create_mirror_view() because v_fad is unmanaged
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

// This test checks we can allocate a view
// with SFad without specifying the fad size in the constructor
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, SFadNoSizeArg, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v("view", num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j) = f;
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Copy back
  Kokkos::deep_copy(h_v, v);

  // Check
  success = true;
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(v), fad_size+1, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(h_v), fad_size+1, out, success);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      success = success && checkFads(f, h_v(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Partition, FadType, Layout, Device )
{
#if !defined(SACADO_VIEW_CUDA_HIERARCHICAL) && !defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ViewType ("view", num_rows, num_cols);
#else
  v = ViewType ("view", num_rows, num_cols, fad_size+1);
#endif
  host_view_type h_v = Kokkos::create_mirror_view(v);

  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j) = f;
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Copy back
  Kokkos::deep_copy(h_v, v);

  // Partition derivative array of h_v into 2, first one starting at index 0,
  // the second at 1
  const size_type stride = 2;
  auto h_v1 = Kokkos::partition<2>(h_v, 0, stride);
  auto h_v2 = Kokkos::partition<2>(h_v, 1, stride);

  // Check
  const size_type fad_size_1 = (fad_size + stride - 0 - 1) / stride;
  const size_type fad_size_2 = (fad_size + stride - 1 - 1) / stride;
  success = true;
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(h_v1), fad_size_1+1, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(h_v2), fad_size_2+1, out, success);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      Sacado::Fad::DFad<double> f1( fad_size_1, f.val() );
      Sacado::Fad::DFad<double> f2( fad_size_2, f.val() );
      for (unsigned int k=0; k<fad_size_1; ++k)
        if (2*k < fad_size) f1.fastAccessDx(k) = f.dx(2*k);
      for (unsigned int k=0; k<fad_size_2; ++k)
        if (2*k+1 < fad_size) f2.fastAccessDx(k) = f.dx(2*k+1);
      success = success && checkFads(f1, h_v1(i,j), out);
      success = success && checkFads(f2, h_v2(i,j), out);
    }
  }
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AssignLayoutContiguousToLayoutStride, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Kokkos::LayoutContiguous<Layout>,Device> ContViewType;
  typedef Kokkos::View<FadType**,Kokkos::LayoutStride,Device> StrideViewType;
  typedef typename ContViewType::size_type size_type;
  typedef typename ContViewType::HostMirror cont_host_view_type;
  typedef typename StrideViewType::HostMirror stride_host_view_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create and fill view
  ContViewType v;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v = ContViewType ("view", num_rows, num_cols);
#else
  v = ContViewType ("view", num_rows, num_cols, fad_size+1);
#endif
  cont_host_view_type h_v = Kokkos::create_mirror_view(v);

  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      h_v(i,j) = f;
    }
  }
  Kokkos::deep_copy(v, h_v);

  // Assign to LayoutStride view
  StrideViewType vs = v;

  // Copy back
  // Note:  don't create h_vs through create_mirror_view and deep_copy
  // since Kokkos doesn't support deep_copy of non-contiguous views
  //stride_host_view_type h_vs = Kokkos::create_mirror_view(vs);
  //Kokkos::deep_copy(h_vs, vs);
  stride_host_view_type h_vs = h_v;

  // Check
  success = true;
  TEUCHOS_TEST_EQUALITY(h_vs.extent(0), num_rows, out, success);
  TEUCHOS_TEST_EQUALITY(h_vs.extent(1), num_cols, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(h_vs), fad_size+1, out, success);
  for (size_type i=0; i<num_rows; ++i) {
    for (size_type j=0; j<num_cols; ++j) {
      FadType f = generate_fad<FadType>(num_rows, num_cols, fad_size, i, j);
      success = success && checkFads(f, h_vs(i,j), out);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, CommonViewAllocMixedSpec, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Kokkos::LayoutContiguous<Layout>,Device> ContViewType;
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
  typedef typename ContViewType::size_type size_type;

  const size_type num_rows = global_num_rows;
  const size_type num_cols = global_num_cols;
  const size_type fad_size = global_fad_size;

  // Create contiguous view
  ContViewType v1;
#if defined (SACADO_DISABLE_FAD_VIEW_SPEC)
  v1 = ContViewType ("view", num_rows, num_cols);
#else
  v1 = ContViewType ("view", num_rows, num_cols, fad_size+1);
#endif

  // Create non-contiguous view using commen_view_alloc_prop
  auto cprop = Kokkos::common_view_alloc_prop(v1);
  ViewType v2(Kokkos::view_alloc("v2",cprop), num_rows, num_cols);

  // Check dimensions are correct for v2
  success = true;
  TEUCHOS_TEST_EQUALITY(v2.extent(0), num_rows, out, success);
  TEUCHOS_TEST_EQUALITY(v2.extent(1), num_cols, out, success);
  TEUCHOS_TEST_EQUALITY(Kokkos::dimension_scalar(v2), fad_size+1, out, success);
}

#else

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, ShmemSize, FadType, Layout, Device )
{
  typedef Kokkos::View<FadType**,Layout,Device> ViewType;
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
    (( sizeof(FadType) * global_num_rows * global_num_cols + mask ) & ~mask) + sizeof(typename ViewType::traits::value_type);
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

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, SFadNoSizeArg, FadType, Layout, Device ) {}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, Partition, FadType, Layout, Device ) {}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, AssignLayoutContiguousToLayoutStride, FadType, Layout, Device ) {}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(
  Kokkos_View_Fad, CommonViewAllocMixedSpec, FadType, Layout, Device ) {}

#endif

#define VIEW_FAD_TESTS_FLD( F, L, D )                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Size, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy_ConstantScalar, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy_ConstantZero, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy_ConstantFad, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DeepCopy_ConstantFadFull, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, LocalDeepCopy, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, LocalDeepCopyTeam, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, ScalarAssign, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, ValueAssign, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Resize, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Unmanaged, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Unmanaged2, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, UnmanagedConst, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, UnmanagedConst2, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Multiply, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, MultiplyUpdate, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, MultiplyConst, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, MultiplyMixed, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Rank8, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Roger, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AtomicAdd, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AssignDifferentStrides, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, ScalarValue, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DynRankDimensionScalar, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DynRankAssignStatic0, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DynRankAssignStatic1, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DynRankAssignStatic2, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, DynRankMultiply, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, SubdynrankviewCol, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, SubdynrankviewRow, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, SubdynrankviewScalar, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Subview, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Subview2, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, ShmemSize, F, L, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, ConstViewAssign, F, L, D )

#define VIEW_FAD_TESTS_SFLD( F, L, D )                                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, SFadNoSizeArg, F, L, D )

#define VIEW_FAD_TESTS_FDI( F, D )                                       \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_FAD_TESTS_FLD( F, LayoutLeft, D )                                \
  VIEW_FAD_TESTS_FLD( F, LayoutRight, D )                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AssignLayoutContiguousToLayoutStride, F, LayoutLeft, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, AssignLayoutContiguousToLayoutStride, F, LayoutRight, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, CommonViewAllocMixedSpec, F, LayoutLeft, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, CommonViewAllocMixedSpec, F, LayoutRight, D )

#define VIEW_FAD_TESTS_SFDI( F, D )                                     \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_FAD_TESTS_SFLD( F, LayoutLeft, D )                               \
  VIEW_FAD_TESTS_SFLD( F, LayoutRight, D )

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)
typedef Kokkos::LayoutContiguous<Kokkos::LayoutLeft> LeftContiguous;
typedef Kokkos::LayoutContiguous<Kokkos::LayoutRight> RightContiguous;
#define VIEW_FAD_TESTS_FDC( F, D )                                      \
  VIEW_FAD_TESTS_FLD( F, LeftContiguous, D )                            \
  VIEW_FAD_TESTS_FLD( F, RightContiguous, D )                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Partition, F, LeftContiguous, D ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Kokkos_View_Fad, Partition, F, RightContiguous, D )

#define VIEW_FAD_TESTS_SFDC( F, D )                                     \
  VIEW_FAD_TESTS_SFLD( F, LeftContiguous, D )                           \
  VIEW_FAD_TESTS_SFLD( F, RightContiguous, D )
#else
#define VIEW_FAD_TESTS_FDC( F, D ) /* */
#define VIEW_FAD_TESTS_SFDC( F, D ) /* */
#endif

#define VIEW_FAD_TESTS_FD( F, D )                                       \
  VIEW_FAD_TESTS_FDI( F, D )                                            \
  VIEW_FAD_TESTS_FDC( F, D )

#define VIEW_FAD_TESTS_SFD( F, D )                                      \
  VIEW_FAD_TESTS_SFDI( F, D )                                           \
  VIEW_FAD_TESTS_SFDC( F, D )

// We've unified the implementation for the different Fad variants, so
// there is no reason to test ELRFad, CacheFad, and ELRCacheFad.
typedef Sacado::Fad::DFad<double> DFadType;
typedef Sacado::Fad::SLFad<double,2*global_fad_size> SLFadType;
typedef Sacado::Fad::SFad<double,global_fad_size> SFadType;

/*
typedef Sacado::ELRFad::DFad<double> ELRDFadType;
typedef Sacado::ELRFad::SLFad<double,2*global_fad_size> ELRSLFadType;
typedef Sacado::ELRFad::SFad<double,global_fad_size> ELRSFadType;

typedef Sacado::CacheFad::DFad<double> CacheDFadType;
typedef Sacado::CacheFad::SLFad<double,2*global_fad_size> CacheSLFadType;
typedef Sacado::CacheFad::SFad<double,global_fad_size> CacheSFadType;

typedef Sacado::ELRCacheFad::DFad<double> ELRCacheDFadType;
typedef Sacado::ELRCacheFad::SLFad<double,2*global_fad_size> ELRCacheSLFadType;
typedef Sacado::ELRCacheFad::SFad<double,global_fad_size> ELRCacheSFadType;
*/

// We can't use DFad unless we use the View specialization
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC) && SACADO_TEST_DFAD
#define VIEW_FAD_TESTS_D( D )                            \
  VIEW_FAD_TESTS_FD( SFadType, D )                       \
  VIEW_FAD_TESTS_FD( SLFadType, D )                      \
  VIEW_FAD_TESTS_FD( DFadType, D )                       \
  VIEW_FAD_TESTS_SFD( SFadType, D )

#if 0
  VIEW_FAD_TESTS_FD( ELRSFadType, D )                    \
  VIEW_FAD_TESTS_FD( ELRSLFadType, D )                   \
  VIEW_FAD_TESTS_FD( ELRDFadType, D )                    \
  VIEW_FAD_TESTS_FD( CacheSFadType, D )                  \
  VIEW_FAD_TESTS_FD( CacheSLFadType, D )                 \
  VIEW_FAD_TESTS_FD( CacheDFadType, D )                  \
  VIEW_FAD_TESTS_FD( ELRCacheSFadType, D )               \
  VIEW_FAD_TESTS_FD( ELRCacheSLFadType, D )              \
  VIEW_FAD_TESTS_FD( ELRCacheDFadType, D )               \
  VIEW_FAD_TESTS_SFD( SFadType, D )                      \
  VIEW_FAD_TESTS_SFD( ELRSFadType, D )                   \
  VIEW_FAD_TESTS_SFD( CacheSFadType, D )                 \
  VIEW_FAD_TESTS_SFD( ELRCacheSFadType, D )
#endif

#else

#define VIEW_FAD_TESTS_D( D )                        \
  VIEW_FAD_TESTS_FD( SFadType, D )                   \
  VIEW_FAD_TESTS_FD( SLFadType, D )                  \
  VIEW_FAD_TESTS_SFD( SFadType, D )

#if 0
  VIEW_FAD_TESTS_FD( ELRSFadType, D )                \
  VIEW_FAD_TESTS_FD( ELRSLFadType, D )               \
  VIEW_FAD_TESTS_FD( CacheSFadType, D )              \
  VIEW_FAD_TESTS_FD( CacheSLFadType, D )             \
  VIEW_FAD_TESTS_FD( ELRCacheSFadType, D )           \
  VIEW_FAD_TESTS_FD( ELRCacheSLFadType, D )          \
  VIEW_FAD_TESTS_SFD( SFadType, D )                  \
  VIEW_FAD_TESTS_SFD( ELRSFadType, D )               \
  VIEW_FAD_TESTS_SFD( CacheSFadType, D )             \
  VIEW_FAD_TESTS_SFD( ELRCacheSFadType, D )
#endif

#endif
