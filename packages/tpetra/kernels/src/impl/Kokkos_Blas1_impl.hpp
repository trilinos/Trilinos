/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOS_BLAS1_IMPL_HPP_
#define KOKKOS_BLAS1_IMPL_HPP_

#include <TpetraKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_InnerProductSpaceTraits.hpp>

namespace KokkosBlas {
namespace Impl {

/// \brief Functor that implements the single-vector, two-argument
///   version of KokkosBlas::dot (dot product of two vectors).
///
/// \tparam XVector Type of the first vector x; 1-D View
/// \tparam YVector Type of the second vector y; 1-D View
/// \tparam SizeType Type of the row index used in the dot product.
///   For best performance, use int instead of size_t here.
template<class XVector, class YVector, typename SizeType = typename XVector::size_type>
struct DotFunctor
{
  typedef typename XVector::execution_space        execution_space;
  typedef SizeType                                 size_type;
  typedef typename XVector::non_const_value_type   xvalue_type;
  typedef Kokkos::Details::InnerProductSpaceTraits<xvalue_type>  IPT;
  typedef typename IPT::dot_type                   value_type;

  XVector  m_x;
  YVector  m_y;

  DotFunctor (const XVector& x, const YVector& y) : m_x (x), m_y (y) {}

  // Prefer const size_type& to const size_type or size_type,
  // since the compiler has an easier time inlining the former.
  KOKKOS_FORCEINLINE_FUNCTION void
  operator() (const size_type &i, value_type& sum) const
  {
    sum += IPT::dot (m_x(i), m_y(i));  // m_x(i) * m_y(i)
  }

  KOKKOS_INLINE_FUNCTION void
  init (volatile value_type& update) const
  {
    update = Kokkos::Details::ArithTraits<value_type>::zero ();
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& update,
        const volatile value_type& source) const
  {
    update += source ;
  }
};

/// \brief Struct that implements the single-vector, two-argument
///   version of KokkosBlas::dot (dot product of two vectors).
///
/// The first five template parameters correspond to the template
/// parameters of XVector in DotFunctor (see above).  The last five
/// template parameters correspond to the template parameters of
/// YVector in DotFunctor.
///
/// This is the "generic" implementation of the dot product.  We
/// declare full specializations in this header file.  Definitions of
/// those specializations live in one or more .cpp files in this
/// directory.
template<class XT, class XL, class XD, class XM, class XS,
         class YT, class YL, class YD, class YM, class YS>
struct Dot {
  typedef Kokkos::View<XT,XL,XD,XM,XS> XVector;
  typedef Kokkos::View<YT,YL,YD,YM,YS> YVector;
  typedef typename Kokkos::Details::InnerProductSpaceTraits<typename XVector::non_const_value_type>::dot_type dot_type;

  //! Return the dot product of x and y.
  static dot_type
  dot (const XVector& x, const YVector& y)
  {
    dot_type result;

    // With Intel 15, Scalar = double, Kokkos::Serial execution space,
    // with 100000 entries for 1000 trials, using int instead of
    // size_t is 2x faster.
    if (x.dimension_0 () < typename XVector::size_type (INT_MAX)) {
      Kokkos::RangePolicy<typename XVector::execution_space, int> policy (0, x.dimension_0 ());
      typedef DotFunctor<XVector, YVector, int> functor_type;
      Kokkos::parallel_reduce (policy, functor_type (x, y), result);
    }
    else {
      typedef typename XVector::size_type size_type;
      Kokkos::RangePolicy<typename XVector::execution_space, size_type> policy (0, x.dimension_0 ());
      typedef DotFunctor<XVector, YVector, size_type> functor_type;
      Kokkos::parallel_reduce (policy, functor_type (x, y), result);
    }
    return result;
  }

  /// \brief Return the dot product of x and y.
  ///
  /// \param x [in] First vector.
  /// \param y [in] Second vector.
  /// \param label [in] This is passed into Kokkos::parallel_for, as a
  ///   label for profiling.
  static dot_type
  dot (const std::string& label, const XVector& x, const YVector& y)
  {
    dot_type result;
    if (x.dimension_0 () < typename XVector::size_type (INT_MAX)) {
      Kokkos::RangePolicy<typename XVector::execution_space, int> policy (0, x.dimension_0 ());
      typedef DotFunctor<XVector, YVector, int> functor_type;
      Kokkos::parallel_reduce (label, policy, functor_type (x, y), result);
    }
    else {
      typedef typename XVector::size_type size_type;
      Kokkos::RangePolicy<typename XVector::execution_space, size_type> policy (0, x.dimension_0 ());
      typedef DotFunctor<XVector, YVector, size_type> functor_type;
      Kokkos::parallel_reduce (label, policy, functor_type (x, y), result);
    }
    return result;
  }
};

//
// Macro that declares a full specialization of the Dot struct.
//
#define KOKKOSBLAS_IMPL_V_DOT_DECL( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
template<> \
struct Dot<const SCALAR*, \
           LAYOUT, \
           Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
           Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
           Kokkos::Impl::ViewDefault, \
           const SCALAR*, \
           LAYOUT, \
           Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
           Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
           Kokkos::Impl::ViewDefault> \
{ \
  typedef const SCALAR* XT; \
  typedef LAYOUT XL; \
  typedef Kokkos::Device<EXEC_SPACE, MEM_SPACE> XD; \
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> XM; \
  typedef Kokkos::Impl::ViewDefault XS; \
  typedef Kokkos::View<XT,XL,XD,XM,XS> XVector; \
  typedef const SCALAR* YT; \
  typedef LAYOUT YL; \
  typedef Kokkos::Device<EXEC_SPACE, MEM_SPACE> YD; \
  typedef Kokkos::MemoryTraits<Kokkos::Unmanaged> YM; \
  typedef Kokkos::Impl::ViewDefault YS; \
  typedef Kokkos::View<YT,YL,YD,YM,YS> YVector; \
  typedef Kokkos::Details::InnerProductSpaceTraits<XVector::non_const_value_type>::dot_type dot_type; \
 \
  static dot_type \
  dot (const XVector& x, const XVector& y); \
 \
  static dot_type \
  dot (const std::string& label, const XVector& x, const XVector& y); \
};

//
// Declare full specializations of the Dot struct.
//

#ifdef KOKKOS_HAVE_SERIAL

  KOKKOSBLAS_IMPL_V_DOT_DECL( double, Kokkos::LayoutLeft, Kokkos::Serial, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP

  KOKKOSBLAS_IMPL_V_DOT_DECL( double, Kokkos::LayoutLeft, Kokkos::OpenMP, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD

  KOKKOSBLAS_IMPL_V_DOT_DECL( double, Kokkos::LayoutLeft, Kokkos::Threads, Kokkos::HostSpace )

#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA

  KOKKOSBLAS_IMPL_V_DOT_DECL( double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaSpace )

#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_CUDA

  KOKKOSBLAS_IMPL_V_DOT_DECL( double, Kokkos::LayoutLeft, Kokkos::Cuda, Kokkos::CudaUVMSpace )

#endif // KOKKOS_HAVE_CUDA

//
// Macro that defines a full specialization of the Dot struct.
// This macro is invoked in one or more .cpp files in this directory.
//
#define KOKKOSBLAS_IMPL_V_DOT_DEF( SCALAR, LAYOUT, EXEC_SPACE, MEM_SPACE ) \
Dot<const SCALAR*, \
    LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
    Kokkos::Impl::ViewDefault, \
    const SCALAR*, \
    LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
    Kokkos::Impl::ViewDefault>::dot_type \
Dot<const SCALAR*, \
    LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
    Kokkos::Impl::ViewDefault, \
    const SCALAR*, \
    LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
    Kokkos::Impl::ViewDefault>:: \
dot (const XVector& x, const XVector& y) \
{ \
  dot_type result; \
  if (x.dimension_0 () < XVector::size_type (INT_MAX)) { \
    Kokkos::RangePolicy<XVector::execution_space, int> policy (0, x.dimension_0 ()); \
    typedef DotFunctor<XVector, YVector, int> functor_type; \
    Kokkos::parallel_reduce (policy, functor_type (x, y), result); \
  } \
  else { \
    typedef XVector::size_type size_type; \
    Kokkos::RangePolicy<XVector::execution_space, size_type> policy (0, x.dimension_0 ()); \
    typedef DotFunctor<XVector, YVector, size_type> functor_type; \
    Kokkos::parallel_reduce (policy, functor_type (x, y), result); \
  } \
  return result; \
} \
 \
Dot<const SCALAR*, \
    LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
    Kokkos::Impl::ViewDefault, \
    const SCALAR*, \
    LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
    Kokkos::Impl::ViewDefault>::dot_type \
Dot<const SCALAR*, \
    LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
    Kokkos::Impl::ViewDefault, \
    const SCALAR*, \
    LAYOUT, \
    Kokkos::Device<EXEC_SPACE, MEM_SPACE>, \
    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
    Kokkos::Impl::ViewDefault>:: \
dot (const std::string& label, const XVector& x, const XVector& y) \
{ \
  dot_type result; \
  if (x.dimension_0 () < XVector::size_type (INT_MAX)) { \
    Kokkos::RangePolicy<XVector::execution_space, int> policy (0, x.dimension_0 ()); \
    typedef DotFunctor<XVector, YVector, int> functor_type; \
    Kokkos::parallel_reduce (label, policy, functor_type (x, y), result); \
  } \
  else { \
    typedef XVector::size_type size_type; \
    Kokkos::RangePolicy<XVector::execution_space, size_type> policy (0, x.dimension_0 ()); \
    typedef DotFunctor<XVector, YVector, size_type> functor_type; \
    Kokkos::parallel_reduce (label, policy, functor_type (x, y), result); \
  } \
  return result; \
}

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_BLAS1_IMPL_HPP_
