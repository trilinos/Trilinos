/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOSSPARSE_IMPL_SPMV_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPMV_SPEC_HPP_

#include <KokkosKernels_config.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <cxxabi.h>

#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Controls.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY 
#include <KokkosSparse_spmv_impl.hpp>
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM>
struct spmv_eti_spec_avail {
  enum : bool { value = false };
};
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM,
         const bool integerScalarType =
           std::is_integral<typename std::decay<AT>::type>::value>
struct spmv_mv_eti_spec_avail {
  enum : bool { value = false };
};

}
}


#define KOKKOSSPARSE_SPMV_ETI_SPEC_AVAIL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template<> \
    struct spmv_eti_spec_avail<const SCALAR_TYPE, \
                  const ORDINAL_TYPE, \
                  Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                  const OFFSET_TYPE, \
                  SCALAR_TYPE const*, \
                  LAYOUT_TYPE, \
                  Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
                  SCALAR_TYPE*, \
                  LAYOUT_TYPE, \
                  Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                  Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
    { enum : bool { value = true }; };


#define KOKKOSSPARSE_SPMV_MV_ETI_SPEC_AVAIL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
    template<> \
    struct spmv_mv_eti_spec_avail <const SCALAR_TYPE, \
                                       const ORDINAL_TYPE, \
                                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                                       const OFFSET_TYPE, \
                                       SCALAR_TYPE const**, \
                                       LAYOUT_TYPE, \
                                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
                                       SCALAR_TYPE**, \
                                       LAYOUT_TYPE, \
                                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
      { enum : bool { value = true }; };


// Include the actual specialization declarations
#include<KokkosSparse_spmv_tpl_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosSparse_spmv_eti_spec_avail.hpp>
#include<generated_specializations_hpp/KokkosSparse_spmv_mv_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::spmv (sparse matrix - dense
///   vector multiply) for single vectors (1-D Views).
///
/// The first 5 template parameters are the same as those of
/// KokkosSparse::CrsMatrix.  In particular:
///
/// AT: type of each entry of the sparse matrix
/// AO: ordinal type (type of column indices) of the sparse matrix
/// AS: offset type (type of row offsets) of the sparse matrix
///
/// The next 4 template parameters (that start with X) correspond to
/// the input Kokkos::View.  The last 4 template parameters (that start
/// with Y) correspond to the output Kokkos::View.
///
/// For the implementation of KokkosSparse::spmv for multivectors (2-D
/// Views), see the SPMV_MV struct below.
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM,
         bool tpl_spec_avail =
             spmv_tpl_spec_avail< AT, AO, AD, AM, AS,
                                  XT, XL, XD, XM,
                                  YT, YL, YD, YM>::value,
         bool eti_spec_avail =
             spmv_eti_spec_avail< AT, AO, AD, AM, AS,
                                  XT, XL, XD, XM,
                                  YT, YL, YD, YM>::value >
struct SPMV{
  typedef CrsMatrix<AT,AO,AD,AM,AS> AMatrix;
  typedef Kokkos::View<XT,XL,XD,XM> XVector;
  typedef Kokkos::View<YT,YL,YD,YM> YVector;

  typedef typename YVector::non_const_value_type coefficient_type;

  static void spmv (const KokkosKernels::Experimental::Controls& controls,
		    const char mode[],
		    const coefficient_type& alpha,
		    const AMatrix& A,
		    const XVector& x,
		    const coefficient_type& beta,
		    const YVector& y);
};

// Unification layer
/// \brief Implementation of KokkosBlas::spmv (sparse matrix - dense
///   vector multiply) for multiple vectors at a time (multivectors)
///   and possibly multiple coefficients at a time.
///
/// This struct implements the following operations:
///
///   1. Y(:,j) := beta(j) * Y(:,j) + alpha(j) * Op(A) * X(:,j)
///   2. Y(:,j) := beta(j) * Y(:,j) + alpha * Op(A) * X(:,j)
///   3. Y(:,j) := beta * Y(:,j) + alpha(j) * Op(A) * X(:,j)
///   4. Y(:,j) := beta * Y(:,j) + alpha * Op(A) * X(:,j)
///
/// In #1 and #2 above, beta is a 1-D View of coefficients, one for
/// each column of Y.  In #1 and #3 above, alpha is a 1-D View of
/// coefficients, one for each column of X.  Otherwise, alpha
/// resp. beta are each a single coefficient.  In all of these
/// operations, X and Y are 2-D Views ("multivectors").  A is a sparse
/// matrix, and Op(A) is either A itself, its transpose, or its
/// conjugate transpose, depending on the 'mode' argument.
///
/// The first 5 template parameters are the template parameters of the
/// input 1-D View of coefficients 'alpha'.  The next 5 template
/// parameters are the same as those of KokkosSparse::CrsMatrix.  In
/// particular:
///
/// AT: type of each entry of the sparse matrix
/// AO: ordinal type (type of column indices) of the sparse matrix
/// AS: offset type (type of row offsets) of the sparse matrix
///
/// The next 4 template parameters (that start with X) correspond to
/// the input Kokkos::View.  The 4 template parameters after that
/// (that start with lower-case b) are the template parameters of the
/// input 1-D View of coefficients 'beta'.  Next, the 5 template
/// parameters that start with Y correspond to the output
/// Kokkos::View.  The last template parameter indicates whether the
/// matrix's entries have integer type.  Per Github Issue #700, we
/// don't optimize as heavily for that case, in order to reduce build
/// times and library sizes.
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM,
         const bool integerScalarType =
           std::is_integral<typename std::decay<AT>::type>::value,
         bool tpl_spec_avail =
             spmv_mv_tpl_spec_avail< AT, AO, AD, AM, AS,
                                  XT, XL, XD, XM,
                                  YT, YL, YD, YM>::value,
         bool eti_spec_avail =
             spmv_mv_eti_spec_avail< AT, AO, AD, AM, AS,
                                  XT, XL, XD, XM,
                                  YT, YL, YD, YM>::value >
struct SPMV_MV{
  typedef CrsMatrix<AT,AO,AD,AM,AS> AMatrix;
  typedef Kokkos::View<XT,XL,XD,XM> XVector;
  typedef Kokkos::View<YT,YL,YD,YM> YVector;
  typedef typename YVector::non_const_value_type coefficient_type;

  static void
  spmv_mv (const char mode[],
           const coefficient_type& alpha,
           const AMatrix& A,
           const XVector& x,
           const coefficient_type& beta,
           const YVector& y);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
//! Full specialization of spmv for single vectors (1-D Views).
// Unification layer
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM>
struct SPMV < AT, AO, AD, AM, AS,
              XT, XL, XD, XM,
              YT, YL, YD, YM, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{

  typedef CrsMatrix<AT,AO,AD,AM,AS> AMatrix;
  typedef Kokkos::View<XT,XL,XD,XM> XVector;
  typedef Kokkos::View<YT,YL,YD,YM> YVector;
  typedef typename YVector::non_const_value_type coefficient_type;

  static void
  spmv (const KokkosKernels::Experimental::Controls& controls,
	const char mode[],
	const coefficient_type& alpha,
	const AMatrix& A,
	const XVector& x,
	const coefficient_type& beta,
	const YVector& y)
  {
    typedef Kokkos::Details::ArithTraits<coefficient_type> KAT;

    typedef Kokkos::Details::ArithTraits<coefficient_type> KAT;

    if (alpha == KAT::zero ()) {
      if (beta != KAT::one ()) {
        KokkosBlas::scal (y, beta, y);
      }
      return;
    }

    if (beta == KAT::zero ()) {
      spmv_beta<AMatrix, XVector, YVector, 0> (controls, mode, alpha, A, x, beta, y);
    }
    else if (beta == KAT::one ()) {
      spmv_beta<AMatrix, XVector, YVector, 1> (controls, mode, alpha, A, x, beta, y);
    }
    else if (beta == -KAT::one ()) {
      spmv_beta<AMatrix, XVector, YVector, -1> (controls, mode, alpha, A, x, beta, y);
    }
    else {
      spmv_beta<AMatrix, XVector, YVector, 2> (controls, mode, alpha, A, x, beta, y);
    }
  }
};

//! Full specialization of spmv_mv for single vectors (2-D Views).
// Unification layer
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM>
struct SPMV_MV<AT, AO, AD, AM, AS,
               XT, XL, XD, XM,
               YT, YL, YD, YM,
               false,
               false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{
  typedef CrsMatrix<AT,AO,AD,AM,AS> AMatrix;
  typedef Kokkos::View<XT,XL,XD,XM> XVector;
  typedef Kokkos::View<YT,YL,YD,YM> YVector;
  typedef typename YVector::non_const_value_type coefficient_type;

  static void
  spmv_mv (const char mode[],
           const coefficient_type& alpha,
           const AMatrix& A,
           const XVector& x,
           const coefficient_type& beta,
           const YVector& y)
  {
    typedef Kokkos::Details::ArithTraits<coefficient_type> KAT;

    if (alpha == KAT::zero ()) {
      spmv_alpha_mv<AMatrix, XVector, YVector, 0> (mode, alpha, A, x, beta, y);
    }
    else if (alpha == KAT::one ()) {
      spmv_alpha_mv<AMatrix, XVector, YVector, 1> (mode, alpha, A, x, beta, y);
    }
    else if (alpha == -KAT::one ()) {
      spmv_alpha_mv<AMatrix, XVector, YVector, -1> (mode, alpha, A, x, beta, y);
    }
    else {
      spmv_alpha_mv<AMatrix, XVector, YVector, 2> (mode, alpha, A, x, beta, y);
    }
  }
};

template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM>
struct SPMV_MV<AT, AO, AD, AM, AS,
               XT, XL, XD, XM,
               YT, YL, YD, YM,
               true,
               false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY>{
  typedef CrsMatrix<AT,AO,AD,AM,AS> AMatrix;
  typedef Kokkos::View<XT,XL,XD,XM> XVector;
  typedef Kokkos::View<YT,YL,YD,YM> YVector;
  typedef typename YVector::non_const_value_type coefficient_type;

  static void
  spmv_mv (const char mode[],
           const coefficient_type& alpha,
           const AMatrix& A,
           const XVector& x,
           const coefficient_type& beta,
           const YVector& y)
  {
    static_assert (std::is_integral<AT>::value,
                   "This implementation is only for integer Scalar types.");
    typedef SPMV<AT, AO, AD, AM, AS,
      typename XVector::value_type*, XL, XD, XM,
      typename YVector::value_type*, YL, YD, YM> impl_type;
    //Create a default Controls object (the impl_type::spmv below is for rank-1, which requires it)
    KokkosKernels::Experimental::Controls defaultControls;
    for (typename AMatrix::non_const_size_type j = 0; j < x.extent(1); ++j) {
      auto x_j = Kokkos::subview (x, Kokkos::ALL (), j);
      auto y_j = Kokkos::subview (y, Kokkos::ALL (), j);
      impl_type::spmv (defaultControls, mode, alpha, A, x_j, beta, y_j);
    }
  }
};
#endif



}
}

//
// Macro for declaration of full specialization of
// KokkosSparse::Impl::SpMV.  This is NOT for users!!!  All
// the declarations of full specializations go in this header file.
// We may spread out definitions (see _DEF macro below) across one or
// more .cpp files.
//
#define KOKKOSSPARSE_SPMV_ETI_SPEC_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
    extern template struct  \
    SPMV<const SCALAR_TYPE, \
         const ORDINAL_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
         const OFFSET_TYPE, \
         SCALAR_TYPE const*, \
         LAYOUT_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
         SCALAR_TYPE*, \
         LAYOUT_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged>, false, true >;


#define KOKKOSSPARSE_SPMV_ETI_SPEC_INST( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template struct  \
    SPMV<const SCALAR_TYPE, \
         const ORDINAL_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
         const OFFSET_TYPE, \
         SCALAR_TYPE const*, \
         LAYOUT_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
         SCALAR_TYPE*, \
         LAYOUT_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged>, false, true >;

#define KOKKOSSPARSE_SPMV_MV_ETI_SPEC_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
    extern template struct  \
    SPMV_MV<const SCALAR_TYPE, \
         const ORDINAL_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
         const OFFSET_TYPE, \
         SCALAR_TYPE const**, \
         LAYOUT_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
         SCALAR_TYPE**, \
         LAYOUT_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged>, std::is_integral<typename std::decay<SCALAR_TYPE>::type>::value, false, true >;


#define KOKKOSSPARSE_SPMV_MV_ETI_SPEC_INST( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE) \
    template struct  \
    SPMV_MV<const SCALAR_TYPE, \
         const ORDINAL_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
         const OFFSET_TYPE, \
         SCALAR_TYPE const**, \
         LAYOUT_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
         SCALAR_TYPE**, \
         LAYOUT_TYPE, \
         Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
         Kokkos::MemoryTraits<Kokkos::Unmanaged>, std::is_integral<typename std::decay<SCALAR_TYPE>::type>::value, false, true >;

#include<KokkosSparse_spmv_tpl_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosSparse_spmv_eti_spec_decl.hpp>
#include<generated_specializations_hpp/KokkosSparse_spmv_mv_eti_spec_decl.hpp>

#endif // KOKKOSSPARSE_IMPL_SPMV_SPEC_HPP_
