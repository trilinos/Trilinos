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

#ifndef KOKKOS_SPARSE_IMPL_SPMV_DECL_HPP_
#define KOKKOS_SPARSE_IMPL_SPMV_DECL_HPP_

#include "TpetraKernels_config.h"
#include "Kokkos_Sparse_CrsMatrix.hpp"
#include "TpetraKernels_ETIHelperMacros.h"

namespace KokkosSparse {
namespace Impl {

template<class InputType, class DeviceType>
struct GetCoeffView {
  typedef Kokkos::View<InputType*,Kokkos::LayoutLeft,DeviceType> view_type;
  typedef Kokkos::View<typename view_type::non_const_value_type*,
                       Kokkos::LayoutLeft,DeviceType> non_const_view_type;
  static non_const_view_type get_view(const InputType in, const int size) {
    non_const_view_type aview("CoeffView",size);
    if(size>0)
      Kokkos::deep_copy(aview,in);
    return aview;
  }
};

template<class IT, class IL, class ID, class IM, class IS, class DeviceType>
struct GetCoeffView<Kokkos::View<IT*,IL,ID,IM,IS>,DeviceType> {
  typedef Kokkos::View<IT*,IL,ID,IM,IS> view_type;
  static Kokkos::View<IT*,IL,ID,IM,IS> get_view(const Kokkos::View<IT*,IL,ID,IM,IS>& in, int size) {
    return in;
  }
};

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
/// The next 5 template parameters (that start with X) correspond to
/// the input Kokkos::View.  The last 5 template parameters (that start
/// with Y) correspond to the output Kokkos::View.
///
/// For the implementation of KokkosSparse::spmv for multivectors (2-D
/// Views), see the SPMV_MV struct below.
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM>
struct SPMV
{
  typedef CrsMatrix<AT,AO,AD,AM,AS> AMatrix;
  typedef Kokkos::View<XT,XL,XD,XM> XVector;
  typedef Kokkos::View<YT,YL,YD,YM> YVector;
  typedef typename YVector::non_const_value_type coefficient_type;

  static void
  spmv (const char mode[],
        const coefficient_type& alpha,
        const AMatrix& A,
        const XVector& x,
        const coefficient_type& beta,
        const YVector& y);
};

//
// Macro for declaring a full specialization of the SPMV struct, which
// implements KokkosSparse::spmv for single vectors (1-D Views).  We
// use this macro below.
//
// SCALAR_TYPE: The type of each entry in the sparse matrix
// ORDINAL_TYPE: The type of each column index in the sparse matrix
// OFFSET_TYPE: The type of each row offset in the sparse matrix
// LAYOUT_TYPE: The layout of the Kokkos::View vector arguments
//   of the sparse matrix-vector multiply
// EXEC_SPACE_TYPE: The execution space type
// MEM_SPACE_TYPE: The memory space type
//
#define KOKKOSSPARSE_IMPL_SPMV_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
template<> \
struct SPMV<const SCALAR_TYPE, \
            ORDINAL_TYPE, \
            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
            Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
            OFFSET_TYPE, \
            const SCALAR_TYPE*, \
            LAYOUT_TYPE, \
            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
            Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
            SCALAR_TYPE*, \
            LAYOUT_TYPE, \
            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
            Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
{ \
  typedef CrsMatrix<const SCALAR_TYPE, \
                    ORDINAL_TYPE, \
                    Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                    OFFSET_TYPE> AMatrix; \
  typedef Kokkos::View<const SCALAR_TYPE*, \
                       LAYOUT_TYPE, \
                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector; \
  typedef Kokkos::View<SCALAR_TYPE*, \
                       LAYOUT_TYPE, \
                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector; \
  typedef typename YVector::non_const_value_type coefficient_type; \
  \
  static void \
  spmv (const char mode[], \
        const coefficient_type& alpha, \
        const AMatrix& A, \
        const XVector& x, \
        const coefficient_type& beta, \
        const YVector& y); \
};

//
// Macro for declaring a full specialization of the SPMV struct, which
// implements KokkosSparse::spmv for single vectors (1-D Views).  This
// version of the macro uses the default OFFSET_TYPE, instead of
// letting users specify it (as with the above macro).  It also uses
// LAYOUT_TYPE = Kokkos::LayoutLeft (Tpetra's layout).
//
// We need to redefine this macro in full, rather than calling the one
// above, because macros don't allow arguments with commas in them.
// The correct OFFSET_TYPE default would otherwise (as of 18 Mar 2016;
// see Tpetra::CrsGraph public typedef 'local_graph_type') be
// Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft,
// EXEC_SPACE_TYPE>::size_type.
//
// SCALAR_TYPE: The type of each entry in the sparse matrix
// ORDINAL_TYPE: The type of each column index in the sparse matrix
// EXEC_SPACE_TYPE: The execution space type
// MEM_SPACE_TYPE: The memory space type
//
#define KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL( SCALAR_TYPE, ORDINAL_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
template<> \
struct SPMV<const SCALAR_TYPE, \
            ORDINAL_TYPE, \
            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
            Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
            Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type, \
            const SCALAR_TYPE*, \
            Kokkos::LayoutLeft, \
            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
            Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
            SCALAR_TYPE*, \
            Kokkos::LayoutLeft, \
            Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
            Kokkos::MemoryTraits<Kokkos::Unmanaged> > \
{ \
  typedef CrsMatrix<const SCALAR_TYPE, \
                    ORDINAL_TYPE, \
                    Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                    Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type> AMatrix; \
  typedef Kokkos::View<const SCALAR_TYPE*, \
                       Kokkos::LayoutLeft, \
                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector; \
  typedef Kokkos::View<SCALAR_TYPE*, \
                       Kokkos::LayoutLeft, \
                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector; \
  typedef typename YVector::non_const_value_type coefficient_type; \
  \
  static void \
  spmv (const char mode[], \
        const coefficient_type& alpha, \
        const AMatrix& A, \
        const XVector& x, \
        const coefficient_type& beta, \
        const YVector& y); \
};

//
// Declarations of full specializations of the SPMV struct.
// Definitions go in various .cpp file(s) in this directory.
//

TPETRAKERNELS_ETI_MANGLING_TYPEDEFS()

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_SERIAL
#define KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_SERIAL( SCALAR, LO ) \
  KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL( SCALAR, LO, Kokkos::Serial, Kokkos::HostSpace )

TPETRAKERNELS_INSTANTIATE_SL( KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_SERIAL )

#undef KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_SERIAL
#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_SERIAL

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_OPENMP
#define KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_OPENMP( SCALAR, LO ) \
  KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL( SCALAR, LO, Kokkos::OpenMP, Kokkos::HostSpace )

TPETRAKERNELS_INSTANTIATE_SL( KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_OPENMP )

#undef KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_OPENMP
#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_OPENMP

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_PTHREAD
#define KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_PTHREAD( SCALAR, LO ) \
  KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL( SCALAR, LO, Kokkos::Threads, Kokkos::HostSpace )

TPETRAKERNELS_INSTANTIATE_SL( KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_PTHREAD )

#undef KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_PTHREAD
#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_PTHREAD

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_CUDA
#define KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_CUDA( SCALAR, LO ) \
  KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL( SCALAR, LO, Kokkos::Cuda, Kokkos::CudaUVMSpace )

TPETRAKERNELS_INSTANTIATE_SL( KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_CUDA )

#undef KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DECL_CUDA
#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_CUDA


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
/// The next 5 template parameters (that start with X) correspond to
/// the input Kokkos::View.  The 5 template parameters after that
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
           std::is_integral<typename std::decay<AT>::type>::value>
struct SPMV_MV
{
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

//
// Macro for declaring a full specialization of the SPMV_MV struct.
// We use this macro below.
//
// SCALAR_TYPE: The type of each entry in the sparse matrix
// ORDINAL_TYPE: The type of each column index in the sparse matrix
// OFFSET_TYPE: The type of each row offset in the sparse matrix
// LAYOUT_TYPE: The layout of the Kokkos::View vector arguments
//   of the sparse matrix-vector multiply
// EXEC_SPACE_TYPE: The execution space type
// MEM_SPACE_TYPE: The memory space type
//
#define KOKKOSSPARSE_IMPL_SPMV_MV_DECL( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
template<> \
struct SPMV_MV<const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        OFFSET_TYPE, \
        const SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
        SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        std::is_integral<SCALAR_TYPE>::value> \
{ \
  typedef CrsMatrix<const SCALAR_TYPE, \
                    ORDINAL_TYPE, \
                    Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                    OFFSET_TYPE> AMatrix; \
  typedef Kokkos::View<const SCALAR_TYPE**, \
                       LAYOUT_TYPE, \
                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector; \
  typedef Kokkos::View<SCALAR_TYPE**, \
                       LAYOUT_TYPE, \
                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector; \
  typedef YVector::non_const_value_type coefficient_type; \
  \
  static void \
  spmv_mv (const char mode[], \
           const coefficient_type& alpha, \
           const AMatrix& A, \
           const XVector& x, \
           const coefficient_type& beta, \
           const YVector& y); \
};

//
// Macro for declaring a full specialization of the SPMV_MV struct,
// but filling in defaults for LAYOUT_TYPE and OFFSET_TYPE as follows:
//
//  - LAYOUT_TYPE = Kokkos::LayoutLeft
//  - OFFSET_TYPE = the graph's default offset type
//
// We use this macro below.
//
// SCALAR_TYPE: The type of each entry in the sparse matrix
// ORDINAL_TYPE: The type of each column index in the sparse matrix
// EXEC_SPACE_TYPE: The execution space type
// MEM_SPACE_TYPE: The memory space type
//
#define KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL( SCALAR_TYPE, ORDINAL_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
template<> \
struct SPMV_MV<const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type, \
        const SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
        SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        std::is_integral<SCALAR_TYPE>::value> \
{ \
  typedef CrsMatrix<const SCALAR_TYPE, \
                    ORDINAL_TYPE, \
                    Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                    Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
                    Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type> AMatrix; \
  typedef Kokkos::View<const SCALAR_TYPE**, \
                       Kokkos::LayoutLeft, \
                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> > XVector; \
  typedef Kokkos::View<SCALAR_TYPE**, \
                       Kokkos::LayoutLeft, \
                       Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > YVector; \
  typedef YVector::non_const_value_type coefficient_type; \
  \
  static void \
  spmv_mv (const char mode[], \
           const coefficient_type& alpha, \
           const AMatrix& A, \
           const XVector& x, \
           const coefficient_type& beta, \
           const YVector& y); \
};

//
// Declarations of full specializations of the SPMV_MV struct.
// Definitions go in various .cpp file(s) in this directory.
//

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_SERIAL
#define KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_SERIAL( SCALAR, LO ) \
  KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL( SCALAR, LO, Kokkos::Serial, Kokkos::HostSpace )

TPETRAKERNELS_INSTANTIATE_SL( KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_SERIAL )

#undef KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_SERIAL
#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_SERIAL

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_OPENMP
#define KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_OPENMP( SCALAR, LO ) \
  KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL( SCALAR, LO, Kokkos::OpenMP, Kokkos::HostSpace )

TPETRAKERNELS_INSTANTIATE_SL( KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_OPENMP )

#undef KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_OPENMP
#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_OPENMP

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_PTHREAD
#define KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_PTHREAD( SCALAR, LO ) \
  KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL( SCALAR, LO, Kokkos::Threads, Kokkos::HostSpace )

TPETRAKERNELS_INSTANTIATE_SL( KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_PTHREAD )

#undef KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_PTHREAD
#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_PTHREAD

#ifdef TPETRAKERNELS_BUILD_EXECUTION_SPACE_CUDA
#define KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_CUDA( SCALAR, LO ) \
  KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL( SCALAR, LO, Kokkos::Cuda, Kokkos::CudaUVMSpace )

TPETRAKERNELS_INSTANTIATE_SL( KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_CUDA )

#undef KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DECL_CUDA
#endif // TPETRAKERNELS_BUILD_EXECUTION_SPACE_CUDA

} // namespace Impl
} // namespace KokkosSparse

#endif // KOKKOS_SPARSE_IMPL_SPMV_DECL_HPP_
