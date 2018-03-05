/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_BLOCKTRIDICONTAINER_DEF_HPP
#define IFPACK2_BLOCKTRIDICONTAINER_DEF_HPP

#include <Teuchos_Details_MpiTypeTraits.hpp>

#include <Tpetra_Distributor.hpp>
#include <Tpetra_BlockMultiVector.hpp>

#include <Kokkos_ArithTraits.hpp>
#include <KokkosBatched_Util.hpp>
#include <KokkosBatched_Vector.hpp>
#if defined(__KOKKOSBATCHED_PROMOTION__)
#include <KokkosBatched_AddRadial_Decl.hpp>
#include <KokkosBatched_AddRadial_Impl.hpp>
#endif // __KOKKOSBATCHED_PROMOTION__
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBatched_Gemm_Serial_Impl.hpp>
#include <KokkosBatched_Gemv_Decl.hpp>
#include <KokkosBatched_Gemv_Serial_Impl.hpp>
#include <KokkosBatched_Trsm_Decl.hpp>
#include <KokkosBatched_Trsm_Serial_Impl.hpp>
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBatched_Trsv_Serial_Impl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>

#include "Ifpack2_BlockTriDiContainer_decl.hpp"

#include <memory>

//todo Move to KokkosKernels.
#if !defined(__KOKKOSBATCHED_PROMOTION__)
namespace KokkosBatched {
namespace Details {
namespace Todo {
using KokkosBatched::Experimental::Vector;
using KokkosBatched::Experimental::VectorTag;
using KokkosBatched::Experimental::AVX;

template <typename ST1, typename ST2> KOKKOS_INLINE_FUNCTION ST1
add_radial (const ST1& a, const ST2& r) { return a >= 0 ? a+r : a-r; }

#ifdef HAVE_TEUCHOS_COMPLEX
template <typename ST1, typename ST2> KOKKOS_INLINE_FUNCTION Kokkos::complex<ST1>
add_radial (const Kokkos::complex<ST1>& a, const ST2& r) {
  const auto a_abs = Kokkos::abs(a);
  const auto sign = a_abs == Kokkos::ArithTraits<ST1>::zero() ? 1 : a/a_abs;
  return sign*(a_abs + r);
}
#endif

#if defined(__AVX__) || defined(__AVX2__) || defined(__AVX512F__)
inline static __m256d add_radial_impl (const __m256d& a, const double r) {
  // Make the mask 10...0.
  const __m256d mask = _mm256_set1_pd(-0.);
  // Use it to get the sign.
  const __m256d sign = _mm256_and_pd(mask, a);
  // Load r.
  const __m256d pr = _mm256_set1_pd(r);
  // Now use the opposite mask, 01...1, to get b = abs(a).
  __m256d pb = _mm256_andnot_pd(mask, a);
  // Compute abs(a) + r.
  pb = _mm256_add_pd(pb, pr);
  // Bit-or the sign back in.
  pb = _mm256_or_pd(pb, sign);
  return pb;
}
#endif

#if defined(__AVX__) || defined(__AVX2__)
template<typename SpT, typename ScalarType>
inline
static Vector<VectorTag<AVX<double,SpT>,4> >
add_radial (Vector<VectorTag<AVX<double,SpT>,4> > const& a, const ScalarType& r) {
  return add_radial_impl(a, r);
}
#endif

#if defined(__AVX512F__)
template<typename SpT, typename ScalarType>
inline
static Vector<VectorTag<AVX<double,SpT>,8> >
add_radial (Vector<VectorTag<AVX<double,SpT>,8> > const& a, const ScalarType& r) {
  // AVX512F doesn't support _mm512_{and,or,andnot}_pd, so do two 256-bit
  // ops. Might do slightly better by doing the add instruction with a 512-bit
  // operation. It's also possible we could use _mm512_{and,or,andnot}_epi64 by
  // reinterpreting __m512d as __mm512i, but I decided not to do that, as I'm
  // not sure how legal that is.
  const __m256d* pa = (const __m256d*) &a;
  __m512d out;
  __m256d* pout = (__m256d*) &out;
  pout[0] = add_radial_impl(pa[0], r);
  pout[1] = add_radial_impl(pa[1], r);
  return out;
}
#endif

template <typename ValueType, typename ScalarType>
KOKKOS_INLINE_FUNCTION int
Serial_LU_Internal_Unblocked_invoke(
  const int m, const int n, ValueType *__restrict__ A, const int as0, const int as1,
  const ScalarType &diag_safety)
{
  typedef ValueType value_type;
  const int k = (m < n ? m : n);
  if (k <= 0) return 0;

  for (int p=0;p<k;++p) {
    const int iend = m-p-1, jend = n-p-1;

    const int idx = p*(as0+as1);
    A[idx] = add_radial(A[idx], diag_safety);

    const value_type
      alpha11 = A[idx],
      *__restrict__ a12t = A+(p  )*as0+(p+1)*as1;

    value_type
      *__restrict__ a21  = A+(p+1)*as0+(p  )*as1,
      *__restrict__ A22  = A+(p+1)*as0+(p+1)*as1;

    for (int i=0;i<iend;++i) {
      // a21[i*as0] *= inv_alpha11;
      a21[i*as0] /= alpha11;
      for (int j=0;j<jend;++j)
        A22[i*as0+j*as1] -= a21[i*as0] * a12t[j*as1];
    }
  }
  return 0;
}

} // namespace Todo
} // namespace Details
} // namespace KokkosBatched

#endif // __KOKKOSBATCHED_PROMOTION__

namespace Ifpack2 {
#if !defined(__KOKKOSBATCHED_PROMOTION__)
namespace Details {
namespace Todo {

template <typename Object, typename Vector>
int test_add_radial_impl (Object& o, Vector& v, const int len) {
  using KokkosBatched::Details::Todo::add_radial;
  int nerr = 0;
  const double r = 1.7;
  for (const double a : {-1.2, 4.2}) {
    const double t = (a >= 0 ? 1 : -1)*(std::abs(a) + r);
    for (int i = 0; i < len; ++i) v[i] = a;
    o = add_radial(o, r);
    for (int i = 0; i < len; ++i)
      if (v[i] != t)
        ++nerr;
  }
  return nerr;
}

inline int test_add_radial () {
  using KokkosBatched::Experimental::Vector;
  using KokkosBatched::Experimental::VectorTag;
  using KokkosBatched::Experimental::AVX;

  int nerr = 0;
  {
    double v;
    double* v_vector = &v;
    nerr += test_add_radial_impl(v, v_vector, 1);
  }
#if defined(__AVX__) || defined(__AVX2__)
  {
    Vector<VectorTag<AVX<double>, 4> > v;
    nerr += test_add_radial_impl<>(v, v, v.vector_length);
  }
#endif
#if defined(__AVX512F__)
  {
    Vector<VectorTag<AVX<double>, 8> > v;
    nerr += test_add_radial_impl<>(v, v, v.vector_length);
  }
#endif
  return nerr;
}
}
}
#endif // __KOKKOSBATCHED_PROMOTION__

namespace Details {
namespace Batched {
#if defined(__KOKKOSBATCHED_PROMOTION__)
using KokkosBatched::Experimental::Vector;
using KokkosBatched::Experimental::SIMD;
using KokkosBatched::Experimental::DefaultVectorLength;
#else // __KOKKOSBATCHED_PROMOTION__
using KokkosBatched::Experimental::Vector;
using KokkosBatched::Experimental::VectorTag;
using KokkosBatched::Experimental::AVX;
#endif // __KOKKOSBATCHED_PROMOTION__
template <typename ExeSpace>
struct DeviceType {
   typedef Kokkos::Device<typename ExeSpace::execution_space,
                          typename ExeSpace::memory_space> type;
};

typedef Kokkos::Device<Kokkos::DefaultExecutionSpace::execution_space,
                       Kokkos::DefaultExecutionSpace::memory_space> DefaultDeviceType;

template <typename ExeSpace, typename Scalar>
struct BatchedNonVector {
  typedef Scalar value_type;
  typedef ExeSpace exec_space;
};

template <typename ExeSpace, typename Scalar>
struct DefaultVectorizationMethod {
  typedef BatchedNonVector<ExeSpace, Scalar> type;
};

template <typename VectorizationMethod>
struct VectorizationTraits {
  typedef typename VectorizationMethod::value_type value_type;
  typedef typename VectorizationMethod::exec_space exec_space;
  typedef value_type vector_type;
  enum : int { vector_length = 1 };
};

template <typename VectorType>
struct VectorLength {
  enum : int { value = 1 };
};

template <typename VectorType>
struct ScalarType {
  typedef VectorType type;
};

#if defined(__KOKKOSBATCHED_PROMOTION__)
#if defined(KOKKOS_ENABLE_SERIAL)
template <>
struct DefaultVectorizationMethod<Kokkos::Serial, double> {
  typedef SIMD<double> type;
};
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template <>
struct DefaultVectorizationMethod<Kokkos::OpenMP, double> {
  typedef SIMD<double> type;
};
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template <>
struct DefaultVectorizationMethod<Kokkos::Cuda, double> {
  typedef BatchedNonVector<Kokkos::Cuda, double> type;
};
#endif
template <>
struct VectorizationTraits<SIMD<double> > {
  typedef double value_type;
  typedef Kokkos::DefaultHostExecutionSpace exec_space; 
  enum : int { vector_length = DefaultVectorLength<double,typename exec_space::memory_space>::value };
  typedef Vector<SIMD<double>, vector_length> vector_type;
};
template <>
struct VectorLength<VectorizationTraits<SIMD<double> >::vector_type> {
  enum : int { value = VectorizationTraits<SIMD<double> >::vector_length };
};
template <>
struct ScalarType<VectorizationTraits<SIMD<double> >::vector_type> {
  typedef VectorizationTraits<SIMD<double> >::vector_type::value_type type;
};
#else //__KOKKOSBATCHED_PROMOTION__
#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
template <typename ExeSpace>
struct DefaultVectorizationMethod<ExeSpace, double> {
  typedef AVX<double> type;
};
#endif

#if defined(__AVX512F__)
template <>
struct VectorizationTraits<AVX<double> > {
  typedef typename AVX<double>::value_type value_type;
  typedef typename AVX<double>::exec_space exec_space;
  enum : int { vector_length = 8 };
  typedef Vector<VectorTag<AVX<double>, vector_length> > vector_type;
};
#elif defined(__AVX2__) || defined(__AVX__)
template <>
struct VectorizationTraits<AVX<double> > {
  typedef typename AVX<double>::value_type value_type;
  typedef typename AVX<double>::exec_space exec_space;
  enum : int { vector_length = 4 };
  typedef Vector<VectorTag<AVX<double>, vector_length> > vector_type;
};
#endif

#if defined(__AVX512F__) || defined(__AVX2__) || defined(__AVX__)
template <>
struct VectorLength<VectorizationTraits<AVX<double> >::vector_type> {
  enum : int { value = VectorizationTraits<AVX<double> >::vector_length };
};
template <>
struct ScalarType<VectorizationTraits<AVX<double> >::vector_type> {
  typedef VectorizationTraits<AVX<double> >::vector_type::value_type type;
};
#endif
#endif //__KOKKOSBATCHED_PROMOTION__


template <typename Int1, typename Int2, typename View, typename... Indices>
KOKKOS_INLINE_FUNCTION
typename std::enable_if<VectorLength<typename View::value_type>::value == 1,
                        typename View::reference_type>::type
idx (const View& v, const Kokkos::pair<Int1,Int2>& i0, Indices... indices) {
  return v(i0.first, indices...);
}

template <typename Int1, typename Int2, typename View, typename... Indices>
KOKKOS_INLINE_FUNCTION
typename std::enable_if<(VectorLength<typename View::value_type>::value > 1),
                        typename View::value_type::value_type&>::type
idx (const View& v, const Kokkos::pair<Int1,Int2>& i0, Indices... indices) {
  return v(i0.first, indices...)[i0.second % VectorLength<typename View::value_type>::value];
}

// No %, but you have to do the % yourself.
template <typename Int1, typename Int2, typename View, typename... Indices>
KOKKOS_FORCEINLINE_FUNCTION
typename std::enable_if<VectorLength<typename View::value_type>::value == 1,
                        typename View::reference_type>::type
fastidx (const View& v, const Kokkos::pair<Int1,Int2>& i0, Indices... indices) {
  return v(i0.first, indices...);
}

template <typename Int1, typename Int2, typename View, typename... Indices>
KOKKOS_FORCEINLINE_FUNCTION
typename std::enable_if<(VectorLength<typename View::value_type>::value > 1),
                        typename View::value_type::value_type&>::type
fastidx (const View& v, const Kokkos::pair<Int1,Int2>& i0, Indices... indices) {
  return v(i0.first, indices...)[i0.second];
}

} // namespace Batched
} // namespace Details

namespace BlockTriDiContainerDetails {
// Implementations of the operation Y = B - R*X, where A = D + R. There are
// three implementations, each described below.

template <typename Scalar> KOKKOS_INLINE_FUNCTION
void copy (const int n, const Scalar* b, const int ldb, Scalar* y, const int ldy, const int nvec) {
  int vi = 0;
  for (;;) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#   pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#   pragma unroll
#endif
    for (int i = 0; i < n; ++i)
      y[i] = b[i];
    ++vi;
    if (vi == nvec) break;
    b += ldb;
    y += ldy;
  }
}

template <int BSZ, typename Scalar> KOKKOS_INLINE_FUNCTION
void copy (const int bsz, const Scalar* const x, const int ldx, Scalar* const y, const int ldy,
           const int nvec) {
  if (BSZ == 0)
    copy(bsz, x, ldx, y, ldy, nvec);
  else
    copy(BSZ, x, ldx, y, ldy, nvec);
}

template <int sign, typename Scalar> KOKKOS_INLINE_FUNCTION
void gemv (const Scalar* const A, const int n, const Scalar* x, const int ldx, Scalar* y,
           const int ldy, const int nvec) {
  int vi = 0;
  for (;;) {
    for (int r = 0; r < n; ++r) {
      Scalar accum = 0;
      const Scalar* const Ar = A + r*n;
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#     pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#     pragma unroll
#endif
      for (int c = 0; c < n; ++c)
        accum += Ar[c] * x[c];
      y[r] += static_cast<Scalar>(sign) * accum;
    }
    ++vi;
    if (vi == nvec) break;
    x += ldx;
    y += ldy;
  }
}

template <int BSZ, int sign, typename Scalar> KOKKOS_INLINE_FUNCTION
void gemv (const Scalar* const A, const int bsz, const Scalar* const x, const int ldx,
           Scalar* const y, const int ldy, const int nvec) {
  if (BSZ == 0)
    gemv<sign>(A, bsz, x, ldx, y, ldy, nvec);
  else
    gemv<sign>(A, BSZ, x, ldx, y, ldy, nvec);
}

template <typename ExeSpace, int BSZ, typename Scalar, typename Int, typename Size, typename Colidx>
class ComputeBMinusRX {
  const Scalar* const b;
  const Int ldb;
  const Scalar* const x;
  const Int ldx;
  Scalar* const y;
  const Int ldy;
  const Int nvec;
  const int block_size;
  const Int nrows;
  const size_t* const rowptr;
  const Colidx* const colindsub;
  const Size* const A_rowptr;
  const Int* const A_colind;
  const Scalar* const A_values;
  const Int bs2;

public:
  ComputeBMinusRX (
    const Scalar* const b_, const Int ldb_, const Scalar* const x_, const Int ldx_,
    Scalar* const y_, const Int ldy_, const Int nvec_, const int block_size_,
    const Int nrows_, const size_t* const rowptr_, const Colidx* const colindsub_,
    const Size* const A_rowptr_, const Int* const A_colind_, const Scalar* const A_values_) :
    b(b_), ldb(ldb_), x(x_), ldx(ldx_), y(y_), ldy(ldy_), nvec(nvec_), block_size(block_size_),
    nrows(nrows_), rowptr(rowptr_), colindsub(colindsub_), A_rowptr(A_rowptr_), A_colind(A_colind_),
    A_values(A_values_), bs2(block_size*block_size)
  {
    Kokkos::parallel_for(Kokkos::RangePolicy<ExeSpace>(0, nrows), *this);
  }

  KOKKOS_INLINE_FUNCTION void operator() (const Int& lr) const {
    const Int os = lr*block_size;
    const Scalar* const br = b + os;
    Scalar* const yr = y + os;
    // y := b
    copy<BSZ>(block_size, br, ldb, yr, ldy, nvec);
    // y -= R x
    const size_t A_k0 = A_rowptr[lr];
    for (size_t k = rowptr[lr]; k < rowptr[lr+1]; ++k) {
      const size_t j = A_k0 + colindsub[k];
      gemv<BSZ, -1>(A_values + j*bs2, block_size, x + A_colind[j]*block_size, ldx, yr, ldy, nvec);
    }
  }
};

template <typename ExeSpace, int BSZ, typename Scalar, typename Int, typename Size, typename Colidx>
void compute_b_minus_Rx (const Scalar* const b, const Int ldb, const Scalar* const x, const Int ldx,
                         Scalar* const y, const Int ldy, const Int nvec, const int block_size,
                         const Int nrows, const size_t* const rowptr, const Colidx* const colindsub,
                         const Size* const A_rowptr, const Int* const A_colind,
                         const Scalar* const A_values) {
  ComputeBMinusRX<ExeSpace, BSZ, Scalar, Int, Size, Colidx>(
    b, ldb, x, ldx, y, ldy, nvec, block_size, nrows, rowptr, colindsub, A_rowptr, A_colind, A_values);
}

// This first version is for the seq_method_ = true code path. This is the
// straightforward operation y = b - R*x, where A = (D + R) and D is the
// preconditioner.
template <typename MvView,
          typename RowptrView, typename ColindView,
          typename ARowptrView, typename ValuesView>
void compute_b_minus_Rx (MvView& B, MvView& X, MvView& Y, const int block_size,
                         const RowptrView& rowptrv, const ColindView& colindsubv,
                         const ARowptrView& A_rowptrv, const ColindView& A_colindv,
                         const ValuesView& valuesv) {
  //todo Replace this placeholder impl with Kokkos-based placeholders, to be
  // replaced later by a good impl of matvec.
  typedef typename ColindView::value_type LO;
  const LO nrows = rowptrv.size() - 1;
  const LO nvec = X.extent(1);
  const auto* rowptr = rowptrv.data();
  const auto* colindsub = colindsubv.data();
  const auto* A_rowptr = A_rowptrv.data();
  const auto* A_colind = A_colindv.data();
  const auto* values = valuesv.data();
  const LO ldb = B.extent(0), ldx = X.extent(0), ldy = Y.extent(0);
#define BLOCKTRIDICONTAINERDETAILS_CBMR(N)                              \
  compute_b_minus_Rx<typename MvView::execution_space, N>(              \
    B.data(), ldb, X.data(), ldx, Y.data(), ldy, nvec, block_size,      \
    nrows, rowptr, colindsub, A_rowptr, A_colind, values);              \
  break
  switch (block_size) {
  case 5 : BLOCKTRIDICONTAINERDETAILS_CBMR(5);
  case 9 : BLOCKTRIDICONTAINERDETAILS_CBMR(9);
  default: BLOCKTRIDICONTAINERDETAILS_CBMR(0);
  }
#undef BLOCKTRIDICONTAINERDETAILS_CBMR
}

template <typename Int, typename Scalar, typename ConstMvView> KOKKOS_INLINE_FUNCTION
void copy (const Int& bs, const ConstMvView& X, const Int& xos, const Int& col, Scalar* const y) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
# pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
# pragma unroll
#endif
  for (Int i = 0; i < bs; ++i)
    y[i] = X(xos+i, col);
}

template <int BSZ, typename Int, typename Scalar, typename ConstMvView> KOKKOS_INLINE_FUNCTION
void copy (const Int& bs, const ConstMvView& X, const Int& xos, const Int& col, Scalar* const y) {
  if (BSZ == 0)
    copy(bs,  X, xos, col, y);
  else
    copy(BSZ, X, xos, col, y);
}

template <int sign, typename Int, typename Scalar> KOKKOS_INLINE_FUNCTION
void gemv (const Int& bs, const Int& col, const Scalar* const A, const Scalar* const x,
           Scalar* const y) {
  for (Int r = 0; r < bs; ++r) {
    const Scalar* const Ar = A + r*bs;
    Scalar accum = 0;
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#   pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#   pragma unroll
#endif
    for (Int c = 0; c < bs; ++c)
      accum += Ar[c] * x[c];
    y[r] += static_cast<Scalar>(sign) * accum;
  }
}

template <int BSZ, int sign, typename Int, typename Scalar> KOKKOS_INLINE_FUNCTION
void gemv (const Int& bs, const Int& col, const Scalar* const A, const Scalar* const x,
           Scalar* const y) {
  if (BSZ == 0)
    gemv<sign>(bs,  col, A, x, y);
  else
    gemv<sign>(BSZ, col, A, x, y);
}

template <int BSZ, typename Int, typename Size, typename Scalar,
          typename MvView, typename PackedView, typename ConstLOList>
class ComputeBMinusRXOverlap {
  const size_t* const rowptr;
  const Int* const colindsub;
  const Size* const A_rowptr;
  const Int* const A_colind;
  const Scalar* const A_values;
  const MvView B;
  const Scalar* const x;
  const Int ldx;
  const PackedView Y;
  const ConstLOList part2rowidx0, part2packrowidx0, lclrow;
  const Int* const dm2cm;
  const bool first_apply;
  const Int bs, bs2, nvec, nparts;
  const Int nlclrow;

  static constexpr Int bufsz = BSZ > 0 ? BSZ : 64;

public:
  ComputeBMinusRXOverlap (
    const size_t* const rowptr_, const Int* const colindsub_,
    const Size* const A_rowptr_, const Int* const A_colind_, const Scalar* const A_values_,
    const MvView& B_, const Scalar* const x_, const Int& ldx_, const PackedView& Y_,
    const ConstLOList& part2rowidx0_, const ConstLOList& part2packrowidx0_,
    const ConstLOList& lclrow_, const Int* const dm2cm_, const bool first_apply_) :
    rowptr(rowptr_), colindsub(colindsub_), A_rowptr(A_rowptr_), A_colind(A_colind_), A_values(A_values_),
    B(B_), x(x_), ldx(ldx_), Y(Y_), part2rowidx0(part2rowidx0_), part2packrowidx0(part2packrowidx0_),
    lclrow(lclrow_), dm2cm(dm2cm_), first_apply(first_apply_), bs(Y.extent(1)), bs2(bs*bs),
    nvec(Y.extent(2)), nparts(part2rowidx0.size() - 1), nlclrow(lclrow.size())
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(bs > bufsz, "bufsz is too small; contact the developer.");
    Kokkos::parallel_for(Kokkos::RangePolicy<typename MvView::execution_space>(0, nparts), *this);
  }

  KOKKOS_INLINE_FUNCTION void operator() (const Int& partidx) const {
    Scalar ybuf[bufsz] = {0};
    const Int pri0 = part2packrowidx0(partidx);
    const Int ri0 = part2rowidx0(partidx);
    const Int nrows = part2rowidx0(partidx+1) - ri0;
    // Get the position within the pack.
    const Int vi = partidx % Details::Batched::VectorLength<typename PackedView::value_type>::value;
    for (Int i = 0; i < nrows; ++i) {
      Kokkos::pair<Int,Int> yp(pri0 + i, vi);
      const Int ri = ri0 + i;
      const Int lr = lclrow(ri);
      for (Int col = 0; col < nvec; ++col) {
        if (first_apply) {
          // y := b
          copy<BSZ>(bs, B, lr*bs, col, ybuf);
          // y -= R x_owned
          const size_t A_k0 = A_rowptr[lr];
          for (size_t k = rowptr[lr]; k < rowptr[lr+1]; ++k) {
            const size_t j = A_k0 + colindsub[k];
            const auto Acj = A_colind[j];
            const auto xos = dm2cm ? dm2cm[Acj] : Acj;
            gemv<BSZ, -1>(bs, col, A_values + j*bs2, x + col*ldx + xos*bs, ybuf);
          }
          for (Int r = 0; r < bs; ++r) Details::Batched::fastidx(Y, yp, r, col)  = ybuf[r];
        } else {
          for (Int r = 0; r < bs; ++r) ybuf[r] = 0;
          const size_t A_k0 = A_rowptr[lr];
          // y -= R x_remote
          for (size_t k = rowptr[lr]; k < rowptr[lr+1]; ++k) {
            const size_t j = A_k0 + colindsub[k];
            const auto Acj = A_colind[j] - nlclrow;
            gemv<BSZ, -1>(bs, col, A_values + j*bs2, x + col*ldx + Acj*bs, ybuf);
          }
          for (Int r = 0; r < bs; ++r)
            Details::Batched::fastidx(Y, yp, r, col) += ybuf[r];
        }
      }
    }
  }
};

template <int BSZ, typename Int, typename Size, typename Scalar,
          typename MvView, typename PackedView, typename ConstLOList>
void compute_b_minus_Rx (
  const size_t* const rowptr, const Int* const colindsub,
  const Size* const A_rowptr, const Int* const A_colind, const Scalar* const A_values,
  const MvView& B, const Scalar* const x, const Int& ldx, const PackedView& Y,
  const ConstLOList& part2rowidx0, const ConstLOList& part2packrowidx0, const ConstLOList& lclrow,
  const Int* const dm2cm, const bool first_apply)
{
  ComputeBMinusRXOverlap<BSZ, Int, Size, Scalar, MvView, PackedView, ConstLOList>(
    rowptr, colindsub, A_rowptr, A_colind, A_values, B, x, ldx, Y, part2rowidx0, part2packrowidx0,
    lclrow, dm2cm, first_apply);
}

// This second version is to overlap communication and computation.
//   If first_apply, then y = b - R_owned*x_owned is performed, where _owned
// means that the data are for the owned GIDs.
//   If ! apply_first, then y -= R_remote*x_remote is performed. x_remote was
// being communicated while x_owned was being used and is now available for
// computation.
//   dm2cm is the permutation of owned IDs from domain map to column map. If it
// is the identity permutation, then it can be empty.
template <typename RowptrView, typename ColindView, typename ARowptrView, typename ValuesView,
          typename MvView, typename PackedView, typename ConstLOList>
void compute_b_minus_Rx (
  const RowptrView& rowptrv, const ColindView& colindsubv,
  const ARowptrView& A_rowptrv, const ColindView& A_colindv, const ValuesView& valuesv,
  const MvView& B, const MvView& X, const PackedView& Y,
  const ConstLOList& part2rowidx0, const ConstLOList& part2packrowidx0, const ConstLOList& lclrow,
  const ConstLOList& dm2cmv, const bool first_apply)
{
  typedef typename ColindView::value_type LO;
  const LO block_size = Y.extent(1);
  const auto* rowptr = rowptrv.data();
  const auto* colindsub = colindsubv.data();
  const auto* A_rowptr = A_rowptrv.data();
  const auto* A_colind = A_colindv.data();
  const auto* values = valuesv.data();
  const auto* dm2cm = dm2cmv.data();
  const LO ldx = X.extent(0);
#define BLOCKTRIDICONTAINERDETAILS_CBMR(N)                              \
  compute_b_minus_Rx<N>(rowptr, colindsub, A_rowptr, A_colind, values, B, X.data(), ldx, Y, \
                        part2rowidx0, part2packrowidx0, lclrow, dm2cm, first_apply); \
  break
  switch (block_size) {
  case 5 : BLOCKTRIDICONTAINERDETAILS_CBMR(5);
  case 9 : BLOCKTRIDICONTAINERDETAILS_CBMR(9);
  default: BLOCKTRIDICONTAINERDETAILS_CBMR(0);
  }
#undef BLOCKTRIDICONTAINERDETAILS_CBMR
}

template <int BSZ, typename Int, typename Size, typename Scalar,
          typename MvView, typename PackedView, typename ConstLOList>
class ComputeBMinusRXOwnedRemote {
  const size_t* const rowptr;
  const Int* const colindsub;
  const Size* const A_rowptr;
  const Int* const A_colind;
  const Scalar* const A_values;
  const MvView B;
  const Scalar* const xo;
  const Int ldxo;
  const Scalar* const xr;
  const Int ldxr;
  const PackedView Y;
  const ConstLOList part2rowidx0, part2packrowidx0, lclrow;
  const Int* const dm2cm;
  const Int bs, bs2, nvec, nparts;
  const Int nlclrow;

  static constexpr Int bufsz = BSZ > 0 ? BSZ : 64;

public:
  ComputeBMinusRXOwnedRemote (
    const size_t* const rowptr_, const Int* const colindsub_,
    const Size* const A_rowptr_, const Int* const A_colind_, const Scalar* const A_values_,
    const MvView& B_, const Scalar* const xo_, const Int& ldxo_, const Scalar* const xr_, const Int& ldxr_,
    const PackedView& Y_, const ConstLOList& part2rowidx0_, const ConstLOList& part2packrowidx0_,
    const ConstLOList& lclrow_, const Int* const dm2cm_) :
    rowptr(rowptr_), colindsub(colindsub_), A_rowptr(A_rowptr_), A_colind(A_colind_), A_values(A_values_),
    B(B_), xo(xo_), ldxo(ldxo_), xr(xr_), ldxr(ldxr_), Y(Y_), part2rowidx0(part2rowidx0_),
    part2packrowidx0(part2packrowidx0_), lclrow(lclrow_), dm2cm(dm2cm_), bs(Y.extent(1)), bs2(bs*bs),
    nvec(Y.extent(2)), nparts(part2rowidx0.size() - 1), nlclrow(lclrow.size())
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(bs > bufsz, "bufsz is too small; contact the developer.");
    Kokkos::parallel_for(Kokkos::RangePolicy<typename MvView::execution_space>(0, nparts), *this);
  }

  KOKKOS_INLINE_FUNCTION void operator() (const Int& partidx) const {
    Scalar ybuf[bufsz];
    const Int pri0 = part2packrowidx0(partidx);
    const Int ri0 = part2rowidx0(partidx);
    const Int nrows = part2rowidx0(partidx+1) - ri0;
    // Get the position within the pack.
    const Int vi = partidx % Details::Batched::VectorLength<typename PackedView::value_type>::value;
    for (Int i = 0; i < nrows; ++i) {
      Kokkos::pair<Int,Int> yp(pri0 + i, vi);
      const Int ri = ri0 + i;
      const Int lr = lclrow(ri);
      for (Int col = 0; col < nvec; ++col) {
        // y := b
        copy<BSZ>(bs, B, lr*bs, col, ybuf);
        // y -= R x_owned
        const size_t A_k0 = A_rowptr[lr];
        for (size_t k = rowptr[lr]; k < rowptr[lr+1]; ++k) {
          const size_t j = A_k0 + colindsub[k];
          const auto Acj = A_colind[j];
          if (Acj < nlclrow) {
            const auto xos = dm2cm ? dm2cm[Acj] : Acj;
            gemv<BSZ, -1>(bs, col, A_values + j*bs2, xo + col*ldxo + xos*bs, ybuf);
          } else {
            const auto xos = Acj - nlclrow;
            gemv<BSZ, -1>(bs, col, A_values + j*bs2, xr + col*ldxr + xos*bs, ybuf);
          }
        }
        for (Int r = 0; r < bs; ++r)
          Details::Batched::fastidx(Y, yp, r, col) = ybuf[r];
      }
    }
  }
};

template <int BSZ, typename Int, typename Size, typename Scalar,
          typename MvView, typename PackedView, typename ConstLOList>
void compute_b_minus_Rx (
  const size_t* const rowptr, const Int* const colindsub,
  const Size* const A_rowptr, const Int* const A_colind, const Scalar* const A_values,
  const MvView& B, const Scalar* const xo, const Int& ldxo, const Scalar* const xr, const Int& ldxr,
  const PackedView& Y, const ConstLOList& part2rowidx0, const ConstLOList& part2packrowidx0,
  const ConstLOList& lclrow, const Int* const dm2cm)
{
  ComputeBMinusRXOwnedRemote<BSZ, Int, Size, Scalar, MvView, PackedView, ConstLOList>(
    rowptr, colindsub, A_rowptr, A_colind, A_values, B, xo, ldxo, xr, ldxr, Y, part2rowidx0,
    part2packrowidx0, lclrow, dm2cm);
}

// This third version is similar to the second in that it distinguishes between
// _owned and _remote data. But y = b - R*x is done all at once, just with x in
// two separate arrays.
template <typename RowptrView, typename ColindView, typename ARowptrView, typename ValuesView,
          typename MvView, typename PackedView, typename ConstLOList>
void compute_b_minus_Rx (
  const RowptrView& rowptrv, const ColindView& colindsubv,
  const ARowptrView& A_rowptrv, const ColindView& A_colindv, const ValuesView& valuesv,
  const MvView& B, const MvView& X_owned, const MvView& X_remote, const PackedView& Y,
  const ConstLOList& part2rowidx0, const ConstLOList& part2packrowidx0, const ConstLOList& lclrow,
  const ConstLOList& dm2cmv)
{
  typedef typename ColindView::value_type LO;
  const LO block_size = Y.extent(1);
  const auto* rowptr = rowptrv.data();
  const auto* colindsub = colindsubv.data();
  const auto* A_rowptr = A_rowptrv.data();
  const auto* A_colind = A_colindv.data();
  const auto* values = valuesv.data();
  const auto* dm2cm = dm2cmv.data();
  const LO ldxo = X_owned.extent(0);
  const LO ldxr = X_remote.extent(0);
#define BLOCKTRIDICONTAINERDETAILS_CBMR(N)                              \
  compute_b_minus_Rx<N>(rowptr, colindsub, A_rowptr, A_colind, values, B, \
                        X_owned.data(), ldxo, X_remote.data(), ldxr,    \
                        Y, part2rowidx0, part2packrowidx0, lclrow, dm2cm); \
  break
  // This compute_b_minus_Rx implementation is the most important, so
  // instantiate several more sizes requested by others.
  switch (block_size) {
  case 5 : BLOCKTRIDICONTAINERDETAILS_CBMR(5);
  case 6 : BLOCKTRIDICONTAINERDETAILS_CBMR(6);
  case 7 : BLOCKTRIDICONTAINERDETAILS_CBMR(7);
  case 9 : BLOCKTRIDICONTAINERDETAILS_CBMR(9);
  case 10: BLOCKTRIDICONTAINERDETAILS_CBMR(10);
  case 11: BLOCKTRIDICONTAINERDETAILS_CBMR(11);
  case 16: BLOCKTRIDICONTAINERDETAILS_CBMR(16);
  case 17: BLOCKTRIDICONTAINERDETAILS_CBMR(17);
  case 18: BLOCKTRIDICONTAINERDETAILS_CBMR(18);
  default: BLOCKTRIDICONTAINERDETAILS_CBMR(0);
  }
#undef BLOCKTRIDICONTAINERDETAILS_CBMR
}

template <typename T> KOKKOS_INLINE_FUNCTION
typename Kokkos::ArithTraits<T>::mag_type abs2 (T& v) {
  const typename Kokkos::ArithTraits<T>::mag_type av = Kokkos::ArithTraits<T>::abs(v);
  return av*av;
}
template <> KOKKOS_INLINE_FUNCTION double abs2 (double& v) { return v*v; }
template <> KOKKOS_INLINE_FUNCTION float abs2 (float& v) { return v*v; }

template <typename MatrixType, typename LocalScalarType, typename ExeSpace>
class Impl {
private:
  typedef typename MatrixType::scalar_type scalar_type;
  // See discussion in Tpetra_MultiVector_decl.hpp. In short, this is to map
  // std::complex to Kokkos::complex for the Cuda build.
  typedef typename Kokkos::Details::ArithTraits<scalar_type>::val_type impl_scalar_type;
  typedef typename Kokkos::ArithTraits<impl_scalar_type>::mag_type magnitude_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> mv_type;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;

  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;
  typedef Tpetra::BlockCrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>
          block_crs_matrix_type;

  static_assert(std::is_same<scalar_type, LocalScalarType>::value,
                "BlockTriDiContainer supports only LocalScalarType same as scalar_type.");

  typedef typename MatrixType::node_type::device_type tpetra_device_type;

  typedef local_ordinal_type LO;
  typedef global_ordinal_type GO;
  typedef size_t Size;

  typedef typename block_crs_matrix_type::crs_graph_type::local_graph_type::row_map_type
    CrsGraphRowMap;
  typedef Kokkos::View<impl_scalar_type*, tpetra_device_type> BlockCrsValues;

  typedef typename Kokkos::View<double***>::array_layout default_layout_type;
  typedef default_layout_type layout_type;

  typedef typename Details::Batched::DeviceType<typename tpetra_device_type::execution_space>::type
    device_type;
#ifdef HAVE_IFPACK2_CUDA
  enum : bool {
    device_is_cuda = std::is_same<typename device_type::execution_space, Kokkos::Cuda>::value
  };
#else
  enum : bool { device_is_cuda = false };
#endif

  typedef typename Details::Batched::DefaultVectorizationMethod<
    typename tpetra_device_type::execution_space, impl_scalar_type>::type vectorization_method;
  typedef Details::Batched::VectorizationTraits<vectorization_method> vectorization_traits;
  typedef typename vectorization_traits::vector_type vector_type;

  // Some Kokkos::View decorators for unmanaged and const memory.
  template <typename MemoryTraitsType, Kokkos::MemoryTraitsFlags flag>
  using MemoryTraits = Kokkos::MemoryTraits<MemoryTraitsType::Unmanaged |
                                            MemoryTraitsType::RandomAccess |
                                            MemoryTraitsType::Atomic |
                                            flag>;

  template <typename View>
  using Unmanaged = Kokkos::View<typename View::data_type,
                                 typename View::array_layout,
                                 typename View::device_type,
                                 MemoryTraits<typename View::memory_traits,
                                              Kokkos::Unmanaged> >;
  template <typename View>
  using Const = Kokkos::View<typename View::const_data_type, typename View::array_layout,
                             typename View::device_type, typename View::memory_traits>;
  template <typename View>
  using ConstUnmanaged = Const<Unmanaged<View> >;

  typedef Kokkos::View<Size*, Kokkos::LayoutLeft, device_type> SizeList;
  typedef Const<SizeList> ConstSizeList;
  typedef Kokkos::View<LO*, Kokkos::LayoutLeft, device_type> LOList;
  typedef Const<LOList> ConstLOList;

private:

  // The block tridiagonal matrices extracted from A, i.e., D in the splitting A
  // = D + R.
  struct Tridiags {
    typedef Kokkos::View<vector_type***, layout_type, device_type> Values;
    typedef Kokkos::View<vector_type***, layout_type, device_type> PackedMultiVector;

    // flat_td_ptr(i) is the index into flat-array values of the start of the
    // i'th tridiag. pack_td_ptr is the same, but for packs. If vector_length ==
    // 1, pack_td_ptr is the same as flat_td_ptr; if vector_length > 1, then i %
    // vector_length is the position in the pack.
    SizeList flat_td_ptr, pack_td_ptr;
    // List of local column indices into A from which to grab
    // data. flat_td_ptr(i) points to the start of the i'th tridiag's data.
    LOList A_colindsub;
    // Tridiag block values. pack_td_ptr(i) points to the start of the i'th
    // tridiag's pack, and i % vector_length gives the position in the pack.
    Values values;

    // Index into row-major block of a tridiag.
    template <typename Int> static KOKKOS_INLINE_FUNCTION
    Int ind2row (const Int& ind) { return (ind + 1) / 3; }
    // Given a row of a row-major tridiag, return the index of the first block
    // in that row.
    template <typename Int> static KOKKOS_INLINE_FUNCTION
    Int row2ind (const Int& row) { return row > 0 ? 3*row - 1 : 0; }
    // Number of blocks in a tridiag having a given number of rows.
    template <typename Int> static KOKKOS_INLINE_FUNCTION
    Int nblks (const Int& nrows) { return nrows > 0 ? 3*nrows - 2 : 0; }
  };

  // A - Tridiags(A), i.e., R in the splitting A = D + R.
  struct AmD {
    // rowptr points to the start of each row of A_colindsub.
    SizeList rowptr, rowptr_remote;
    // Indices into A's rows giving the blocks to extract. rowptr(i) points to
    // the i'th row. Thus, g.entries(A_colindsub(rowptr(row) : rowptr(row+1))),
    // where g is A's graph, are the columns AmD uses. If seq_method_, then
    // A_colindsub contains all the LIDs and A_colindsub_remote is empty. If !
    // seq_method_, then A_colindsub contains owned LIDs and A_colindsub_remote
    // contains the remote ones.
    LOList A_colindsub, A_colindsub_remote;

    // Currently always true.
    bool is_tpetra_block_crs;
    // If is_tpetra_block_crs, then this is a pointer to A_'s value data.
    BlockCrsValues tpetra_values;
  };

public: // for Cuda

  class AsyncableImport;
  class NormManager;

private:
  BlockTriDiContainer<MatrixType, LocalScalarType>& container_;
  Teuchos::RCP<const block_crs_matrix_type> A_bcrs_;
  // Use the slow, straightforward import(dm->cm), matvec, solve pattern?
  bool seq_method_;

  // Some terms:
  //   The matrix A is split as A = D + R, where D is the matrix of tridiag
  // blocks and R is the remainder.
  //   A part is roughly a synonym for a tridiag. The distinction is that a part
  // is the set of rows belonging to one tridiag and, equivalently, the off-diag
  // rows in R associated with that tridiag. In contrast, the term tridiag is
  // used to refer specifically to tridiag data, such as the pointer into the
  // tridiag data array.
  //   Local (lcl) row arge the LIDs. lclrow lists the LIDs belonging to each
  // tridiag, and partptr points to the beginning of each tridiag. This is the
  // LID space.
  //   Row index (idx) is the ordinal in the tridiag ordering. lclrow is indexed
  // by this ordinal. This is the 'index' space.
  //   A flat index is the mathematical index into an array. A pack index
  // accounts for SIMD packing.

  // Local row LIDs. Permutation from caller's index space to tridiag index
  // space.
  LOList lclrow_;
  // partptr_ is the pointer array into lclrow_.
  LOList partptr_; // np+1
  // packptr_(i), for i the pack index, indexes partptr_. partptr_(packptr_(i))
  // is the start of the i'th pack.
  LOList packptr_; // npack+1
  // part2rowidx0_(i) is the flat row index of the start of the i'th part. It's
  // an alias of partptr_ in the case of no overlap.
  LOList part2rowidx0_; // np+1
  // part2packrowidx0_(i) is the packed row index. If vector_length is 1, then
  // it's the same as part2rowidx0_; if it's > 1, then the value is combined
  // with i % vector_length to get the location in the packed data.
  LOList part2packrowidx0_; // np+1
  LO part2packrowidx0_back; // So we don't need to grab the array from the GPU.
  // rowidx2part_ maps the row index to the part index.
  LOList rowidx2part_; // nr
  // If needed, permute owned domain main LIDs to owned column map LIDs.
  LOList dm2cm_;
  // True if lcl{row|col} is at most a constant away from row{idx|col}. In
  // practice, this knowledge is not particularly useful, as packing for batched
  // processing is done at the same time as the permutation from LID to index
  // space. But it's easy to detect, so it's recorded in case an optimization
  // can be made based on it.
  bool row_contiguous_, col_contiguous_;

  Tridiags btdm_;
  AmD amd_;

  bool overlap_comm_;
  Teuchos::RCP<const import_type> importer_;
  std::shared_ptr<AsyncableImport> async_import_;
  mutable Teuchos::RCP<mv_type> Y_remote_;
  mutable typename Tridiags::PackedMultiVector pmv_;

  mutable bool have_norms_;
  mutable typename NormManager::Ptr norm_mgr_;

  bool validate_;

public: // Intended to be private, but public for Cuda.

  // Set the tridiags to be I to the full pack block size. That way, if a
  // tridiag within a pack is shorter than the longest one, the extra blocks are
  // processed in a safe way. Similarly, in the solve phase, if the extra blocks
  // in the packed multvector are 0, and the tridiag LU reflects the extra I
  // blocks, then the solve proceeds as though the extra blocks aren't
  // present. Since this extra work is part of the SIMD calls, it's not actually
  // extra work. Instead, it means we don't have to put checks or masks in, or
  // quiet NaNs. This functor has to be called just once, in the symbolic phase,
  // since the numeric phase fills in only the used entries, leaving these I
  // blocks intact.
  class SetTridiagsToI {
    Tridiags btdm;
    LOList packptr;

  public:
    SetTridiagsToI (const Tridiags& btdm_, const LOList& packptr_)
      : btdm(btdm_), packptr(packptr_)
    {}

    KOKKOS_INLINE_FUNCTION void operator() (const LO& k) const {
      const Size is = btdm.pack_td_ptr(packptr(k)), ie = btdm.pack_td_ptr(packptr(k+1));
      const LO bs = btdm.values.extent(1);
      for (Size i = is; i < ie; i += 3)
        for (LO j = 0; j < bs; ++j)
          btdm.values(i,j,j) = 1;
    }

    void run () {
      Kokkos::RangePolicy<typename device_type::execution_space> rp(0, packptr.size() - 1);
      Kokkos::parallel_for(rp, *this);
    }
  };

  // Partial replacement for forward-mode MultiVector::doImport. Permits
  // overlapped communication and computation, but also supports sync'ed. I'm
  // finding that overlapped comm/comp can give quite poor performance on some
  // platforms, so we can't just use it straightforwardly always.
  class AsyncableImport {
    typedef Kokkos::View<impl_scalar_type*, layout_type, device_type> Buffer;
    typedef typename mv_type::dual_view_type::t_dev MvView;

    // Need to use raw MPI so I can use waitany on a vector of MPI_Requests.
#ifndef HAVE_MPI
    typedef int MPI_Request;
    typedef int MPI_Comm;
#endif

    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    MvView X_remote_;

    // (s)end and (r)eceive data:
    // MPI ranks.
    LO bsz_;
    Teuchos::Array<int> spids_, rpids_;
    std::vector<MPI_Request> sreqs_, rreqs_;
    // Offsets into lists and buffers.
    std::vector<size_t> sos_, ros_;
    // Lists of LIDs.
    LOList slids_, rlids_;
    // Buffer for MPI.
    Buffer sbuf_, rbuf_;

  private:

    void init_mem (const LO& nvec) {
      if (X_remote_.extent(0) != rlids_.size()*bsz_ ||
          static_cast<LO>(X_remote_.extent(1)) != nvec) {
        X_remote_ = MvView("X_remote", rlids_.size()*bsz_, nvec);
        sbuf_ = Buffer("sbuf", sos_.back()*bsz_*nvec);
        rbuf_ = Buffer("rbuf", ros_.back()*bsz_*nvec);
      }
    }

    void init (const Teuchos::ArrayView<const size_t>& lens, std::vector<size_t>& os) {
      os.resize(lens.size() + 1);
      os[0] = 0;
      for (int i = 0; i < lens.size(); ++i)
        os[i+1] = os[i] + lens[i];
    }

    void init (const Teuchos::RCP<const map_type>& src_map,
               const Teuchos::RCP<const map_type>& tgt_map, const LO bsz) {
      bsz_ = bsz;
      comm_ = tgt_map->getComm();

      const import_type import(src_map, tgt_map);

      { // Send and recv buffers by PID.
        Tpetra::Distributor& d = import.getDistributor();
        init(d.getLengthsTo(), sos_);
        init(d.getLengthsFrom(), ros_);
        rpids_ = d.getProcsFrom();
        rreqs_.resize(rpids_.size());
        spids_ = d.getProcsTo();
        sreqs_.resize(spids_.size());
      }

      { // For each remote PID, the list of LIDs to receive.
        const auto remote_lids = import.getRemoteLIDs();
        rlids_ = LOList("rlids", remote_lids.size());
        const auto rlids = Kokkos::create_mirror_view(rlids_);
        std::copy(remote_lids.begin(), remote_lids.end(), rlids.data());
        Kokkos::deep_copy(rlids_, rlids);
      }

      { // For each export PID, the list of LIDs to send.
        auto epids = import.getExportPIDs();
        auto elids = import.getExportLIDs();
        TEUCHOS_ASSERT(epids.size() == elids.size());
        slids_ = LOList("slids", elids.size());
        const auto slids = Kokkos::create_mirror_view(slids_);
        int cnt = 0;
        for (int i = 0; i < spids_.size(); ++i) {
          for (int j = 0; j < epids.size(); ++j)
            if (epids[j] == spids_[i])
              slids[cnt++] = elids[j];
          TEUCHOS_ASSERT(static_cast<size_t>(cnt) == sos_[i+1]);
        }
        Kokkos::deep_copy(slids_, slids);
      }

      // Optimistically set up for 1-vector case.
      init_mem(1);
    }

    template <typename T>
    static int isend (const MPI_Comm comm, const T* buf, int count, int dest, int tag,
                      MPI_Request* ireq) {
#ifdef HAVE_MPI
      const auto dt = Teuchos::Details::MpiTypeTraits<scalar_type>::getType();
      MPI_Request ureq;
      MPI_Request* req = ireq ? ireq : &ureq;
      int ret = MPI_Isend(const_cast<T*>(buf), count, dt, dest, tag, comm, req);
      if ( ! ireq) MPI_Request_free(req);
      return ret;
#else
      return 0;
#endif
    }

    template <typename T>
    static int irecv (const MPI_Comm comm, T* buf, int count, int src, int tag,
                      MPI_Request* ireq) {
#ifdef HAVE_MPI
      const auto dt = Teuchos::Details::MpiTypeTraits<scalar_type>::getType();
      MPI_Request ureq;
      MPI_Request* req = ireq ? ireq : &ureq;
      int ret = MPI_Irecv(buf, count, dt, src, tag, comm, req);
      if ( ! ireq) MPI_Request_free(req);
      return ret;
#else
      return 0;
#endif
    }

    static int waitany (int count, MPI_Request* reqs, int* index) {
#ifdef HAVE_MPI
      return MPI_Waitany(count, reqs, index, MPI_STATUS_IGNORE);
#else
      return 0;
#endif
    }

    static int waitall (int count, MPI_Request* reqs) {
#ifdef HAVE_MPI
      return MPI_Waitall(count, reqs, MPI_STATUS_IGNORE);
#else
      return 0;
#endif
    }

  public: // (for Cuda)

    class UnPack {
    public:
      enum Direction { pack, unpack };

    private:
      size_t is_, ie_, n_;
      LO bsz_, bos0_;
      Direction dir_;
      ConstUnmanaged<LOList> lids_;
      Unmanaged<Buffer> buf_;
      Unmanaged<MvView> X_;

      template <typename Tag> KOKKOS_FORCEINLINE_FUNCTION void compute (const LO& k) const {
        const auto col  = k / n_;
        const auto i    = k % n_;
        const auto bos  = bos0_ + bsz_*k;
        const auto Xos  = bsz_*lids_(is_ + i);
        // if-else should compile out.
        if (std::is_same<Tag, PackTag>::value)
          for (int j = 0; j < bsz_; ++j)
            buf_(bos + j) = X_(Xos + j, col);
        else
          for (int j = 0; j < bsz_; ++j)
            X_(Xos + j, col) = buf_(bos + j);
      }

      template <typename Tag> KOKKOS_FORCEINLINE_FUNCTION void computeone (const LO& k) const {
        const auto i = is_ + k;
        const auto bos = bsz_*i;
        const auto Xos = bsz_*lids_(i);
        if (std::is_same<Tag, PackTag>::value)
          for (int j = 0; j < bsz_; ++j)
            buf_(bos + j) = X_(Xos + j, 0);
        else
          for (int j = 0; j < bsz_; ++j)
            X_(Xos + j, 0) = buf_(bos + j);
      }

    public:
      struct PackTag {};
      struct UnpackTag {};
      struct PackOneTag {};
      struct UnpackOneTag {};

      UnPack (const LOList& lids, const Buffer& buf, const size_t& is, const size_t& ie,
              const LO& bsz, const Direction& dir, const MvView& X)
        : is_(is), ie_(ie), bsz_(bsz), dir_(dir), lids_(lids), buf_(buf), X_(X)
      {}

      KOKKOS_INLINE_FUNCTION void operator() (const PackTag  &, const LO& k) const
      { compute<PackTag  >(k); }
      KOKKOS_INLINE_FUNCTION void operator() (const UnpackTag&, const LO& k) const
      { compute<UnpackTag>(k); }
      KOKKOS_INLINE_FUNCTION void operator() (const PackOneTag  &, const LO& k) const
      { computeone<PackTag  >(k); }
      KOKKOS_INLINE_FUNCTION void operator() (const UnpackOneTag&, const LO& k) const
      { computeone<UnpackTag>(k); }

      void run () {
        namespace ko = Kokkos;
        using es = typename device_type::execution_space;
        if (X_.extent(1) == 1) {
          if (dir_ == pack)
            ko::parallel_for(ko::RangePolicy<PackOneTag  , es>(0, ie_ - is_), *this);
          else
            ko::parallel_for(ko::RangePolicy<UnpackOneTag, es>(0, ie_ - is_), *this);
        } else {
          n_ = ie_ - is_;
          bos0_ = bsz_*is_*X_.extent(1);
          if (dir_ == pack)
            ko::parallel_for(ko::RangePolicy<PackTag  , es>(0, n_*X_.extent(1)), *this);
          else
            ko::parallel_for(ko::RangePolicy<UnpackTag, es>(0, n_*X_.extent(1)), *this);
        }
      }
    };

  public:
    AsyncableImport (const Teuchos::RCP<const map_type>& src_map,
                     const Teuchos::RCP<const map_type>& tgt_map, const LO block_size) {
      init(src_map, tgt_map, block_size);
    }

    void async_setup (const mv_type& X) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Async_Setup::async_setup");
#endif
#ifdef HAVE_MPI
      const auto nvec = X.getNumVectors();
      const auto bsz_nvec = bsz_*nvec;
      init_mem(nvec);

      const auto mpi_comm =
        Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm_);
      TEUCHOS_ASSERT( ! mpi_comm.is_null());
      const MPI_Comm comm = *mpi_comm->getRawMpiComm();

      // Send async.
      for (int i = 0; i < spids_.size(); ++i) {
        UnPack(slids_, sbuf_, sos_[i], sos_[i+1], bsz_, UnPack::pack,
               X.template getLocalView<tpetra_device_type>()).run();
        isend(comm, sbuf_.data() + sos_[i]*bsz_nvec, (sos_[i+1] - sos_[i])*bsz_nvec, spids_[i], 42,
              &sreqs_[i]);
      }

      // Init receives async.
      for (int i = 0; i < rpids_.size(); ++i)
        irecv(comm, rbuf_.data() + ros_[i]*bsz_nvec, (ros_[i+1] - ros_[i])*bsz_nvec, rpids_[i], 42,
              &rreqs_[i]);

      // I find that issuing an Iprobe seems to nudge some MPIs into action,
      // which helps with overlapped comm/comp performance.
      for (int i = 0; i < rpids_.size(); ++i) {
        int flag;
        MPI_Status stat;
        MPI_Iprobe(rpids_[i], 42, comm, &flag, &stat);
      }
#endif
    }

    void cancel () {
#ifdef HAVE_MPI
      for (int i = 0; i < rpids_.size(); ++i)
        MPI_Cancel(&rreqs_[i]);
#endif
    }

    void sync_receive (const bool need_waitany = true) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Async_Setup::sync_receive");
#endif
#ifdef HAVE_MPI
      // Receive async.
      for (int ir = 0; ir < rpids_.size(); ++ir) {
        int i = ir;
        if (need_waitany)
          waitany(rpids_.size(), rreqs_.data(), &i);
        UnPack(rlids_, rbuf_, ros_[i], ros_[i+1], bsz_, UnPack::unpack, X_remote_).run();
      }
      // Wait on the sends to match all Isends with a cleanup operation.
      waitall(sreqs_.size(), sreqs_.data());
#endif
    }

    void sync_exchange (const mv_type& X) {
      async_setup(X);
      sync_receive();
    }

    const Teuchos::RCP<const Teuchos::Comm<int> >& get_comm () const { return comm_; }

    const MvView& get_mv_remote () const { return X_remote_; }
  };

  // Manage the distributed part of the computation of residual norms.
  class NormManager {
    bool collective_;
    int sweep_step_, bsz_, nvec_;
    std::vector<magnitude_type> wrk_;
    magnitude_type n0_;
#ifdef HAVE_MPI
    MPI_Request mpi_request_;
    MPI_Comm comm_;
#endif

  public:
    typedef std::shared_ptr<NormManager> Ptr;

    NormManager (const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
      sweep_step_ = 1;
      n0_ = 0;
      collective_ = comm->getSize() > 1;
      if (collective_) {
#ifdef HAVE_MPI
        const auto mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
        TEUCHOS_ASSERT( ! mpi_comm.is_null());
        comm_ = *mpi_comm->getRawMpiComm();
#endif
      }
    }

    // Resize the buffer to accommodate nvec vectors in the multivector, for a
    // matrix having block size block_size.
    void resize (const int& block_size, const int& nvec) {
      bsz_ = block_size;
      nvec_ = nvec;
      wrk_.resize((2*block_size + 1)*nvec);
    }

    // Check the norm every sweep_step sweeps.
    void set_check_frequency (const int& sweep_step) {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(sweep_step < 1, "sweep step must be >= 1");
      sweep_step_ = sweep_step;
    }

    // Get the buffer into which to store rank-local squared norms.
    magnitude_type* get_buffer () { return wrk_.data(); }

    // Call MPI_Iallreduce to find the global squared norms.
    void ireduce (const int& sweep, const bool force = false) {
      if ( ! force && sweep % sweep_step_) return;
      const int n = bsz_*nvec_;
      if (collective_) {
        std::copy(wrk_.begin(), wrk_.begin() + n, wrk_.begin() + n);
#ifdef HAVE_MPI
# if MPI_VERSION >= 3
        MPI_Iallreduce(wrk_.data() + n, wrk_.data(), n,
                       Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                       MPI_SUM, comm_, &mpi_request_);
# else
        MPI_Allreduce (wrk_.data() + n, wrk_.data(), n,
                       Teuchos::Details::MpiTypeTraits<magnitude_type>::getType(),
                       MPI_SUM, comm_);
# endif
#endif
      }
    }

    // Check if the norm-based termination criterion is met. tol2 is the
    // tolerance squared. Sweep is the sweep index. If not every iteration is
    // being checked, this function immediately returns false. If a check must
    // be done at this iteration, it waits for the reduction triggered by
    // ireduce to complete, then checks the global norm against the tolerance.
    bool check_done (const int& sweep, const magnitude_type& tol2, const bool force = false) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::NormManager::check_done");
#endif
      TEUCHOS_ASSERT(sweep >= 1);
      if ( ! force && (sweep - 1) % sweep_step_) return false;
      if (collective_) {
#ifdef HAVE_MPI
# if MPI_VERSION >= 3
        MPI_Wait(&mpi_request_, MPI_STATUS_IGNORE);
# else
        // Do nothing.
# endif
#endif
      }
      const auto n = bsz_*nvec_;
      if (sweep == 1) {
        magnitude_type* const n0 = wrk_.data() + 2*n;
        for (int v = 0; v < nvec_; ++v) {
          const magnitude_type* const dn0 = wrk_.data() + v*bsz_;
          magnitude_type mdn0 = 0;
          for (int i = 0; i < bsz_; ++i)
            mdn0 = std::max(mdn0, dn0[i]);
          n0[v] = mdn0;
        }
        return false;
      } else {
        const auto n0 = wrk_.data() + 2*n;
        bool done = true;
        for (int v = 0; v < nvec_; ++v) {
          const magnitude_type* const dnf = wrk_.data() + v*bsz_;
          magnitude_type mdnf = 0;
          for (int i = 0; i < bsz_; ++i)
            mdnf = std::max(mdnf, dnf[i]);
          if (mdnf > tol2*n0[v]) {
            done = false;
            break;
          }
        }
        return done;
      }
    }

    // After termination has occurred, finalize the norms for use in
    // get_norms{0,final}.
    void finalize () {
      for (int v = 0; v < nvec_; ++v) {
        const magnitude_type* const dnf = wrk_.data() + v*bsz_;
        magnitude_type mdnf = 0;
        for (int i = 0; i < bsz_; ++i)
          mdnf = std::max(mdnf, dnf[i]);
        // This overwrites the receive buffer, but that's OK; at the time of
        // this write, we no longer need the data in this slot.
        wrk_[v] = mdnf;
      }
      for (int i = 0; i < nvec_; ++i)
        wrk_[i] = std::sqrt(wrk_[i]);
      magnitude_type* const nf = wrk_.data() + 2*bsz_*nvec_;
      for (int v = 0; v < nvec_; ++v)
        nf[v] = std::sqrt(nf[v]);
    }

    // Report norms to the caller.
    const magnitude_type* get_norms0 () const { return wrk_.data() + 2*bsz_*nvec_; }
    const magnitude_type* get_norms_final () const { return wrk_.data(); }
  };

  // Map SIMD-packed multivectors to Tpetra::MultiVectors, accounting for
  // preconditioner partitions. At the same time (for efficiency), compute the
  // rank-local squared norms, if requested.
  class PermuteAndRepack {
  private:
    enum Enum { left, right };

    typename mv_type::dual_view_type::t_dev mv;
    typename Tridiags::PackedMultiVector pmv;
    Enum side;
    impl_scalar_type damping_factor;
    ConstUnmanaged<LOList> part2rowidx0, part2packrowidx0, lclrow, packptr;
    bool y_is_zero;
    magnitude_type* norm2sqr;

    void init (const mv_type& mv_, const typename Tridiags::PackedMultiVector& pmv_,
               const Enum& side_, const impl_scalar_type& damping_factor_, const LOList& part2rowidx0_,
               const LOList& part2packrowidx0_, const LOList& lclrow_, const LOList& packptr_,
               const bool y_is_zero_ = false, magnitude_type* norm2sqr_ = nullptr) {
      mv = mv_.template getLocalView<tpetra_device_type>();
      pmv = pmv_;
      value_count = pmv.extent(1)*pmv.extent(2);
      side = side_;
      damping_factor = damping_factor_;
      part2rowidx0 = part2rowidx0_;
      part2packrowidx0 = part2packrowidx0_;
      lclrow = lclrow_;
      packptr = packptr_;
      y_is_zero = y_is_zero_;
      norm2sqr = norm2sqr_;
      TEUCHOS_ASSERT(lclrow.extent(0)*pmv.extent(1) == mv.extent(0));
    }

    template <typename Tag>
    void run () {
      const Kokkos::RangePolicy<Tag, typename device_type::execution_space>
        rp(0, packptr.size() - 1);
      if (norm2sqr) {
        for (int i = 0; i < value_count; ++i)
          norm2sqr[i] = Kokkos::ArithTraits<magnitude_type>::zero();
        Kokkos::parallel_reduce(rp, *this, norm2sqr);
      } else
        Kokkos::parallel_for(rp, *this);
    }

  public:
    struct LeftTag {};
    struct LeftDFTag {};
    struct RightTag {};

    // For || reduce.
    typedef magnitude_type value_type[];
    int value_count;

    // Packed multivector <- Tpetra::MultiVector.
    PermuteAndRepack (const mv_type& mv_, typename Tridiags::PackedMultiVector& pmv_,
                      const LOList& part2rowidx0_, const LOList& part2packrowidx0_,
                      const LOList& lclrow_, const LOList& packptr_)
    {
      const auto one = Kokkos::ArithTraits<magnitude_type>::one();
      init(mv_, pmv_, right, one, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_);
    }

    // Tpetra::MultiVector <- packed multivector.
    PermuteAndRepack (const typename Tridiags::PackedMultiVector& pmv_, mv_type& mv_,
                      const impl_scalar_type& damping_factor_, const LOList& part2rowidx0_,
                      const LOList& part2packrowidx0_, const LOList& lclrow_, const LOList& packptr_,
                      const bool y_is_zero_, magnitude_type* norm2sqr_ = nullptr)
    {
      init(mv_, pmv_, left, damping_factor_, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_,
           y_is_zero_, norm2sqr_);
    }

    KOKKOS_INLINE_FUNCTION void init (magnitude_type* dst) const {
      for (int i = 0; i < value_count; ++i) dst[i] = 0;
    }

    KOKKOS_INLINE_FUNCTION
    void join (volatile magnitude_type* dst, const volatile magnitude_type* src) const {
      for (int i = 0; i < value_count; ++i) dst[i] += src[i];
    }

    template <typename Tag> KOKKOS_INLINE_FUNCTION
    void operator() (const Tag&, const LO& packidx, magnitude_type* const nr = nullptr) const {
      constexpr auto vl = vectorization_traits::vector_length;
      const LO bs = pmv.extent(1), nrhs = pmv.extent(2);
      LO partidx = packptr(packidx);
      LO npacks = packptr(packidx+1) - partidx;
      const LO pri0 = part2packrowidx0(partidx);
      LO ri0[vl] = {0}, nrows[vl] = {0};
      for (LO vi = 0; vi < npacks; ++vi, ++partidx) {
        ri0[vi] = part2rowidx0(partidx);
        nrows[vi] = part2rowidx0(partidx+1) - ri0[vi];
      }
      for (LO j = 0; j < nrows[0]; ++j) {
        for (LO vi = 1; vi < npacks; ++vi)
          if (j == nrows[vi]) {
            npacks = vi;
            break;
          }
        const LO pri = pri0 + j;
        for (LO col = 0; col < nrhs; ++col) {
          // In both cases, the written memory is traversed as contiguously as
          // possible. The if statements should compile away since for a given
          // Tag, the condition is always true or false at compile time.
          if (std::is_same<Tag, RightTag>::value) {
            for (LO i = 0; i < bs; ++i)
              for (LO vi = 0; vi < npacks; ++vi)
                Details::Batched::fastidx(pmv, Kokkos::make_pair(pri, vi), i, col) =
                  mv(bs*lclrow(ri0[vi] + j) + i, col);
          } else {
            if (nr) {
              for (LO vi = 0; vi < npacks; ++vi) {
                const LO lr0 = bs*lclrow(ri0[vi] + j);
                const auto pair = Kokkos::make_pair(pri, vi);
                for (LO i = 0; i < bs; ++i) {
                  auto& y = mv(lr0 + i, col);
                  const auto& yc = Details::Batched::fastidx(pmv, pair, i, col);
                  auto d = yc;
                  if ( ! y_is_zero) d -= y;
                  nr[bs*col + i] += BlockTriDiContainerDetails::abs2(d);
                  if (std::is_same<Tag, LeftTag>::value)
                    y = yc;
                  else if (std::is_same<Tag, LeftDFTag>::value) {
                    if (y_is_zero)
                      y = damping_factor * yc;
                    else
                      y += damping_factor * d;
                  }
                }
              }
            } else {
              for (LO vi = 0; vi < npacks; ++vi) {
                const LO lr0 = bs*lclrow(ri0[vi] + j);
                const auto pair = Kokkos::make_pair(pri, vi);
                if (std::is_same<Tag, LeftTag>::value) {
                  for (LO i = 0; i < bs; ++i)
                    mv(lr0 + i, col) = Details::Batched::fastidx(pmv, pair, i, col);
                } else {
                  for (LO i = 0; i < bs; ++i) {
                    auto& y = mv(lr0 + i, col);
                    const auto& yc = Details::Batched::fastidx(pmv, pair, i, col);
                    if (y_is_zero)
                      y = damping_factor * yc;
                    else
                      y += damping_factor * (yc - y);
                  }
                }
              }
            }
          }
        }
      }
    }

    void run () {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::PermuteAndRepack");
#endif
      if (side == right)
        run<RightTag>();
      else if (damping_factor == Kokkos::ArithTraits<impl_scalar_type>::one())
        run<LeftTag>();
      else
        run<LeftDFTag>();
    }
  };

  // Repack from the BlockCrsMatrix to the tridiagonal storage, and factorize
  // the data in the latter.
  class ExtractAndFactorizeTridiags {
    Tridiags btdm;
    ConstUnmanaged<LOList> partptr, lclrow, packptr;
    ConstUnmanaged<CrsGraphRowMap> A_rowptr;
    ConstUnmanaged<BlockCrsValues> A_values;
    ConstUnmanaged<SizeList> pack_td_ptr;
    Unmanaged<typename Tridiags::Values> values;
    magnitude_type diag_safety;
    bool add_to_diag;

    template <LO vl> KOKKOS_FORCEINLINE_FUNCTION void
    pack (const impl_scalar_type* const block[vl], const LO& bs, const LO& pi, const LO& npacks,
          typename std::enable_if<(vl > 1)>::type* = 0) const {
      for (LO bi = 0; bi < bs; ++bi)
        for (LO bj = 0; bj < bs; ++bj) {
          const LO idx = bi*bs + bj;
          auto& v = values(pi, bi, bj);
          for (LO vi = 0; vi < npacks; ++vi)
            v[vi] = block[vi][idx];
        }
    }

    template <LO vl> KOKKOS_FORCEINLINE_FUNCTION void
    pack (const impl_scalar_type* const block[vl], const LO& bs, const LO& pi, const LO& npacks,
          typename std::enable_if<vl == 1>::type* = 0) const {
      for (LO bi = 0; bi < bs; ++bi)
        for (LO bj = 0; bj < bs; ++bj)
          for (LO vi = 0; vi < npacks; ++vi)
            values(pi, bi, bj) = block[0][bi*bs + bj];
    }

    KOKKOS_INLINE_FUNCTION void extract (const LO& packidx) const {
      constexpr auto vl = vectorization_traits::vector_length;
      const LO bs = values.extent(1), bs2 = bs*bs;
      LO partidx = packptr(packidx);
      // Recall problems are packed from longest to shortest. We exploit this
      // monotonicity property here in how we deal with different-size
      // tridiags within a pack.
      LO npacks = packptr(packidx+1) - partidx;
      const Size kps = btdm.pack_td_ptr(partidx);
      LO kfs[vl] = {0}, ri0[vl] = {0}, nrows[vl] = {0};
      for (LO vi = 0; vi < npacks; ++vi, ++partidx) {
        kfs[vi] = btdm.flat_td_ptr(partidx);
        ri0[vi] = partptr(partidx);
        nrows[vi] = partptr(partidx+1) - ri0[vi];
      }
      for (LO tr = 0, j = 0; tr < nrows[0]; ++tr) {
        for (int e = 0; e < 3; ++e) {
          const impl_scalar_type* block[vl] = {0};
          for (LO vi = 0; vi < npacks; ++vi) {
            const Size Aj = A_rowptr(lclrow(ri0[vi] + tr)) + btdm.A_colindsub(kfs[vi] + j);
            block[vi] = &A_values(Aj*bs2);
          }
          const Size pi = kps + j;
          ++j;
          pack<vl>(block, bs, pi, npacks);
          if (nrows[0] == 1) break;
          if (e == 1 && (tr == 0 || tr+1 == nrows[0])) break;
          for (LO vi = 1; vi < npacks; ++vi)
            if ((e == 0 && nrows[vi] == 1) || (e == 1 && tr+1 == nrows[vi])) {
              npacks = vi;
              break;
            }
        }
      }
    }

    template <typename AViewType>
    KOKKOS_INLINE_FUNCTION void LU (const AViewType& A) const {
#if defined(__KOKKOSBATCHED_PROMOTION__)
      {
        namespace kbe = KokkosBatched::Experimental;
        if (add_to_diag) {
          kbe::SerialAddRadial::invoke(diag_safety, A);
        }
        kbe::SerialLU<kbe::Algo::LU::Blocked>::invoke(A);
      }
#else // __KOKKOSBATCHED_PROMOTION__
      if (add_to_diag)
        KokkosBatched::Details::Todo::Serial_LU_Internal_Unblocked_invoke(
          A.dimension_0(), A.dimension_1(), A.data(), A.stride_0(), A.stride_1(), diag_safety);
      else {
        namespace kbe = KokkosBatched::Experimental;
        kbe::SerialLU<kbe::Algo::LU::Blocked>::invoke(A);
      }
#endif // __KOKKOSBATCHED_PROMOTION__
    }

    KOKKOS_INLINE_FUNCTION void factorize (const LO& packidx) const {
      using Kokkos::subview;
      using Kokkos::ALL;
      namespace kbe = KokkosBatched::Experimental;
      const auto one = Kokkos::ArithTraits<magnitude_type>::one();
      Size i0 = pack_td_ptr(packptr(packidx));
      const LO nrows = Tridiags::ind2row(pack_td_ptr(packptr(packidx+1)) - i0 - 1) + 1;
      auto A = subview(values, i0, ALL(), ALL());
      LU(A);
      for (LO i = 1; i < nrows; ++i) {
        /* Factorization of a 2x2 set of blocks from a block tridiag
           matrix. Given
             E = [(al,au) b; c d],
           form
             L = [al 0; L21 L22]
             U = [au U12; 0 U22],
           yielding the identities
           b = al U12
             c = L21 au => au' L21' = c'
             d = L21 U12 + L22 U22.
           Algorithm:
           Factor (EL,EU) <- E in place:
             1a. Solve al b = b.
             1b. Solve au' c' = c'.
             2.  d := d - c*b.
             3.  Factor d.
        */
        const auto B = subview(values, i0+1, ALL(), ALL());
        kbe::SerialTrsm<kbe::Side::Left, kbe::Uplo::Lower, kbe::Trans::NoTranspose, kbe::Diag::Unit,
                        kbe::Algo::Trsm::Blocked>
          ::invoke(one, A, B);
        const auto C = subview(values, i0+2, ALL(), ALL());
        kbe::SerialTrsm<kbe::Side::Right, kbe::Uplo::Upper, kbe::Trans::NoTranspose, kbe::Diag::NonUnit,
                        kbe::Algo::Trsm::Blocked>
          ::invoke(one, A, C);
        A = subview(values, i0+3, ALL(), ALL());
        kbe::SerialGemm<kbe::Trans::NoTranspose, kbe::Trans::NoTranspose, kbe::Algo::Gemm::Blocked>
          ::invoke(-one, C, B, one, A);
        LU(A);
        i0 += 3;
      }
    }

  public:
    ExtractAndFactorizeTridiags (const Tridiags& btdm_, const LOList& partptr_, const LOList& lclrow_,
                                 const LOList& packptr_, const CrsGraphRowMap& A_rowptr_,
                                 const BlockCrsValues& A_values_, const magnitude_type& diag_safety_)
      : btdm(btdm_), partptr(partptr_), lclrow(lclrow_), packptr(packptr_),
        A_rowptr(A_rowptr_), A_values(A_values_),
        pack_td_ptr(btdm.pack_td_ptr), values(btdm.values),
        diag_safety(diag_safety_), add_to_diag(diag_safety != Kokkos::ArithTraits<magnitude_type>::zero())
    {}

    KOKKOS_INLINE_FUNCTION void operator() (const LO& packidx) const {
      extract(packidx);
      factorize(packidx);
    }

    void run () {
      Kokkos::RangePolicy<typename device_type::execution_space> rp(0, packptr.size() - 1);
      Kokkos::parallel_for(rp, *this);
    }
  };

  // Solve D X = Y, where A + D + R. This general implementation uses subviews
  // and works on any platform.
  template <typename Layout, typename PartialSpecialization>
  class Solve {
  public: // for Cuda
    // Whether RHS is a vector or a multivector.
    struct VecTag {};
    struct MatTag {};

  private:
    ConstUnmanaged<SizeList> pack_td_ptr;
    ConstUnmanaged<typename Tridiags::Values> values;
    ConstUnmanaged<LOList> packptr, part2packrowidx0;
    Unmanaged<typename Tridiags::PackedMultiVector> X;

    template <typename Tag> void run () {
      Kokkos::RangePolicy<Tag, typename device_type::execution_space> rp(0, packptr.size() - 1);
      Kokkos::parallel_for(rp, *this);
    }

    // X := a T \ X
    template <typename Tag, typename UploType, typename DiagType,
              typename Scalar, typename MatA, typename MatX>
    KOKKOS_FORCEINLINE_FUNCTION static void
    trsmv (const Scalar& a, const MatA& A, MatX& X,
           typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialTrsv<UploType, kbe::Trans::NoTranspose, DiagType, kbe::Algo::Trsv::Unblocked>
        ::invoke(a, A, X);
    }

    template <typename Tag, typename UploType, typename DiagType,
              typename Scalar, typename MatA, typename MatX>
    KOKKOS_FORCEINLINE_FUNCTION static void
    trsmv (const Scalar& a, const MatA& A, MatX& X,
           typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialTrsm<kbe::Side::Left, UploType, kbe::Trans::NoTranspose, DiagType, kbe::Algo::Trsm::Blocked>
        ::invoke(a, A, X);
    }

    // C := b C + a A B
    template <typename Tag, typename Scalar, typename MatA, typename MatB, typename MatC>
    KOKKOS_FORCEINLINE_FUNCTION static void
    gemmv (const Scalar& a, const MatA& A, const MatB& B, const Scalar& b, MatC& C,
           typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialGemv<kbe::Trans::NoTranspose, kbe::Algo::Gemv::Unblocked>::invoke(a, A, B, b, C);
    }

    template <typename Tag, typename Scalar, typename MatA, typename MatB, typename MatC>
    KOKKOS_FORCEINLINE_FUNCTION static void
    gemmv (const Scalar& a, const MatA& A, const MatB& B, const Scalar& b, MatC& C,
           typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialGemm<kbe::Trans::NoTranspose, kbe::Trans::NoTranspose, kbe::Algo::Gemm::Blocked>
        ::invoke(a, A, B, b, C);
    }

  public:
    Solve (const Tridiags& btdm, const LOList& packptr_, const LOList& part2packrowidx0_,
           const typename Tridiags::PackedMultiVector& X_)
      : pack_td_ptr(btdm.pack_td_ptr), values(btdm.values),
        packptr(packptr_), part2packrowidx0(part2packrowidx0_),
        X(X_)
    {}

    template <typename Tag>
    KOKKOS_INLINE_FUNCTION void operator() (const Tag&, const LO& packidx) const {
      using Kokkos::subview;
      using Kokkos::ALL;
      namespace kbe = KokkosBatched::Experimental;

      const auto one = Kokkos::ArithTraits<magnitude_type>::one();

      const auto partidx = packptr(packidx);
      Size i0 = pack_td_ptr(partidx);
      LO r0 = part2packrowidx0(partidx);
      const LO nrows = part2packrowidx0(packptr(packidx+1)) - r0;

      // Solve L x = x.
      auto A = subview(values, i0, ALL(), ALL());
      auto X1 = subview(X, r0, ALL(), ALL());
      trsmv<Tag, kbe::Uplo::Lower, kbe::Diag::Unit>(one, A, X1);
      for (LO i = 1; i < nrows; ++i) {
        const auto B = subview(values, i0+2, ALL(), ALL());
        r0++;
        const auto X2 = subview(X, r0, ALL(), ALL());
        gemmv<Tag>(-one, B, X1, one, X2);
        i0 += 3;
        A = subview(values, i0, ALL(), ALL());
        trsmv<Tag, kbe::Uplo::Lower, kbe::Diag::Unit>(one, A, X2);
        X1 = X2;
      }

      // Solve U x = x.
      trsmv<Tag, kbe::Uplo::Upper, kbe::Diag::NonUnit>(one, A, X1);
      for (LO i = nrows; i > 1; --i) {
        i0 -= 3;
        const auto B = subview(values, i0+1, ALL(), ALL());
        r0--;
        const auto X2 = subview(X, r0, ALL(), ALL());
        gemmv<Tag>(-one, B, X1, one, X2);
        A = subview(values, i0, ALL(), ALL());
        trsmv<Tag, kbe::Uplo::Upper, kbe::Diag::NonUnit>(one, A, X2);
        X1 = X2;
      }
    }

    void run () {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Solve");
#endif
      if (X.extent(2) == 1)
        run<VecTag>();
      else
        run<MatTag>();
    }
  };

  // This specialization is for CPU and KNL. On KNL, it speeds up the solve by
  // ~20% for block size 5, more for smaller, less for bigger. The importance of
  // the apply phase of a preconditioner justifies the extra complexity of this
  // specialization.
  template <typename PartialSpecialization>
  class Solve<Kokkos::LayoutRight, PartialSpecialization> {
  public: // for Cuda
    // Whether RHS is a vector or a multivector.
    struct VecTag {};
    struct MatTag {};

  private:
    ConstUnmanaged<SizeList> pack_td_ptr;
    ConstUnmanaged<typename Tridiags::Values> values;
    ConstUnmanaged<LOList> packptr, part2packrowidx0;
    Unmanaged<typename Tridiags::PackedMultiVector> X;
    const int bs2, xsz;

    template <typename Tag> void run () {
      Kokkos::RangePolicy<Tag, typename device_type::execution_space> rp(0, packptr.size() - 1);
      Kokkos::parallel_for(rp, *this);
    }

    // For speed in this important, O(nnz) operation, we don't use subviews, but
    // instead go directly to the raw-array interface.

    // X := a T \ X

    template <typename Tag>
    KOKKOS_FORCEINLINE_FUNCTION void
    trsmvlo (const magnitude_type& a, const LO& i0, const LO& r0,
             typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialTrsvInternalLower<kbe::Algo::Trsv::Unblocked>::invoke(
        true,
        values.dimension_1(),
        a,
        values.data() + i0*bs2, values.stride_1(), values.stride_2(),
        X.data() + r0*xsz, X.stride_1());
    }

    template <typename Tag>
    KOKKOS_FORCEINLINE_FUNCTION void
    trsmvlo (const magnitude_type& a, const LO& i0, const LO& r0,
             typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialTrsmInternalLeftLower<kbe::Algo::Trsm::Blocked>::invoke(
        true,
        values.dimension_1(), X.dimension_2(),
        a,
        values.data() + i0*bs2, values.stride_1(), values.stride_2(),
        X.data() + r0*xsz, X.stride_1(), X.stride_2());
    }

    template <typename Tag>
    KOKKOS_FORCEINLINE_FUNCTION void
    trsmvup (const magnitude_type& a, const LO& i0, const LO& r0,
             typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialTrsvInternalUpper<kbe::Algo::Trsv::Unblocked>::invoke(
        false,
        values.dimension_1(),
        a,
        values.data() + i0*bs2, values.stride_1(), values.stride_2(),
        X.data() + r0*xsz, X.stride_1());
    }

    template <typename Tag>
    KOKKOS_FORCEINLINE_FUNCTION void
    trsmvup (const magnitude_type& a, const LO& i0, const LO& r0,
             typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialTrsmInternalLeftUpper<kbe::Algo::Trsm::Blocked>::invoke(
        false,
        values.dimension_1(), X.dimension_2(),
        a,
        values.data() + i0*bs2, values.stride_1(), values.stride_2(),
        X.data() + r0*xsz, X.stride_1(), X.stride_2());
    }

    // C := b C + a A B

    template <typename Tag>
    KOKKOS_FORCEINLINE_FUNCTION void
    gemmv (const magnitude_type& a, const magnitude_type& b, const LO& a0, const LO& b0, const LO& c0,
           typename std::enable_if<std::is_same<Tag, VecTag>::value>::type* = 0) const {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialGemvInternal<kbe::Algo::Gemv::Unblocked>::invoke(
        values.dimension_1(), values.dimension_2(),
        a,
        values.data() + a0*bs2, values.stride_1(), values.stride_2(),
        X.data() + b0*xsz, X.stride_1(),
        b,
        X.data() + c0*xsz, X.stride_1());
    }

    template <typename Tag>
    KOKKOS_FORCEINLINE_FUNCTION void
    gemmv (const magnitude_type& a, const magnitude_type& b, const LO& a0, const LO& b0, const LO& c0,
           typename std::enable_if<std::is_same<Tag, MatTag>::value>::type* = 0) const {
      namespace kbe = KokkosBatched::Experimental;
      kbe::SerialGemmInternal<kbe::Algo::Gemm::Blocked>::invoke(
        X.dimension_1(), X.dimension_2(), values.dimension_2(),
        a,
        values.data() + a0*bs2, values.stride_1(), values.stride_2(),
        X.data() + b0*xsz, X.stride_1(), X.stride_2(),
        b,
        X.data() + c0*xsz, X.stride_1(), X.stride_2());
    }

  public:
    Solve (const Tridiags& btdm, const LOList& packptr_, const LOList& part2packrowidx0_,
           const typename Tridiags::PackedMultiVector& X_)
      : pack_td_ptr(btdm.pack_td_ptr), values(btdm.values),
        packptr(packptr_), part2packrowidx0(part2packrowidx0_),
        X(X_), bs2(values.dimension_1()*values.dimension_1()),
        xsz(values.dimension_1()*X.dimension_2())
    {}

    template <typename Tag>
    KOKKOS_INLINE_FUNCTION void operator() (const Tag&, const LO& packidx) const {
      using Kokkos::subview;
      using Kokkos::ALL;
      namespace kbe = KokkosBatched::Experimental;

      const auto one = Kokkos::ArithTraits<magnitude_type>::one();

      const auto partidx = packptr(packidx);
      Size i0 = pack_td_ptr(partidx);
      LO r0 = part2packrowidx0(partidx);
      const LO nrows = part2packrowidx0(packptr(packidx+1)) - r0;

      // Solve L x = x.
      trsmvlo<Tag>(one, i0, r0);
      for (LO i = 1; i < nrows; ++i) {
        r0++;
        gemmv<Tag>(-one, one, i0+2, r0-1, r0);
        i0 += 3;
        trsmvlo<Tag>(one, i0, r0);
      }

      // Solve U x = x.
      trsmvup<Tag>(one, i0, r0);
      for (LO i = nrows; i > 1; --i) {
        i0 -= 3;
        r0--;
        gemmv<Tag>(-one, one, i0+1, r0+1, r0);
        trsmvup<Tag>(one, i0, r0);
      }
    }

    void run () {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
      TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::Solve");
#endif
      if (X.extent(2) == 1)
        run<VecTag>();
      else
        run<MatTag>();
    }
  };

private:

  template <typename T> static KOKKOS_INLINE_FUNCTION constexpr T square (const T& x) { return x*x; }

  // Consolidate error output prefix.
  std::string get_msg_prefix () const {
    const auto& comm = A_bcrs_->getRowMap()->getComm();
    const auto rank = comm->getRank();
    const auto nranks = comm->getSize();
    std::stringstream ss;
    ss << "Rank " << rank << " of " << nranks << ": ";
    return ss.str();
  }

  template <typename View>
  typename View::HostMirror tohost (const View& v) {
    const auto hv = Kokkos::create_mirror_view(v);
    Kokkos::deep_copy(hv, v);
    return hv;
  }

  // Top-level Impl object initialization.
  void init (const Teuchos::RCP<const row_matrix_type>& matrix,
             const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
             const Teuchos::RCP<const import_type>& importer, int overlapLevel,
             bool overlapCommAndComp, bool useSeqMethod) {
    seq_method_ = useSeqMethod;
    overlap_comm_ = overlapCommAndComp;
    validate_ = true;

    TEUCHOS_TEST_FOR_EXCEPT_MSG(
      overlapLevel != 0,
      "BlockTriDiContainer does not curently support OverlapLevel != 0; user provided "
      << overlapLevel);

    A_bcrs_ = Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(matrix);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
      A_bcrs_.is_null(), "BlockTriDiContainer currently supports Tpetra::BlockCrsMatrix only.");

    init_importer(importer);
    init_parts(partitions);
    init_btdm();
  }

  // If possible, create a Tpetra::Map that is for remotes only, instead of for
  // remotes and owned IDs. This map is used in the two versions of
  // compute_b_minus_Rx that distinguish between R_owned and R_remote.
  static Teuchos::RCP<const map_type>
  create_parsimonious_map (const map_type& dm, const map_type& cm, bool& need_owned_permutation) {
    need_owned_permutation = false;
    std::vector<GO> gids;
    bool separate_remotes = true, found_first = false;
    for (size_t i = 0; i < cm.getNodeNumElements(); ++i) {
      const GO gid = cm.getGlobalElement(i);
      if ( ! dm.isNodeGlobalElement(gid)) {
        found_first = true;
        gids.push_back(gid);
      } else if (found_first) {
        separate_remotes = false;
        break;
      }
      if ( ! need_owned_permutation && dm.getLocalElement(gid) != static_cast<LO>(i)) {
        // The owned part of the domain and column maps are different
        // orderings. We *could* do a super efficient impl of this case in the
        // num_sweeps > 1 case by adding complexity to PermuteAndRepack. But,
        // really, if a caller cares about speed, they wouldn't make different
        // local permutations like this. So we punt on the best impl and go for
        // a pretty good one: the permutation is done in place in
        // compute_b_minus_Rx for the pure-owned part of the MVP. The only cost
        // is the presumably worse memory access pattern of the input vector.
        need_owned_permutation = true;
      }
    }
    if ( ! separate_remotes) {
      // The overlapped comm/comp method doesn't handle the case of remote LIDs
      // mixed with owned LIDs. I strongly suspect this is not supported in
      // Trilinos in general, but I'm not positive. In any case, return null,
      // which triggers switching to seq_method_ = true.
      return Teuchos::null;
    }
    const auto invalid = Teuchos::OrdinalTraits<GO>::invalid();
    return Teuchos::rcp(new map_type(invalid, gids.data(), gids.size(), 0, dm.getComm()));
  }

  void init_importer (const Teuchos::RCP<const import_type>& importer) {
    if (A_bcrs_.is_null()) {
      importer_ = importer;
    } else {
      typedef Tpetra::BlockMultiVector<scalar_type, LO, GO, node_type> BMV;
      const auto blocksz = A_bcrs_->getBlockSize();
      const auto g = A_bcrs_->getCrsGraph();
      bool need_owned_permutation;
      const auto parsimonious_col_map = seq_method_ ? Teuchos::null :
        create_parsimonious_map(*g.getDomainMap(), *g.getColMap(), need_owned_permutation);
      if (parsimonious_col_map.is_null()) {
        seq_method_ = true;
        overlap_comm_ = false;
      } else {
        // Make the importer only if needed.
        if (parsimonious_col_map->getGlobalNumElements() > 0) {
          try {
            async_import_ = std::make_shared<AsyncableImport>(
              g.getDomainMap(), parsimonious_col_map, blocksz);
          } catch (...) {
            // Only if import is needed and it fails to build do we need to
            // switch to the slow method.
            seq_method_ = true;
          }
        }
        if (need_owned_permutation) {
          const auto& dm = *g.getDomainMap();
          const auto& cm = *g.getColMap();
          dm2cm_ = LOList("dm2cm_", dm.getNodeNumElements());
          const auto dm2cm = Kokkos::create_mirror_view(dm2cm_);
          for (size_t i = 0; i < dm.getNodeNumElements(); ++i)
            dm2cm(i) = dm.getLocalElement(cm.getGlobalElement(i));
          Kokkos::deep_copy(dm2cm_, dm2cm);
        }
      }
      if (seq_method_) {
        const auto cpm = Teuchos::rcp(
          new typename BMV::map_type(
            BMV::makePointMap(seq_method_ ? *g.getColMap() : *parsimonious_col_map, blocksz)));
        const auto dpm = Teuchos::rcp(
          new typename BMV::map_type(BMV::makePointMap(*g.getDomainMap(), blocksz)));
        importer_ = Teuchos::rcp(new import_type(dpm, cpm));
      }
    }
  }

  // Reorder parts to maximize SIMD packing efficiency.
  static void reorder_parts (const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& parts,
                             std::vector<LO>& p) {
    const LO nparts = parts.size();
    p.resize(nparts);

    if (vectorization_traits::vector_length == 1) {
      for (LO i = 0; i < nparts; ++i)
        p[i] = i;
    } else {
      typedef std::pair<LO,LO> SzIdx;
      std::vector<SzIdx> partsz(nparts);
      for (LO i = 0; i < nparts; ++i)
        partsz[i] = SzIdx(parts[i].size(), i);
      std::sort(partsz.begin(), partsz.end(),
                [=] (const SzIdx& x, const SzIdx& y) { return x.first > y.first; });
      for (LO i = 0; i < nparts; ++i)
        p[i] = partsz[i].second;
    }
  }

  // Initialize the parts using the Container partitions array.
  void init_parts (const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions) {
    const bool jacobi = partitions.size() == 0;
    const LO A_n_lclrows = A_bcrs_->getNodeNumRows();
    const LO nparts = jacobi ? A_n_lclrows : partitions.size();

    LO nrows = 0;
    if (jacobi)
      nrows = nparts;
    else
      for (LO i = 0; i < nparts; ++i)
        nrows += partitions[i].size();
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
      nrows != A_n_lclrows, get_msg_prefix() << "The #rows implied by the local partition is not "
      << "the same as getNodeNumRows: " << nrows << " vs " << A_n_lclrows);

    std::vector<LO> pp;
    if ( ! jacobi)
      reorder_parts(partitions, pp);

    partptr_ = LOList("partptr_", nparts + 1);
    lclrow_ = LOList("lclrow_", nrows);
    part2rowidx0_ = LOList("part2rowidx0_", nparts + 1);
    part2packrowidx0_ = LOList("part2packrowidx0_", nparts + 1);
    rowidx2part_ = LOList("rowidx2part_", nrows);
    const auto partptr = Kokkos::create_mirror_view(partptr_);
    const auto lclrow = Kokkos::create_mirror_view(lclrow_);
    const auto part2rowidx0 = Kokkos::create_mirror_view(part2rowidx0_);
    const auto part2packrowidx0 = Kokkos::create_mirror_view(part2packrowidx0_);
    const auto rowidx2part = Kokkos::create_mirror_view(rowidx2part_);

    // Determine parts.
    row_contiguous_ = true;
    partptr(0) = 0;
    part2rowidx0(0) = 0;
    part2packrowidx0(0) = 0;
    LO pack_nrows = 0;
    for (LO ip = 0; ip < nparts; ++ip) {
      const auto* part = jacobi ? nullptr : &partitions[pp[ip]];
      const LO ipnrows = jacobi ? 1 : part->size();
      TEUCHOS_ASSERT(ip == 0 || vectorization_traits::vector_length == 1 ||
                     (jacobi || ipnrows <= static_cast<LO>(partitions[pp[ip-1]].size())));
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
        ipnrows == 0, get_msg_prefix() << "partition " << pp[ip]
        << " is empty, which is not allowed.");
      //assume No overlap.
      part2rowidx0(ip+1) = part2rowidx0(ip) + ipnrows;
      // Since parts are ordered in nonincreasing size, the size of the first
      // part in a pack is the size for all parts in the pack.
      if (ip % vectorization_traits::vector_length == 0)
        pack_nrows = ipnrows;
      part2packrowidx0(ip+1) = part2packrowidx0(ip) +
        ((ip+1) % vectorization_traits::vector_length == 0 || ip+1 == nparts ? pack_nrows : 0);
      const LO os = partptr(ip);
      for (LO i = 0; i < ipnrows; ++i) {
        const auto lcl_row = jacobi ? ip : (*part)[i];
        TEUCHOS_TEST_FOR_EXCEPT_MSG(
          lcl_row < 0 || lcl_row >= A_n_lclrows, get_msg_prefix() << "partitions[" << pp[ip] << "]["
          << i << "] = " << lcl_row << " but input matrix implies limits of [0, " << A_n_lclrows-1
          << "].");
        lclrow(os+i) = lcl_row;
        rowidx2part(os+i) = ip;
        if (row_contiguous_ && os+i > 0 && lclrow((os+i)-1) + 1 != lcl_row)
          row_contiguous_ = false;
      }
      partptr(ip+1) = os + ipnrows;
    }
    TEUCHOS_ASSERT(partptr(nparts) == nrows);
    if (lclrow(0) != 0) row_contiguous_ = false;

    Kokkos::deep_copy(partptr_, partptr);
    Kokkos::deep_copy(lclrow_, lclrow);
    //assume No overlap. Thus:
    part2rowidx0_ = partptr_;
    if (vectorization_traits::vector_length > 1)
      Kokkos::deep_copy(part2packrowidx0_, part2packrowidx0);
    else
      part2packrowidx0_ = part2rowidx0_;
    part2packrowidx0_back = part2packrowidx0(part2packrowidx0.size() - 1);
    Kokkos::deep_copy(rowidx2part_, rowidx2part);

    { // Fill packptr.
      LO npacks = 0;
      for (LO ip = 1; ip <= nparts; ++ip)
        if (part2packrowidx0(ip) != part2packrowidx0(ip-1))
          ++npacks;
      packptr_ = LOList("packptr_", npacks + 1);
      const auto packptr = Kokkos::create_mirror_view(packptr_);
      packptr(0) = 0;
      for (LO ip = 1, k = 1; ip <= nparts; ++ip)
        if (part2packrowidx0(ip) != part2packrowidx0(ip-1))
          packptr(k++) = ip;
      Kokkos::deep_copy(packptr_, packptr);
    }
  }

  // Initialize the block tridiagonal matrices data structure.
  void init_btdm () {
    const auto partptr = tohost(partptr_);
    const LO ntridiags = partptr.size() - 1;

    { // Construct the flat index pointers into the tridiag values array.
      btdm_.flat_td_ptr = SizeList("btdm.flat_td_ptr", ntridiags + 1);
      const auto flat_td_ptr = Kokkos::create_mirror_view(btdm_.flat_td_ptr);
      flat_td_ptr(0) = 0;
      for (LO ti = 1; ti <= ntridiags; ++ti) {
        const LO nrows = partptr(ti) - partptr(ti-1);
        flat_td_ptr(ti) = flat_td_ptr(ti-1) + (3*nrows - 2);
      }
      const Size nnz = 3*partptr(partptr.size() - 1) - 2*ntridiags;
      TEUCHOS_ASSERT(flat_td_ptr(ntridiags) == nnz);
      Kokkos::deep_copy(btdm_.flat_td_ptr, flat_td_ptr);
    }

    // And the packed index pointers.
    constexpr auto vl = vectorization_traits::vector_length;
    if (vl == 1) {
      btdm_.pack_td_ptr = btdm_.flat_td_ptr;
    } else {
      const auto packptr = tohost(packptr_);
      const LO npacks = packptr.size() - 1;
      btdm_.pack_td_ptr = SizeList("btdm.pack_td_ptr", ntridiags + 1);
      const auto pack_td_ptr = Kokkos::create_mirror_view(btdm_.pack_td_ptr);
      Size nblks = 0;
      for (LO pki = 0; pki < npacks; ++pki) {
        const LO parti = packptr(pki);
        for (LO pti = parti; pti < packptr(pki+1); ++pti)
          pack_td_ptr(pti) = nblks;
        nblks += Tridiags::nblks(partptr(parti+1) - partptr(parti));
      }
      pack_td_ptr(ntridiags) = nblks;
      Kokkos::deep_copy(btdm_.pack_td_ptr, pack_td_ptr);
    }
  }

  void find_col2row (std::vector<LO>& col2row) {
    Teuchos::RCP<const map_type> rowmap, colmap, dommap; {
      const auto& g = A_bcrs_->getCrsGraph();
      rowmap = g.getRowMap();
      colmap = g.getColMap();
      dommap = g.getDomainMap();
      TEUCHOS_ASSERT( ! (rowmap.is_null() || colmap.is_null(), dommap.is_null()));
    }
    const LO nrows = partptr_(partptr_.size() - 1);
    col2row.resize(A_bcrs_->getNodeNumCols(), Teuchos::OrdinalTraits<LO>::invalid());
    for (LO lr = 0; lr < nrows; ++lr) {
      const GO gid = rowmap->getGlobalElement(lr);
      TEUCHOS_ASSERT(gid != Teuchos::OrdinalTraits<GO>::invalid());
      if ( ! dommap->isNodeGlobalElement(gid)) continue;
      const LO lc = colmap->getLocalElement(gid);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(
        lc == Teuchos::OrdinalTraits<LO>::invalid(), get_msg_prefix() << "GID " << gid
        << " gives an invalid local column.");
      col2row[lc] = lr;
    }
  }

  // Wrappers to the compute_b_minus_Rx implementations.
  // Y := B - R X, where A = D + R.
  void compute_b_minus_Rx (const mv_type& B, const mv_type& X, mv_type& Y) const {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
#endif
    const auto Bv = B.template getLocalView<tpetra_device_type>();
    const auto Xv = X.template getLocalView<tpetra_device_type>();
    const auto Yv = Y.template getLocalView<tpetra_device_type>();
    if (amd_.is_tpetra_block_crs) {
      const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
      BlockTriDiContainerDetails::compute_b_minus_Rx(Bv, Xv, Yv, A_bcrs_->getBlockSize(),
                                                     amd_.rowptr, amd_.A_colindsub,
                                                     g.row_map, g.entries, amd_.tpetra_values);
    }
  }

  template <typename MvView>
  void compute_b_minus_Rx (const mv_type& B, const MvView& Xv,
                           typename Tridiags::PackedMultiVector& Y,
                           const bool first_apply) const {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
#endif
    const auto Bv = B.template getLocalView<tpetra_device_type>();
    if (amd_.is_tpetra_block_crs) {
      const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
      BlockTriDiContainerDetails::compute_b_minus_Rx(
        first_apply ? amd_.rowptr : amd_.rowptr_remote,
        first_apply ? amd_.A_colindsub : amd_.A_colindsub_remote,
        g.row_map, g.entries, amd_.tpetra_values,
        Bv, Xv, Y, part2rowidx0_, part2packrowidx0_, lclrow_, dm2cm_, first_apply);
    }
  }

  template <typename MvView>
  void compute_b_minus_Rx (const mv_type& B, const MvView& X_owned, const MvView& X_remote,
                           typename Tridiags::PackedMultiVector& Y) const {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::compute_b_minus_Rx");
#endif
    const auto Bv = B.template getLocalView<tpetra_device_type>();
    if (amd_.is_tpetra_block_crs) {
      const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
      BlockTriDiContainerDetails::compute_b_minus_Rx(
        amd_.rowptr, amd_.A_colindsub,
        g.row_map, g.entries, amd_.tpetra_values,
        Bv, X_owned, X_remote, Y, part2rowidx0_, part2packrowidx0_, lclrow_, dm2cm_);
    }
  }

public:
  Impl (BlockTriDiContainer<MatrixType, LocalScalarType>& container,
        const Teuchos::RCP<const row_matrix_type>& matrix,
        const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
        const Teuchos::RCP<const import_type>& importer, int overlapLevel,
        bool overlapCommAndComp = false, bool useSeqMethod = false)
    : container_(container)
  {
    init(matrix, partitions, importer, overlapLevel, overlapCommAndComp, useSeqMethod);
  }

  std::string describe () const {
    std::stringstream ss;
    ss << "seq_method " << seq_method_
       << " overlap_comm " << overlap_comm_
       << " dm2cm " << (dm2cm_.data() ? true : false);
    return ss.str();
  }

  // Top-level symbolic phase.
  void symbolic () {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::symbolic");
#endif

    const LO nrows = partptr_(partptr_.size() - 1);
    const auto lclrow = tohost(lclrow_);

    std::vector<LO> col2row;
    find_col2row(col2row);

    col_contiguous_ = col2row[0] == 0;
    if (col_contiguous_) {
      for (LO lr = 1; lr < nrows; ++lr)
        if (lclrow(col2row[lr-1]) + 1 != lclrow(col2row[lr])) {
          col_contiguous_ = false;
          break;
        }
    }

    { // Construct the D and R graphs in A = D + R.
      const auto& g = A_bcrs_->getCrsGraph().getLocalGraph();
      const auto& g_rowptr = g.row_map;
      TEUCHOS_ASSERT(g_rowptr.size() == static_cast<size_t>(nrows + 1));
      const auto& g_colind = g.entries;

      const auto rowidx2part = tohost(rowidx2part_);
      const auto part2rowidx0 = tohost(part2rowidx0_);

      //assume No overlap.
      std::vector<LO> lclrow2idx(nrows);
      for (LO i = 0; i < nrows; ++i)
        lclrow2idx[lclrow(i)] = i;

      // Count (block) nnzs in D and R.
      Size D_nnz = 0, R_nnz_owned = 0, R_nnz_remote = 0;
      for (LO lr = 0; lr < nrows; ++lr) {
        // LID -> index.
        const LO ri0 = lclrow2idx[lr];
        const LO pi0 = rowidx2part(ri0);
        for (Size j = g_rowptr(lr); j < g_rowptr(lr+1); ++j) {
          const LO lc = g_colind(j);
          const LO lc2r = col2row[lc];
          bool incr_R = false;
          do { // breakable
            if (lc2r == Teuchos::OrdinalTraits<LO>::invalid()) {
              incr_R = true;
              break;
            }
            const LO ri = lclrow2idx[lc2r];
            const LO pi = rowidx2part(ri);
            if (pi != pi0) {
              incr_R = true;
              break;
            }
            // Test for being in the tridiag. This is done in index space. In
            // LID space, tridiag LIDs in a row are not necessarily related by
            // {-1, 0, 1}.
            if (ri0 + 1 >= ri && ri0 <= ri + 1)
              ++D_nnz;
            else
              incr_R = true;
          } while (0);
          if (incr_R) {
            if (lc < nrows) ++R_nnz_owned;
            else ++R_nnz_remote;
          }
        }
      }
      if (! overlap_comm_) {
        R_nnz_owned += R_nnz_remote;
        R_nnz_remote = 0;
      }

      { // Construct the D graph.
        btdm_.A_colindsub = LOList("btdm.A_colindsub", D_nnz);
        const auto D_A_colindsub = Kokkos::create_mirror_view(btdm_.A_colindsub);
        if (validate_)
          Kokkos::deep_copy(D_A_colindsub, Teuchos::OrdinalTraits<LO>::invalid());
        const auto partptr = tohost(partptr_);
        const auto flat_td_ptr = tohost(btdm_.flat_td_ptr);
        const LO nparts = partptr.size() - 1;
        for (LO pi0 = 0; pi0 < nparts; ++pi0) {
          const LO part_ri0 = part2rowidx0_(pi0);
          LO offset = 0;
          for (LO ri0 = partptr(pi0); ri0 < partptr(pi0+1); ++ri0) {
            const LO td_row_os = Tridiags::row2ind(ri0 - part_ri0) + offset;
            offset = 1;
            const LO lr0 = lclrow(ri0);
            const Size j0 = g_rowptr(lr0);
            for (Size j = j0; j < g_rowptr(lr0+1); ++j) {
              const LO lc = g_colind(j);
              const LO lc2r = col2row[lc];
              if (lc2r == Teuchos::OrdinalTraits<LO>::invalid()) continue;
              const LO ri = lclrow2idx[lc2r];
              const LO pi = rowidx2part(ri);
              if (pi != pi0) continue;
              if (ri + 1 < ri0 || ri > ri0 + 1) continue;
              const LO row_entry = j - j0;
              D_A_colindsub(flat_td_ptr(pi0) + ((td_row_os + ri) - ri0)) = row_entry;
            }
          }
        }
        if (validate_)
          for (size_t i = 0; i < D_A_colindsub.size(); ++i)
            TEUCHOS_ASSERT(D_A_colindsub(i) != Teuchos::OrdinalTraits<LO>::invalid());
        Kokkos::deep_copy(btdm_.A_colindsub, D_A_colindsub);
        { // Allocate values.
          const auto packptr = tohost(packptr_);
          const LO npacks = packptr.size() - 1;
          LO nblks = 0; // Number of tridiag blocks, accounting for packing.
          for (LO pai = 0; pai < npacks; ++pai) {
            const LO pti = packptr(pai);
            const LO inrows = partptr(pti+1) - partptr(pti);
            nblks += Tridiags::row2ind(inrows);
          }
          const auto bs = A_bcrs_->getBlockSize();
          btdm_.values = typename Tridiags::Values("btdm.values", nblks, bs, bs);
          if (vectorization_traits::vector_length > 1)
            SetTridiagsToI(btdm_, packptr_).run();
        }
      }

      { // Construct the R graph.
        amd_.is_tpetra_block_crs = true;
        amd_.rowptr = SizeList("amd.rowptr", nrows + 1);
        amd_.A_colindsub = LOList("amd.A_colindsub", R_nnz_owned);
        const auto R_rowptr = Kokkos::create_mirror_view(amd_.rowptr);
        const auto R_A_colindsub = Kokkos::create_mirror_view(amd_.A_colindsub);
        R_rowptr(0) = 0;
        if (overlap_comm_) {
          amd_.rowptr_remote = SizeList("amd.rowptr_remote", nrows + 1);
          amd_.A_colindsub_remote = LOList("amd.A_colindsub_remote", R_nnz_remote);
        }
        const auto R_rowptr_remote = Kokkos::create_mirror_view(amd_.rowptr_remote);
        const auto R_A_colindsub_remote = Kokkos::create_mirror_view(amd_.A_colindsub_remote);
        if (overlap_comm_) R_rowptr_remote(0) = 0;
        for (LO lr = 0; lr < nrows; ++lr) {
          const LO ri0 = lclrow2idx[lr];
          const LO pi0 = rowidx2part(ri0);
          R_rowptr(lr+1) = R_rowptr(lr);
          if (overlap_comm_) R_rowptr_remote(lr+1) = R_rowptr_remote(lr);
          const Size j0 = g_rowptr(lr);
          for (Size j = j0; j < g_rowptr(lr+1); ++j) {
            const LO lc = g_colind(j);
            const LO lc2r = col2row[lc];
            if (lc2r != Teuchos::OrdinalTraits<LO>::invalid()) {
              const LO ri = lclrow2idx[lc2r];
              const LO pi = rowidx2part(ri);
              if (pi == pi0 && ri + 1 >= ri0 && ri <= ri0 + 1)
                continue;
            }
            const LO row_entry = j - j0;
            if ( ! overlap_comm_ || lc < nrows) {
              R_A_colindsub(R_rowptr(lr+1)) = row_entry;
              ++R_rowptr(lr+1);
            } else {
              R_A_colindsub_remote(R_rowptr_remote(lr+1)) = row_entry;
              ++R_rowptr_remote(lr+1);
            }
          }
        }
        TEUCHOS_ASSERT(R_rowptr(nrows) == R_nnz_owned);
        Kokkos::deep_copy(amd_.rowptr, R_rowptr);
        Kokkos::deep_copy(amd_.A_colindsub, R_A_colindsub);
        if (overlap_comm_) {
          TEUCHOS_ASSERT(R_rowptr_remote(nrows) == R_nnz_remote);
          Kokkos::deep_copy(amd_.rowptr_remote, R_rowptr_remote);
          Kokkos::deep_copy(amd_.A_colindsub_remote, R_A_colindsub_remote);
        }
        // Allocate or view values.
        if (amd_.is_tpetra_block_crs)
          amd_.tpetra_values = (const_cast<block_crs_matrix_type*>(A_bcrs_.get())->
                                template getValues<typename tpetra_device_type::memory_space>());
      }
    }
  }

  static void debug_print (const Tridiags& t) {
    const auto& v = t.values;
    std::stringstream ss;
    ss << "v = [";
    for (size_t pi = 0; pi < t.pack_td_ptr.size() - 1; ++pi)
      for (LO i1 = 0; i1 < vectorization_traits::vector_length; ++i1)
        for (Size ind = t.pack_td_ptr(pi); ind < t.pack_td_ptr(pi+1); ++ind) {
          const auto i = Kokkos::make_pair(ind,i1);
          for (LO j = 0; j < v.extent_int(1); ++j)
            for (LO k = 0; k < v.extent_int(2); ++k)
              ss << " " << Details::Batched::idx(v,i,j,k);
        }
    ss << "]\n";
    std::cout << ss.str();
  }

  // Top-level numeric phase.
  void numeric (const magnitude_type add_to_diag = 0) {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::numeric");
#endif
    TEUCHOS_ASSERT( ! A_bcrs_.is_null());
    const auto A_rowptr = A_bcrs_->getCrsGraph().getLocalGraph().row_map;
    ExtractAndFactorizeTridiags(btdm_, partptr_, lclrow_, packptr_, A_rowptr,
                                amd_.tpetra_values, add_to_diag).run();
    TEUCHOS_ASSERT(amd_.is_tpetra_block_crs);
  }

  // Top-level apply phase.
  int applyInverseJacobi (const mv_type& X, mv_type& Y, const impl_scalar_type& damping_factor,
                          bool y_is_zero, const int max_num_sweeps, magnitude_type tolerance,
                          const int check_tol_every) const {
#ifdef HAVE_IFPACK2_BLOCKTRIDICONTAINER_TIMERS
    TEUCHOS_FUNC_TIME_MONITOR("BlockTriDi::applyInverseJacobi");
#endif

    const bool do_norm = tolerance > Kokkos::ArithTraits<magnitude_type>::zero();
    if (do_norm) tolerance *= tolerance;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
      do_norm && seq_method_,
      "The seq method for applyInverseJacobi, which in any case is for developer use only, " <<
      "does not support norm-based termination.");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(max_num_sweeps <= 0, "Maximum number of sweeps must be >= 1.");

    // Set up work space if needed.
    if (pmv_.size() == 0 || Y.getNumVectors() != pmv_.extent(2)) {
      //todo Instead of reallocating, could get a subview if nvec is smaller
      // than the capacity.
      const auto nvec = Y.getNumVectors();
      const auto bs = A_bcrs_->getBlockSize();
      if (seq_method_)
        Y_remote_ = Teuchos::rcp(new mv_type(importer_->getTargetMap(), nvec));
      const auto nblks = part2packrowidx0_back;
      pmv_ = typename Tridiags::PackedMultiVector("pmv_" , nblks, bs, nvec);
    }

    have_norms_ = do_norm;
    if (do_norm) {
      if ( ! norm_mgr_) norm_mgr_ = std::make_shared<NormManager>(Y.getMap()->getComm());
      norm_mgr_->resize(A_bcrs_->getBlockSize(), Y.getNumVectors());
      norm_mgr_->set_check_frequency(check_tol_every);
    }

    int sweep;
    for (sweep = 0; sweep < max_num_sweeps; ++sweep) {
      if (y_is_zero) {
        // pmv := x(lclrow)
        PermuteAndRepack(X, pmv_, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_).run();

      } else {
        if (seq_method_) {
          // y := x - R y.
          Y_remote_->doImport(Y, *importer_, Tpetra::REPLACE);
          compute_b_minus_Rx(X, *Y_remote_, Y);
          // pmv := y(lclrow).
          PermuteAndRepack(Y, pmv_, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_).run();

        } else {

          // Fused y := x - R y and pmv := y(lclrow).
          if (overlap_comm_ || ! async_import_) {
            // Do R y_owned followed by R y_remote.
            if (async_import_)
              async_import_->async_setup(Y);
            compute_b_minus_Rx(X, Y.template getLocalView<tpetra_device_type>(), pmv_, true);
            if (do_norm && sweep > 0 && norm_mgr_->check_done(sweep, tolerance)) {
              if (async_import_) async_import_->cancel();
              break;
            }
            if (async_import_) {
              async_import_->sync_receive();
              compute_b_minus_Rx(X, async_import_->get_mv_remote(), pmv_, false);
            }
          } else {
            // Use our custom import, but don't overlap comm with R y_owned.
            async_import_->sync_exchange(Y);
            if (do_norm && sweep > 0 && norm_mgr_->check_done(sweep, tolerance)) break;
            compute_b_minus_Rx(X, Y.template getLocalView<tpetra_device_type>(),
                               async_import_->get_mv_remote(), pmv_);
          }

        }
      }

      // pmv := inv(D) pmv.
      Solve<layout_type, void>(btdm_, packptr_, part2packrowidx0_, pmv_).run();

      // y(lclrow) = (b - a) y(lclrow) + a pmv, with b = 1 always.
      PermuteAndRepack(pmv_, Y, damping_factor, part2rowidx0_, part2packrowidx0_, lclrow_, packptr_,
                       y_is_zero, do_norm ? norm_mgr_->get_buffer() : nullptr).run();

      if (do_norm) {
        if (sweep + 1 == max_num_sweeps) {
          norm_mgr_->ireduce(sweep, true);
          norm_mgr_->check_done(sweep + 1, tolerance, true);
        } else {
          norm_mgr_->ireduce(sweep);
        }
      }

      y_is_zero = false;
    }

    // Sqrt the norms for the caller's use.
    if (do_norm) norm_mgr_->finalize();

    return sweep;
  }

  // Report norms to caller.
  const magnitude_type* get_norms0 () const {
    if ( ! have_norms_) return nullptr;
    return norm_mgr_->get_norms0();
  }

  const magnitude_type* get_norms_final () const {
    if ( ! have_norms_) return nullptr;
    return norm_mgr_->get_norms_final();
  }

  int get_block_size () const { return pmv_.extent(1); }
  int get_nvec () const { return pmv_.extent(2); }
};

// Base class for any unimplemented typename combinations.
template <typename MatrixType, typename LocalScalarType, typename ExeSpace>
struct UnImpl {
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> mv_type;
  typedef Tpetra::Import<local_ordinal_type, global_ordinal_type, node_type> import_type;
  typedef Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> row_matrix_type;

  std::string describe () const { return ""; }
  int get_block_size () const { return 0; }
  int get_nvec () const { return 0; }
  const magnitude_type* get_norms0 () const { return nullptr; }
  const magnitude_type* get_norms_final () const { return nullptr; }
  void symbolic () {}
  void numeric (const magnitude_type add_to_diag = 0) {}
  int applyInverseJacobi (const mv_type& X, mv_type& Y, const scalar_type damping_factor,
                          bool y_is_zero, const int max_num_sweeps, magnitude_type tolerance,
                          const int check_tol_every) const
  { return 0; }
};

#if defined HAVE_STOKHOS_IFPACK2
// BlockTriDiContainer is built on KokkosBatch's SIMD type, which doesn't work
// directly with Stokhos composite types. Hence, catch these with UnImpl.
#define IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Type)                \
  template <typename MatrixType, typename T, typename ExeSpace>         \
  struct Impl<MatrixType, Type<T>, ExeSpace> : UnImpl<MatrixType, Type<T>, ExeSpace> { \
    typedef UnImpl<MatrixType, Type<T>, ExeSpace> ui;                   \
    Impl (BlockTriDiContainer<MatrixType, typename ui::scalar_type>& container, \
          const Teuchos::RCP<const typename ui::row_matrix_type>& matrix, \
          const Teuchos::Array<Teuchos::Array<typename ui::local_ordinal_type> >& partitions, \
          const Teuchos::RCP<const typename ui::import_type>& importer, int overlapLevel, \
          bool overlapCommAndComp = false, bool useSeqMethod = false) { \
      TEUCHOS_TEST_FOR_EXCEPT_MSG(                                      \
        true, "BlockTriDiContainer is not currently supported for Stokhos composite types."); \
    }                                                                   \
  };
IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Sacado::MP::Vector)
IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS(Sacado::UQ::PCE)
#undef IFPACK2_BLOCKTRIDICONTAINER_UNIMPL_STOKHOS
#endif

} // namespace BlockTriDiContainerDetails

template <typename MatrixType, typename LocalScalarType>
BlockTriDiContainer<MatrixType, LocalScalarType>::
BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                     const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                     const Teuchos::RCP<const import_type>& importer,
                     int OverlapLevel, scalar_type DampingFactor)
  : Container<MatrixType>(matrix, partitions, importer, OverlapLevel, DampingFactor)
{
  impl_ = Teuchos::rcp(new ImplWithES(*this, matrix, partitions, importer, OverlapLevel, false));
}

template <typename MatrixType, typename LocalScalarType>
BlockTriDiContainer<MatrixType, LocalScalarType>::
BlockTriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                     const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                     const bool overlapCommAndComp, const bool useSeqMethod)
  : Container<MatrixType>(matrix, partitions, Teuchos::null, 0,
                          Kokkos::ArithTraits<magnitude_type>::one())
{
  impl_ = Teuchos::rcp(new ImplWithES(*this, matrix, partitions, Teuchos::null, 0, overlapCommAndComp,
                                      useSeqMethod));
}

template <typename MatrixType, typename LocalScalarType>
BlockTriDiContainer<MatrixType, LocalScalarType>::~BlockTriDiContainer ()
{
}

template <typename MatrixType, typename LocalScalarType>
bool BlockTriDiContainer<MatrixType, LocalScalarType>::isInitialized () const
{
  return IsInitialized_;
}

template <typename MatrixType, typename LocalScalarType>
bool BlockTriDiContainer<MatrixType, LocalScalarType>::isComputed () const
{
  return IsComputed_;
}

template <typename MatrixType, typename LocalScalarType>
void BlockTriDiContainer<MatrixType, LocalScalarType>::
setParameters (const Teuchos::ParameterList& /* List */)
{
  // the solver doesn't currently take any parameters
}

template <typename MatrixType, typename LocalScalarType>
void BlockTriDiContainer<MatrixType, LocalScalarType>::initialize ()
{
  IsInitialized_ = true;
  // We assume that if you called this method, you intend to recompute
  // everything.
  IsComputed_ = false;
  TEUCHOS_ASSERT( ! impl_.is_null());
  impl_->symbolic();
}

template <typename MatrixType, typename LocalScalarType>
void BlockTriDiContainer<MatrixType, LocalScalarType>::compute ()
{
  IsComputed_ = false;
  if ( ! this->isInitialized())
    this->initialize();
  impl_->numeric();
  IsComputed_ = true;
}

template <typename MatrixType, typename LocalScalarType>
void BlockTriDiContainer<MatrixType, LocalScalarType>::clearBlocks ()
{
  impl_ = Teuchos::null;
  IsInitialized_ = false;
  IsComputed_ = false;
}

template <typename MatrixType, typename LocalScalarType>
void BlockTriDiContainer<MatrixType, LocalScalarType>
::applyInverseJacobi (const mv_type& X, mv_type& Y, bool zeroStartingSolution,
                      int numSweeps) const
{
  impl_->applyInverseJacobi(X, Y, this->DampingFactor_, zeroStartingSolution, numSweeps,
                            Kokkos::ArithTraits<magnitude_type>::zero(), 1);
}

template <typename MatrixType, typename LocalScalarType>
BlockTriDiContainer<MatrixType, LocalScalarType>::ComputeParameters::ComputeParameters ()
  : addRadiallyToDiagonal(Kokkos::ArithTraits<magnitude_type>::zero())
{}

template <typename MatrixType, typename LocalScalarType>
typename BlockTriDiContainer<MatrixType, LocalScalarType>::ComputeParameters
BlockTriDiContainer<MatrixType, LocalScalarType>::createDefaultComputeParameters () const
{
  return ComputeParameters();
}

template <typename MatrixType, typename LocalScalarType>
void BlockTriDiContainer<MatrixType, LocalScalarType>::compute (const ComputeParameters& in)
{
  IsComputed_ = false;
  if ( ! this->isInitialized())
    this->initialize();
  impl_->numeric(in.addRadiallyToDiagonal);
  IsComputed_ = true;
}

template <typename MatrixType, typename LocalScalarType>
BlockTriDiContainer<MatrixType, LocalScalarType>::ApplyParameters::ApplyParameters ()
  : zeroStartingSolution(false), dampingFactor(Kokkos::ArithTraits<magnitude_type>::one()),
    maxNumSweeps(1), tolerance(Kokkos::ArithTraits<magnitude_type>::zero()),
    checkToleranceEvery(1)
{}

template <typename MatrixType, typename LocalScalarType>
typename BlockTriDiContainer<MatrixType, LocalScalarType>::ApplyParameters
BlockTriDiContainer<MatrixType, LocalScalarType>::createDefaultApplyParameters () const
{
  ApplyParameters in;
  in.dampingFactor = this->DampingFactor_;
  return in;
}

template <typename MatrixType, typename LocalScalarType>
int BlockTriDiContainer<MatrixType, LocalScalarType>
::applyInverseJacobi (const mv_type& X, mv_type& Y, const ApplyParameters& in) const
{
  return impl_->applyInverseJacobi(X, Y, in.dampingFactor, in.zeroStartingSolution, in.maxNumSweeps,
                                   in.tolerance, in.checkToleranceEvery);
}

template <typename MatrixType, typename LocalScalarType>
const Teuchos::ArrayRCP<const typename BlockTriDiContainer<MatrixType, LocalScalarType>::magnitude_type>
BlockTriDiContainer<MatrixType, LocalScalarType>::getNorms0 () const {
  const auto p = impl_->get_norms0();
  return Teuchos::arcp(p, 0, p ? impl_->get_nvec() : 0, false);
}

template <typename MatrixType, typename LocalScalarType>
const Teuchos::ArrayRCP<const typename BlockTriDiContainer<MatrixType, LocalScalarType>::magnitude_type>
BlockTriDiContainer<MatrixType, LocalScalarType>::getNormsFinal () const {
  const auto p = impl_->get_norms_final();
  return Teuchos::arcp(p, 0, p ? impl_->get_nvec() : 0, false);
}

template <typename MatrixType, typename LocalScalarType>
void BlockTriDiContainer<MatrixType, LocalScalarType>
::apply (HostView& X, HostView& Y, int blockIndex, int stride, Teuchos::ETransp mode,
         scalar_type alpha, scalar_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(
    true, "BlockTriDiContainer::apply is not implemented. You may have reached this message "
    << "because you want to use this container's performance-portable Jacobi iteration. In "
    << "that case, set \"relaxation: type\" to \"MT Split Jacobi\" rather than \"Jacobi\".");
}

template <typename MatrixType, typename LocalScalarType>
void BlockTriDiContainer<MatrixType, LocalScalarType>
::weightedApply (HostView& X, HostView& Y, HostView& D, int blockIndex, int stride,
                 Teuchos::ETransp mode, scalar_type alpha, scalar_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "BlockTriDiContainer::weightedApply is not implemented.");
}

template <typename MatrixType, typename LocalScalarType>
std::ostream& BlockTriDiContainer<MatrixType, LocalScalarType>::print (std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return os;
}

template <typename MatrixType, typename LocalScalarType>
std::string BlockTriDiContainer<MatrixType, LocalScalarType>::description () const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }

  oss << "}";
  return oss.str();
}

template <typename MatrixType, typename LocalScalarType>
void
BlockTriDiContainer<MatrixType, LocalScalarType>::
describe (Teuchos::FancyOStream& os,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << endl
     << "Ifpack2::BlockTriDiContainer" << endl
     << "Number of blocks        = " << this->numBlocks_ << endl
     << "isInitialized()         = " << IsInitialized_ << endl
     << "isComputed()            = " << IsComputed_ << endl
     << "impl                    = " << impl_->describe() << endl
     << "================================================================================" << endl
     << endl;
}

template <typename MatrixType, typename LocalScalarType>
std::string BlockTriDiContainer<MatrixType, LocalScalarType>::getName() { return "BlockTriDi"; }

#define IFPACK2_BLOCKTRIDICONTAINER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::BlockTriDiContainer< Tpetra::RowMatrix<S, LO, GO, N>, S >;

} // namespace Ifpack2

#endif
