// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef INCLUDE_SHYLU_HTS_IMPL_DEF_HPP
#define INCLUDE_SHYLU_HTS_IMPL_DEF_HPP

#ifdef _OPENMP
# include <omp.h>
#else
namespace Experimental {
namespace htsimpl {
// Protect against a header that has #define'd replacements for OpenMP
// functions.
#ifndef omp_get_max_threads
inline int omp_get_max_threads () { return 1; }
#endif
#ifndef omp_set_num_threads
inline void omp_set_num_threads (const int&) {}
#endif
#ifndef omp_get_num_threads
inline int omp_get_num_threads () { return 1; }
#endif
#ifndef omp_get_thread_num
inline int omp_get_thread_num () { return 0; }
#endif
}
}
#endif

#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>

#ifdef HAVE_SHYLU_NODEHTS_COMPLEX
# include <complex>
#endif
#ifdef HAVE_SHYLU_NODEHTS_MKL
# include <mkl.h>
# ifdef HAVE_SHYLU_NODEHTS_COMPLEX
#  include <type_traits> // std::is_same, etc.
# endif // HAVE_SHYLU_NODEHTS_COMPLEX
#endif

#ifdef HAVE_SHYLU_NODEHTS_KOKKOSKERNELS
# include <Kokkos_ArithTraits.hpp>
#endif

#include "shylu_hts_impl.hpp"

namespace Experimental {
namespace htsimpl {

static const int parfor_static_size = 20;

#if defined(HAVE_SHYLU_NODEHTS_BLAS) || defined(HAVE_SHYLU_NODEHTS_MKL)
typedef int blas_int;

template<typename T> void gemm(
  char transa, char transb, blas_int m, blas_int nrhs, blas_int n, T alpha,
  const T* a, blas_int lda, const T* b, blas_int ldb, T beta,
  T* c, blas_int ldc);

extern "C" {
  void F77_BLAS_MANGLE(sgemm,SGEMM)(
    const char*, const char*, const blas_int*, const blas_int*, const blas_int*,
    const float*, const float*, const blas_int*, const float*, const blas_int*,
    const float*, float*, const blas_int*);
  void F77_BLAS_MANGLE(dgemm,DGEMM)(
    const char*, const char*, const blas_int*, const blas_int*, const blas_int*,
    const double*, const double*, const blas_int*, const double*, const blas_int*,
    const double*, double*, const blas_int*);
#ifdef HAVE_SHYLU_NODEHTS_COMPLEX
  void F77_BLAS_MANGLE(cgemm,CGEMM)(
    const char*, const char*, const blas_int*, const blas_int*, const blas_int*,
    const std::complex<float>*, const std::complex<float>*, const blas_int*,
    const std::complex<float>*, const blas_int*, const std::complex<float>*,
    std::complex<float>*, const blas_int*);
  void F77_BLAS_MANGLE(zgemm,ZGEMM)(
    const char*, const char*, const blas_int*, const blas_int*, const blas_int*,
    const std::complex<double>*, const std::complex<double>*, const blas_int*,
    const std::complex<double>*, const blas_int*, const std::complex<double>*,
    std::complex<double>*, const blas_int*);
#endif
}

template<> inline void gemm<float> (
  char transa, char transb, blas_int m, blas_int nrhs, blas_int n, float alpha,
  const float* a, blas_int lda, const float* b, blas_int ldb, float beta,
  float* c, blas_int ldc)
{
  F77_BLAS_MANGLE(sgemm,SGEMM)(
    &transa, &transb, &m, &nrhs, &n, &alpha, a, &lda,
    b, &ldb, &beta, c, &ldc);
}

template<> inline void gemm<double> (
  char transa, char transb, blas_int m, blas_int nrhs, blas_int n, double alpha,
  const double* a, blas_int lda, const double* b, blas_int ldb, double beta,
  double* c, blas_int ldc)
{
  F77_BLAS_MANGLE(dgemm,DGEMM)(
    &transa, &transb, &m, &nrhs, &n, &alpha, a, &lda,
    b, &ldb, &beta, c, &ldc);
}

#ifdef HAVE_SHYLU_NODEHTS_COMPLEX
template<> inline void gemm<std::complex<float> > (
  char transa, char transb, blas_int m, blas_int nrhs, blas_int n,
  std::complex<float> alpha, const std::complex<float>* a, blas_int lda,
  const std::complex<float>* b, blas_int ldb, std::complex<float> beta,
  std::complex<float>* c, blas_int ldc)
{
  F77_BLAS_MANGLE(cgemm,CGEMM)(
    &transa, &transb, &m, &nrhs, &n, &alpha,
    a, &lda,
    b, &ldb, &beta,
    c, &ldc);
}

template<> inline void gemm<std::complex<double> > (
  char transa, char transb, blas_int m, blas_int nrhs, blas_int n,
  std::complex<double> alpha, const std::complex<double>* a, blas_int lda,
  const std::complex<double>* b, blas_int ldb, std::complex<double> beta,
  std::complex<double>* c, blas_int ldc)
{
  F77_BLAS_MANGLE(zgemm,ZGEMM)(
    &transa, &transb, &m, &nrhs, &n, &alpha,
    a, &lda,
    b, &ldb, &beta,
    c, &ldc);
}
#endif
#endif

#ifdef HAVE_SHYLU_NODEHTS_MKL
// sparse A * dense x
template<typename T> void hts_mkl_csrmm(
  const bool transp, const MKL_INT m, const MKL_INT n, const T* d,
  const MKL_INT* ir, const MKL_INT* jc, const T* x, const int ldx,
  T* y, const int ldy, const MKL_INT nrhs);

template<> inline void hts_mkl_csrmm<float> (
  const bool transp, const MKL_INT m, const MKL_INT n, const float* d,
  const MKL_INT* ir, const MKL_INT* jc, const float* x, const int ldx,
  float* y, const int ldy, const MKL_INT nrhs)
{
  char transa = transp ? 'T' : 'N';
  static const char A_descr[6] = {'G', '*', '*', 'C', '*', '*'};
  float alpha = -1, beta = 1;
  for (int k = 0; k < nrhs; ++k)
    mkl_scsrmv(
      &transa, const_cast<MKL_INT*>(&m), const_cast<MKL_INT*>(&n),
      &alpha, const_cast<char*>(A_descr), const_cast<float*>(d),
      const_cast<MKL_INT*>(jc), const_cast<MKL_INT*>(ir),
      const_cast<MKL_INT*>(ir+1), const_cast<float*>(x + k*ldx), &beta,
      y + k*ldy);
}

template<> inline void hts_mkl_csrmm<double> (
  const bool transp, const MKL_INT m, const MKL_INT n, const double* d,
  const MKL_INT* ir, const MKL_INT* jc, const double* x, const int ldx,
  double* y, const int ldy, const MKL_INT nrhs)
{
  char transa = transp ? 'T' : 'N';
  static const char A_descr[6] = {'G', '*', '*', 'C', '*', '*'};
  double alpha = -1, beta = 1;
  for (int k = 0; k < nrhs; ++k)
    mkl_dcsrmv(
      &transa, const_cast<MKL_INT*>(&m), const_cast<MKL_INT*>(&n),
      &alpha, const_cast<char*>(A_descr), const_cast<double*>(d),
      const_cast<MKL_INT*>(jc), const_cast<MKL_INT*>(ir),
      const_cast<MKL_INT*>(ir+1), const_cast<double*>(x + k*ldx), &beta,
      y + k*ldy);
}

#ifdef HAVE_SHYLU_NODEHTS_COMPLEX
template<> inline void hts_mkl_csrmm<std::complex<float> > (
  const bool transp, const MKL_INT m, const MKL_INT n,
  const std::complex<float>* d, const MKL_INT* ir, const MKL_INT* jc,
  const std::complex<float>* x, const int ldx, std::complex<float>* y,
  const int ldy, const MKL_INT nrhs)
{
  char transa = transp ? 'T' : 'N';
  static const char A_descr[6] = {'G', '*', '*', 'C', '*', '*'};
  std::complex<float> alpha(-1, 0), beta(1, 0);

  static_assert (sizeof (std::complex<float>) == sizeof (MKL_Complex8),
                 "This code assumes that sizeof(std::complex<float>) == "
                 "sizeof(MKL_Complex8).");

  const MKL_Complex8* alpha_ptr = const_cast<const MKL_Complex8*> (reinterpret_cast<MKL_Complex8*> (&alpha));
  const MKL_Complex8* beta_ptr = const_cast<const MKL_Complex8*> (reinterpret_cast<MKL_Complex8*> (&beta));
  const MKL_Complex8* val = reinterpret_cast<const MKL_Complex8*> (d);

  for (int k = 0; k < nrhs; ++k) {
    const MKL_Complex8* x_k = reinterpret_cast<const MKL_Complex8*> (x + k*ldx);
    MKL_Complex8* y_k = reinterpret_cast<MKL_Complex8*>(y + k*ldy);

    mkl_ccsrmv(
      &transa, &m, &n,
      alpha_ptr, const_cast<const char*>(A_descr), val,
      jc, ir,
      ir + 1, x_k,
      beta_ptr, y_k);
  }
}

template<> inline void hts_mkl_csrmm<std::complex<double> > (
  const bool transp, const MKL_INT m, const MKL_INT n,
  const std::complex<double>* d, const MKL_INT* ir, const MKL_INT* jc,
  const std::complex<double>* x, const int ldx, std::complex<double>* y,
  const int ldy, const MKL_INT nrhs)
{
  char transa = transp ? 'T' : 'N';
  static const char A_descr[6] = {'G', '*', '*', 'C', '*', '*'};
  std::complex<double> alpha(-1, 0), beta(1, 0);

  static_assert (sizeof (std::complex<double>) == sizeof (MKL_Complex16),
                 "This code assumes that sizeof(std::complex<double>) == "
                 "sizeof(MKL_Complex16).");

  const MKL_Complex16* alpha_ptr = const_cast<const MKL_Complex16*> (reinterpret_cast<MKL_Complex16*> (&alpha));
  const MKL_Complex16* beta_ptr = const_cast<const MKL_Complex16*> (reinterpret_cast<MKL_Complex16*> (&beta));
  const MKL_Complex16* val = reinterpret_cast<const MKL_Complex16*> (d);

  for (int k = 0; k < nrhs; ++k) {
    const MKL_Complex16* x_k = reinterpret_cast<const MKL_Complex16*> (x + k*ldx);
    MKL_Complex16* y_k = reinterpret_cast<MKL_Complex16*>(y + k*ldy);

    mkl_zcsrmv(
      &transa, &m, &n,
      alpha_ptr, const_cast<const char*>(A_descr), val,
      jc, ir,
      ir+1, x_k,
      beta_ptr, y_k);
  }
}
#endif
#endif

template<typename T> inline void touch (T* const p, const size_t n,
                                        const T& init = T()) {
  // 1 KB should be a safe lower bound on page size. Touch enough to touch every
  // page; I don't think there's any need to touch more memory than that.
#if ! defined __MIC__
  for (size_t i = 0; i < n; i += 1024 / sizeof(T))
    p[i] = init;
  // Make sure the last part is touched.
  if (n) p[n-1] = init;
#endif
}

template<typename T> inline T*
allocn (const size_t n, const bool first_touch = false) {
  if ( ! n) return 0;
  T* p = new T[n];
  if (first_touch) touch(p, n);
  return p;
}
template<typename T> inline void deln (T*& p) {
  if (p) delete[] p;
  p = 0;
}
template<typename T> inline void deln_const (const T* p) {
  if (p) delete[] p;
}
template<typename T> inline void del (T*& p) {
  if (p) delete p;
  p = 0;
}

// For exception safety when allocating.
template<typename T> class Allocnator {
  T* p_;
  bool dealloc_;
public:
  Allocnator () : p_(0), dealloc_(true) {}
  // Try to allocate memory.
  Allocnator (const size_t n, const char* msg = "",
              const bool first_touch = false)
    : p_(0), dealloc_(true)
  { resize(n, msg, first_touch); }
  void resize (const size_t n, const char* msg = "",
               const bool first_touch = false) {
    if (p_ && dealloc_) deln<T>(p_);
    p_ = 0;
    dealloc_ = true;
    if ( ! n) return;
    try { p_ = allocn<T>(n, first_touch); }
    catch (...) {
      throw hts::Exception(std::string(msg) + ": failed to allocate.");
    }
  }
  T* get () { return p_; }
  // Release the pointer to the user and subsequently don't dealloc.
  T* release () { dealloc_ = false; return p_; }
  // Dealloc only if the user hasn't released the pointer yet.
  ~Allocnator () { if (p_ && dealloc_) deln<T>(p_); }
};

template<typename T>
inline void Array<T>::init () {
  n_ = cap_ = 0;
  p_ = 0;
}

template<typename T>
inline Array<T>::Array (std::size_t n)
  : p_(0), n_(0), cap_(0)
{ optclear_and_resize(n); }

template<typename T>
inline Array<T>::Array (std::size_t n, const T& val)
  : p_(0), n_(0), cap_(0)
{ optclear_and_resize(n, val); }

template<typename T>
inline void Array<T>::clear () {
  n_ = cap_ = 0;
  deln(p_);
}

template<typename T>
inline void Array<T>::optclear_and_reserve (std::size_t n) {
  n_ = 0;
  if (n <= cap_) return;
  clear();
  p_ = allocn<T>(n);
  cap_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_reserve_ft (std::size_t n) {
  n_ = 0;
  if (n <= cap_) return;
  clear();
  p_ = allocn<T>(n, true);
  cap_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize (std::size_t n) {
  if (n <= cap_) {
    n_ = n;
    return;
  }
  optclear_and_reserve(n);
  n_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize_ft (std::size_t n) {
  if (n <= cap_) {
    n_ = n;
    return;
  }
  optclear_and_reserve_ft(n);
  n_ = n;
}

template<typename T>
inline void Array<T>::optclear_and_resize (std::size_t n, const T& val) {
  optclear_and_resize(n);
  for (std::size_t i = 0; i < n_; ++i)
    p_[i] = val;
}

template<typename T>
inline void Array<T>::unsafe_push_back (const T& e) {
  assert(n_ < cap_);
  p_[n_++] = e;
}

template<typename T> inline T square (const T& x) { return x*x; }

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
set_options (const typename ihts::Options& os, Options& od) {
  od.min_block_size = os.min_block_size;
  od.min_parallel_rows = os.min_parallel_rows;
  od.pp_min_block_size = os.pp_min_block_size;
  od.min_dense_density = os.min_dense_density;
  od.ls_blk_sz = os.levelset_block_size;
  od.lset_min_size = os.min_lset_size;
  od.lset_min_size_scale_with_nthreads = os.lset_min_size_scale_with_nthreads;
  od.lset_max_bad_fraction = os.lset_max_bad_fraction;
  od.profile = os.profile;
  od.printlvl = os.print_level;
}

template<typename Int, typename Size, typename Sclr> Impl<Int, Size, Sclr>::
Options::Options () {
  set_options(typename HTS<Int, Size, Sclr>::Options(), *this);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::Options::print (std::ostream& os) const {
  os << " hts min_block_size " << min_block_size
     << " min_parallel_rows " << min_parallel_rows
     << " pp_min_block_size " << pp_min_block_size
     << " min_dense_density " << min_dense_density
     << " ls_blk_sz " << ls_blk_sz
     << " lset_min_size " << lset_min_size
     << " lset_min_size_scale_with_nthreads "
     << lset_min_size_scale_with_nthreads
     << " lset_max_bad_fraction " << lset_max_bad_fraction
     << " profile " << profile;
}

static inline void print_compiletime_options(std::ostream& os) {
#ifdef HAVE_SHYLU_NODEHTS_BLAS
  os << " HAVE_SHYLU_NODEHTS_BLAS";
#endif
#ifdef HAVE_SHYLU_NODEHTS_MKL
  os << " HAVE_SHYLU_NODEHTS_MKL";
#endif
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::print_options (const Options& o, std::ostream& os) {
  print_compiletime_options(os);
  os << std::endl;
  o.print(os);
  os << std::endl;
}

template<typename Int, typename Size, typename Sclr> Impl<Int, Size, Sclr>::
ConstCrsMatrix::~ConstCrsMatrix () {
  if ( ! deallocate_) return;
  assert( ! deallocator);
  deln_const(ir); deln_const(jc); deln_const(d);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::ConstCrsMatrix::deallocate () const {
  if ( ! deallocator) return;
  deallocator->counter--;
  if (deallocator->counter > 0) return;
  assert(deallocator->counter == 0);
  deallocator->free_CrsMatrix_data();
}

template<typename Int, typename Size, typename Sclr> Impl<Int, Size, Sclr>::
CrsMatrix::~CrsMatrix () {
  deln_const(ir); deln_const(jc); deln_const(d);
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::Partition::alloc_d () {
  assert( ! cm->d);
  cm->d = Allocnator<Sclr>(cm->ir[cm->m], "Partition::alloc_d").release();
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::Partition::alloc_A_idxs (const Size innz) {
  assert( ! A_idxs);
  this->nnz = innz;
  A_idxs = Allocnator<Size>(nnz, "Partition::alloc_A_idxs").release();
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::Partition::clear () {
  del(cm);
  deln(A_idxs);
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::Partition::clear_d () { deln(cm->d); }

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
partition_n_uniformly (const Int n, const Int nparts, Array<Int>& p) {
  p.optclear_and_resize(nparts + 1);
  const Int base = n / nparts;
  Int rem = n - base*nparts;
  Int extra = rem > 0 ? 1 : 0;
  p[0] = 0;
  for (Int i = 1; i <= nparts; ++i) {
    p[i] = p[i-1] + base + extra;
    if (rem > 0) {
      --rem;
      if (rem == 0) extra = 0;
    }
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
SparseData::init (const Int m, const Size nnz, const char* fail_msg,
                  const bool touch) {
  free();
  dealloc_ = true;
  try {
    ir = allocn<Size>(m+1, touch);
    if (nnz > 0) {
      jc = allocn<Int>(nnz, touch);
      d = allocn<Sclr>(nnz, touch);
    }
    ir[0] = 0;
    ir[m] = nnz;
  } catch (...) {
    free();
    std::stringstream ss;
    ss << fail_msg << ": SparseData failed to allocate. m = "
       << m << " nnz = " << nnz;
    throw hts::Exception(ss.str());
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::SparseData::free () {
  deln(ir);
  deln(jc);
  deln(d);
}

struct NumThreads {
  int omp, mkl;
  bool mkl_dynamic;
};

inline void set_num_threads (const int nthreads, NumThreads& save) {
  save.omp = omp_get_max_threads();
#ifdef HAVE_SHYLU_NODEHTS_MKL
  save.mkl = mkl_get_max_threads();
  save.mkl_dynamic = mkl_get_dynamic();
#endif
  omp_set_num_threads(nthreads);
#ifdef HAVE_SHYLU_NODEHTS_MKL
  // We never use MKL threading.
  mkl_set_dynamic(0);
  mkl_set_num_threads(1);
#endif
}

inline void restore_num_threads (const NumThreads& save) {
#ifdef HAVE_SHYLU_NODEHTS_MKL
  mkl_set_dynamic(save.mkl_dynamic);
#endif

  // This ruins performance on especially the KNC if the number of threads is
  // actually switching. Need to think about this more. Right now, the interface
  // does not promise OMP state will remain the same.
  return;
//   omp_set_num_threads(save.omp);
// #ifdef HAVE_SHYLU_NODEHTS_MKL
//   mkl_set_num_threads(save.mkl);
// #endif
}

inline bool check_nthreads (const int nt_requested, const int nt_rcvd,
                            std::string& msg) {
  if (nt_requested != nt_rcvd) {
    std::stringstream ss;
    ss << "set_num_threads: " << nt_requested << " requested, but "
       << nt_rcvd << " obtained.";
    msg = ss.str();
    return false;
  }
  return true;
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
throw_if_nthreads_not_ok (const int nthreads) {
  int nr;
#ifdef _OPENMP
# pragma omp parallel
#endif
  nr = omp_get_num_threads();
  std::string msg;
  if ( ! check_nthreads(nthreads, nr, msg))
    throw hts::Exception(msg);
}

// Return i such that ir[i] is the first value >= c. If no such i exists,
// return n.
template<typename Int>
inline Int find_first (const Int* jc, const Int n, const Int c) {
  return c == 0 ? 0 : std::lower_bound(jc, jc+n, c) - jc;
}

// Return the number of nonzeros in row r that are in [c_first, c_last). The
// corresponding indices, relative to the start of the row, are i_first:i_last.
template<typename Int, typename Size, typename Sclr>
inline Int Impl<Int, Size, Sclr>::
find_first_and_last (const Size* const ir, const Int r, const Int* const jc,
                     const Int c_first, const Int c_last,
                     Int& i_first, Int& i_last) {
  assert(c_last >= c_first);
  const Size
    iri = ir[r],
    irip1 = ir[r+1];
  const Int n = static_cast<Int>(irip1 - iri);
  i_first = find_first(jc + iri, n, c_first);
  if (i_first == n) {
    i_last = i_first;
    return 0;
  }
  i_last = i_first + find_first<Int>(jc + iri + i_first, n - i_first, c_last);
  // A return value of n - i_first is OK.
  return i_last - i_first;
}

// Crop the submatrix A(b) such that A(cb) has no 0 border.
template<typename Int, typename Size, typename Sclr>
Size Impl<Int, Size, Sclr>::
crop_matrix (const CrsMatrix& T, const Box& b, Box& cb) {
  cb.r0 = -1;
  Int r1 = -1;
  cb.c0 = b.c0 + b.nc;
  Int c1 = b.c0;
  Size nnz = 0;
  for (Int r = b.r0, lrow = 0; r < b.r0 + b.nr; ++r, ++lrow) {
    Int i_first, i_last;
    const Int cnt = find_first_and_last(T.ir, r, T.jc, b.c0, b.c0 + b.nc,
                                        i_first, i_last);
    if (cnt) {
      nnz += cnt;
      const Size irr = T.ir[r];
      cb.c0 = std::min(cb.c0, T.jc[irr + i_first]);
      c1 = std::max(c1, T.jc[irr + i_last - 1]);
      r1 = r;
      if (cb.r0 == -1) cb.r0 = r;
    }
  }
  if (cb.r0 == -1) {
    cb.r0 = b.r0;
    cb.nr = 0;
  } else cb.nr = r1 - cb.r0 + 1;
  if (cb.c0 > c1) {
    cb.c0 = b.c0;
    cb.nc = 0;
  } else cb.nc = c1 - cb.c0 + 1;
  return nnz;
}

// Decide how many level sets to keep.
template<typename Int, typename Size, typename Sclr>
Int Impl<Int, Size, Sclr>::
decide_level_set_max_index (const Array<Int>& N, const Int size_thr,
                            const Options& o) {
  Int N_end = static_cast<Int>(N.size());
  while (N_end > 0 && N[N_end-1] < size_thr) --N_end;
  Int nrows_total = 0, nrows_under = 0;
  for (Int i = 0; i < N_end; ++i) {
    const Int Ni = N[i];
    nrows_total += Ni;
    if (Ni < size_thr) nrows_under += Ni;
  }
  Int i;
  for (i = N_end - 1; i >= 0; --i) {
    const Int Ni = N[i];
    if (Ni >= size_thr && nrows_under <= o.lset_max_bad_fraction * nrows_total)
      break;
    nrows_total -= Ni;
    nrows_under -= Ni;
  }
  return i;
}

// Allocate lsets.
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
alloc_lsets (
  const Int lsmi, const Int sns, const Array<Int>& level, const Array<Int>& n,
  typename LevelSetter::LevelSets& lsets)
{
  if (lsmi < 0) return;
  const Int Lm_sr = static_cast<Int>(level.size());
  lsets.optclear_and_resize(lsmi+1);
  for (Int i = 0; i <= lsmi; ++i) {
    lsets[i].init();
    lsets[i].optclear_and_reserve(sns * n[i]);
  }
  // Load.
  for (Int i = 0; i < Lm_sr; ++i) {
    const Int ilev = level[i];
    if (ilev <= lsmi)
      for (Int j = 0; j < sns; ++j)
        lsets[ilev].unsafe_push_back(i * sns + j);
  }
}

template<typename Int, typename Size, typename Sclr>
Int Impl<Int, Size, Sclr>::
locrsrow_schedule_serial (const ConstCrsMatrix& L, const Int sns,
                          Array<Int>& w) {
  // Eq. 18 in Y. Saad's 1989 SIAM J Sci Stat Comput paper.
  Int max_level = -1;
  if (sns == 1) {
    w.optclear_and_resize(L.m, -1);
    for (Int r = 0; r < L.m; ++r) {
      Int level = -1;
      for (Size j = L.ir[r]; j < L.ir[r+1]; ++j)
        level = std::max(level, w[L.jc[j]]);
      ++level;
      w[r] = level;
      max_level = std::max(max_level, level);
    }
  } else {
    // Implement it for a blocked matrix, where the block size is sns >= 1.
    const Int Lm_sr = L.m / sns;
    w.optclear_and_resize(Lm_sr, -1);
    for (Int sr = 0, r = 0; sr < Lm_sr; ++sr) {
      Int level = -1;
      for (Int i = 0; i < sns; ++i, ++r) {
        // For each row in the block row:
        for (Size j = L.ir[r]; j < L.ir[r+1]; ++j) {
          const Int sc = L.jc[j] / sns;
          level = std::max(level, w[sc]);
        }
      }
      ++level;
      w[sr] = level;
      max_level = std::max(max_level, level);
    }
  }
  return max_level;
}

template<typename Int, typename Size, typename Sclr>
Int Impl<Int, Size, Sclr>::
locrsrow_schedule_sns1 (const ConstCrsMatrix& L, Array<Int>& w,
                        const Options& o) {
  const Int
    nthreads = omp_get_max_threads(),
    blksz = nthreads*((o.pp_min_block_size + nthreads)/nthreads),
    rows_per_thread = std::max(1, blksz / nthreads);
  if (blksz > L.m)
    return locrsrow_schedule_serial(L, 1, w);
  Array<Size> frontier(blksz);
  for (Int i = 0; i < blksz; ++i) frontier[i] = L.ir[i];
  w.optclear_and_resize(L.m);
  Int max_level = -1;
  volatile Int done = -1;
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
#ifdef _OPENMP
#   pragma omp for
#endif
    for (Int i = 0; i < L.m; ++i) w[i] = -1;
    const Size* const ir = L.ir;
    const Int* const jc = L.jc;
    for (Int c = 0; c < L.m; c += blksz) {
      const Int tlim = std::min<Int>(c + blksz, L.m);
      // On-diag serial triangle.
#ifdef _OPENMP
#     pragma omp single nowait
#endif
      {
        for (Int r = c; r < tlim; ++r) {
          // w[r] contains the max level seen so far.
          Int level = w[r];
          for (Size j = frontier[r - c], jlim = ir[r+1]; j < jlim; ++j)
            level = std::max(level, w[jc[j]]);
          ++level;
          w[r] = level;
          max_level = std::max(max_level, level);
        }
        done = c;
      }
      if (tlim == L.m) break;
      // Off-diag parallel block row.
      const Int rlim = std::min<Int>(tlim + blksz, L.m);
      while (done != c)
        ;
#ifdef _OPENMP
#     pragma omp for schedule(static, rows_per_thread)
#endif
      for (Int r = tlim; r < rlim; ++r) {
        Int level = -1;
        const Size jlim = ir[r+1];
        for (Size j = ir[r]; j < jlim; ++j) {
          const Int col = jc[j];
          if (col >= tlim) {
            frontier[r - tlim] = j;
            w[r] = level;
            break;
          }
          level = std::max(level, w[col]);
        }
      }
      // Implied barrier from parfor.
    }
  }
  return max_level;
}

template<typename Int, typename Size, typename Sclr>
Int Impl<Int, Size, Sclr>::
locrsrow_schedule (const ConstCrsMatrix& L, const Int sns, Array<Int>& w,
                   const Options& o) {
  assert(L.m > 0);
  if (sns == 1) return locrsrow_schedule_sns1(L, w, o);
  const Int
    Lm_sr = L.m / sns,
    blksz = (o.pp_min_block_size + sns) / sns,
    bnr = sns*blksz,
    nthreads = omp_get_max_threads(),
    rows_per_thread = std::max(1, (blksz + nthreads) / nthreads);
  if (blksz > Lm_sr)
    return locrsrow_schedule_serial(L, sns, w);
  Array<Size> frontier(bnr);
  for (Int i = 0; i < bnr; ++i) frontier[i] = L.ir[i];
  w.optclear_and_resize(Lm_sr);
  Int max_level = -1;
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
#ifdef _OPENMP
#   pragma omp for
#endif
    for (Int i = 0; i < Lm_sr; ++i) w[i] = -1;
    const Size* const ir = L.ir;
    const Int* const jc = L.jc;
    for (Int sc = 0; sc < Lm_sr; sc += blksz) {
      const Int stlim = std::min<Int>(sc + blksz, Lm_sr);
      // On-diag serial triangle.
#ifdef _OPENMP
#     pragma omp single nowait
#endif
      {
        const Int c = sns*sc;
        for (Int sr = sc, r = c; sr < stlim; ++sr) {
          Int level = w[sr];
          for (Int i = 0; i < sns; ++i, ++r)
            for (Size j = frontier[r - c], jlim = ir[r+1]; j < jlim; ++j) {
              const Int jsc = jc[j] / sns;
              assert(jsc >= sc);
              level = std::max(level, w[jsc]);
            }
          ++level;
          w[sr] = level;
          max_level = std::max(max_level, level);
        }
      }
      if (stlim == Lm_sr) break;
      // Off-diag parallel block row.
#ifdef _OPENMP
#     pragma omp barrier
#endif
      const Int
        srlim = std::min<Int>(stlim + blksz, Lm_sr),
        tlim = sns*stlim;
#ifdef _OPENMP
#     pragma omp for schedule(static, rows_per_thread)
#endif
      for (Int sr = stlim; sr < srlim; ++sr) {
        Int level = -1;
        for (Int i = 0, r = sns*sr; i < sns; ++i, ++r) {
          const Size jlim = ir[r+1];
          for (Size j = ir[r]; j < jlim; ++j) {
            const Int c = jc[j];
            if (c >= tlim) {
              frontier[r - tlim] = j;
              break;
            }
            const Int scol = c / sns;
            level = std::max(level, w[scol]);
          }
        }
        w[sr] = level;
      }
      // Implied barrier from parfor.
    }
  }
  return max_level;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
find_row_level_sets_Lcrs (const ConstCrsMatrix& L, const Int sns,
                          Int size_thr, typename LevelSetter::LevelSets& lsets,
                          const Options& o) {
  assert(L.m % sns == 0);

  Array<Int> w;
#ifdef __MIC__
  // || is working pretty well on MIC, but not on CPU.
  const Int max_level = locrsrow_schedule(L, sns, w, o);
#else
  const Int max_level = locrsrow_schedule_serial(L, sns, w);
#endif

  // Count level set sizes.
  Array<Int> n(max_level+1, 0);
  for (size_t i = 0; i < w.size(); ++i)
    ++n[w[i]];
  // Cutoff.
  const Int lsmi = decide_level_set_max_index(n, size_thr, o);
  // Fill lsets.
  alloc_lsets(lsmi, sns, w, n, lsets);
}

// Upper tri, CRS, col (not row) level sets. Equivalent to lower tri, CCS, row
// level sets.
template<typename Int, typename Size, typename Sclr>
Int Impl<Int, Size, Sclr>::
upcrscol_schedule_serial (const ConstCrsMatrix& U, const Int sns,
                          Array<Int>& w) {
  Int max_level = -1;
  if (sns == 1) {
    w.optclear_and_resize(U.m, -1);
    for (Int r = 0; r < U.m; ++r) {
      ++w[r];
      const Int level = w[r];
      max_level = std::max(level, max_level);
      for (Size j = U.ir[r]; j < U.ir[r+1]; ++j) {
        const Int c = U.jc[j];
        w[c] = std::max(w[c], level);
      }
    }
  } else {
    const Int Um_sr = U.m / sns;
    w.optclear_and_resize(Um_sr, -1);
    for (Int sr = 0, r = 0; sr < Um_sr; ++sr) {
      ++w[sr];
      const Int level = w[sr];
      max_level = std::max(level, max_level);
      for (Int i = 0; i < sns; ++i, ++r)
        for (Size j = U.ir[r]; j < U.ir[r+1]; ++j) {
          const Int sc = U.jc[j] / sns;
          w[sc] = std::max(w[sc], level);
        }
    }
  }
  return max_level;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
find_col_level_sets_Ucrs (const ConstCrsMatrix& U, const Int sns,
                          Int size_thr, typename LevelSetter::LevelSets& lsets,
                          const Options& o) {
  assert(U.m % sns == 0);

  Array<Int> w;
  const Int max_level = upcrscol_schedule_serial(U, sns, w);

  // Count level set sizes.
  Array<Int> n(max_level+1, 0);
  for (size_t i = 0; i < w.size(); ++i)
    ++n[w[i]];
  // Cutoff.
  const Int lsmi = decide_level_set_max_index(n, size_thr, o);
  // Fill lsets.
  alloc_lsets(lsmi, sns, w, n, lsets);
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
find_level_sets (
  const ConstCrsMatrix& T, const Int sns, const Int size_thr, const bool is_lo,
  typename LevelSetter::LevelSets& lsets, const Options& o)
{
  if (is_lo)
    find_row_level_sets_Lcrs(T, sns, size_thr, lsets, o);
  else
    find_col_level_sets_Ucrs(T, sns, size_thr, lsets, o);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
LevelSetter::init (const ConstCrsMatrix& T, const Int size_thr,
                   const bool is_lo, const Options& o) {
  lsets_.clear();
  is_lo_ = is_lo;
  // Guard against an invalid setting.
  ls_blk_sz_ = T.m % o.ls_blk_sz == 0 ? o.ls_blk_sz : 1;
  find_level_sets(T, ls_blk_sz_, size_thr, is_lo_, lsets_, o);
}

template<typename Int, typename Size, typename Sclr>
const Array<Int>& Impl<Int, Size, Sclr>::
LevelSetter::lset (const size_t i) const {
  return is_lo_ ? lsets_[i] : lsets_[lsets_.size() - i - 1];
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
LevelSetter::reverse_variable_order (Int n) {
  --n;
  for (size_t i = 0; i < lsets_.size(); ++i) {
    Array<Int>& ls = lsets_[i];
    for (size_t j = 0; j < ls.size(); ++j)
      ls[j] = n - ls[j];
    std::reverse(ls.begin(), ls.end());
  }
}

template<typename Int, typename Size, typename Sclr>
typename Impl<Int, Size, Sclr>::CrsMatrix* Impl<Int, Size, Sclr>::
get_matrix_p (const CrsMatrix& A, const Array<Int>& p,
              const bool set_diag_reciprocal) {
  const Int n = static_cast<Int>(p.size());
  Size nnz = 0;
  for (Int i = 0; i < n; ++i) {
    const Int r = p[i];
    nnz += A.ir[r+1] - A.ir[r];
  }

  SparseData sd(n, nnz, "get_matrix_p");
  for (Int i = 0; i < n; ++i) {
    const Int r = p[i];
    const Size nc = A.ir[r+1] - A.ir[r];
    sd.ir[i+1] = sd.ir[i] + nc;
    memcpy(sd.jc + sd.ir[i], A.jc + A.ir[r], nc*sizeof(*sd.jc));
    Sclr* const d_start = sd.d + sd.ir[i];
    memcpy(d_start, A.d + A.ir[r], nc*sizeof(*sd.d));
    if (set_diag_reciprocal)
      d_start[nc-1] = Sclr{1.0}/d_start[nc-1];
  }

  CrsMatrix* cm;
  try { cm = new CrsMatrix(n, A.n, sd.ir, sd.jc, sd.d); }
  catch (...) { throw hts::Exception("get_matrix_p failed to alloc."); }
  sd.release();
  return cm;
}

template<typename Int, typename Size, typename Sclr>
typename Impl<Int, Size, Sclr>::ConstCrsMatrix* Impl<Int, Size, Sclr>::
permute_to_other_tri (const ConstCrsMatrix& U) {
  const Int n = U.m;
  const Size nnz = U.ir[n];
  SparseData sd(n, nnz, "permute_to_other_tri");
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
#ifdef _OPENMP
#   pragma omp for schedule(static)
#endif
    for (Int k = 1; k <= n; ++k)
      sd.ir[k] = nnz - U.ir[n-k];
#ifdef _OPENMP
#   pragma omp for schedule(static)
#endif
    for (Size k = 0; k < nnz; ++k) {
      const Size i = nnz - k - 1;
      sd.jc[k] = n - U.jc[i] - 1;
      sd.d[k] = U.d[i];
    }
  }
  ConstCrsMatrix* ccm = 0;
  try {
    ccm = new ConstCrsMatrix(n, n, sd.ir, sd.jc, sd.d, U.dir, U.conj,
                             true, U.unitdiag, ! U.is_lo);
  } catch (...) { throw hts::Exception("permute_to_other_tri"); }
  sd.release();
  return ccm;
}

// Partition 1:n into lsis, the set of level scheduled rows, and dpis, the set
// of data-|| rows.
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
get_idxs (const Int n, const LevelSetter& lstr, Array<Int>& lsis,
          Array<Int>& dpis) {
  std::size_t lsis_sz = 0;
  for (Int i = 0; i < lstr.size(); ++i) {
    const Array<Int>& lset = lstr.lset(i);
    lsis_sz += lset.size();
  }
  lsis.optclear_and_resize(lsis_sz);

  Array<char> dpisb(n, 1);
  for (Int i = 0, il = 0; i < lstr.size(); ++i) {
    const Array<Int>& lset = lstr.lset(i);
    for (std::size_t j = 0; j < lset.size(); ++j) {
      const Int lj = lset[j];
      lsis[il++] = lj;
      dpisb[lj] = 0;
    }
  }

  dpis.optclear_and_resize(n - lsis.size());
  for (size_t i = 0, dk = 0; i < dpisb.size(); ++i)
    if (dpisb[i]) dpis[dk++] = i;
}

template<typename Int, typename Size, typename Sclr>
typename Impl<Int, Size, Sclr>::Shape Impl<Int, Size, Sclr>::
determine_shape (const ConstCrsMatrix& A) {
  int red_is_lower = 0, red_is_upper = 0, red_has_full_diag = 0,
    red_has_no_diag = 0, red_nthreads = 0;
#ifdef _OPENMP
# pragma omp parallel \
         reduction(+: red_is_lower, red_is_upper, red_has_full_diag, \
                      red_has_no_diag, red_nthreads)
#endif
  {
    bool tid_used = false, has_full_diag = true, has_no_diag = true,
      is_lower = false, is_upper = false;
#ifdef _OPENMP
#   pragma omp for schedule(static, parfor_static_size)
#endif
    for (Int r = 0; r < A.m; ++r) {
      tid_used = true;
      bool diag_fnd = false;
      for (Size j = A.ir[r]; j < A.ir[r+1]; ++j) {
        const Int c = A.jc[j];
        if (c != r) {
          if (c < r) is_lower = true;
          else is_upper = true;
        } else
          diag_fnd = true;
      }
      if ( ! diag_fnd)
        has_full_diag = false;
      else
        has_no_diag = false;
    }
    if (tid_used) {
      ++red_nthreads;
      if (has_full_diag) ++red_has_full_diag;
      if (has_no_diag) ++red_has_no_diag;
      if (is_upper) ++red_is_upper;
      if (is_lower) ++red_is_lower;
    }
  }
  const bool
    // Each thread saw a triangle having the same orientation.
    is_tri = ! (red_is_lower > 0 && red_is_upper > 0),
    // Every thread saw a full diagonal.
    has_full_diag = red_has_full_diag == red_nthreads,
    // No thread saw any diagonal.
    has_no_diag = red_has_no_diag == red_nthreads,
    // Valid only if is_tri.
    is_lower = red_is_upper == 0;
  // If ! tri_determined, then T must be a diag matrix. Can treat as lower,
  // which is is_lower's default value.
  return Shape(is_lower, is_tri, has_full_diag, has_no_diag);
}

template<typename Int, typename Size>
Int partition_ir (const Int n, const Size* const ir, const Int nparts,
                  Array<Int>& start) {
  const Int nseg = std::min<Int>(nparts, n);
  Int i0 = 1, j0 = 1;
  start.optclear_and_resize(nseg + 1);
  start[0] = 0;
  const Size nnz = ir[n];
  for (Int i = i0; i < nseg; ++i) {
    const double d = ((double) i / nseg)*nnz;
    Int j = static_cast<Int>(std::upper_bound(ir + j0, ir + n, d) - ir);
    if (d - ir[j-1] > ir[j] - d) ++j;
    j = std::min<Int>(j, n - nseg + i);
    if (j < j0) {
      // If not all the threads will have work, let the earlier threads have
      // the work.
      j = std::min<Int>(j0, n);
      assert(j == start[i-1] + 1);
    }
    start[i] = j;
    j0 = j + 1;
  }
  start[nseg] = n;

  return nseg;
}

template <typename T> inline T& conjugate (T& v) {
#ifdef HAVE_SHYLU_NODEHTS_KOKKOSKERNELS
  v = Kokkos::ArithTraits<T>::conj(v);
#else
  v = T(v.real(), -v.imag());
#endif
  return v;
}
template <> inline double& conjugate (double& v) { return v; }
template <> inline float& conjugate (float& v) { return v; }

/* The following routines extract the three types of matrices (level set, big
 * MVP, data parallel) in the two cases of L and U.
 *   The extracted matrices' rows do not need to be sorted. In the level set
 * matrix, the whole row is extracted and not further subdivided. The same is
 * true of the big scatter block. The data parallel block's permutation is
 * already sorted; therefore, even though it is broken up subsequently, it's
 * already sorted.
 *   However, sorting tends to improve the solver's performance, probably
 * because it increases data locality when accessing x.
 */

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
get_matrix_common_with_covers_all (
  const ConstCrsMatrix& A, const PermVec& pv, const PermVec& qv, Partition& p,
  const bool get_A_idxs, const bool pp)
{
  // Count nnz.
  Size nnz = 0;
#ifdef _OPENMP
# pragma omp parallel for reduction(+:nnz)
#endif
  for (size_t i = 0; i < pv.size(); ++i) {
    const Int r = pv.get(i);
    nnz += A.ir[r+1] - A.ir[r];
  }
  if (A.unitdiag) nnz += pv.size();

  SparseData sd(pv.size(), nnz, "get_matrix_pp_with_covers_all");
  Int extra = A.unitdiag ? 1 : 0;
  for (size_t ipv = 0; ipv < pv.size(); ++ipv) {
    const Int i = pv.get(ipv);
    const Size nc = A.ir[i+1] - A.ir[i] + extra;
    sd.ir[ipv+1] = sd.ir[ipv] + nc;
  }
  assert(sd.ir[pv.size()] == nnz);
  if (get_A_idxs) p.alloc_A_idxs(nnz);

  const int nthreads = omp_get_max_threads();
  Array<Int> start;
  const Int nseg = partition_ir<Int, Size>(pv.size(), sd.ir, nthreads, start);
#ifdef _OPENMP
# pragma omp parallel
#endif
  do {
    const int tid = omp_get_thread_num();
    if (tid >= nseg) break;
    for (Int ipv = start[tid], ipvlim = start[tid+1]; ipv < ipvlim; ++ipv) {
      const Int i = pv.get(ipv);
      Size Anc = A.ir[i+1] - A.ir[i], Aj0 = A.ir[i];
      Size Bj0 = sd.ir[ipv];
      if (pp) {
        for (Size j = 0; j < Anc; ++j) {
          const Size Aj = Aj0 + j, Bj = Bj0 + j;
          sd.jc[Bj] = pv.to_block(A.jc[Aj]);
          sd.d[Bj] = A.d[Aj];
          if (A.conj) conjugate(sd.d[Bj]);
          if (get_A_idxs) p.A_idxs[Bj] = Aj;
        }
      } else {
        for (Size j = 0; j < Anc; ++j) {
          const Size
            Aj = Aj0 + j, Bj = Bj0 + j,
            Ac = A.jc[Aj];
          sd.jc[Bj] = qv.has(Ac) ? qv.to_block(Ac) : qv.size() + pv.to_block(Ac);
          sd.d[Bj] = A.d[Aj];
          if (A.conj) conjugate(sd.d[Bj]);
          if (get_A_idxs) p.A_idxs[Bj] = Aj;
        }
      }
      assert(A.is_lo);
      if (A.unitdiag) {
        const Size Bj = Bj0 + Anc;
        sd.jc[Bj] = ipv + (pp ? 0 : qv.size());
        sd.d[Bj] = static_cast<Sclr>(1);
        p.A_invalidate(Bj);
      }
    }
  } while (0);

  try { p.cm = new CrsMatrix(pv.size(), A.n, sd.ir, sd.jc, sd.d); }
  catch (...)
  { throw hts::Exception("get_matrix_pp_with_covers_all failed to alloc."); }
  sd.release();

  sort(p, start);
}

// Extract A(p,p) given that A(p,~p) is empty.
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
get_matrix_pp_with_covers_all (const ConstCrsMatrix& A, const PermVec& pv, Partition& p,
                               const bool get_A_idxs) {
  get_matrix_common_with_covers_all(A, pv, pv, p, get_A_idxs, true);
}

// Extract B = A(p,s), s = [q p], given that A(p,~s) is empty.
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
get_matrix_p_qp_with_covers_all (const ConstCrsMatrix& A, const PermVec& pv,
                                 const PermVec& qv, Partition& p,
                                 const bool get_A_idxs) {
  get_matrix_common_with_covers_all(A, pv, qv, p, get_A_idxs, false);
}

template<typename Int, typename Size, typename T> struct SortEntry {
  Size A_idx;
  Int i;
  T d;
  SortEntry () {}
  bool operator< (const SortEntry& se) const { return i < se.i; }
};

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::sort (
  Partition& p,
  // Use the nonzero partitioning the caller already obtained to ||ize sort's
  // loop. The O(nc log nc) complexity here relative to the O(nc) complexity in
  // the caller probably doesn't invalidate the partitioning.
  const Array<Int>& start)
{
  CrsMatrix& A = *p.cm;
  assert(A.m == start.back());
  bool ok = true;
#ifdef _OPENMP
# pragma omp parallel
#endif
  do {
    const int tid = omp_get_thread_num();
    if (tid + 1 >= static_cast<int>(start.size())) break;
    Size max_nc = 0;
    for (Int r = start[tid], rlim = start[tid+1]; r < rlim; ++r)
      max_nc = std::max(max_nc, A.ir[r+1] - A.ir[r]);
    Array<SortEntry<Int, Size, Sclr> > ses;
    try { ses.optclear_and_resize(max_nc); }
    catch (...) {
      ok = false;
      break;
    }
    for (Int r = start[tid], rlim = start[tid+1]; r < rlim; ++r) {
      const Size irr = A.ir[r], irrp1 = A.ir[r+1], nc = irrp1 - irr;
      for (Size j = 0; j < nc; ++j) {
        const Size Aj = irr + j;
        ses[j].i = A.jc[Aj];
        ses[j].d = A.d[Aj];
        if (p.A_idxs) ses[j].A_idx = p.A_idxs[Aj];
      }
      std::sort(ses.begin(), ses.begin() + nc);
      for (Size j = 0; j < nc; ++j) {
        const Size Aj = irr + j;
        A.jc[Aj] = ses[j].i;
        A.d[Aj] = ses[j].d;
        if (p.A_idxs) p.A_idxs[Aj] = ses[j].A_idx;
      }
    }
  } while (0);
  if ( ! ok) throw hts::Exception("sort: Could not allocate ses.");
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
reverse_A_idxs (const Size nnz, Partition& p) {
#ifdef _OPENMP
# pragma omp for schedule(static)
#endif
  for (Size i = 0; i < p.nnz; ++i)
    if (p.A_valid(i))
      p.A_idxs[i] = nnz - p.A_idxs[i] - 1;
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
copy_partition (const ConstCrsMatrix& A, Partition& p) {
  const Size ilim = p.nnz;
  if (A.unitdiag) {
    if (A.conj) {
#ifdef _OPENMP
#     pragma omp parallel for
#endif
      for (Size i = 0; i < ilim; ++i) {
        if (p.A_valid(i)) {
          p.cm->d[i] = A.d[p.A_idxs[i]];
          conjugate(p.cm->d[i]);
        } else {
          p.cm->d[i] = static_cast<Sclr>(1);
        }
      }
    } else {
#ifdef _OPENMP
#     pragma omp parallel for
#endif
      for (Size i = 0; i < ilim; ++i)
        p.cm->d[i] = p.A_valid(i) ? A.d[p.A_idxs[i]] : static_cast<Sclr>(1);
    }
  } else {
    if (A.conj) {
#ifdef _OPENMP
#     pragma omp parallel for
#endif
      for (Size i = 0; i < ilim; ++i) {
        p.cm->d[i] = A.d[p.A_idxs[i]];
        conjugate(p.cm->d[i]);
      }
    } else {
#ifdef _OPENMP
#     pragma omp parallel for
#endif
      for (Size i = 0; i < ilim; ++i)
        p.cm->d[i] = A.d[p.A_idxs[i]];
    }
  }
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
repartition_into_2_blocks (Partition* const p, const ConstCrsMatrix& A) {
  for (int i = 0; i < 2; ++i)
    if (p[i].cm)
      copy_partition(A, p[i]);
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
CrsSegmenter::count_nnz_by_row_loop (const Int i, Array<Int>& rcnt) {
  Int i_first, i_last;
  rcnt[i] = find_first_and_last(A_.ir, r0_ + i, A_.jc, c0_, c0_ + nc_,
                                i_first, i_last);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
CrsSegmenter::count_nnz_by_row (Array<Int>& rcnt) {
  rcnt.optclear_and_resize(nr_);
  // Don't allow nested ||ism.
  if (omp_get_num_threads() == 1) {
#ifdef _OPENMP
#   pragma omp parallel for schedule(guided)
#endif
    for (Int i = 0; i < nr_; ++i)
      count_nnz_by_row_loop(i, rcnt);
  } else
    for (Int i = 0; i < nr_; ++i)
      count_nnz_by_row_loop(i, rcnt);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
CrsSegmenter::init_nnz (const Array<Int>& rcnt) {
  const Array<Int>& p = this->p_;
  const Int nseg = p.size() - 1;
  this->nnz_.optclear_and_resize(nthreads_, 0);
  for (Int i = 0; i < nseg; ++i)
    for (Int j = p[i] - p[0]; j < p[i+1] - p[0]; ++j)
      this->nnz_[i] += rcnt[j];
}

// p partitions rows. It may not be a perfect partitioning in the sense of
// optimal given that each thread must handle a mutually exclusive set of
// rows. Attempt to remove spikes in #rows.
template<typename Int, typename Size>
inline void
smooth_spikes (const Array<Int>& rcnt, Array<Int>& p, Array<Size>& nnz,
               const bool ignore_0) {
  const Int n = static_cast<Int>(p.size()) - 1;
  if (n <= 1) return;
  Int first = ignore_0 ? 1 : 0;
  for (int it = 0; ; ++it) {
    // Find a spike.
    Size s_val = nnz[first];
    Int s_i = 0;
    for (Int i = first + 1; i < n; ++i)
      if (nnz[i] > s_val) {
        s_i = i;
        s_val = nnz[i];
      }

    // See if it decreases the spikiness to move a row from one set to another.
    Size rnnz, s;
    Int d, idx;
    if (s_i == 0) {
      idx = p[1] - p[0] - 1;
      rnnz = rcnt[idx];
      d = 1;
      s = std::max(nnz[1] + rnnz, nnz[0] - rnnz);
      if (s >= s_val) break;
      --p[1];
    } else if (s_i == n-1) {
      idx = p[s_i] - p[0];
      rnnz = rcnt[idx];
      d = -1;
      s = std::max(nnz[s_i-1] + rnnz, nnz[s_i] - rnnz);
      if (s >= s_val) break;
      ++p[s_i];
    } else {
      const Int
        idx_l = p[s_i] - p[0],
        idx_r = p[s_i+1] - 1 - p[0];
      const Size
        s_l = std::max(nnz[s_i-1] + rcnt[idx_l], nnz[s_i] - rcnt[idx_l]),
        s_r = std::max(nnz[s_i+1] + rcnt[idx_r], nnz[s_i] - rcnt[idx_r]);
      s = std::min(s_l, s_r);
      if (s >= s_val) break;
      if (s_l < s_r) {
        idx = idx_l;
        d = -1;
        ++p[s_i];
      } else {
        idx = idx_r;
        d = 1;
        --p[s_i+1];
      }
      rnnz = rcnt[idx];
    }

    nnz[s_i+d] += rnnz;
    nnz[s_i] -= rnnz;
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::CrsSegmenter::segment () {
  assert(nr_ % ls_blk_sz_ == 0);
  Array<Int>& p = this->p_;

  const Int nseg = std::min<Int>(nthreads_, nr_ / ls_blk_sz_);
  if (nseg == 0) {
    assert(nr_ == 0);
    p.optclear_and_resize(1, r0_);
    return;
  }

  Array<Int> rcnt;
  count_nnz_by_row(rcnt);

  // cumsum(rcnt).
  Array<Int> cs_rcnt(rcnt.size() / ls_blk_sz_);
  if (nseg > 0) {
    cs_rcnt[0] = 0;
    for (int j = 0; j < ls_blk_sz_; ++j) cs_rcnt[0] += rcnt[j];
    for (size_t i = 1; i < cs_rcnt.size(); ++i) {
      cs_rcnt[i] = cs_rcnt[i-1];
      for (int j = 0; j < ls_blk_sz_; ++j)
        cs_rcnt[i] += rcnt[i * ls_blk_sz_ + j];
    }
  }

  bool tid_empty = block_0_nnz_os_ && nseg > 1 &&
    block_0_nnz_os_ >= (Real) (cs_rcnt.back() + block_0_nnz_os_) / nseg;
  if ( ! tid_empty) {
    for (size_t i = 0; i < cs_rcnt.size(); ++i)
      cs_rcnt[i] += block_0_nnz_os_;
    rcnt[0] += block_0_nnz_os_;
  }

  p.optclear_and_resize(nseg + 1, r0_);
  if (nseg > 0) {
    if (cs_rcnt.back() == 0) {
      for (Int i = 1; i < nseg; ++i)
        p[i] = r0_ + i;
    } else {
      Int i0 = 1, j0 = 1, nparts = nseg, start = 0;
      if (tid_empty) {
        // Adjust so that the first thread gets nothing and the nonzeros are
        // balanced among the rest.
        p[1] = r0_;
        ++i0;
        --nparts;
        ++start;
      }
      for (Int i = i0; i < nseg; ++i) {
        const Real d = ((Real) (i - start) / nparts)*cs_rcnt.back();
        Int j = (std::upper_bound(cs_rcnt.begin() + j0, cs_rcnt.end(), d)
                 - cs_rcnt.begin());
        if (d - cs_rcnt[j-1] > cs_rcnt[j] - d) ++j;
        j = std::min<Int>(j, cs_rcnt.size() - nseg + i);
        if (j < j0) {
          // If not all the threads will have work, let the earlier threads have
          // the work.
          j = std::min<Int>(j0, cs_rcnt.size());
          assert(r0_ + j*ls_blk_sz_ == p[i-1] + 1);
        }
        p[i] = r0_ + j*ls_blk_sz_;
        j0 = j + 1;
      }
    }
  }
  p[nseg] = r0_ + nr_;

  init_nnz(rcnt);
  if (ls_blk_sz_ == 1) smooth_spikes(rcnt, p, this->nnz_, tid_empty);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TMatrix::clear () {
  bs_.clear();
  ros_.clear();
  is_empty_ = true;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TMatrix::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
               const InitInfo& in, const Int block_0_nnz_os,
               const int tid_offset) {
  init_metadata(A, r0, c0, nr, nc, in, block_0_nnz_os, tid_offset);
  init_memory(in);
#ifdef _OPENMP
# pragma omp parallel
#endif
  init_numeric(A, omp_get_thread_num());
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TMatrix::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
               const InitInfo& in, const CrsSegmenter& seg) {
  init_metadata(A, r0, c0, nr, nc, in, seg);
  init_memory(in);
#ifdef _OPENMP
# pragma omp parallel
#endif
  init_numeric(A, omp_get_thread_num());
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TMatrix::init_metadata (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                        const InitInfo& in, const Int block_0_nnz_os,
                        const int tid_offset) {
  clear();
  nr_ = nr;
  tid_os_ = tid_offset;

  Box b;
  const Int nnz = crop_matrix(A, Box(r0, c0, nr, nc), b);
  const Int roff = b.r0 - r0;
  coff_ = b.c0 - c0;
  r0 = b.r0; c0 = b.c0; nr = b.nr; nc = b.nc;

  // p2p sync'ing favors fewer threads. So don't run a bunch of threads on small
  // blocks.
  const Int nt = nnz / square(in.min_parallel_rows) + 1;
  const Int nseg = std::max<Int>(1, std::min(in.nthreads, nt));
  CrsSegmenter seg(A, r0, c0, nr, nc, nseg, 1, block_0_nnz_os);

  if (nnz)
    init_metadata_with_seg(A, r0, c0, nr, nc, roff, in, seg);
  else
    is_parallel_ = false;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TMatrix::init_metadata (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                        const InitInfo& in, const CrsSegmenter& seg) {
  clear();
  nr_ = nr;
  init_metadata_with_seg(A, r0, c0, nr, nc, 0, in, seg);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::TMatrix::
init_metadata_with_seg (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                        const Int roff, const InitInfo& in,
                        const CrsSegmenter& seg) {
  // Serial if block is too small.
  const Int nseg = seg.get_p().size() - 1;
  is_parallel_ = nseg > 1 || tid_os_ > 0;
  if (nseg == 1) {
    bs_.optclear_and_reserve(1);
    bs_.unsafe_push_back(SerialBlock());
    bs_.back().init_metadata(A, r0, c0, nr, nc, in);
    ros_.optclear_and_resize(1, roff);
  } else {
    const Array<Int>& p = seg.get_p();
    const Int n = static_cast<Int>(p.size()) - 1;
    ros_.optclear_and_resize(n, 0);
    bs_.optclear_and_resize(n, SerialBlock());
    for (Int id = 0; id < n; ++id) {
      ros_[id] = p[id] - r0 + roff;
      const Int nri = p[id+1] - p[id];
      bs_[id].init_metadata(A, p[id], c0, nri, nc, in);
    }
  }

  is_empty_ = true;
  for (size_t i = 0; i < bs_.size(); ++i)
    if (bs_[i].get_nnz() != 0) {
      is_empty_ = false;
      break;
    }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TMatrix::init_memory (const InitInfo& in) {
#ifdef _OPENMP
# pragma omp parallel
#endif
  { const int tid = omp_get_thread_num(), id = tid - tid_os_;
    if (id >= 0 && id < (int) bs_.size())
      bs_[id].init_memory(in);
  }
}

// possibly inside || {}
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TMatrix::init_numeric (const CrsMatrix& A, const int tid) {
  const int id = tid - tid_os_;
  if (id >= 0 && id < (int) bs_.size())
    bs_[id].init_numeric(A);
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
TMatrix::reinit_numeric (const CrsMatrix& A, const int tid) {
  const int id = tid - tid_os_;
  if (id >= 0 && id < (int) bs_.size())
    bs_[id].reinit_numeric(A);
}

template<typename Int, typename Size, typename Sclr>
inline Int Impl<Int, Size, Sclr>::
TMatrix::block_r0 (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return bs_[id].get_r0();
}

template<typename Int, typename Size, typename Sclr>
inline Int Impl<Int, Size, Sclr>::
TMatrix::block_nr (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return bs_[id].get_nr();
}

template<typename Int, typename Size, typename Sclr>
inline const typename Impl<Int, Size, Sclr>::SerialBlock*
Impl<Int, Size, Sclr>::TMatrix::block (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return &bs_[id];
}

template<typename Int, typename Size, typename Sclr>
inline typename Impl<Int, Size, Sclr>::SerialBlock*
Impl<Int, Size, Sclr>::TMatrix::block (const int tid) {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return &bs_[id];
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
SerialBlock::clear () {
  if (deallocate_) {
    deln(ir_);
    deln(jc_);
    deln(d_);
  }
  ir_ = 0;
  jc_ = 0;
  roff_ = coff_ = 0;
  d_ = 0;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
SerialBlock::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                   const InitInfo& in) {
  init_metadata(A, r0, c0, nr, nc, in);
  init_memory(in);
  init_numeric(A);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
SerialBlock::init_metadata (const CrsMatrix& A, Int r0, Int c0, Int nr,
                            Int nc, const InitInfo& in) {
  clear();

  Box b;
  nnz_ = crop_matrix(A, Box(r0, c0, nr, nc), b);
  roff_ = b.r0 - r0;
  r0 = b.r0; nr = b.nr;
  coff_ = b.c0 - c0;
  c0 = b.c0; nc = b.nc;

  r0_ = r0; c0_ = c0; nr_ = nr; nc_ = nc;
  if (nr_ == 0 || nc_ == 0) return;
  assert(nnz_ > 0);
  assert(nnz_ <= (1.0*nr_)*nc_);
  is_dense_ = nnz_ >= (in.min_dense_density*nr_)*nc_;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
SerialBlock::init_memory (const InitInfo& in) {
  if (nr_ == 0 || nc_ == 0) return;
  if (is_dense_) {
    // First touch is measurably required.
    d_ = Allocnator<Sclr>(nr_*nc_, "SerialBlock::init dense d", true).release();
  } else {
    Allocnator<Size> air(nr_+1, "SerialBlock::init sparse ir", true);
    if (nnz_ == 0) {
      ir_ = air.release();
      ir_[nr_] = 0;
      return;
    }
    Allocnator<Int> ajc(nnz_, "SerialBlock::init sparse jc", true);
    Allocnator<Sclr> ad(nnz_, "SerialBlock::init sparse d", true);
    ir_ = air.release();
    jc_ = ajc.release();
    d_ = ad.release();
  }
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
SerialBlock::init_numeric (const CrsMatrix& A) {
  if (is_dense_)
  {
    for (Int i = 0; i < nr_ * nc_; i++)
      d_[i] = Sclr(0);
  }
  reinit_numeric(A);
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
SerialBlock::reinit_numeric (const CrsMatrix& A) {
  if (ir_) reinit_numeric_spars(A);
  else reinit_numeric_dense(A);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
SerialBlock::reinit_numeric_dense (const CrsMatrix& A) {
  const Int ilim = r0_ + nr_;
  for (Int i = r0_; i < ilim; ++i) {
    const Int lrow = i - r0_;
    const Size
      iri = A.ir[i],
      irip1 = A.ir[i+1];
    for (Size j = iri + find_first<Int>(A.jc + iri, irip1 - iri, c0_);
         j < irip1; ++j) {
      const Int lcol = A.jc[j] - c0_;
      if (lcol >= nc_) break;
      const Size k = nc_*lrow + lcol;
      d_[k] = A.d[j];
    }
  }
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
SerialBlock::reinit_numeric_spars (const CrsMatrix& A) {
  Size k = 0;
  ir_[0] = 0;
  const Int ilim = r0_ + nr_;
  for (Int i = r0_; i < ilim; ++i) {
    const Size iri = A.ir[i], irip1 = A.ir[i+1];
    for (Size g = iri + find_first<Int>(A.jc + iri, irip1 - iri, c0_);
         g < irip1; ++g) {
      const Int jcg = A.jc[g];
      if (jcg >= c0_ + nc_) break;
      jc_[k] = jcg - c0_;
      d_[k] = A.d[g];
      ++k;
    }
    ir_[i - r0_ + 1] = k;
  }
  assert(k == ir_[nr_]);
}

template<typename Int> inline Int ntri (const int n) { return (n*(n + 1))/2; }

template<typename Int, typename Size, typename Sclr>
inline Int Impl<Int, Size, Sclr>::
count_nnz_lotri (const CrsMatrix& T, const Int r0, const Int c0,
                 const Int n) {
  Int nnz = 0;
  for (Int r = r0, lrow = 0; r < r0 + n; ++r, ++lrow) {
    Int i_first, i_last;
    nnz += find_first_and_last(T.ir, r, T.jc, c0, c0 + lrow + 1,
                               i_first, i_last);
  }
  return nnz;
}

template<typename Int, typename Size, typename Sclr>
Impl<Int, Size, Sclr>::
OnDiagTri::OnDiagTri ()
  : c0_(0), nnz_(0), d_(0), m_(0), dense_(true)
{}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::clear () {
  deln(d_);
  del(m_);
  t_.clear();
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::init (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                 const InitInfo& in) {
  init_metadata(T, r0, c0, n, in);
  init_memory(in);
  init_numeric(T);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::init (const Int r0, const Int c0, const Int n) {
  clear();
  this->n_ = n; this->r0_ = r0; c0_ = c0;
}

template<typename Int, typename Size, typename Sclr>
inline bool Impl<Int, Size, Sclr>::
is_dense_tri (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
              const InitInfo& in, Size* innz) {
  const Size nnz = count_nnz_lotri(T, r0, c0, n);
  if (innz) *innz = nnz;
  return nnz >= 0.5*in.min_dense_density*ntri<Int>(n);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::init_metadata (const CrsMatrix& T, const Int r0, const Int c0,
                          const Int n, const InitInfo& in) {
  clear();
  this->n_ = n; this->r0_ = r0; c0_ = c0;
  if ( ! this->n_) {
    nnz_ = 0;
    return;
  }
  dense_ = is_dense_tri(T, r0, c0, n, in, &nnz_);
  if (dense_) inv_init_metadata(in);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::init_memory (const InitInfo& in) {
  if ( ! nnz_) return;
  if (dense_) {
    Allocnator<Sclr> ad(ntri<Int>(this->n_), "OnDiagTri::init dense", true);
    if ( ! t_.empty()) inv_init_memory();
    d_ = ad.release();
  } else {
    SparseData sd(this->n_, nnz_, "OnDiagTri::init_memory");
    try { m_ = new CrsMatrix(this->n_, this->n_, sd.ir, sd.jc, sd.d); }
    catch (...)
    { throw hts::Exception("OnDiagTri::init_memory failed to alloc."); }
    sd.release();
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::init_numeric (const CrsMatrix& T, const bool invert) {
  if (dense_) {
    reinit_numeric(T, invert);
  } else {
    Size* const ir = m_->ir;
    Int* const jc = m_->jc;
    Sclr* const d = m_->d;
    for (Int grow = this->r0_; grow < this->r0_ + this->n_; ++grow) {
      const Int lrow = grow - this->r0_;
      const Size
        irg = T.ir[grow],
        irgp1 = T.ir[grow+1];
      ir[lrow+1] = ir[lrow];
      for (Size k = irg + find_first<Int>(T.jc + irg, irgp1 - irg, c0_);
           k < irgp1; ++k) {
        const Int lcol = T.jc[k] - c0_;
        if (lcol >= this->n_) break;
        Size& i = ir[lrow+1];
        jc[i] = lcol;
        d[i] = lrow == lcol ? Sclr{1.0}/T.d[k] : T.d[k];
        ++i;
      }
    }
    assert(ir[this->n_] == nnz_);
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::reinit_numeric (const CrsMatrix& T, const bool invert) {
  if (d_) {
    {
      Int nd = ntri<Int>(this->n_);
      for (Int i = 0; i < nd; i++)
        d_[i] = Sclr(0);
    }
    Size nnz = 0;
    for (Int grow = this->r0_; grow < this->r0_ + this->n_; ++grow) {
      const Int lrow = grow - this->r0_;
      const Size
        irg = T.ir[grow],
        irgp1 = T.ir[grow+1];
      for (Size k = irg + find_first<Int>(T.jc + irg, irgp1 - irg, c0_);
           k < irgp1; ++k) {
        const Int lcol = T.jc[k] - c0_;
        if (lcol >= this->n_) break;
        const Size di = (lrow*(lrow + 1))/2 + lcol;
        const Sclr dv = lrow == lcol ? Sclr{1.0}/T.d[k] : T.d[k];
        // Compressed dense triangle.
        d_[di] = dv;
        ++nnz;
      }
    }
    if (invert && ! t_.empty())
      inv_reinit_numeric();
  } else if (m_) {
    Sclr* const d = m_->d;
    Int i = 0;
    for (Int grow = this->r0_; grow < this->r0_ + this->n_; ++grow) {
      const Int lrow = grow - this->r0_;
      const Size
        irg = T.ir[grow],
        irgp1 = T.ir[grow+1];
      for (Size k = irg + find_first<Int>(T.jc + irg, irgp1 - irg, c0_);
           k < irgp1; ++k) {
        const Int lcol = T.jc[k] - c0_;
        if (lcol >= this->n_) break;
        d[i] = lrow == lcol ? Sclr{1.0}/T.d[k] : T.d[k];
        ++i;
      }
    }
  }
}

template<typename Int, typename Size, typename Sclr>
inline Int Impl<Int, Size, Sclr>::
OnDiagTri::nthreads () const {
  return std::max<Int>(1, static_cast<Int>(t_.size()));
}

template<typename Int, typename Size, typename Sclr>
inline Int Impl<Int, Size, Sclr>::
OnDiagTri::block_row_start (const int tid) const {
  return t_.empty() ? 0 : t_[tid].r0;
}

template<typename Int, typename Size, typename Sclr>
inline Int Impl<Int, Size, Sclr>::
OnDiagTri::block_nr (const int tid) const {
  return t_.empty() ? this->n_ : t_[tid].nr;
}

template<typename Int, typename Size, typename Sclr>
Impl<Int, Size, Sclr>::
OnDiagTri::Thread::~Thread () { deln(d); }

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::inv_init_metadata (const InitInfo& in) {
  const Int
    nnz = ntri<Int>(this->n_),
    nt = std::min((int) in.nthreads, 1 + nnz / square(in.min_parallel_rows));
  if (nt <= 1) return;

  t_.optclear_and_resize(nt);
  const Real nnz_per_thread = static_cast<Real>(nnz) / nt;
  Int r0 = 0;
  for (Int tid = 0; tid < nt; ++tid) {
    t_[tid].r0 = r0;
    t_[tid].d = 0;
    const Int n_max = this->n_ - r0;
    if (tid+1 == nt)
      t_[tid].nr = n_max;
    else {
      // Solve for n in
      //   ((r0 + n) (r0 + n + 1))/2 - (r0 (r0 + 1))/2 = nnz_per_thread.
      const Int b = 1 + 2*r0;
      t_[tid].nr = std::min<Real>(
        static_cast<Real>(n_max),
        round(0.5*(std::sqrt(b*b + 8*nnz_per_thread) - b)));
    }
    r0 += t_[tid].nr;
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
OnDiagTri::inv_init_memory () {
  bool ok = true;
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
    const int tid = omp_get_thread_num();
    if (tid < nthreads()) {
      Thread& t = t_[tid];
      const Int nnz = ntri<Int>(t.r0 + t.nr) - ntri<Int>(t.r0);
      try {
        t.d = Allocnator<Sclr>(nnz, "OnDiagTri::inv_init_memory", true).release();
      } catch (...) {
        ok = false;
      }
    }
  }
  if ( ! ok) {
    for (std::size_t i = 0; i < t_.size(); ++i)
      deln(t_[i].d);
    throw hts::Exception("OnDiagTri::inv_init_memory failed to allocate.");
  }
}

// T is in row-major compressed dense tri format. The diag of T is the
// reciprocal.
template<typename Int, typename Sclr>
inline void invert (Sclr* T, const Int n, Sclr* w) {
  for (Int c = 0; c < n; ++c) {
    // Solve for column c. That involves only the (n-c)x(n-c) lower-right
    // subblock of T. Store column in w.
    w[0] = T[0];
    Sclr* Tp = T + c + 1;
    for (Int r = 1; r < n - c; ++r) {
      w[r] = 0;
      for (Int k = 0; k < r; ++k)
        w[r] -= w[k]*Tp[k];
      w[r] *= Tp[r];
      Tp += r + c + 1;
    }
    // Copy w to column c.
    Tp = T;
    for (Int r = 0; r < n - c; ++r) {
      *Tp = w[r];
      Tp += r + c + 1;
    }
    T += c + 2;
  }
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::OnDiagTri::inv_reinit_numeric () {
  {
    Array<Sclr> w(this->n_);
    invert(d_, this->n_, w.data());
  }
  inv_copy();
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::OnDiagTri::inv_copy () {
  for (int tid = 0; tid < nthreads(); ++tid) {
    Thread& t = t_[tid];
    const Int nr0 = ntri<Int>(t.r0);
    memcpy(t.d, d_ + nr0, (ntri<Int>(t.r0 + t.nr) - nr0)*sizeof(*t.d));
  }
}

template<typename Int>
inline Int split (const Int n, const Int nthreads) {
  const Int s = n/2;
  if (n <= nthreads) return s;
  const Int smod = s % nthreads;
  if ( ! smod) return s;
  return s + (nthreads - smod);
}

template<typename Int>
inline bool empty (const Int* const col, const Int ncol,
                   // Empty inside [rs, rei]?
                   const Int rs, const Int rei) {
  for (Int i = ncol - 1; i >= 0; --i) {
    const Int c = col[i];
    if (c < rs) break;
    if (c > rei) continue;
    return false;
  }
  return true;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
find_split_rows (const CrsMatrix& T, const Int r0, const Int c0,
                 const Int n, const InitInfo& in,
                 std::vector<Int>& split_rows) {
  const Int sz = in.min_blksz / 2;
  bool ok = true;
#ifdef _OPENMP
# pragma omp parallel for schedule(static, parfor_static_size)
#endif
  for (Int r = r0 + 1; r < r0 + n; ++r) {
    bool found = true;
    for (Int r1 = r, r1_lim = std::min<Int>(T.m, r + sz); r1 < r1_lim; ++r1)
      if ( ! empty<Int>(T.jc + T.ir[r1], T.ir[r1+1] - T.ir[r1],
                        std::max<Int>(c0, c0 + r - sz), c0 + r - 1)) {
        found = false;
        break;
      }
    if (found) {
#ifdef _OPENMP
#     pragma omp critical (find_split_rows)
#endif
      {
        try { split_rows.push_back(r); }
        catch (...) { ok = false; }
      }
    }
    if ( ! ok) continue;
  }
  if ( ! ok) throw hts::Exception("find_split_rows could not allocate.");
  std::sort(split_rows.begin(), split_rows.end());
}

template<typename Int>
inline double calc_centeredness (const Int i0, const Int n, const Int i) {
  return std::abs(i - (2*i0 + n - 1)/2.0);
}

template<typename Int>
inline Int find_split_row_if_available (
  const std::vector<Int>& split_rows, const Int r, const Int n,
  Int& split_row, double& centeredness)
{
  split_row = -1;
  centeredness = 0;
  Int pos = -1;
  const Int nsr = static_cast<Int>(split_rows.size());
  for (Int i = find_first<Int>(split_rows.data(), nsr, r); i < nsr; ++i) {
    const Int i_split_row = split_rows[i];
    if (i_split_row > r && i_split_row < r + n) {
      const double i_centeredness = calc_centeredness(r, n, i_split_row);
      if (split_row < 0 || i_centeredness < centeredness) {
        split_row = i_split_row;
        centeredness = i_centeredness;
        pos = i;
      }
    }
  }
  return pos;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::build_recursive_tri_r (
  const CrsMatrix& T, const Int r, const Int c, const Int n,
  const InitInfo& in, std::vector<Int>& split_rows,
  std::list<Box>& b)
{
  Int n1 = -1;
  {
    // Find a split if one exists and is suitable.
    Int split_row;
    double centeredness;
    const Int pos =
      find_split_row_if_available(split_rows, r, n, split_row, centeredness);
    if (split_row >= 0 &&
         // Constrain the split to be not too far from the center ...
        (centeredness <= 0.25*n ||
         // ... unless we're at or nearly at a leaf triangle.
         n <= 2*in.min_blksz)) {
      // Accept the split.
      n1 = split_row - r;
      // Remove it from the arry to reduce subsequent search time.
      split_rows.erase(split_rows.begin() + pos);
    }
  }

  if (n1 < 0) {
    // A split was not available.

    // Check if this is a leaf triangle.
    bool checked_density = false, is_dense = false;
    if (n <= in.min_blksz) {
      bool is_leaf = n <= in.min_parallel_rows;
      if ( ! is_leaf) {
        is_dense = is_dense_tri(T, r, c, n, in);
        checked_density = true;
        is_leaf = is_dense;
      }
      if (is_leaf) {
        // Leaf triangle.
        b.push_back(Box(r, c, n, n));
        return;
      }
    }

    // Not a leaf and no split is available, so bisect based on size.
    const Int blksz = n > in.min_blksz ? in.min_blksz : in.min_parallel_rows;
    const Int
      n_blks = (n + blksz - 1) / blksz,
      n1_blks = n_blks / 2;
    if (n_blks == 2 && n < 2*blksz) {
      // At the leaves associated with the remainder of the triangle, don't let
      // a block get tiny; merge into bigger blocks, instead.
      bool is_leaf = blksz == in.min_parallel_rows;
      if ( ! is_leaf) {
        if ( ! checked_density)
          is_dense = is_dense_tri(T, r, c, n, in);
        is_leaf = is_dense;
      }
      if (is_leaf)
        n1 = n;
      else
        n1 = n/2;
    } else
      n1 = n1_blks*blksz;
  }

  // Adjustments based on size have made this a leaf triangle, after all.
  if (n1 == n) {
    b.push_back(Box(r, c, n, n));
    return;
  }

  // Recurse.
  const Int
    n2 = n - n1,
    r1 = r + n1;
  build_recursive_tri_r(T, r, c, n1, in, split_rows, b);
  b.push_back(Box(r1, c, n2, n1));
  build_recursive_tri_r(T, r1, c + n1, n2, in, split_rows, b);
}

/* Decompose a shape
 *      _
 *     | |\
 *     |1|2\
 *     |_|__\
 *
 * where block 1 has mvp_block_nc columns, the top-left corner of block 1 is at
 * (r,c), and block 2 is nxn. Block 2 is recursively partitioned.
 *   b is returned in solution order.
 */
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
build_recursive_tri (const CrsMatrix& T, const Int r, const Int c, const Int n,
                     const Int mvp_block_nc, const InitInfo& in,
                     Array<Box>& bv) {
  // Find clear divisions in this triangle that we can exploit to pack the data
  // efficiently.
  std::vector<Int> split_rows;
  find_split_rows(T, r, c + mvp_block_nc, n, in, split_rows);

  std::list<Box> bl;
  // This is the large MVP block that scatters the LS part to the RB part. It is
  // present only in the T = L case; in the T = U case, it is in the LS data
  // structure.
  if (mvp_block_nc)
    bl.push_back(Box(r, c, n, mvp_block_nc));

  build_recursive_tri_r(T, r, c + mvp_block_nc, n, in, split_rows, bl);

  // list -> vector.
  bv.optclear_and_resize(bl.size());
  Int i = 0;
  for (typename std::list<Box>::const_iterator it = bl.begin();
       it != bl.end(); ++it)
    bv[i++] = *it;
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
RecursiveTri::init_invert_ondiag_tris_separately () {
  Int nt = 0;
  for (size_t i = 0; i < nd_.t.size(); ++i)
    if (nd_.t[i].is_dense_inverted())
      ++nt;
  invert_separately_ = nthreads_ > nt;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
RecursiveTri::init (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                    const InitInfo& in, const Int mvp_block_nc) {
  clear();
  this->n_ = n; this->r0_ = r0; nthreads_ = in.nthreads;

  Array<Box> bs;
  build_recursive_tri(T, r0, c0, n, mvp_block_nc, in, bs);

  const Int
    nblock = bs.size(),
    ntri = ((mvp_block_nc ? 2 : 1) + nblock) / 2,
    nmvp = ntri - 1;
  nd_.t.optclear_and_resize(ntri, OnDiagTri());
  nd_.s.optclear_and_resize(nmvp, TMatrix());
  nd_.os.optclear_and_resize(ntri);

  if (mvp_block_nc)
    nd_.t[0].init(0, 0, mvp_block_nc);
  const Int i0 = mvp_block_nc ? 1 : 0;
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
#ifdef _OPENMP
#   pragma omp for schedule(static)
#endif
    for (Int i = i0; i < ntri; ++i) {
      const Box& b = bs[i0 + 2*(i - i0)];
      nd_.t[i].init_metadata(T, b.r0, b.c0, b.nr, in);
    }
#ifdef _OPENMP
#   pragma omp for schedule(static)
#endif
    for (Int i = 0; i < nmvp; ++i) {
      const Box& b = bs[i0 + 2*(i - i0) + 1];
      nd_.s[i].init_metadata(T, b.r0, b.c0, b.nr, b.nc, in);
      nd_.os[i] = b.c0;
    }
  }
  // Initialize memory in order.
  nd_.os.back() = 0;
  for (Int i = 0; i < nmvp; ++i) {
    nd_.t[i].init_memory(in);
    nd_.s[i].init_memory(in);
  }
  nd_.t.back().init_memory(in);
  // Initialize numerical data.
  init_invert_ondiag_tris_separately();
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
#ifdef _OPENMP
#   pragma omp for
#endif
    for (Int i = i0; i < ntri; ++i)
      nd_.t[i].init_numeric(T, ! invert_separately_);
    const Size ilim = nmvp*nthreads_;
#ifdef _OPENMP
#   pragma omp for schedule(static, parfor_static_size)
#endif
    for (Size i = 0; i < ilim; ++i) {
      const Int si = i / nthreads_, tid = i % nthreads_;
      nd_.s[si].init_numeric(T, tid);
    }
  }
  if (invert_separately_) invert_ondiag_tris();

  { // Initialize threading data for the inverse of the on-diag tri.
    max_diag_tri_ = 0;
    Int max_nthreads = 0;
    for (size_t i = 0; i < nd_.t.size(); ++i) {
      max_diag_tri_ = std::max(max_diag_tri_, nd_.t[i].get_n());
      max_nthreads = std::max(max_nthreads, nd_.t[i].nthreads());
    }
    wrk_.optclear_and_resize(max_diag_tri_ * in.max_nrhs);
    nd_.inv_tri_done.optclear_and_resize(max_nthreads);
  }

  p2p_init();
}

template<typename Int, typename Size, typename Sclr>
inline Impl<Int, Size, Sclr>::
DenseTrisInverter::DenseTrisInverter (Array<Tri>& tris)
  : tris_(tris)
{
  nthreads_ = omp_get_max_threads();
  Int n = 0;
  for (size_t i = 0; i < tris_.size(); ++i)
    n = std::max(n, tris_[i].n);
  w_.optclear_and_resize(ntri<Int>(n));
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
DenseTrisInverter::invert (Sclr* const T, const Int n, Sclr* const w) {
#ifdef _OPENMP
# pragma omp for schedule(static, 1)
#endif
  for (Int c = 0; c < n; ++c) {
    // Solve for column c. That involves only the (n-c)x(n-c) lower-right
    // subblock of T. T is row major. w is col major.
    const Int nt = ntri<Int>(c - 1); // 0 - 1 = -1 is ok.
    Sclr* const Tc = T + 2*c + nt;
    Sclr* const wc = w + n*c - nt;
    wc[0] = Tc[0];
    Sclr* Tp = Tc + c + 1;
    for (Int r = 1; r < n - c; ++r) {
      Sclr acc = 0;
      for (Int k = 0; k < r; ++k)
        acc -= wc[k]*Tp[k];
      wc[r] = acc*Tp[r];
      Tp += r + c + 1;
    }
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
DenseTrisInverter::copy (Sclr* const T, const Int n, Sclr* const w) {
#ifdef _OPENMP
# pragma omp for schedule(static, 1)
#endif
  for (Int c = 0; c < n; ++c) {
    const Int nt = ntri<Int>(c - 1);
    Sclr* Tc = T + 2*c + nt;
    Sclr* const wc = w + n*c - nt;
    for (Int r = 0; r < n - c; ++r) {
      *Tc = wc[r];
      Tc += r + c + 1;
    }
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
DenseTrisInverter::compute () {
  for (size_t i = 0; i < tris_.size(); ++i) {
    Tri& t = tris_[i];
    invert(t.d, t.n, w_.data());
    copy(t.d, t.n, w_.data());
#ifdef _OPENMP
#   pragma omp barrier
#endif
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
RecursiveTri::invert_ondiag_tris () {
  typedef typename DenseTrisInverter::Tri DtiTri;
  Int ntri = 0;
  for (std::size_t i = 0; i < nd_.t.size(); ++i) {
    OnDiagTri& t = nd_.t[i];
    if (t.is_dense_inverted()) ++ntri;
  }
  Array<DtiTri> dtris(ntri);
  for (std::size_t i = 0, k = 0; i < nd_.t.size(); ++i) {
    OnDiagTri& t = nd_.t[i];
    if (t.is_dense_inverted()) {
      DtiTri& dti = dtris[k++];
      dti.idx = i;
      dti.d = t.dense_tri();
      dti.n = t.get_n();
    }
  }
  DenseTrisInverter dti(dtris);
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
    dti.compute();
#ifdef _OPENMP
#   pragma omp for
#endif
    for (size_t i = 0; i < dtris.size(); ++i)
      nd_.t[dtris[i].idx].inv_copy();
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
RecursiveTri::clear () {
  nd_.t.clear();
  nd_.s.clear();
  nd_.os.clear();
  nd_.s_done.clear();
  nd_.t_ids.clear();
  nd_.t_idx.clear();
  nd_.s_ids.clear();
  nd_.s_idx.clear();
  wrk_.clear();
  nd_.inv_tri_done.clear();
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
RecursiveTri::reinit_numeric (const CrsMatrix& T) {
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
    const Int n = static_cast<Int>(nd_.t.size());
#ifdef _OPENMP
#   pragma omp for
#endif
    for (Int i = 0; i < n; ++i)
      nd_.t[i].reinit_numeric(T, ! invert_separately_);
    const Int nmvp = static_cast<Int>(nd_.s.size());
    const Size ilim = nmvp*nthreads_;
#ifdef _OPENMP
#   pragma omp for
#endif
    for (Size i = 0; i < ilim; ++i) {
      const Int si = i / nthreads_;
      if (nd_.s[si].empty()) continue;
      const Int tid = i % nthreads_;
      nd_.s[si].reinit_numeric(T, tid);
    }
  }
  if (invert_separately_) invert_ondiag_tris();
}

// Form the lists of p2p dependencies. A dependency in this algorithm is of two
// types. One is the usual dependency: a variable has to be solved for before
// the next. The second imposes an ordering to assure there is no write race
// condition when computing a row's dot product in pieces.
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
RecursiveTri::p2p_init () {
  // For L, in the case where the LS block has rows > 0 and so there is an MVP
  // block, we rely on the fact that n_ doesn't count the rows associated with
  // the space-holding triangle.
  const Int
    tn = static_cast<Int>(nd_.t.size()),
    sn = static_cast<Int>(nd_.s.size());
  assert(sn+1 == tn);

  // w[r] is the current MVP block that has precedence for row r.
  const Size marker = tn*nthreads_;
  Array<Size> w(this->n_, marker);
  // Fast way to get a set of unique elements.
  Array<Int> w_cnted(sn * nthreads_, 0);
  Int w_cnted_symbol = 1;

  nd_.s_done.optclear_and_resize(nthreads_*sn, 0);
  nd_.t_idx.optclear_and_resize(tn);
  if (tn) nd_.t_idx[0] = 0;
  const Size nst = sn * nthreads_;
  Array<std::vector<Size> > s_ids(nst, std::vector<Size>());
  std::vector<Size> t_ids;
  Array<Int> s_ids_cnt(nthreads_, 0);
  for (Int ti = 0; ti < tn; ++ti) {
    // Consider each tri or MVP block in solution order.
    if (ti > 0) { // Tri block. First tri has no dependencies.
      const Tri& t = nd_.t[ti];
      const Int r0 = t.get_r0(), nr = t.get_n();
      Int k = nd_.t_idx[ti-1];
      for (Int r = r0, rlim = r0 + nr; r < rlim; ++r) {
        assert(r < this->n_);
        const Size wr = w[r];
        if (wr != marker && w_cnted[wr] != w_cnted_symbol) {
          const Int sid = rb_p2p_ind2sid(wr), tid = rb_p2p_ind2tid(wr);
          t_ids.push_back(rb_p2p_sub2ind(sid, tid));
          ++k;
          w_cnted[wr] = w_cnted_symbol;
        }
      }
      nd_.t_idx[ti] = k;
      ++w_cnted_symbol;
    }

    if (ti+1 < tn) { // MVP block.
      const TMatrix& s = nd_.s[ti];
      if (s.empty()) continue;
      for (Int bi = 0, bn = s.nblocks(); bi < bn; ++bi) {
        const Int r0 = s.block_r0(bi), nr = s.block_nr(bi);
        if (nr == 0) continue;
        const Size ind = rb_p2p_sub2ind(ti, bi);
        for (Int r = r0, rlim = r0 + nr; r < rlim; ++r) {
          const Size wr = w[r];
          // If row depends on an MVP block, and that block has not yet been
          // recorded in my dependency list, record it.
          if (wr != marker && w_cnted[wr] != w_cnted_symbol) {
            // If tids are the same, program order takes care of the dep.
            const Int tid = rb_p2p_ind2tid(wr);
            if (tid != bi) {
              s_ids[ind].push_back(wr);
              ++s_ids_cnt[bi];
            }
            w_cnted[wr] = w_cnted_symbol;
          }
          // I now have precedence.
          w[r] = ind;
        }
        ++w_cnted_symbol;
      }
    }
  }

  nd_.s_ids.optclear_and_resize(nthreads_);
  nd_.s_idx.optclear_and_resize(nthreads_);
  bool ok = true;
#ifdef _OPENMP
# pragma omp parallel
#endif
  do {
    const int tid = omp_get_thread_num();
    try {
      nd_.s_ids[tid].init();
      nd_.s_ids[tid].optclear_and_reserve_ft(s_ids_cnt[tid]);
      nd_.s_idx[tid].init();
      nd_.s_idx[tid].optclear_and_reserve_ft(sn + 1);
      nd_.s_idx[tid].unsafe_push_back(0);
    } catch (...) { ok = false; }
  } while (0);
  if ( ! ok) throw hts::Exception("RecursiveTri::p2p_init failed to allocate");
  for (Size i = 0; i < nst; ++i) {
    const Int tid = rb_p2p_ind2tid(i);
    const Int ndep = static_cast<Int>(s_ids[i].size());
    nd_.s_idx[tid].unsafe_push_back(nd_.s_idx[tid].back() + ndep);
    for (Int j = 0; j < ndep; ++j)
      nd_.s_ids[tid].unsafe_push_back(s_ids[i][j]);
  }
  s_ids.clear();
  nd_.t_ids.optclear_and_resize(t_ids.size());
  memcpy(nd_.t_ids.data(), t_ids.data(), t_ids.size()*sizeof(Size));
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::RecursiveTri::p2p_reset () const {
  nd_.t_barrier = -1;
  for (size_t i = 0; i < nd_.inv_tri_done.size(); ++i)
    nd_.inv_tri_done[i] = -1;
  ++nd_.done_symbol;
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
RecursiveTri::reset_max_nrhs (const Int max_nrhs) {
  wrk_.optclear_and_resize(max_diag_tri_ * max_nrhs);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
LevelSetTri::init_lsets (const LevelSetter& lstr,
                         const bool save_for_reprocess) {
  ls_blk_sz_ = lstr.ls_blk_sz();
  save_for_reprocess_ = save_for_reprocess;
  lsp_.optclear_and_resize(lstr.size() + 1);
  lsp_[0] = 0;
  for (Int i = 0; i < lstr.size(); ++i)
    lsp_[i+1] = lsp_[i] + lstr.lset(i).size();
  nlvls_ = static_cast<Int>(lsp_.size()) - 1;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
LevelSetTri::init (const CrsMatrix& T, const Int, const Int,
                   const Int in, const InitInfo& info) {
  n_ = in;
  t_.optclear_and_resize(info.nthreads, Thread());

  ps_.optclear_and_resize(info.nthreads);
  bool ok = true;
#ifdef _OPENMP
# pragma omp parallel
#endif
  do {
    const int tid = omp_get_thread_num();
    Thread& t = t_[tid];

    // Allocate and first touch.
    try {
      t.p.optclear_and_resize_ft(nlvls_);
      t.lsp.optclear_and_resize_ft(nlvls_ + 1);
    } catch (...) {
      ok = false;
      break;
    }

    // Get a thread's chunk of each level set.
    t.lsp[0] = 0;
#ifdef _OPENMP
#   pragma omp barrier
#   pragma omp for
#endif
    for (Int ils = 0; ils < nlvls_; ++ils) {
      const Int
        r0 = lsp_[ils],
        c0 = 0,
        nr = lsp_[ils+1] - r0,
        nc = lsp_[ils+1] + mvp_block_nc_;

      CrsSegmenter seg(T, r0, c0, nr, nc, info.nthreads, ls_blk_sz_);
      const Array<Int>& p = seg.get_p();
      assert(p[1] >= p[0]);
      const Int nseg = static_cast<Int>(p.size()) - 1;

      for (int tid1 = 0; tid1 < info.nthreads; ++tid1) {
        Thread& t1 = t_[tid1];
        if (tid1 < nseg) {
          t1.p[ils] = p[tid1];
          t1.lsp[ils+1] = p[tid1+1] - p[tid1];
        } else {
          t1.p[ils] = -1;
          t1.lsp[ils+1] = 0;
        }
      }
    }

    // Cumsum lsp.
    for (Int ils = 1; ils <= nlvls_; ++ils)
      t.lsp[ils] += t.lsp[ils-1];

    // Fill ps_.
    Array<Int>& p = ps_[tid];
    try {
      p.init();
      p.optclear_and_resize(t.lsp[nlvls_]);
    } catch (...) {
      ok = false;
      break;
    }
    Int k = 0;
    for (Int ils = 0; ils < nlvls_; ++ils) {
      const Int is = t.p[ils];
      if (is == -1) continue;
      const Int n = t.lsp[ils+1] - t.lsp[ils];
      for (Int i = 0; i < n; ++i)
        p[k++] = is + i;
    }
    assert(k == t.lsp[nlvls_]);
  } while (0);
  if ( ! ok) throw hts::Exception("LevelSetTri::init failed to allocate.");

  init_numeric(T);
  if ( ! save_for_reprocess_) ps_.clear();
}

// The original ordering of the LS block is by level, like this (l1t0 = lev 1,
// thread 0):
//     [l1t0, l1t1, ..., l2t0, l2t1, ...].
// The problem with this ordering is that it spreads out the variables relevant
// to a particular thread, causing a lot of paging. The solution is to reorder
// like this:
//     [l1t0, l2t0, ..., l1t1, l2t1, ...].
// The result matrix is in fact only psychologically triangular, but by then we
// already have fully analyzed it. Moreover, this reordering only requires two
// data changes: lsis is permuted, and t.m->jc is renumbered accordingly. In
// particular, t.m->d and t.m->ir are unchanged.
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
LevelSetTri::update_permutation (Array<Int>& lsis, const Partition& p) {
  Array<Int> old_lsis(lsis.size());
  Array<Int> q(p.cm->n);
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
#ifdef _OPENMP
#   pragma omp for
#endif
    for (std::size_t i = 0; i < lsis.size(); ++i)
      old_lsis[i] = lsis[i];
    const int tid = omp_get_thread_num();
    Thread& t = t_[tid];
    int start = 0;
    for (int i = 0; i < tid; ++i)
      start += t_[i].m->m;
    for (Int ils = 0, k = 0; ils < nlvls_; ++ils) {
      const Int p_ils = t.p[ils];
      t.p[ils] = start + k;
      const Int n = t.lsp[ils+1] - t.lsp[ils];
      for (Int i = 0; i < n; ++i, ++k) {
        lsis[start + k] = old_lsis[p_ils + i];
        q[p_ils + i] = start + k;
      }
    }
#ifdef _OPENMP
#   pragma omp barrier
#endif
    for (Size i = 0, nnz = t.m->ir[t.m->m]; i < nnz; ++i) {
      Int& jci = t.m->jc[i];
      if (jci < mvp_block_nc_) continue;
      jci = mvp_block_nc_ + q[jci - mvp_block_nc_];
    }
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::LevelSetTri::
find_task_responsible_for_variable (p2p_Pair* const pairs) {
  const int tid = omp_get_thread_num();
  // Find the level and tid responsible for variable i.
  Thread& t = t_[tid];
  const Array<Int>& lsp = t.lsp;
  const Array<Int>& p = t.p;
  for (Int ils = 0, ils_lim = static_cast<Int>(lsp.size()) - 1;
       ils < ils_lim; ++ils) {
    const Int r0 = p[ils];
    for (Int j = 0, jlim = lsp[ils+1] - lsp[ils]; j < jlim; ++j) {
      const Int i = r0 + j;
      pairs[i].tid = tid;
      pairs[i].lvl = ils;
    }
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
Int Impl<Int, Size, Sclr>::LevelSetTri::
fill_graph (const p2p_Pair* const pairs, Size* const g, Size* const gp,
            Size* const wrk) {
  // O(nnz/nthreads) time and O(nnz) space. Reduce the entry-wise graph to the
  // (lvl, tid)-wise dependency graph.
  const int tid = omp_get_thread_num();
  const Size n = nlvls_*t_.size();
  Thread& t = t_[tid];
  const Array<Int>& lsp = t.lsp;
  // Max value of e + 1, used as a marker.
  const Size& me = n;
  // Build the graph g. g[e] is the list of dependencies for (level, tid) e.
  const CrsMatrix* const cm = t.m;
  const Size* const ir = cm->ir;
  const Int* const jc = cm->jc;

  // Get an unused piece of the workspace.
  Size* Te = wrk;
  for (int i = 0; i < tid; ++i)
    if (t_[i].m) {
      const CrsMatrix* const m = t_[i].m;
      Te += m->ir[m->m];
    }

  // Count entries in g.
  gp[ls_p2p_sub2ind(1, tid)] = 0;
  for (Int ils = 1, ils_lim = nlvls_; ils < ils_lim; ++ils) {
    const Size e = ls_p2p_sub2ind(ils, tid);
    const Size Ten = ir[lsp[ils+1]] - ir[lsp[ils]];
    // Record all e = (lvl, tid) dependencies in a CRS matrix structurally
    // identical to cm.
    for (Int r = lsp[ils], rlim = lsp[ils+1]; r < rlim; ++r)
      for (Size j = ir[r], jlim = ir[r+1]; j < jlim; ++j) {
        const Int c = jc[j] - mvp_block_nc_;
        if (c < 0) {
          Te[j] = me;
          continue;
        }
        const Int r_tid = pairs[c].tid, r_lvl = pairs[c].lvl;
        const Size ed = ls_p2p_sub2ind(r_lvl, r_tid);
        Te[j] = ed;
      }
    // Sort the values ...
    Size* const Te0 = Te + ir[lsp[ils]];
    std::sort(Te0, Te0 + Ten);
    // ... so that it's fast to count the number of unique entries.
    Int k = 0;
    Size prev = me;
    for (Size i = 0; i < Ten; ++i) {
      if (Te0[i] == me) {
        // All me values sort to the end, so we can break.
        break;
      }
      if (Te0[i] != prev) {
        prev = Te0[i];
        const Int r_lvl = ls_p2p_ind2lvl(prev);
        if (r_lvl != ils)
          ++k;
      }
    }
    // Now we know how much space to allocate.
    gp[e+1] = k;
  }

  // Cumsum to get row pointers.
  Int max_gelen = 0;
#ifdef _OPENMP
# pragma omp barrier
# pragma omp master
#endif
  for (Size i = 1; i <= n; ++i) {
    if (static_cast<Int>(gp[i]) > max_gelen) max_gelen = gp[i];
    gp[i] += gp[i-1];
  }
#ifdef _OPENMP
# pragma omp barrier
#endif

  // Fill g. Can reuse the data in Te. Everything in this loop is identical to
  // the previous except that we can now record the unique values.
  for (Int ils = 1, ils_lim = static_cast<Int>(lsp.size()) - 1;
       ils < ils_lim; ++ils) {
    const Size e = ls_p2p_sub2ind(ils, tid);
    Size* const Te0 = Te + ir[lsp[ils]];
    const Size Ten = ir[lsp[ils+1]] - ir[lsp[ils]];
    Size* const ge_g = &g[gp[e]];
    Int k = 0;
    Size prev = me;
    for (Size i = 0; i < Ten; ++i) {
      if (Te0[i] == me) break;
      if (Te0[i] != prev) {
        prev = Te0[i];
        const Int r_lvl = ls_p2p_ind2lvl(prev);
        if (r_lvl != ils)
          ge_g[k++] = prev;
      }
    }
    assert(k == static_cast<Int>(gp[e+1] - gp[e]));
    // Sort for the prune_graph phase.
    std::sort(ge_g, ge_g + k);
  }

  return max_gelen;
}

//pre a is sorted.
//pre b is sorted.
//pre len(mark) == an.
//pre mark can == a.
// O(an log bn)
template<typename Int, typename Size> inline Int
mark_intersection (const Size* a, const Int an, const Size* b, Int bn,
                   Size* mark /* len(mark) == an */, const Size marker) {
  Int nmarked = 0;
  for (Int i = 0; i < an; ++i) {
    if (a[i] == marker) continue;
    const Int n = static_cast<Int>(std::lower_bound(b, b + bn, a[i]) - b);
    if (n < bn && b[n] == a[i]) {
      mark[i] = marker;
      ++nmarked;
      b += n;
      bn -= n;
    }
  }
  return nmarked;
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::LevelSetTri::
prune_graph (const Size* const gc, const Size* const gp, Size* const g,
             Int* const gsz, Size* const wrk, const Int max_gelen) {
  // Space is O(#threads max |g(e)|).
  const int tid = omp_get_thread_num();
  const Size n = nlvls_*t_.size();
  const Size& me = n; // marker
  Size* const mark = &wrk[tid*max_gelen];
  // I find that it's more efficient to use a parfor than to work on a thread's
  // own e's.
#ifdef _OPENMP
# pragma omp for schedule(static,1)
#endif
  for (Size e = 0; e < n; ++e) {
    const Int gcelen = static_cast<Int>(gp[e+1] - gp[e]);
    gsz[e] = 0;
    if (gcelen == 0) continue;
    // e's dependencies.
    const Size* const gce = &gc[gp[e]];
    // Record only the deepest task in a tid group. The deps are sorted first by
    // tid group, then by lvl, which we exploit here.
    Int mark_len = 0;
    Int prev_lvl_len = 0;
    {
      const Int elvl = ls_p2p_ind2lvl(e), etid = ls_p2p_ind2tid(e);
      Size prev_e = me;
      Int prev_tid = static_cast<Int>(t_.size());
      for (Int i = 0; i <= gcelen; ++i) {
        // A little bit of messiness so we can grab the last one.
        Int edtid = 0;
        if (i < gcelen)
          edtid = ls_p2p_ind2tid(gce[i]);
        if (i == gcelen || edtid != prev_tid) {
          if (prev_e != me) {
            const Int edlvl = ls_p2p_ind2lvl(prev_e);
            if (edlvl+1 == elvl) {
              // If a dep is only one level shallower, it cannot be
              // pruned. Record it for later insertion unless its tid is e's
              // tid, in which case it can and will be pruned.
              if (prev_tid != etid)
                mark[gcelen - prev_lvl_len++ - 1] = prev_e;
            } else {
              // The dep is prunable.
              mark[mark_len++] = prev_e;
            }
          }
          prev_tid = edtid;
        }
        if (i == gcelen) break;
        prev_e = gce[i];
      }
    }
    assert(mark_len + prev_lvl_len <= gcelen);
    Int nmarked = 0;
    for (Int ied = 0; ied < gcelen; ++ied) { // For each of e's deps:
      const Size ed = gce[ied];
      assert(ed >= 0 && ed < n);
      if (ls_p2p_ind2lvl(ed) == 0) continue; // No parent deps to check.
      // ed's dependencies.
      const Size* const gced = &gc[gp[ed]];
      const Int gcedlen = static_cast<Int>(gp[ed+1] - gp[ed]);
      nmarked += mark_intersection(mark, mark_len, gced, gcedlen, mark, me);
      if (nmarked == mark_len) break;
    }
    // Insert the pruned set of dependencies.
    Int k = 0;
    Size* const ge = &g[gp[e]];
    const Int etid = ls_p2p_ind2tid(e);
    for (Int i = 0; i < mark_len; ++i) {
      const Size ed = mark[i];
      if (ed == me) continue;
      const Int edtid = ls_p2p_ind2tid(ed);
      if (edtid != etid) ge[k++] = ed;
    }
    for (Int i = 0; i < prev_lvl_len; ++i)
      ge[k++] = mark[gcelen - i - 1];
    gsz[e] = k;
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::LevelSetTri::
fill_dependencies (const Size* const g, const Size* const gp,
                   const Int* const gsz) {
  const int tid = omp_get_thread_num();
  Thread& t = t_[tid];
  const Array<Int>& lsp = t.lsp;
  // Allocate inside this thread up front. I could probably do this even more
  // efficiently, but fill_dependencies is negligible compared with fill_graph
  // and prune_graph.
  t.p2p_depends_p.optclear_and_reserve(nlvls_ + 1);
  Int sz = 0;
  for (Int ils = 1, ils_lim = static_cast<Int>(lsp.size()) - 1;
       ils < ils_lim; ++ils)
    sz += gsz[ls_p2p_sub2ind(ils, tid)];
  t.p2p_depends.optclear_and_reserve(sz);
  // Make the final lists of dependencies.
  t.p2p_depends_p.unsafe_push_back(0);
  for (Int ils = 1, ils_lim = static_cast<Int>(lsp.size()) - 1;
       ils < ils_lim; ++ils) {
    const Size e = ls_p2p_sub2ind(ils, tid);
    const Size* const ed = &g[gp[e]];
    const Int edsz = gsz[e];
    // Insert.
    t.p2p_depends_p.unsafe_push_back(t.p2p_depends_p[ils-1]);
    for (Int ied = 0; ied < edsz; ++ied) {
      t.p2p_depends.unsafe_push_back(ed[ied]);
      ++t.p2p_depends_p.back();
    }
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::LevelSetTri::p2p_init () {
  if (t_[0].lsp.empty()) return;

  if (p2p_done_.empty()) {
    p2p_done_value_ = 0;
    p2p_done_.optclear_and_resize(t_[0].lsp.size()*t_.size(), 1);
  }
  if (t_[0].lsp.size() <= 1) return;

  Size nnz = 0;
  for (size_t i = 0, n = t_.size(); i < n; ++i)
    if (t_[i].m) {
      const CrsMatrix* const cm = t_[i].m;
      nnz += cm->ir[cm->m];
    }
  const Size n = t_.size() * nlvls_;

  // g is a graph. nnz is an upper bound on the memory needed.
  //   g(e) is the set of dependencies (incoming edges) for node e. gp is the
  // pointer to g(e). So g(e) is written g[gp[e]]. gsz is the size of g(e).
  //   A node e represents a pair (level, thread id) that is a task.
  Array<Size> as;
  Size* g, * gc, * gp;
  Array<Int> gsz;
  // Thread and level responsible for a variable.
  Array<p2p_Pair> pairs;
  Array<Size> wrk;
  try {
    wrk.optclear_and_resize(nnz);
    gsz.optclear_and_resize(n);
    pairs.optclear_and_resize(n_);
    as.optclear_and_resize(2*nnz + n + 1);
    g = as.data();
    gp = g + nnz;
    gp[0] = 0;
    gc = gp + n + 1;
  } catch (...) {
    std::stringstream ss;
    ss << "p2p_init failed to resize: n = " << n << " nnz = " << nnz;
    throw hts::Exception(ss.str());
  }

  Int max_gelen;
  bool ok = true;
#ifdef _OPENMP
# pragma omp parallel
#endif
  do {
    find_task_responsible_for_variable(pairs.data());
#ifdef _OPENMP
#   pragma omp barrier
#endif
    const Int max_gelen_t = fill_graph(pairs.data(), g, gp, wrk.data());
#ifdef _OPENMP
#   pragma omp barrier
#   pragma omp master
#endif
    {
      pairs.clear();
      // Tell all threads; only master's max_gelen_t is valid.
      max_gelen = max_gelen_t;
      // In the unlikely case that max_gelen * #threads > nnz, allocate more
      // workspace.
      const Size space = max_gelen * t_.size();
      if (space > nnz) {
        try { wrk.optclear_and_resize(space); }
        catch (...) { ok = false; }
      }
    }
    // Keep the original graph.
#ifdef _OPENMP
#   pragma omp for
#endif
    for (Size i = 0; i < nnz; ++i)
      gc[i] = g[i];
    if ( ! ok) break;
    prune_graph(gc, gp, g, gsz.data(), wrk.data(), max_gelen);
#ifdef _OPENMP
#   pragma omp barrier
#endif
    try { fill_dependencies(g, gp, gsz.data()); }
    catch (...) { ok = false; }
  } while (0);
  if ( ! ok) throw hts::Exception("p2p_init failed to reallocate wrk.");
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::LevelSetTri::p2p_reset () const {
  p2p_done_value_ = ! p2p_done_.empty() ? (p2p_done_[0] + 1) % 2 : 0;
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
LevelSetTri::init_numeric (const CrsMatrix& T) {
  bool ok = true;
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
    const int tid = omp_get_thread_num();
    del(t_[tid].m);
    try { t_[tid].m = get_matrix_p(T, ps_[tid], true); }
    catch (...) { ok = false; }
  }
  if ( ! ok)
    throw hts::Exception("LevelSetTri::init_numeric failed to get_matrix_p.");
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
LevelSetTri::reinit_numeric (const CrsMatrix& T) {
  bool nthreads_ok;
  std::string msg;
#ifdef _OPENMP
# pragma omp parallel
#endif
  do {
    nthreads_ok = check_nthreads(t_.size(), omp_get_num_threads(), msg);
    if ( ! nthreads_ok) break;
    const int tid = omp_get_thread_num();
    Thread& t = t_[tid];
    assert(t.m);
    const Array<Int>& p = ps_[tid];
    for (size_t i = 0; i < p.size(); ++i) {
      const Int r = p[i], nc = static_cast<Int>(T.ir[r+1] - T.ir[r]);
      assert(nc == static_cast<Int>(t.m->ir[i+1] - t.m->ir[i]));
      Sclr* const d_start = t.m->d + t.m->ir[i];
      memcpy(d_start, T.d + T.ir[r], nc*sizeof(*t.m->d));
      d_start[nc-1] = Sclr{1.0}/d_start[nc-1];
    }
  } while (0);
  if ( ! nthreads_ok) throw hts::Exception(msg);
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::Permuter::clear () {
  if (q_ && q_ != p_) { deln(q_); }
  deln(p_);
  deln(scale_);
  deln(px_);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::Permuter::init (
  const Int n, const bool is_lo, const Array<Int>& lsis, const Array<Int>& dpis,
  const Int nthreads, const Int max_nrhs, const Int* p, const Int* q,
  const Real* scale)
{
  clear();
  n_ = n;
  max_nrhs_ = max_nrhs;
  is_lo_ = is_lo;
  p_ = q_ = 0;
  px_ = 0;
  scale_ = 0;
  try {
    p_ = allocn<Int>(n_);
    px_ = allocn<Sclr>(n_*max_nrhs);
    q_ = p || q ? allocn<Int>(n_) : p_;
    if (scale) {
      scale_ = allocn<Real>(n_);
      reinit_numeric(scale);
    }
  } catch (...) {
    clear();
    throw hts::Exception("Permuter::init failed to allocate.");
  }

  partition_n_uniformly(n_, nthreads, part_);

  // Load p_ and q_ with internal and possibly user permutations.
  if (p || q) {
    const Int dpis_sz = static_cast<Int>(dpis.size());
    // Incorporate user's permutations.
    for (Int i = 0, ilim = static_cast<Int>(lsis.size()),
           k = is_lo_ ? 0 : dpis_sz;
         i < ilim; ++i) {
      Int li = lsis[i];
      if ( ! is_lo_) li = n_ - li - 1;
      p_[k + i] = p ? p[li] : li;
      q_[k + i] = q ? q[li] : li;
    }
    for (Int i = 0, ilim = dpis_sz, k = is_lo_ ? n - dpis_sz : 0;
         i < ilim; ++i) {
      Int di = dpis[i];
      if ( ! is_lo_) di = n_ - di - 1;
      p_[k + i] = p ? p[di] : di;
      q_[k + i] = q ? q[di] : di;
    }
  } else {
    // Just the internal (level set and data parallel) permutations.
    if ( ! lsis.empty())
      memcpy(p_ + (is_lo_ ? 0 : dpis.size()), lsis.data(),
             lsis.size()*sizeof(*p_));
    if ( ! dpis.empty())
      memcpy(p_ + (is_lo_ ? n - dpis.size() : 0), dpis.data(),
             dpis.size()*sizeof(*p_));
    if ( ! is_lo_)
      for (Int i = 0; i < n_; ++i) p_[i] = n_ - p_[i] - 1;
  }
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
Permuter::reinit_numeric (const Real* scale) {
  if (scale)
    memcpy(scale_, scale, n_*sizeof(*scale_));
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
Permuter::reset_max_nrhs (const Int max_nrhs) {
  if (max_nrhs_ == max_nrhs) return;
  max_nrhs_ = max_nrhs;
  deln(px_);
  px_ = Allocnator<Sclr>(n_*max_nrhs_, "Permuter::reset_max_nrhs").release();
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
Permuter::check_nrhs (const Int nrhs) const {
  if (nrhs > max_nrhs_)
    throw hts::Exception("nrhs is > max_nrhs.");
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::TriSolver::clear () {}

// len(irt) == n-1 because of how the transpose algorithm offsets irt. Partition
// by nonzeros. At exit, start[i] is the starting column for the i'th thread.
template<typename Int, typename Size>
void partition_irt (const Int n, const Size* const irt, const Size nnz,
                    const Int nparts, Array<Int>& start) {
  const Int nseg = std::min<Int>(nparts, n);
  Int i0 = 1, j0 = 1;
  start.optclear_and_resize(nseg);
  start[0] = 0;

  for (Int i = i0; i < nseg; ++i) {
    const double d = ((double) i / nseg)*nnz;
    Int j = static_cast<Int>(std::upper_bound(irt + j0, irt + n - 1, d) - irt);
    if (d - irt[j-1] > irt[j] - d) ++j;
    j = std::min<Int>(j, n - nseg + i);
    if (j < j0) {
      // If not all the threads will have work, let the earlier threads have
      // the work.
      j = std::min<Int>(j0, n);
      assert(j == start[i-1] + 1);
    }
    start[i] = j;
    j0 = j + 1;
  }
}

template<typename Int, typename Size, typename Sclr>
void transpose (
  const Int m, const Size* const ir, const Int* const jc, const Sclr* const d,
  Size* const irt, Int* const jct, Sclr* const dt,
  Array<Size>* transpose_perm = 0)
{
  const Size nnz = ir[m];
  if (transpose_perm) transpose_perm->optclear_and_resize(nnz);
  Array<Int> start;
#ifdef _OPENMP
# pragma omp parallel
#endif
  {
#ifdef _OPENMP
#   pragma omp for schedule(static)
#endif
    for (Int i = 0; i <= m; ++i)
      irt[i] = 0;
#ifdef _OPENMP
#   pragma omp single
#endif
    {
      // 1. Count the number of entries in each col.
      for (Size k = 0; k < nnz; ++k) {
        // Store shifted up by 1. This makes extra workspace unnecessary.
        ++irt[jc[k]+1];
      }
      // 2. Cumsum to get the col pointers.
      // Could do a parallel scan, but this step is extremely negligible.
      Size sum = 0;
      for (Int i = 0; i < m; ++i) {
        const Size incr = irt[i+1];
        irt[i+1] = sum;
        sum += incr;
      }
      assert(sum == nnz);
      // At this point, At.ir[i+1] gives what is normally in At.ir[i].
      // 3. Fill in jc and d.
      partition_irt<Int, Size>(m, irt + 2, nnz, omp_get_num_threads(), start);
    }
    const int tid = omp_get_thread_num();
    if (tid < static_cast<int>(start.size())) {
      const Int
        start_tid = start[tid],
        stop = tid + 1 == static_cast<int>(start.size()) ? m : start[tid+1];
      for (Int r = 0; r < m; ++r) {
        const Size irr = ir[r], jlim = ir[r+1];
        if (irr == jlim || stop <= jc[irr] || start_tid > jc[jlim-1]) continue;
        for (Size j = irr + find_first<Int>(jc + irr, jlim - irr, start_tid);
             j < jlim; ++j) {
          const Int c = jc[j];
          if (c >= stop) break;
          // Points to next entry to be filled.
          // 3a. Increment to next available entry.
          const Size p = irt[c+1]++;
          jct[p] = r;
          dt[p] = d[j];
          if (transpose_perm)
            (*transpose_perm)[p] = j;
        }
      }
    }
  }
  // 4. Step 3a made it such that At.ir[i] is now what we need except for
  // i = 0. At.ir[0] was already set to 0 in step 1, so we're done.
  assert(irt[m] == nnz);
}

template<typename Int, typename Size, typename Sclr>
typename Impl<Int, Size, Sclr>::ConstCrsMatrix* Impl<Int, Size, Sclr>::
transpose (const ConstCrsMatrix& T, Array<Size>* transpose_perm) {
  const Int n = T.m;
  SparseData sd(n, T.ir[n], "transpose");
  htsimpl::transpose(n, T.ir, T.jc, T.d, sd.ir, sd.jc, sd.d, transpose_perm);
  ConstCrsMatrix* ccm = 0;
  try {
    ccm = new ConstCrsMatrix(n, n, sd.ir, sd.jc, sd.d,
                             Direction::opposite(T.dir), T.conj, true,
                             T.unitdiag, ! T.is_lo);
  } catch (...) { throw hts::Exception("transpose failed to alloc."); }
  sd.release();
  return ccm;
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
compose_transpose (const Array<Size>& transpose_perm, Partition& p) {
  Array<Size> A_idxs_copy(p.nnz);
#ifdef _OPENMP
# pragma omp for schedule(static)
#endif
  for (Size i = 0; i < p.nnz; ++i)
    A_idxs_copy[i] = p.A_idxs[i];
#ifdef _OPENMP
# pragma omp for schedule(static)
#endif
  for (Size i = 0; i < p.nnz; ++i) {
    if (p.A_valid(i))
      p.A_idxs[i] = transpose_perm[A_idxs_copy[i]];
    // Otherwise, A_idxs[i] is already invalidated.
  }
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TriSolver::init (const ConstCrsMatrix* T, Int nthreads, const Int max_nrhs,
                 const bool save_for_reprocess, const Int* p, const Int* q,
                 const Real* scale, const Options& o) {
  NumThreads nthreads_state;
  set_num_threads(nthreads, nthreads_state);
  throw_if_nthreads_not_ok(nthreads);
  bool delete_T = false;
  try {
    // Basic size parameters.
    n_ = T->m;
    const Size nnz = T->ir[T->m];
    nthreads_ = nthreads;
    InitInfo in;
    in.nthreads = nthreads_;
    in.min_blksz = o.min_block_size;
    in.min_dense_density = o.min_dense_density;
    in.max_nrhs = max_nrhs;
    in.min_parallel_rows = o.min_parallel_rows;

    const bool make_transpose = T->dir == Direction::transpose;
    Array<Size> transpose_perm;
    if (make_transpose) {
      assert( ! delete_T);
      const ConstCrsMatrix* Tp = T;
      T = transpose(*Tp, save_for_reprocess ? &transpose_perm : 0);
      Tp->deallocate();
      delete_T = true;
    }

    // Determine shape.
    Shape shape = determine_shape(*T);
    if ( ! shape.is_triangular)
      throw hts::NotTriangularException();
    unitdiag_ = T->unitdiag = false;
    if (shape.has_no_diag) {
      // We don't make the user provide this information; inject it now that it
      // is known.
      unitdiag_ = T->unitdiag = true;
    } else if ( ! shape.has_full_diag) {
      throw hts::NotFullDiagonalException();
    }
    is_lo_ = T->is_lo = shape.is_lower;

    // Find level sets.
    LevelSetter lstr;
    const Int lstr_threshold = o.lset_min_size *
      (o.lset_min_size_scale_with_nthreads ? nthreads_ : 1);
    lstr.init(*T, lstr_threshold, is_lo_, o);

    if ( ! is_lo_) {
      const ConstCrsMatrix* Tp = T;
      T = permute_to_other_tri(*Tp);
      if (delete_T) del(Tp); else Tp->deallocate();
      delete_T = true;
      lstr.reverse_variable_order(n_);
    }

    // Separate into three blocks: level set, scatter, and data parallel:
    //     [(1)        0
    //      (scatter) (2)],
    // where if is_lo_, (1) = level sets and (2) = data parallel; if ! is_lo_,
    // then the opposite.
    Array<Int> lsis, dpis;
    get_idxs(n_, lstr, lsis, dpis);
    if (o.printlvl > 0)
      std::cout << "n " << n_ << " |lsis| " << lsis.size() << " |dpis| "
                << dpis.size() << "\n";
    if (is_lo_) {
      // 1. Level-scheduling block.
      {
        PermVec lsis_pv(T->m, lsis);
        get_matrix_pp_with_covers_all(*T, lsis_pv, p_[0], save_for_reprocess);
      }
      lst_.init_lsets(lstr, save_for_reprocess);
      lst_.init(*p_[0].cm, 0, 0, p_[0].cm->m, in);
      lst_.update_permutation(lsis, p_[0]);
      if ( ! save_for_reprocess) p_[0].clear();
      {
        PermVec dpis_pv(T->m, dpis), lsis_pv(T->m, lsis);
        get_matrix_p_qp_with_covers_all(*T, dpis_pv, lsis_pv, p_[1],
                                        save_for_reprocess);
      }
      if (delete_T) del(T); else T->deallocate();
      lst_.p2p_init();
      if (p_[1].cm->m > 0) {
        // 2. No MVP block. It's in the data-parallel block.
        // 3. Data-parallel block (+ MVP block).
        const Int mvp_nc = p_[1].cm->n - p_[1].cm->m;
        t_.init(*p_[1].cm, 0, 0, p_[1].cm->m, in, mvp_nc);
      }
      if ( ! save_for_reprocess) p_[1].clear();
    } else {
      {
        PermVec dpis_pv(T->m, dpis);
        get_matrix_pp_with_covers_all(*T, dpis_pv, p_[1], save_for_reprocess);
        {
          PermVec lsis_pv(T->m, lsis);
          get_matrix_p_qp_with_covers_all(*T, lsis_pv, dpis_pv, p_[0],
                                          save_for_reprocess);
        }
      }
      if (delete_T) del(T); else T->deallocate();
      if (p_[1].cm->m > 0) {
        // 3. Data-parallel block.
        t_.init(*p_[1].cm, 0, 0, p_[1].cm->m, in);
        // 2. No MVP block. It's in the level scheduling block.
      }
      if ( ! save_for_reprocess) p_[1].clear();
      // 1. Level-scheduling block (+ MVP block).
      lst_.init_lsets(lstr, save_for_reprocess);
      lst_.set_mvp_block_nc(dpis.size());
      lst_.init(*p_[0].cm, 0, 0, p_[0].cm->m, in);
      lst_.update_permutation(lsis, p_[0]);
      if ( ! save_for_reprocess) p_[0].clear();
      lst_.p2p_init();
    }
    if (save_for_reprocess) {
      // For 32-bit Int, 64-bit Real, save 2/3 of memory during solves at little
      // extra init_numeric cost.
      for (Int i = 0; i < 2; ++i)
        if (p_[i].cm) {
          p_[i].clear_d();
          if ( ! is_lo_)
            reverse_A_idxs(nnz, p_[i]);
          if (make_transpose)
            compose_transpose(transpose_perm, p_[i]);
        }
    }
    perm_.init(n_, is_lo_, lsis, dpis, nthreads_, max_nrhs, p, q, scale);
  } catch (...) {
    if (delete_T) del(T); else T->deallocate();
    throw;
    // restore_num_threads(nthreads_state); // unreachable
  }
  restore_num_threads(nthreads_state);
}

// Reinitialize numbers, but keep the same structures.
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TriSolver::reinit_numeric (const ConstCrsMatrix* T, const Real* r) {
  T->is_lo = is_lo_;
  T->unitdiag = unitdiag_;
  NumThreads nthreads_state;
  set_num_threads(nthreads_, nthreads_state);
  for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].alloc_d();
  repartition_into_2_blocks(p_, *T);
  lst_.reinit_numeric(*p_[0].cm);
  if (p_[1].cm->m > 0) t_.reinit_numeric(*p_[1].cm);
  for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].clear_d();
  if (r) perm_.reinit_numeric(r);
  restore_num_threads(nthreads_state);
}

template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::
TriSolver::reset_max_nrhs (const Int max_nrhs) {
  perm_.reset_max_nrhs(max_nrhs);
  t_.reset_max_nrhs(max_nrhs);
}

//> Solve code.

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
OnDiagTri::solve (const Sclr* b, const Int ldb, Sclr* x, const Int ldx,
                  const Int nrhs) const {
  if (d_) {
    return t_.empty() ? solve_dense(b, ldb, x, ldx, nrhs) :
      solve_dense_inv(b, ldb, x, ldx, nrhs);
  }
  if (m_) return solve_spars(b, ldb, x, ldx, nrhs);
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
OnDiagTri::solve_dense (const Sclr* b, const Int ldb, Sclr* x, const Int ldx,
                        const Int nrhs) const {
  for (Int irhs = 0; ; ) {
    for (Int j = 0, k = 0; j < this->n_; ++j) {
      Sclr a = b[j];
      const Int ilim = j;
      for (Int i = 0; i < ilim; ++i, ++k)
        a -= x[i]*d_[k];
      x[j] = a * d_[k++];
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldb;
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
void Impl<Int, Size, Sclr>::OnDiagTri::
solve_dense_inv (const Sclr* b, const Int ldb, Sclr* x, const Int ldx,
                 const Int nrhs) const {
  const int tid = omp_get_thread_num();
  if (tid >= nthreads()) return;
  const Thread& t = t_[tid];
  for (Int irhs = 0; ; ) {
    for (Int r = t.r0, rlim = t.r0 + t.nr, k = 0; r < rlim; ++r) {
      Sclr a = 0;
      for (Int c = 0; c < r+1; ++c, ++k)
        a += b[c]*t.d[k];
      x[r] = a;
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldb;
  }
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
OnDiagTri::solve_spars (const Sclr* b, const Int ldb, Sclr* x, const Int ldx,
                        const Int nrhs) const {
  const CrsMatrix& T = *m_;
  const Int m = T.m;
  const Size* const ir = T.ir;
  const Int* const jc = T.jc;
  const Sclr* const d = T.d;
  for (int irhs = 0; ; ) {
    for (int r = 0; r < m; ++r) {
      const Size
        rp_rp1 = ir[r+1],
        jlim = rp_rp1 - 1;
      Sclr a = b[r];
      for (Size j = ir[r]; j < jlim; ++j)
        a -= x[jc[j]] * d[j];
      x[r] = a * d[rp_rp1 - 1];
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldb;
  }
}

template<typename Int, typename Size, typename Sclr>
inline void SerialBlock_n1Axpy_spars (
  const Int nr_, const Int nc_, const Size* const ir_, const Int* const jc_,
  const Sclr* const d_, const Sclr* x, const Int ldx, const Int nrhs,
  Sclr* y, const Int ldy)
{
  for (Int k = 0; ; ) {
    Size iri = ir_[0];
    for (Int i = 0; i < nr_; ++i) {
      const Size irip1 = ir_[i+1];
      const Int N = static_cast<Int>(irip1 - iri);
      if (N == 0) continue;
      Sclr a = 0; {
        const Sclr* const d = d_ + iri;
        const Int* const jc = jc_ + iri;
        for (Int j = 0; j < N; ++j)
          a += d[j] * x[jc[j]];
      }
      y[i] -= a;
      iri = irip1;
    }
    if (++k == nrhs) break;
    x += ldx;
    y += ldy;
  }
}

#ifdef HAVE_SHYLU_NODEHTS_MKL
template<> inline void SerialBlock_n1Axpy_spars<MKL_INT, MKL_INT, float> (
  const MKL_INT nr_, const MKL_INT nc_, const MKL_INT* const ir_,
  const MKL_INT* const jc_, const float* const d_, const float* x,
  const MKL_INT ldx, const MKL_INT nrhs, float* y, const MKL_INT ldy)
{
  hts_mkl_csrmm(false, nr_, nc_, d_, ir_, jc_, x, ldx, y, ldy, nrhs);
}

template<> inline void SerialBlock_n1Axpy_spars<MKL_INT, MKL_INT, double> (
  const MKL_INT nr_, const MKL_INT nc_, const MKL_INT* const ir_,
  const MKL_INT* const jc_, const double* const d_, const double* x,
  const MKL_INT ldx, const MKL_INT nrhs, double* y, const MKL_INT ldy)
{
  hts_mkl_csrmm(false, nr_, nc_, d_, ir_, jc_, x, ldx, y, ldy, nrhs);
}
#endif

template<typename Int, typename Sclr>
inline void SerialBlock_n1Axpy_dense (
  const Int nr, const Int nc, const Sclr* d, const Sclr* x, const Int ldx,
  const Int nrhs, Sclr* y, const Int ldy)
{
  for (Int g = 0; ; ) {
    size_t k = 0;
    for (Int i = 0; i < nr; ++i) {
      Sclr a = 0;
      for (Int j = 0; j < nc; ++j, ++k) a += d[k]*x[j];
      y[i] -= a;
    }
    if (++g == nrhs) break;
    x += ldx;
    y += ldy;
  }
}

#if defined(HAVE_SHYLU_NODEHTS_BLAS) || defined(HAVE_SHYLU_NODEHTS_MKL)
template<> inline void SerialBlock_n1Axpy_dense (
  const blas_int nr, const blas_int nc, const float* d, const float* x,
  const blas_int ldx, const blas_int nrhs, float* y, const blas_int ldy)
{ gemm<float>('t', 'n', nr, nrhs, nc, -1, d, nc, x, ldx, 1, y, ldy); }

template<> inline void SerialBlock_n1Axpy_dense (
  const blas_int nr, const blas_int nc, const double* d, const double* x,
  const blas_int ldx, const blas_int nrhs, double* y, const blas_int ldy)
{ gemm<double>('t', 'n', nr, nrhs, nc, -1, d, nc, x, ldx, 1, y, ldy); }

#ifdef HAVE_SHYLU_NODEHTS_COMPLEX
template<> inline void SerialBlock_n1Axpy_dense (
  const blas_int nr, const blas_int nc, const std::complex<float>* d,
  const std::complex<float>* x, const blas_int ldx, const blas_int nrhs,
  std::complex<float>* y, const blas_int ldy)
{
  gemm<std::complex<float> >('t', 'n', nr, nrhs, nc, -1, d, nc, x, ldx, 1,
                             y, ldy);
}

template<> inline void SerialBlock_n1Axpy_dense (
  const blas_int nr, const blas_int nc, const std::complex<double>* d,
  const std::complex<double>* x, const blas_int ldx, const blas_int nrhs,
  std::complex<double>* y, const blas_int ldy)
{
  gemm<std::complex<double> >('t', 'n', nr, nrhs, nc, -1, d, nc, x, ldx, 1,
                              y, ldy);
}
#endif // HAVE_SHYLU_NODEHTS_COMPLEX
#endif // HAVE_SHYLU_NODEHTS_BLAS || HAVE_SHYLU_NODEHTS_MKL

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
SerialBlock::n1Axpy (const Sclr* x, const Int ldx, const Int nrhs,
                     Sclr* y, const Int ldy) const {
  if ( ! d_) return;
  if (ir_) {
    assert(ir_);
    if (ir_[nr_] == 0) return;
    SerialBlock_n1Axpy_spars(nr_, nc_, ir_, jc_, d_, x + coff_, ldx, nrhs,
                             y + roff_, ldy);
  } else {
    SerialBlock_n1Axpy_dense(nr_, nc_, d_, x + coff_, ldx, nrhs,
                             y + roff_, ldy);
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
TMatrix::n1Axpy (const Sclr* x, const Int ldx, const Int nrhs, Sclr* y,
                 const Int ldy, const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return;
  if ((size_t) id >= bs_.size()) return;
  bs_[id].n1Axpy(x + coff_, ldx, nrhs, y + ros_[id], ldy);
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline Sclr* Impl<Int, Size, Sclr>::
Permuter::from_outside (const Sclr* x, const Int nrhs, Int ldx) const {
  if ( ! ldx) ldx = n_;
  const int tid = omp_get_thread_num();
  const Int i0 = part_[tid], i1 = part_[tid+1];
  Sclr* ppx = px_;
  const Sclr* px = x;
  for (int k = 0; ; ) {
    if (scale_)
      for (Int i = i0; i < i1; ++i) {
        const Int pi = p_[i];
        ppx[i] = px[pi] / scale_[pi];
      }
    else
      for (Int i = i0; i < i1; ++i)
        ppx[i] = px[p_[i]];
    if (++k == nrhs) break;
    ppx += n_;
    px += ldx;
  }
  return px_;
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::Permuter::
to_outside (Sclr* x, const Int nrhs, const Sclr a, const Sclr b,
            Int ldx) const {
  if ( ! ldx) ldx = n_;
  const int tid = omp_get_thread_num();
  const Int i0 = part_[tid], i1 = part_[tid+1];
  Sclr* ppx = px_;
  Sclr* px = x;
  if (b != Sclr{1.0}) {
    if (a != Sclr{0.0}) {
      for (int k = 0; ; ) {
        for (Int i = i0; i < i1; ++i) {
          Sclr* const pxqi = px + q_[i];
          *pxqi = a * *pxqi + b*ppx[i];
        }
        if (++k == nrhs) break;
        ppx += n_;
        px += ldx;
      }
    } else {
      for (int k = 0; ; ) {
        for (Int i = i0; i < i1; ++i)
          px[q_[i]] = b*ppx[i];
        if (++k == nrhs) break;
        ppx += n_;
        px += ldx;
      }
    }
  } else {
    if (a != Sclr{0.0}) {
      for (int k = 0; ; ) {
        for (Int i = i0; i < i1; ++i) {
          Sclr* const pxqi = px + q_[i];
          *pxqi = a * *pxqi + ppx[i];
        }
        if (++k == nrhs) break;
        ppx += n_;
        px += ldx;
      }
    } else {
      for (int k = 0; ; ) {
        for (Int i = i0; i < i1; ++i)
          px[q_[i]] = ppx[i];
        if (++k == nrhs) break;
        ppx += n_;
        px += ldx;
      }
    }
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
LevelSetTri::solve (const Sclr* b, Sclr* x, const Int ldx,
                    const Int nrhs) const {
  if (nlvls_ == 0) return;
  const int tid = omp_get_thread_num();
  const Thread& t = t_[tid];
  const Array<Int>& lsp = t.lsp;
  const CrsMatrix* const T = t.m;
  const Size* const ir = T->ir;
  const Int* const jc = T->jc;
  const Sclr* const d = T->d;
  const Array<Int>& p = t.p;
  const Int lsp_size_m2 = nlvls_ - 1;
  const Array<Int>& p2p_depends_p = t.p2p_depends_p;
  const Array<Size>& p2p_depends = t.p2p_depends;
  p2p_Done p2p_done_value = p2p_done_value_;
  for (Int irhs = 0; irhs < nrhs; ++irhs) {
    Int
      ils = 0,      // level set index
      i = lsp[ils]; // this thread's current row
    Size j = ir[i]; // the usual j index into jc and d
    for ( ; ; ++ils) {
      const Int lsp_ilsp1 = lsp[ils+1];
      if (i == lsp_ilsp1) {
        if (ils == lsp_size_m2) break;
        continue;
      }
      // Synchronize.
      if (ils > 0) { // First level doesn't have dependencies.
        // Range in p2p_depends.
        const Int
          p2p_depends_0 = p2p_depends_p[ils-1],
          p2p_depends_1 = p2p_depends_p[ils];
        for (Int di = p2p_depends_0; di < p2p_depends_1; ++di) {
          volatile p2p_Done* const done = &p2p_done_[p2p_depends[di]];
          while (*done != p2p_done_value) ;
        }
      }
      // Flush (acquire) the global to the local view of x. With OpenMP 5.0, we
      // can specify acquire, release, or both, but for now we don't have 5.0.
#ifdef _OPENMP
#     pragma omp flush
#endif
      for (Int r = mvp_block_nc_ + p[ils];
           i < lsp_ilsp1;
           ++i, ++r) {
        const Size jlim = ir[i+1] - 1;
        Sclr a = b[r];
        for ( ; j < jlim; ++j)
          a -= x[jc[j]] * d[j];
        x[r] = a * d[j++];
      }
      // Flush (release) the local to the global view of x.
#ifdef _OPENMP
#     pragma omp flush
#endif
      if (ils == lsp_size_m2) break;
      // This thread and level is done.
      p2p_done_[ls_p2p_sub2ind(ils, tid)] = p2p_done_value;
    }
    x += ldx;
    b += ldx;
    // Increment the done indicator.
    ++p2p_done_value;
#ifdef _OPENMP
#   pragma omp barrier
#endif
  }
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
rbwait (volatile p2p_Done* const s_done, const Size* s_ids,
        const Int* const s_idx, const Int i, const p2p_Done done_symbol) {
  const Int si = s_idx[i], si1 = s_idx[i+1];
  if (si == si1) {
    // acquire
#ifdef _OPENMP
#   pragma omp flush
#endif
    return;
  }
  const Size* id = s_ids + si;
  const Size* const idn = s_ids + si1;
  while (id != idn) {
    volatile p2p_Done* const d = s_done + *id;
    while (*d != done_symbol) ;
    ++id;
  }
  // acquire
#ifdef _OPENMP
# pragma omp flush
#endif
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::RecursiveTri::
ondiag_solve (const OnDiagTri& t, Sclr* x, const Int ldx, const Int nrhs,
              const int tid, const Int step, volatile Int* const t_barrier,
              volatile Int* const inv_tri_done) const {
  const Int nthreads = t.nthreads();
  if (nthreads == 1) {
    t.solve(x, ldx, x, ldx, nrhs);
    // release
#ifdef _OPENMP
#   pragma omp flush
#endif
    *t_barrier = step;
  } else {
    // Solve T wrk_ = x.
    t.solve(x, ldx, wrk_.data(), t.get_n(), nrhs);
    // release
#ifdef _OPENMP
#   pragma omp flush
#endif
    { // Wait for the block row MVPs to finish.
      const Int done = (step << 1);
      inv_tri_done[tid] = done;
      for (Int i = 0; i < nthreads; ++i)
        while (inv_tri_done[i] < done) ;
    }
    // acquire
#ifdef _OPENMP
#   pragma omp flush
#endif
    // Copy wrk_ to x.
    const Int row_start = t.block_row_start(tid), nr = t.block_nr(tid);
    for (Int irhs = 0; irhs < nrhs; ++irhs)
      memcpy(x + irhs*ldx + row_start,
             wrk_.data() + irhs*t.get_n() + row_start,
             nr*sizeof(Sclr));
    // release
#ifdef _OPENMP
#   pragma omp flush
#endif
    { // Wait for the memcpy's to finish.
      const Int done = (step << 1) + 1;
      inv_tri_done[tid] = done;
      //todo Not every thread necessarily needs this on-diag tri's solution, but
      // our dep graph doesn't encode that yet.
      for (Int i = 0; i < nthreads; ++i)
        while (inv_tri_done[i] < done) ;
    }
    if (tid == 0)
      *t_barrier = step;
  }
}

// inside || {}
template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
RecursiveTri::solve (const Sclr* b, Sclr* x, const Int ldx,
                     const Int nrhs) const {
  if (nd_.t.empty()) return;
  assert(x == b);
  const int tid = omp_get_thread_num();
  Int os = 0;
  const Int sn = static_cast<Int>(nd_.s.size());
  Sclr* x_osi, * x_os;
  volatile Int* const t_barrier = &nd_.t_barrier;
  volatile Int* const inv_tri_done = nd_.inv_tri_done.data();
  volatile p2p_Done* const s_done = nd_.s_done.data();
  const Size* const t_ids = nd_.t_ids.data();
  const Int* const t_idx = nd_.t_idx.data();
  const Size* const s_ids = nd_.s_ids[tid].empty() ? 0 : nd_.s_ids[tid].data();
  const Int* const s_idx = nd_.s_idx[tid].empty() ? 0 : nd_.s_idx[tid].data();
  const p2p_Done done_symbol = nd_.done_symbol;
  {
    const OnDiagTri& t = nd_.t[0];
    if (tid < t.nthreads())
      ondiag_solve(t, x, ldx, nrhs, tid, 0, t_barrier, inv_tri_done);
  }
  if ( ! nd_.os.empty()) {
    os += nd_.t[0].get_n();
    x_osi = x + nd_.os[0];
    x_os = x + os;
  } else {
    if (tid) while (*t_barrier < 0) ;
    return;
  }
  for (Int i = 0; i < sn; ++i) {
    const TMatrix& s = nd_.s[i];
    const OnDiagTri& t = nd_.t[i+1];
    if (tid) while (*t_barrier < i) ;
    if ( ! s.empty() && (s.parallel() || tid == 0)) {
      rbwait(s_done, s_ids, s_idx, i, done_symbol);
      s.n1Axpy(x_osi, ldx, nrhs, x_os, ldx, tid);
      // release
#ifdef _OPENMP
#     pragma omp flush
#endif
      s_done[rb_p2p_sub2ind(i, tid)] = done_symbol;
    }
    if (tid < t.nthreads()) {
      rbwait(s_done, t_ids, t_idx, i, done_symbol);
      ondiag_solve(t, x_os, ldx, nrhs, tid, i+1, t_barrier, inv_tri_done);
    }
    os += t.get_n();
    x_osi = x + nd_.os[i+1];
    x_os = x + os;
  }
}

template<typename Int, typename Size, typename Sclr>
inline void Impl<Int, Size, Sclr>::
TriSolver::solve (const Sclr* b, const Int nrhs, Sclr* x, const Sclr alpha,
                  const Sclr beta, const Int ldb, const Int ldx) const {
  t_.p2p_reset();
  lst_.p2p_reset();
  perm_.check_nrhs(nrhs);
  NumThreads nthreads_save;
  set_num_threads(nthreads_, nthreads_save);
  bool nthreads_ok;
  std::string msg;
#ifdef _OPENMP
# pragma omp parallel
#endif
  do {
    nthreads_ok = check_nthreads(nthreads_, omp_get_num_threads(), msg);
    if ( ! nthreads_ok) break;
    Sclr* px = perm_.from_outside(b, nrhs, ldb);
    if (is_lo_) {
#ifdef _OPENMP
#     pragma omp barrier
#endif
      lst_.solve(px, px, n_, nrhs);
      // No barrier needed here because lst_.solve does it.
      if (t_.get_n()) t_.solve(px, px, n_, nrhs);
#ifdef _OPENMP
#     pragma omp barrier
#endif
      perm_.to_outside(x, nrhs, alpha, beta, ldx);
    } else {
      if (t_.get_n()) {
#ifdef _OPENMP
#       pragma omp barrier
#endif
        t_.solve(px, px, n_, nrhs);
      }
#ifdef _OPENMP
#     pragma omp barrier
#endif
      lst_.solve(px, px, n_, nrhs);
      // No barrier needed because of lst_.solve.
      perm_.to_outside(x, nrhs, alpha, beta, ldx);
    }
  } while (0);
  restore_num_threads(nthreads_save);
  if ( ! nthreads_ok) throw hts::Exception(msg);
}

//< Solve code.

// Solve T x = x with size(x, 2) = nrhs. If is_lo, T is lower tri, else upper.
template<typename Int, typename Size, typename Sclr>
void trisolve_serial (
  const typename HTS<Int, Size, Sclr>::CrsMatrix& T, Sclr* ix, const Int nrhs,
  const bool is_lo, const bool unit_diag, const Int* p, const Int* q,
  const typename hts::ScalarTraits<Sclr>::Real* scale, Sclr* w)
{
  assert(( ! p && ! q) || w);

  const int os = unit_diag ? 0 : 1;
  for (int irhs = 0; irhs < nrhs; ++irhs) {
    // x = x./scale optionally.
    if (scale)
      for (Int r = 0; r < T.m; ++r) ix[r] /= scale[r];

    // x = x(p) optionally.
    Sclr* x;
    if (p) {
      x = w;
      for (Int r = 0; r < T.m; ++r) x[r] = ix[p[r]];
    } else
      x = ix;

    // Solve.
    if (T.dir == Direction::forward) {
      if (T.conj) {
        if (is_lo) {
          for (Int r = 0; r < T.m; ++r) {
            const Size rp_rp1 = T.ir[r+1];
            for (Size j = T.ir[r]; j < rp_rp1 - os; ++j) {
              Sclr dj(T.d[j]);
              conjugate(dj);
              x[r] -= x[T.jc[j]] * dj;
            }
            if ( ! unit_diag) {
              Sclr diag(T.d[rp_rp1 - 1]);
              conjugate(diag);
              x[r] /= diag;
            }
          }
        } else {
          for (Int r = T.m - 1; r >= 0; --r) {
            const Size rp_r = T.ir[r];
            for (Size j = rp_r + os; j < T.ir[r+1]; ++j) {
              Sclr dj(T.d[j]);
              conjugate(dj);
              x[r] -= x[T.jc[j]] * dj;
            }
            if ( ! unit_diag) {
              Sclr diag(T.d[rp_r]);
              conjugate(diag);
              x[r] /= diag;
            }
          }
        }
      } else {
        if (is_lo) {
          for (Int r = 0; r < T.m; ++r) {
            const Size rp_rp1 = T.ir[r+1];
            for (Size j = T.ir[r]; j < rp_rp1 - os; ++j)
              x[r] -= x[T.jc[j]] * T.d[j];
            if ( ! unit_diag)
              x[r] /= T.d[rp_rp1 - 1];
          }
        } else {
          for (Int r = T.m - 1; r >= 0; --r) {
            const Size rp_r = T.ir[r];
            for (Size j = rp_r + os; j < T.ir[r+1]; ++j)
              x[r] -= x[T.jc[j]] * T.d[j];
            if ( ! unit_diag)
              x[r] /= T.d[rp_r];
          }
        }
      }
    } else { // T.dir == Direction::transpose; treat as CCS
      if (T.conj) {
        if (is_lo) {
          for (Int c = T.m - 1; c >= 0; --c) {
            const Size js = T.ir[c], je = T.ir[c+1];
            if ( ! unit_diag) {
              Sclr diag(T.d[je - 1]);
              conjugate(diag);
              x[c] /= diag;
            }
            const Sclr& xc = x[c];
            for (Size j = js; j < je - os; ++j) {
              const Int r = T.jc[j];
              Sclr dj(T.d[j]);
              conjugate(dj);
              x[r] -= dj*xc;
            }
          }
        } else {
          for (Int c = 0; c < T.n; ++c) {
            const Size js = T.ir[c], je = T.ir[c+1];
            if ( ! unit_diag) {
              Sclr diag(T.d[js]);
              conjugate(diag);
              x[c] /= diag;
            }
            const Sclr& xc = x[c];
            for (Size j = js + os; j < je; ++j) {
              const Int r = T.jc[j];
              Sclr dj(T.d[j]);
              conjugate(dj);
              x[r] -= dj*xc;
            }
          }
        }
      } else {
        if (is_lo) {
          for (Int c = T.m - 1; c >= 0; --c) {
            const Size js = T.ir[c], je = T.ir[c+1];
            if ( ! unit_diag) x[c] /= T.d[je - 1];
            const Sclr& xc = x[c];
            for (Size j = js; j < je - os; ++j) {
              const Int r = T.jc[j];
              x[r] -= T.d[j]*xc;
            }
          }
        } else {
          for (Int c = 0; c < T.n; ++c) {
            const Size js = T.ir[c], je = T.ir[c+1];
            if ( ! unit_diag) x[c] /= T.d[js];
            const Sclr& xc = x[c];
            for (Size j = js + os; j < je; ++j) {
              const Int r = T.jc[j];
              x[r] -= T.d[j]*xc;
            }
          }
        }
      }
    }

    // x(q) = x optionally.
    if (q) {
      if ( ! p) memcpy(w, x, T.m*sizeof(*ix));
      for (Int r = 0; r < T.m; ++r) ix[q[r]] = w[r];
    } else if (p)
      memcpy(ix, w, T.m*sizeof(*ix));

    ix += T.m;
  }
}

} // namespace htsimpl

template<typename Int, typename Size, typename Sclr>
typename HTS<Int, Size, Sclr>::CrsMatrix* HTS<Int, Size, Sclr>::
make_CrsMatrix (const Int nrow, const Size* rowptr, const Int* col,
                const Sclr* val, const bool make_transpose,
                const bool make_conj) {
  static const bool Int_is_signed = static_cast<Int>(-1) < 1;
  if ( ! Int_is_signed)
    throw hts::Exception(
      "HTS<Int, Size, Scalar> cannot be used if Int is unsigned.");
  CrsMatrix* m;
  const htsimpl::Direction::Enum dir = (make_transpose ?
                                        htsimpl::Direction::transpose :
                                        htsimpl::Direction::forward);
  try {
    m = new CrsMatrix(nrow, rowptr, col, val, dir, make_conj);
  } catch (...) { throw hts::Exception("make_CrsMatrix failed to alloc."); }
  return m;
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
delete_CrsMatrix (CrsMatrix* T) { htsimpl::del(T); }

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
register_Deallocator (CrsMatrix* T, Deallocator* d) {
  d->counter++;
  T->deallocator = d;
}

template<typename Int, typename Size, typename Sclr>
HTS<Int, Size, Sclr>::
PreprocessArgs::PreprocessArgs ()
  : T(0), max_nrhs(1), nthreads(-1), save_for_reprocess(false), p(0), q(0),
    scale_rhs(0), scale_rhs_by_division(true), scale_solution(0),
    scale_solution_by_division(true), options(0)
{}

template<typename Int, typename Size, typename Sclr>
typename HTS<Int, Size, Sclr>::Impl* HTS<Int, Size, Sclr>::
preprocess (const PreprocessArgs& a) {
  if ( ! a.T)
    throw hts::Exception("T is null.");
  if (a.nthreads <= 0)
    throw hts::Exception("nthreads must be > 0.");
  if (a.scale_solution)
    throw hts::Exception("scale_solution is not impl'ed yet.");
  if (a.scale_rhs && ! a.scale_rhs_by_division)
    throw hts::Exception(
      "scale_rhs && !scale_rhs_by_division is not impl'ed yet");
  Impl* i;
  try { i = new Impl(); }
  catch (...) { throw hts::Exception("preprocess failed to alloc."); }
  if (a.options)
    htsimpl::Impl<Int, Size, Sclr>::set_options(*a.options, i->o);
  try {
    i->ts.init(a.T, a.nthreads, a.max_nrhs, a.save_for_reprocess, a.p, a.q,
               a.scale_rhs, i->o);
  } catch (hts::Exception& e) {
    htsimpl::del(i);
    throw;
  }
  return i;
}

template<typename Int, typename Size, typename Sclr>
typename HTS<Int, Size, Sclr>::Impl* HTS<Int, Size, Sclr>::
preprocess (const CrsMatrix* T, const Int max_nrhs, const Int nthreads,
            const bool save_for_reprocess, const Int* p, const Int* q,
            const Real* r, const Options* options)
{
  PreprocessArgs a;
  a.T = T;
  a.max_nrhs = max_nrhs;
  a.nthreads = nthreads;
  a.save_for_reprocess = save_for_reprocess;
  a.p = p;
  a.q = q;
  a.scale_rhs = r;
  a.scale_rhs_by_division = true;
  a.options = options;
  return preprocess(a);
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
reprocess_numeric (Impl* impl, const CrsMatrix* T, const Real* r) {
  impl->ts.reinit_numeric(T, r);
}

template<typename Int, typename Size, typename Sclr>
bool HTS<Int, Size, Sclr>::
is_lower_tri (const Impl* impl) { return impl->ts.is_lower_tri(); }

template<typename Int, typename Size, typename Sclr>
bool HTS<Int, Size, Sclr>::
has_implicit_unit_diag (const Impl* impl)
{ return impl->ts.has_implicit_unit_diag(); }

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
delete_Impl (Impl* impl) {
  htsimpl::del(impl);
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
reset_max_nrhs (Impl* impl, const Int max_nrhs) {
  impl->ts.reset_max_nrhs(max_nrhs);
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
solve_serial (const CrsMatrix* T, const bool is_lo, const bool unit_diag,
              Sclr* xb, const Int nrhs, const Int* p, const Int* q,
              const typename hts::ScalarTraits<Sclr>::Real* r, Sclr* w) {
  htsimpl::trisolve_serial<Int, Size, Sclr>(*T, xb, nrhs, is_lo, unit_diag,
                                            p, q, r, w);
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
solve_omp (Impl* impl, Sclr* x, const Int nrhs, const Int ldx) {
  impl->ts.solve(x, nrhs, x, 0, 1, ldx, ldx);
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
solve_omp (Impl* impl, const Sclr* b, const Int nrhs, Sclr* x,
           const Int ldb, const Int ldx) {
  impl->ts.solve(b, nrhs, x, 0, 1, ldb, ldx);
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
solve_omp (Impl* impl, const Sclr* b, const Int nrhs, Sclr* x,
           const Sclr alpha, const Sclr beta,
           const Int ldb, const Int ldx) {
  impl->ts.solve(b, nrhs, x, alpha, beta, ldb, ldx);
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
print_options (const Impl* i, std::ostream& os) {
  htsimpl::Impl<Int, Size, Sclr>::print_options(i->o, os);
}

template<typename Int, typename Size, typename Sclr>
void HTS<Int, Size, Sclr>::
set_level_schedule_only (Options& o) {
  o.min_lset_size = 0;
}

template<typename Int, typename Size, typename Sclr>
HTS<Int, Size, Sclr>::Options::Options () {
  min_dense_density = 0.75;
  levelset_block_size = 1;
  lset_min_size_scale_with_nthreads = false;
  profile = false;
  print_level = 0;
  lset_max_bad_fraction = 0.01;
  min_lset_size = lset_min_size_scale_with_nthreads ? 1 : 10;
  min_block_size = 256;
  min_parallel_rows = 64;
  pp_min_block_size = 256;
}

} // namespace Experimental

#endif // INCLUDE_SHYLU_HTS_IMPL_DEF_HPP
