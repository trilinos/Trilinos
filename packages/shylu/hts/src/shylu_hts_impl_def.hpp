#ifndef INCLUDE_SHYLU_HTS_IMPL_DEF_HPP
#define INCLUDE_SHYLU_HTS_IMPL_DEF_HPP

#include <omp.h>
#ifdef HAVE_SHYLUHTS_MKL
# include <mkl.h>
#endif

#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>

#include <Teuchos_BLAS_wrappers.hpp>

#include "shylu_hts.hpp"
#include "shylu_hts_impl.hpp"

#ifdef __MIC__
# include <xmmintrin.h>
#endif

//#define TIME
#ifdef TIME
# define TIMENUM
# include <sys/time.h>
#endif

// For analysis and figures, output some data.
#if ! defined __MIC__
//# define OUTPUT_RECURSIVEBLOCK_PARTITION
#endif

namespace Experimental {
namespace htsimpl {

static const int parfor_static_size = 20;

template<typename T>
inline void compress (std::vector<T>& v) { std::vector<T>(v).swap(v); }

#if defined(HAVE_SHYLUHTS_BLAS) || defined(HAVE_SHYLUHTS_MKL)
typedef int blas_int;

template<typename T> void gemm(
  char transa, char transb, blas_int m, blas_int nrhs, blas_int n, T alpha,
  const T* a, blas_int lda, const T* b, blas_int ldb, T beta,
  const T* c, blas_int ldc);

// Calling Teuchos::BLAS<Int, Real>::GEMM is measurably slower, so use
// Teuchos_BLAS_wrappers directly. I'm copying the following macro from
// Teuchos_BLAS.cpp. Perhaps one should be exposed in Teuchos_BLAS.hpp or
// Teuchos_BLAS_wrappers.hpp.
#ifdef INTEL_CXML
# define SHYLU_HTS_CHAR_MACRO(char_var) &char_var, 1
#else
# define SHYLU_HTS_CHAR_MACRO(char_var) &char_var
#endif

template<> inline void
gemm<float> (
  char transa, char transb, blas_int m, blas_int nrhs, blas_int n, float alpha,
  const float* a, blas_int lda, const float* b, blas_int ldb, float beta,
  const float* c, blas_int ldc)
{
  SGEMM_F77(SHYLU_HTS_CHAR_MACRO(transa), SHYLU_HTS_CHAR_MACRO(transb),
            &m, &nrhs, &n, &alpha, const_cast<float*>(a), &lda,
            const_cast<float*>(b), &ldb, &beta, const_cast<float*>(c), &ldc);
}

template<> inline void
gemm<double> (
  char transa, char transb, blas_int m, blas_int nrhs, blas_int n, double alpha,
  const double* a, blas_int lda, const double* b, blas_int ldb, double beta,
  const double* c, blas_int ldc)
{
  DGEMM_F77(SHYLU_HTS_CHAR_MACRO(transa), SHYLU_HTS_CHAR_MACRO(transb),
            &m, &nrhs, &n, &alpha, const_cast<double*>(a), &lda,
            const_cast<double*>(b), &ldb, &beta, const_cast<double*>(c), &ldc);
}
#endif

#ifdef HAVE_SHYLUHTS_MKL
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
#endif

namespace {
class Timer {
public:
  enum Op { total_pre = 0,
            setup, transpose, tolower, perm,
            lsetfind, lsetinit, lsetinitp2p,
            lsp2p_1, lsp2p_3, lsp2p_6,
            dpb, dpb_getmatrix, dpb_sort, dpb_tinit,
            dpb_tinit_2, dpb_tinit_3,
            for_reprocess,
            NSETUPTIMERS,
            total_re, numthr, numpart, numls, numrbt, numrbm, numperm,
            NNUMTIMERS,
            slvlls, slvlother, slvuls, slvuother,
            NTIMERS };
  static inline void init () {
#ifdef TIME
    for (int i = 0; i < NTIMERS; ++i) et_[i] = 0;
#endif
  }
  static inline void start (const Op op) {
#ifdef TIME
    gettimeofday(&t_start_[op], 0);
#endif
  }
  static inline void stop (const Op op) {
#ifdef TIME
    timeval t2;
    gettimeofday(&t2, 0);
    const timeval& t1 = t_start_[op];
    static const double us = 1.0e6;
    et_[op] += (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
#endif
  }
# define tpr(op) do {                                                   \
    printf("%20s %10.3e %10.1f\n", #op, et_[op], 100*et_[op]/tot);      \
  } while (0)
#ifdef TIME
  static void print_setup () {
    const double tot = et_[total_pre];
    tpr(setup); tpr(transpose); tpr(tolower);
    tpr(lsetfind);
    tpr(perm);
    tpr(lsetinit);
    tpr(lsetinitp2p); tpr(lsp2p_1); tpr(lsp2p_3); tpr(lsp2p_6);
    tpr(dpb); tpr(dpb_getmatrix); tpr(dpb_sort); tpr(dpb_tinit);
    tpr(dpb_tinit_2); tpr(dpb_tinit_3);
    printf("%20s %10.3e %10.1f\n", "total", et_[total_pre], 100.0);
  }
#endif
#ifdef TIMENUM
  static void print_numeric () {
    const double tot = et_[total_re];
    tpr(numthr); tpr(numpart); tpr(numls); tpr(numrbt); tpr(numrbm);
    tpr(numperm);
    printf("%20s %10.3e %10.1f\n", "total", et_[total_re], 100.0);
  }
#endif
#undef tpr
private:
#ifdef TIME
  static timeval t_start_[NTIMERS];
  static double et_[NTIMERS];
#endif
};
#ifdef TIME
timeval Timer::t_start_[Timer::NTIMERS];
double Timer::et_[Timer::NTIMERS];
#endif
} // namespace

template<typename T> inline void touch (T* const p, const size_t n) {
  // 1 KB should be a safe lower bound on page size. Touch enough to touch every
  // page; I don't think there's any need to touch more memory than that. On
  // the KNC, first-touch doesn't matter.
#ifndef __MIC__
  for (size_t i = 0; i < n; i += 1024 / sizeof(T))
    p[i] = 0;
  // Make sure the last part is touched.
  if (n) p[n-1] = 0;
#endif
}

#ifdef __MIC__
template<typename T> inline T*
allocn (const size_t n, const bool first_touch = false) {
  if ( ! n) return 0;
  T* p = (T*) _mm_malloc(n*sizeof(T), 64);
  if ( ! p) throw std::bad_alloc();
  if (first_touch) touch(p, n);
  return p;
}
template<typename T> inline void deln (T*& p) {
  if (p) _mm_free(p);
  p = 0;
}
template<typename T> inline void deln_const (const T* p) {
  if (p) _mm_free(const_cast<T*>(p));
}
#else
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
#endif
template<typename T> inline void del (T*& p) {
  if (p) delete p;
  p = 0;
}

// For exception safety when allocating.
template<typename T> class Allocnator {
  T* p_;
  bool dealloc_;
public:
  // Try to allocate memory.
  Allocnator (const size_t n, const char* msg, const bool first_touch = false) {
    p_ = 0;
    dealloc_ = true;
    try { p_ = allocn<T>(n, first_touch); }
    catch (...) {
      throw hts::Exception(std::string(msg) + ": failed to allocate.");
    }
  }
  // Release the pointer to the user and subsequently don't dealloc.
  T* release () { dealloc_ = false; return p_; }
  // Dealloc only if the user hasn't released the pointer yet.
  ~Allocnator () { if (dealloc_) deln<T>(p_); }
};

template<typename T> inline T square (const T& x) { return x*x; }

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real> Impl<Int, Size, Real>::
Options::Options () {
  set_options(typename HTS<Int, Size, Real>::Options(), *this);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::Options::print (std::ostream& os) const {
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

static void print_compiletime_options(std::ostream& os) {
#ifdef HAVE_SHYLUHTS_BLAS
  os << " HAVE_SHYLUHTS_BLAS";
#endif
#ifdef HAVE_SHYLUHTS_MKL
  os << " HAVE_SHYLUHTS_MKL";
#endif
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::print_options (const Options& o, std::ostream& os) {
  print_compiletime_options(os);
  os << std::endl;
  o.print(os);
  os << std::endl;
}

template<typename Int, typename Size, typename Real> Impl<Int, Size, Real>::
ConstCrsMatrix::~ConstCrsMatrix () {
  if ( ! deallocate_) return;
  deln_const(ir); deln_const(jc); deln_const(d);
}

template<typename Int, typename Size, typename Real> Impl<Int, Size, Real>::
CrsMatrix::~CrsMatrix () {
  deln_const(ir); deln_const(jc); deln_const(d);
}

template<typename Int, typename Size, typename Real>
inline typename Impl<Int, Size, Real>::ConstCrsMatrix::Direction
opposite (const typename Impl<Int, Size, Real>::ConstCrsMatrix::Direction dir) {
  return static_cast<typename Impl<Int, Size, Real>::ConstCrsMatrix::Direction>(
    (dir + 1) % 2);
}

template<typename T> static T* vec2arr (const std::vector<T>& v) {
  if (v.empty()) return 0;
  T* a = Allocnator<T>(v.size(), "vec2arr").release();
  memcpy(a, v.data(), sizeof(T)*v.size());
  return a;
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::Partition::alloc_d () {
  assert( ! cm->d);
  cm->d = Allocnator<Real>(cm->ir[cm->m], "Partition::alloc_d").release();
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::Partition::clear () { del(cm); }

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::Partition::clear_d () { deln(cm->d); }

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
partition_n_uniformly (const Int n, const Int nparts, std::vector<Int>& p) {
  p.resize(nparts + 1);
  const Int base = n / nparts;
  Int rem = n - base*nparts;
  Int extra = rem > 0 ? 1 : 0;
  p[nparts] = 0;
  for (Int i = 1; i <= nparts; ++i) {
    p[i] = p[i-1] + base + extra;
    if (rem > 0) {
      --rem;
      if (rem == 0) extra = 0;
    }
  }
}

template<typename Int, typename Size, typename Real>
Impl<Int, Size, Real>::
SparseData::SparseData (const Int m, const Size nnz, const char* fail_msg,
                        const bool touch) {
  ir = 0;
  jc = 0;
  d = 0;
  dealloc_ = true;
  try {
    ir = allocn<Size>(m+1, touch);
    if (nnz > 0) {
      jc = allocn<Int>(nnz, touch);
      d = allocn<Real>(nnz, touch);
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::SparseData::free () {
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
#ifdef HAVE_SHYLUHTS_MKL
  save.mkl = mkl_get_max_threads();
  save.mkl_dynamic = mkl_get_dynamic();
#endif
  omp_set_num_threads(nthreads);
#ifdef HAVE_SHYLUHTS_MKL
  // We never use MKL threading.
  mkl_set_dynamic(0);
  mkl_set_num_threads(1);
#endif
}

inline void restore_num_threads (const NumThreads& save) {
#ifdef HAVE_SHYLUHTS_MKL
  mkl_set_dynamic(save.mkl_dynamic);
#endif
  omp_set_num_threads(save.omp);
#ifdef HAVE_SHYLUHTS_MKL
  mkl_set_num_threads(save.mkl);
#endif
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

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
throw_if_nthreads_not_ok (const int nthreads) {
  int nr;
# pragma omp parallel
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
template<typename Int, typename Size, typename Real>
inline Int Impl<Int, Size, Real>::
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
template<typename Int, typename Size, typename Real>
Size Impl<Int, Size, Real>::
crop_matrix (const CrsMatrix& T, const Box& b, Box& cb) {
  cb.r0 = -1;
  Int r1 = -1;
  cb.c0 = b.c0 + b.nc;
  Int c1 = b.c0;
  Size nnz = 0;
  //todo ||ize if not nested.
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
template<typename Int, typename Size, typename Real>
Int Impl<Int, Size, Real>::
decide_level_set_max_index (const std::vector<Int>& N, const Int size_thr,
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
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
alloc_lsets (
  const Int lsmi, const Int sns, const std::vector<Int>& level,
  const std::vector<Int>& n, typename LevelSetter::LevelSets& lsets)
{
  if (lsmi < 0) return;
  const Int Lm_sr = static_cast<Int>(level.size());
  lsets.resize(lsmi+1);
  for (Int i = 0; i <= lsmi; ++i)
    lsets[i].reserve(sns * n[i]);
  // Load.
  for (Int i = 0; i < Lm_sr; ++i) {
    const Int ilev = level[i];
    if (ilev <= lsmi)
      for (Int j = 0; j < sns; ++j)
        lsets[ilev].push_back(i * sns + j);
  }
}

template<typename Int, typename Size, typename Real>
Int Impl<Int, Size, Real>::
locrsrow_schedule_serial (const ConstCrsMatrix& L, const Int sns,
                          std::vector<Int>& w) {
  // Eq. 18 in Y. Saad's 1989 SIAM J Sci Stat Comput paper.
  Int max_level = -1;
  if (sns == 1) {
    w.resize(L.m, -1);
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
    w.resize(Lm_sr, -1);
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

template<typename Int, typename Size, typename Real>
Int Impl<Int, Size, Real>::
locrsrow_schedule_sns1 (const ConstCrsMatrix& L, std::vector<Int>& w,
                        const Options& o) {
  const Int
    nthreads = omp_get_max_threads(),
    blksz = nthreads*((o.pp_min_block_size + nthreads)/nthreads),
    rows_per_thread = std::max(1, blksz / nthreads);
  if (blksz > L.m)
    return locrsrow_schedule_serial(L, 1, w);
  std::vector<Size> frontier(blksz);
  for (Int i = 0; i < blksz; ++i) frontier[i] = L.ir[i];
  w.resize(L.m);
  Int max_level = -1;
  volatile Int done = -1;
# pragma omp parallel
  {
#   pragma omp for schedule(static, parfor_static_size)
    for (Int i = 0; i < L.m; ++i) w[i] = -1;
    const Size* const ir = L.ir;
    const Int* const jc = L.jc;
    for (Int c = 0; c < L.m; c += blksz) {
      const Int tlim = std::min<Int>(c + blksz, L.m);
      // On-diag serial triangle.
#     pragma omp single nowait
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
      while (done != c) ;
#     pragma omp for schedule(static, rows_per_thread)
      for (Int r = tlim; r < rlim; ++r) {
        Int level = -1;
        const Size jlim = ir[r+1];
        for (Size j = ir[r]; j < jlim; ++j) {
          const Int c = jc[j];
          if (c >= tlim) {
            frontier[r - tlim] = j;
            w[r] = level;
            break;
          }
          level = std::max(level, w[c]);
        }
      }
      // Implied barrier from parfor.
    }
  }
  return max_level;
}

template<typename Int, typename Size, typename Real>
Int Impl<Int, Size, Real>::
locrsrow_schedule (const ConstCrsMatrix& L, const Int sns,
                   std::vector<Int>& w, const Options& o) {
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
  std::vector<Size> frontier(bnr);
  for (Int i = 0; i < bnr; ++i) frontier[i] = L.ir[i];
  w.resize(Lm_sr);
  Int max_level = -1;
# pragma omp parallel
  {
#   pragma omp for schedule(static, parfor_static_size)
    for (Int i = 0; i < Lm_sr; ++i) w[i] = -1;
    const Size* const ir = L.ir;
    const Int* const jc = L.jc;
    for (Int sc = 0; sc < Lm_sr; sc += blksz) {
      const Int stlim = std::min<Int>(sc + blksz, Lm_sr);
      // On-diag serial triangle.
#     pragma omp single nowait
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
#     pragma omp barrier
      const Int
        srlim = std::min<Int>(stlim + blksz, Lm_sr),
        tlim = sns*stlim;
#     pragma omp for schedule(static, rows_per_thread)
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
            const Int sc = c / sns;
            level = std::max(level, w[sc]);
          }
        }
        w[sr] = level;
      }
      // Implied barrier from parfor.
    }
  }
  return max_level;
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
find_row_level_sets_Lcrs (const ConstCrsMatrix& L, const Int sns,
                          Int size_thr, typename LevelSetter::LevelSets& lsets,
                          const Options& o) {
  assert(L.m % sns == 0);

  std::vector<Int> w;
#ifdef __MIC__
  // || is working pretty well on MIC, but not on CPU.
  const Int max_level = locrsrow_schedule(L, sns, w, o);
#else
  const Int max_level = locrsrow_schedule_serial(L, sns, w);
#endif
 
  // Count level set sizes.
  std::vector<Int> n(max_level+1);
  for (size_t i = 0; i < w.size(); ++i)
    ++n[w[i]];
  // Cutoff.
  const Int lsmi = decide_level_set_max_index(n, size_thr, o);
  // Fill lsets.
  alloc_lsets(lsmi, sns, w, n, lsets);
}

// Upper tri, CRS, col (not row) level sets. Equivalent to lower tri, CCS, row
// level sets.
//todo Need to ||ize, but this is equivalent to ||izing a CSC (unanalyzed) L
// trisolve, which is harder than CSR L trisolve. Not sure yet how to do this
// well.
template<typename Int, typename Size, typename Real>
Int Impl<Int, Size, Real>::
upcrscol_schedule_serial (const ConstCrsMatrix& U, const Int sns,
                          std::vector<Int>& w) {
  Int max_level = -1;
  if (sns == 1) {
    w.resize(U.m, -1);
    for (Int r = 0; r < U.m; ++r) {
      ++w[r];
      const Int level = w[r];
      max_level = std::max(level, max_level);
      for (Size j = U.ir[r] + 1; j < U.ir[r+1]; ++j) {
        const Int c = U.jc[j];
        w[c] = std::max(w[c], level);
      }
    }
  } else {
    const Int Um_sr = U.m / sns;
    w.resize(Um_sr, -1);
    for (Int sr = 0, r = 0; sr < Um_sr; ++sr) {
      ++w[sr];
      const Int level = w[sr];
      max_level = std::max(level, max_level);
      for (Int i = 0; i < sns; ++i, ++r)
        for (Size j = U.ir[r] + 1; j < U.ir[r+1]; ++j) {
          const Int sc = U.jc[j] / sns;
          w[sc] = std::max(w[sc], level);
        }
    }
  }
  return max_level;
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
find_col_level_sets_Ucrs (const ConstCrsMatrix& U, const Int sns,
                          Int size_thr, typename LevelSetter::LevelSets& lsets,
                          const Options& o) {
  assert(U.m % sns == 0);

  std::vector<Int> w;
  const Int max_level = upcrscol_schedule_serial(U, sns, w);
 
  // Count level set sizes.
  std::vector<Int> n(max_level+1);
  for (size_t i = 0; i < w.size(); ++i)
    ++n[w[i]];
  // Cutoff.
  const Int lsmi = decide_level_set_max_index(n, size_thr, o);
  // Fill lsets.
  alloc_lsets(lsmi, sns, w, n, lsets);
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
find_level_sets (
  const ConstCrsMatrix& T, const Int sns, const Int size_thr, const bool is_lo,
  typename LevelSetter::LevelSets& lsets, const Options& o)
{
  if (is_lo)
    find_row_level_sets_Lcrs(T, sns, size_thr, lsets, o);
  else
    find_col_level_sets_Ucrs(T, sns, size_thr, lsets, o);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
LevelSetter::init (const ConstCrsMatrix& T, const Int size_thr,
                   const bool is_lo, const Options& o) {
  lsets_.clear();
  is_lo_ = is_lo;
  // Guard against an invalid setting.
  ls_blk_sz_ = T.m % o.ls_blk_sz == 0 ? o.ls_blk_sz : 1;
  find_level_sets(T, ls_blk_sz_, size_thr, is_lo_, lsets_, o);
}

template<typename Int, typename Size, typename Real>
const std::vector<Int>& Impl<Int, Size, Real>::
LevelSetter::lset (const size_t i) const {
  return is_lo_ ? lsets_[i] : lsets_[lsets_.size() - i - 1];
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
LevelSetter::reverse_variable_order (Int n) {
  --n;
  for (size_t i = 0; i < lsets_.size(); ++i) {
    std::vector<Int>& ls = lsets_[i];
    for (size_t j = 0; j < ls.size(); ++j)
      ls[j] = n - ls[j];
    std::reverse(ls.begin(), ls.end());
  }
}

template<typename Int, typename Size, typename Real>
typename Impl<Int, Size, Real>::CrsMatrix* Impl<Int, Size, Real>::
get_matrix_p (const CrsMatrix& A, const std::vector<Int>& p) {
  const Int n = static_cast<Int>(p.size());
  Size nnz = 0;
  for (size_t i = 0; i < p.size(); ++i) {
    const Int r = p[i];
    nnz += A.ir[r+1] - A.ir[r];
  }

  SparseData sd(n, nnz, "get_matrix_p", true);
  for (size_t i = 0; i < p.size(); ++i) {
    const Int r = p[i];
    const Size nc = A.ir[r+1] - A.ir[r];
    sd.ir[i+1] = sd.ir[i] + nc;
    memcpy(sd.jc + sd.ir[i], A.jc + A.ir[r], nc*sizeof(*sd.jc));
    memcpy(sd.d + sd.ir[i], A.d + A.ir[r], nc*sizeof(*sd.d));
  }

  CrsMatrix* cm;
  try { cm = new CrsMatrix(n, A.n, sd.ir, sd.jc, sd.d); }
  catch (...) { throw hts::Exception("get_matrix_p failed to alloc."); }
  sd.release();
  return cm;
}

template<typename Int, typename Size, typename Real>
typename Impl<Int, Size, Real>::ConstCrsMatrix* Impl<Int, Size, Real>::
permute_to_other_tri (const ConstCrsMatrix& U) {
  const Int n = U.m;
  const Size nnz = U.ir[n];
  SparseData sd(n, nnz, "permute_to_other_tri");
# pragma omp parallel
  {
#   pragma omp for schedule(static)
    for (Int k = 1; k <= n; ++k)
      sd.ir[k] = nnz - U.ir[n-k];
#   pragma omp for schedule(static)
    for (Size k = 0; k < nnz; ++k) {
      const Size i = nnz - k - 1;
      sd.jc[k] = n - U.jc[i] - 1;
      sd.d[k] = U.d[i];
    }
  }
  ConstCrsMatrix* ccm = 0;
  try { ccm = new ConstCrsMatrix(n, n, sd.ir, sd.jc, sd.d, U.dir, true); }
  catch (...) { throw hts::Exception("permute_to_other_tri"); }
  sd.release();
  return ccm;
}

// Partition 1:n into lsis, the set of level scheduled rows, and dpis, the set
// of data-|| rows.
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
get_idxs (const Int n, const LevelSetter& lstr, std::vector<Int>& lsis,
          std::vector<Int>& dpis) {
  std::vector<char> dpisb(n, 1);
  for (Int i = 0; i < lstr.size(); ++i) {
    const std::vector<Int>& lset = lstr.lset(i);
    for (size_t j = 0; j < lset.size(); ++j) {
      const Int lj = lset[j];
      lsis.push_back(lj);
      dpisb[lj] = 0;
    }
  }
  dpis.resize(n - lsis.size());
  for (size_t i = 0, dk = 0; i < dpisb.size(); ++i)
    if (dpisb[i]) dpis[dk++] = i;  
}

template<typename Int, typename Size, typename Real>
typename Impl<Int, Size, Real>::Shape Impl<Int, Size, Real>::
determine_shape (const ConstCrsMatrix& A) {
  int red_is_lower = 0, red_is_tri = 0, red_has_full_diag = 0, red_nthreads = 0;
# pragma omp parallel \
         reduction(+: red_is_lower, red_is_tri, red_has_full_diag, red_nthreads)
  {
    bool tid_used = false, has_full_diag = true, is_lower = false,
      is_upper = false;
#   pragma omp for schedule(static, parfor_static_size)
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
    }
    if (tid_used) {
      ++red_nthreads;
      if (has_full_diag) ++red_has_full_diag;
      if ( ! (is_lower && is_upper)) ++red_is_tri;
      if ( ! is_upper) ++red_is_lower;
    }
  }
  const bool
    is_tri = ((// Each thread saw a triangle.
                red_is_tri == red_nthreads)
              &&
              (// Each thread saw the same orientation.
                red_is_lower == 0 || red_is_lower == red_nthreads)),
    // Every thread saw a full diagonal.
    has_full_diag = red_has_full_diag == red_nthreads,
    // Valid only if is_tri.
    is_lower = red_is_lower;
  // If ! tri_determined, then T must be a diag matrix. Can treat as lower,
  // which is is_lower's default value.
  return Shape(is_lower, is_tri, has_full_diag);
}

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

// Extract A(p,p) given that A(p,~p) is empty.
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
get_matrix_pp_with_covers_all (const ConstCrsMatrix& A, const PermVec& pv,
                               Partition& p) {
  // Count nnz.
  Size nnz = 0;
# pragma omp parallel for reduction(+:nnz)
  for (size_t i = 0; i < pv.size(); ++i) {
    const Int k = pv.get(i);
    nnz += A.ir[k+1] - A.ir[k];
  }
  SparseData sd(pv.size(), nnz, "get_matrix_pp_with_covers_all");
  for (size_t ipv = 0; ipv < pv.size(); ++ipv) {
    const Int i = pv.get(ipv);
    const Size nc = A.ir[i+1] - A.ir[i];
    sd.ir[ipv+1] = sd.ir[ipv] + nc;
  }
  assert(sd.ir[pv.size()] == nnz);
  p.A_idxs.resize(nnz);
  const Int ipv_lim = pv.size();
# pragma omp parallel for schedule(static, parfor_static_size)
  for (Int ipv = 0; ipv < ipv_lim; ++ipv) {
    const Int i = pv.get(ipv);
    const Size
      nc = A.ir[i+1] - A.ir[i],
      Aj0 = A.ir[i],
      Bj0 = sd.ir[ipv];
    for (Size j = 0; j < nc; ++j) {
      const Size Aj = Aj0 + j, Bj = Bj0 + j;
      sd.jc[Bj] = pv.to_block(A.jc[Aj]);
      sd.d[Bj] = A.d[Aj];
      p.A_idxs[Bj] = Aj;
    }
  }
  try { p.cm = new CrsMatrix(pv.size(), pv.size(), sd.ir, sd.jc, sd.d); }
  catch (...)
  { throw hts::Exception("get_matrix_pp_with_covers_all failed to alloc."); }
  sd.release();
}

// Extract B = A(p,s), s = [q p], given that A(p,~s) is empty.
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
get_matrix_p_qp_with_covers_all (
  const ConstCrsMatrix& A, const PermVec& pv, const PermVec& qv, Partition& p)
{
  Size nnz = 0;
# pragma omp parallel for reduction(+:nnz)
  for (size_t i = 0; i < pv.size(); ++i) {
    const Int r = pv.get(i);
    nnz += A.ir[r+1] - A.ir[r];
  }

  SparseData sd(pv.size(), nnz, "get_matrix_p_qp_with_covers_all");
  for (size_t ipv = 0, lim = pv.size(); ipv < lim; ++ipv) {
    const Int r = pv.get(ipv);
    const Size nc = A.ir[r+1] - A.ir[r];
    sd.ir[ipv+1] = sd.ir[ipv] + nc;
  }
  assert(sd.ir[pv.size()] == nnz);

  p.A_idxs.resize(nnz);
  const Int ipv_lim = static_cast<Int>(pv.size());
# pragma omp parallel for schedule(static, parfor_static_size)
  for (Int ipv = 0; ipv < ipv_lim; ++ipv) {
    const Int i = pv.get(ipv);
    const Size
      nc = A.ir[i+1] - A.ir[i],
      Aj0 = A.ir[i],
      Bj0 = sd.ir[ipv];
    for (Size j = 0; j < nc; ++j) {
      const Size
        Aj = Aj0 + j, Bj = Bj0 + j,
        Ac = A.jc[Aj];
      sd.jc[Bj] = qv.has(Ac) ? qv.to_block(Ac) : qv.size() + pv.to_block(Ac);
      sd.d[Bj] = A.d[Aj];
      p.A_idxs[Bj] = Aj;
    }
  }

  try { p.cm = new CrsMatrix(pv.size(), A.n, sd.ir, sd.jc, sd.d); }
  catch (...)
  { throw hts::Exception("get_matrix_p_qp_with_covers_all failed to alloc."); }
  sd.release();
}

template<typename Int, typename Size, typename T> struct SortEntry {
  Size A_idx;
  Int i;
  T d;
  SortEntry () {}
  bool operator< (const SortEntry& se) const { return i < se.i; }
};

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::sort (Partition& p) {
  const int nthreads = omp_get_max_threads();
  CrsMatrix& A = *p.cm;

  //todo ||ize
  Size max_nc = 0;
  for (Int r = 0; r < A.m; ++r)
    max_nc = std::max(max_nc, A.ir[r+1] - A.ir[r]);
  std::vector<SortEntry<Int, Size, Real> > sess(nthreads * max_nc);

# pragma omp parallel for schedule(static, 1)
  for (Int r = 0; r < A.m; ++r) {
    const int tid = omp_get_thread_num();
    const Size irr = A.ir[r], irrp1 = A.ir[r+1], nc = irrp1 - irr;
    SortEntry<Int, Size, Real>* ses = &sess[tid*max_nc];
    for (Size j = 0; j < nc; ++j) {
      const Size Aj = irr + j;
      ses[j].i = A.jc[Aj];
      ses[j].d = A.d[Aj];
      ses[j].A_idx = p.A_idxs[Aj];
    }
    std::sort(ses, ses + nc);
    for (Size j = 0; j < nc; ++j) {
      const Size Aj = irr + j;
      A.jc[Aj] = ses[j].i;
      A.d[Aj] = ses[j].d;
      p.A_idxs[Aj] = ses[j].A_idx;
    }
  }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::partition_into_2_blocks (
  const ConstCrsMatrix& A, const bool is_lo, const std::vector<Int>& lsis,
  const std::vector<Int>& dpis, Partition* p)
{
  PermVec lsis_pv(A.m, lsis), dpis_pv(A.m, dpis);
  if (is_lo) {
    get_matrix_pp_with_covers_all(A, lsis_pv, p[0]);
    get_matrix_p_qp_with_covers_all(A, dpis_pv, lsis_pv, p[1]);
  } else {
    get_matrix_pp_with_covers_all(A, dpis_pv, p[1]);
    get_matrix_p_qp_with_covers_all(A, lsis_pv, dpis_pv, p[0]);
  }
  sort(p[0]);
  sort(p[1]);
}

// inside || {}
template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
reverse_A_idxs (const Size nnz, Partition& p) {
# pragma omp for schedule(static)
  for (size_t i = 0; i < p.A_idxs.size(); ++i)
    p.A_idxs[i] = nnz - p.A_idxs[i] - 1;
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
copy_partition (const ConstCrsMatrix& A, Partition& p) {
  const Size ilim = static_cast<Size>(p.A_idxs.size());
# pragma omp parallel for schedule(static)
  for (Size i = 0; i < ilim; ++i)
    p.cm->d[i] = A.d[p.A_idxs[i]];
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
repartition_into_2_blocks (Partition* const p, const ConstCrsMatrix& A) {
  for (int i = 0; i < 2; ++i)
    if (p[i].cm)
      copy_partition(A, p[i]);
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
CrsSegmenter::count_nnz_by_row_loop (const Int i, std::vector<Int>& rcnt) {
  Int i_first, i_last;
  rcnt[i] = find_first_and_last(A_.ir, r0_ + i, A_.jc, c0_, c0_ + nc_,
                                i_first, i_last);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
CrsSegmenter::count_nnz_by_row (std::vector<Int>& rcnt) {
  rcnt.resize(nr_);
  // Don't allow nested ||ism.
  if (omp_get_num_threads() == 1) {
#   pragma omp parallel for schedule(guided)
    for (Int i = 0; i < nr_; ++i)
      count_nnz_by_row_loop(i, rcnt);
  } else
    for (Int i = 0; i < nr_; ++i)
      count_nnz_by_row_loop(i, rcnt);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
CrsSegmenter::init_nnz (const std::vector<Int>& rcnt) {
  const std::vector<Int>& p = this->p_;
  const Int nseg = p.size() - 1;
  this->nnz_.resize(nthreads_, 0);
  for (Int i = 0; i < nseg; ++i)
    for (Int j = p[i] - p[0]; j < p[i+1] - p[0]; ++j)
      this->nnz_[i] += rcnt[j];
}

// p partitions rows. It may not be a perfect partitioning in the sense of
// optimal given that each thread must handle a mutually exclusive set of
// rows. Attempt to remove spikes in #rows.
template<typename Int, typename Size>
inline void
smooth_spikes (const std::vector<Int>& rcnt, std::vector<Int>& p,
               std::vector<Size>& nnz, const bool ignore_0) {
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::CrsSegmenter::segment () {
  assert(nr_ % ls_blk_sz_ == 0);
  std::vector<Int>& p = this->p_;

  const Int nseg = std::min<Int>(nthreads_, nr_ / ls_blk_sz_);
  if (nseg == 0) {
    assert(nr_ == 0);
    p.resize(1, r0_);
    return;
  }

  std::vector<Int> rcnt;
  count_nnz_by_row(rcnt);

  // cumsum(rcnt).
  std::vector<Int> cs_rcnt(rcnt.size() / ls_blk_sz_);
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

  p.resize(nseg + 1, r0_);
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TMatrix::clear () {
  bs_.clear();
  ros_.clear();
  is_empty_ = true;
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TMatrix::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
               const InitInfo& in, const Int block_0_nnz_os,
               const int tid_offset) {
  init_metadata(A, r0, c0, nr, nc, in, block_0_nnz_os, tid_offset);
  init_memory(in);
# pragma omp parallel
  init_numeric(A, omp_get_thread_num());
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TMatrix::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
               const InitInfo& in, const CrsSegmenter& seg) {
  init_metadata(A, r0, c0, nr, nc, in, seg);
  init_memory(in);
# pragma omp parallel
  init_numeric(A, omp_get_thread_num());
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TMatrix::init_metadata (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                        const InitInfo& in, const CrsSegmenter& seg) {
  clear();
  nr_ = nr;
  init_metadata_with_seg(A, r0, c0, nr, nc, 0, in, seg);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::TMatrix::
init_metadata_with_seg (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                        const Int roff, const InitInfo& in,
                        const CrsSegmenter& seg) {
  // Serial if block is too small.
  const Int nseg = seg.p().size() - 1;
  is_parallel_ = nseg > 1 || tid_os_ > 0;
  if (nseg == 1) {
    bs_.push_back(SerialBlock());
    bs_.back().init_metadata(A, r0, c0, nr, nc, in);
    ros_.push_back(roff);
  } else {
    const std::vector<Int>& p = seg.p();
    const Int n = static_cast<Int>(p.size()) - 1;
    ros_.resize(n, 0);
    bs_.resize(n, SerialBlock());
    for (Int id = 0; id < n; ++id) {
      ros_[id] = p[id] - r0 + roff;
      const Int nri = p[id+1] - p[id];
      bs_[id].init_metadata(A, p[id], c0, nri, nc, in);
    }
  }

  is_empty_ = true;
  for (size_t i = 0; i < bs_.size(); ++i)
    if (bs_[i].nnz() != 0) {
      is_empty_ = false;
      break;
    }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TMatrix::init_memory (const InitInfo& in) {
# pragma omp parallel
  { const int tid = omp_get_thread_num(), id = tid - tid_os_;
    if (id >= 0 && id < (int) bs_.size())
      bs_[id].init_memory(in);
  }
}

// possibly inside || {}
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TMatrix::init_numeric (const CrsMatrix& A, const int tid) {
  const int id = tid - tid_os_;
  if (id >= 0 && id < (int) bs_.size())
    bs_[id].init_numeric(A);
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
TMatrix::reinit_numeric (const CrsMatrix& A) {
  if ( ! is_parallel_) {
    if (bs_.empty()) return;
    bs_[0].reinit_numeric(A);
  } else {
#   pragma omp parallel
    { const int tid = omp_get_thread_num(), id = tid - tid_os_;
      if (id >= 0 && id < (int) bs_.size())
        bs_[id].reinit_numeric(A);
    }
  }
}

template<typename Int, typename Size, typename Real>
inline Int Impl<Int, Size, Real>::
TMatrix::block_r0 (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return bs_[id].r0();
}

template<typename Int, typename Size, typename Real>
inline Int Impl<Int, Size, Real>::
TMatrix::block_nr (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return bs_[id].nr();
}

template<typename Int, typename Size, typename Real>
inline const typename Impl<Int, Size, Real>::SerialBlock*
Impl<Int, Size, Real>::TMatrix::block (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return &bs_[id];
}

template<typename Int, typename Size, typename Real>
inline typename Impl<Int, Size, Real>::SerialBlock*
Impl<Int, Size, Real>::TMatrix::block (const int tid) {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return &bs_[id];
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
SerialBlock::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                   const InitInfo& in) {
  init_metadata(A, r0, c0, nr, nc, in);
  init_memory(in);
  init_numeric(A);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
SerialBlock::init_memory (const InitInfo& in) {
  if (nr_ == 0 || nc_ == 0) return;
  if (is_dense_) {
    // First touch is measurably required.
    d_ = Allocnator<Real>(nr_*nc_, "SerialBlock::init dense d", true).release();
  } else {
    Allocnator<Size> air(nr_+1, "SerialBlock::init sparse ir", true);
    if (nnz_ == 0) {
      ir_ = air.release();
      ir_[nr_] = 0;
      return;
    }
    Allocnator<Int> ajc(nnz_, "SerialBlock::init sparse jc", true);
    Allocnator<Real> ad(nnz_, "SerialBlock::init sparse d", true);
    ir_ = air.release();
    jc_ = ajc.release();
    d_ = ad.release();
  }
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
SerialBlock::init_numeric (const CrsMatrix& A) {
  if (is_dense_) memset(d_, 0, nr_*nc_*sizeof(*d_));
  reinit_numeric(A);
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
SerialBlock::reinit_numeric (const CrsMatrix& A) {
  if (ir_) reinit_numeric_spars(A);
  else reinit_numeric_dense(A);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
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
      const Int k = nc_*lrow + lcol;
      d_[k] = A.d[j];
    }
  }
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real>
inline Int Impl<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real>
Impl<Int, Size, Real>::
OnDiagTri::OnDiagTri ()
  : c0_(0), nnz_(0), d_(0), m_(0), dense_(true)
{}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::clear () {
  deln(d_);
  del(m_);
  t_.clear();
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::init (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                 const InitInfo& in) {
  init_metadata(T, r0, c0, n, in);
  init_memory(in);
  init_numeric(T);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::init (const Int r0, const Int c0, const Int n) {
  clear();
  this->n_ = n; this->r0_ = r0; c0_ = c0;
}

template<typename Int, typename Size, typename Real>
inline bool Impl<Int, Size, Real>::
is_dense_tri (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
              const InitInfo& in, Size* innz) {
  const Size nnz = count_nnz_lotri(T, r0, c0, n);
  if (innz) *innz = nnz;
  return nnz >= 0.5*in.min_dense_density*ntri<Int>(n);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::init_memory (const InitInfo& in) {
  if ( ! nnz_) return;
  if (dense_) {
    d_ = Allocnator<Real>(
      ntri<Int>(this->n_), "OnDiagTri::init dense", true).release();
    if ( ! t_.empty()) inv_init_memory();
  } else {
    SparseData sd(this->n_, nnz_, "OnDiagTri::init_memory");
    try { m_ = new CrsMatrix(this->n_, this->n_, sd.ir, sd.jc, sd.d); }
    catch (...)
    { throw hts::Exception("OnDiagTri::init_memory failed to alloc."); }
    sd.release();
  }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::init_numeric (const CrsMatrix& T) {
  if (dense_) {
    reinit_numeric(T);
  } else {
    Size* const ir = m_->ir;
    Int* const jc = m_->jc;
    Real* const d = m_->d;
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
        d[i] = lrow == lcol ? 1/T.d[k] : T.d[k];
        ++i;
      }
    }
    assert(ir[this->n_] == nnz_);
  }
}

//todo At the cost of doubling the jc-related memory, I could store indices into
// T in init_numeric and then use them here in reinit_numeric so I wouldn't have
// to run find_first. At least record the find_first results.
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::reinit_numeric (const CrsMatrix& T) {
  if (d_) {
    memset(d_, 0, ntri<Int>(this->n_)*sizeof(*d_));
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
        const Real dv = lrow == lcol ? 1/T.d[k] : T.d[k];
        // Compressed dense triangle.
        d_[di] = dv;
        ++nnz;
      }
    }
    if ( ! t_.empty())
      inv_reinit_numeric(T);
  } else if (m_) {
    Real* const d = m_->d;
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
        d[i] = lrow == lcol ? 1/T.d[k] : T.d[k];
        ++i;
      }
    }
  }
}

template<typename Int, typename Size, typename Real>
inline Int Impl<Int, Size, Real>::
OnDiagTri::nthreads () const {
  return std::max<Int>(1, static_cast<Int>(t_.size()));
}

template<typename Int, typename Size, typename Real>
inline Int Impl<Int, Size, Real>::
OnDiagTri::block_row_start (const int tid) const {
  return t_.empty() ? 0 : t_[tid].r0;
}

template<typename Int, typename Size, typename Real>
inline Int Impl<Int, Size, Real>::
OnDiagTri::block_nr (const int tid) const {
  return t_.empty() ? this->n_ : t_[tid].nr;
}

template<typename Int, typename Size, typename Real>
Impl<Int, Size, Real>::
OnDiagTri::Thread::~Thread () { deln(d); }

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::inv_init_metadata (const InitInfo& in) {
  const Int
    nnz = ntri<Int>(this->n_),
    nt = std::min((int) in.nthreads, 1 + nnz / square(in.min_parallel_rows));
  if (nt <= 1) return;

  t_.resize(nt);
  const Real nnz_per_thread = (Real) nnz / nt;
  Int r0 = 0;
  for (Int tid = 0; tid < nt; ++tid) {
    t_[tid].r0 = r0;
    const Int n_max = this->n_ - r0;
    if (tid+1 == nt)
      t_[tid].nr = n_max;
    else {
      // Solve for n in
      //   ((r0 + n) (r0 + n + 1))/2 - (r0 (r0 + 1))/2 = nnz_per_thread.
      const Int b = 1 + 2*r0;
      t_[tid].nr = std::min<Real>(
        (Real) n_max, round(0.5*(std::sqrt(b*b + 8*nnz_per_thread) - b)));
    }
    r0 += t_[tid].nr;
  }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::inv_init_memory () {
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    if (tid < nthreads()) {
      Thread& t = t_[tid];
      const Int nnz = ntri<Int>(t.r0 + t.nr) - ntri<Int>(t.r0);
      t.d = Allocnator<Real>(nnz, "OnDiagTri::inv_init_memory", true).release();
    } }
}

// T is in row-major compressed dense tri format. The diag of T is the
// reciprocal.
template<typename Int, typename Real>
inline void invert (Real* T, const Int n, Real* w) {
  for (Int c = 0; c < n; ++c) {
    // Solve for column c. That involves only the (n-c)x(n-c) lower-right
    // subblock of T. Store column in w.
    w[0] = T[0];
    Real* Tp = T + c + 1;
    for (int r = 1; r < n - c; ++r) {
      w[r] = 0;
      for (int k = 0; k < r; ++k)
        w[r] -= w[k]*Tp[k];
      w[r] *= Tp[r];
      Tp += r + c + 1;
    }
    // Copy w to column c.
    Tp = T;
    for (int r = 0; r < n - c; ++r) {
      *Tp = w[r];
      Tp += r + c + 1;
    }
    T += c + 2;
  }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
OnDiagTri::inv_reinit_numeric (const CrsMatrix& T) {
  { std::vector<Real> w(this->n_);
    invert(d_, this->n_, w.data()); }
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
find_split_rows (const CrsMatrix& T, const Int r0, const Int c0,
                 const Int n, const InitInfo& in,
                 std::vector<Int>& split_rows) {
  const Int sz = in.min_blksz / 2;
# pragma omp parallel for schedule(static, parfor_static_size)
  for (Int r = r0 + 1; r < r0 + n; ++r) {
    bool found = true;
    for (Int r1 = r, r1_lim = std::min<Int>(T.m, r + sz); r1 < r1_lim; ++r1)
      if ( ! empty<Int>(T.jc + T.ir[r1], T.ir[r1+1] - T.ir[r1],
                        std::max<Int>(c0, c0 + r - sz), c0 + r - 1)) {
        found = false;
        break;
      }
    if (found) {
#     pragma omp critical (find_split_rows)
      split_rows.push_back(r);
    }
  }
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::build_recursive_tri_r (
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
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
build_recursive_tri (const CrsMatrix& T, const Int r, const Int c, const Int n,
                     const Int mvp_block_nc, const InitInfo& in,
                     std::vector<Box>& bv) {
  std::list<Box> bl;

  // This is the large MVP block that scatters the LS part to the RB part. It is
  // present only in the T = L case; in the T = U case, it is in the LS data
  // structure.
  if (mvp_block_nc)
    bl.push_back(Box(r, c, n, mvp_block_nc));

  // Find clear divisions in this triangle that we can exploit to pack the data
  // efficiently.
  std::vector<Int> split_rows;
  find_split_rows(T, r, c + mvp_block_nc, n, in, split_rows);

  build_recursive_tri_r(T, r, c + mvp_block_nc, n, in, split_rows, bl);

  // list -> vector.
  bv.resize(bl.size());
  Int i = 0;
  for (typename std::list<Box>::const_iterator it = bl.begin();
       it != bl.end(); ++it)
    bv[i++] = *it;
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
RecursiveTri::init (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                    const InitInfo& in, const Int mvp_block_nc) {
  clear();
  this->n_ = n; this->r0_ = r0; nthreads_ = in.nthreads;

  std::vector<Box> bs;
  build_recursive_tri(T, r0, c0, n, mvp_block_nc, in, bs);

  const Int
    nblock = bs.size(),
    ntri = ((mvp_block_nc ? 2 : 1) + nblock) / 2,
    nmvp = ntri - 1;
  nd_.t.resize(ntri);
  nd_.s.resize(nmvp);
  nd_.os.resize(ntri);

#ifdef __MIC__
  // || is working pretty well on MIC, but not on CPU.
  // Initialize metadata, such as nnz and partitioning information.
  if (mvp_block_nc)
    nd_.t[0].init(0, 0, mvp_block_nc);
  const Int i0 = mvp_block_nc ? 1 : 0;
# pragma omp parallel
  {
#   pragma omp for schedule(static)
    for (Int i = i0; i < ntri; ++i) {
      const Box& b = bs[i0 + 2*(i - i0)];
      nd_.t[i].init_metadata(T, b.r0, b.c0, b.nr, in);
    }
#   pragma omp for schedule(static)
    for (Int i = 0; i < nmvp; ++i) {
      const Box& b = bs[i0 + 2*(i - i0) + 1];
      nd_.s[i].init_metadata(T, b.r0, b.c0, b.nr, b.nc, in);
      nd_.os[i] = b.c0;
    }
  }
  Timer::start(Timer::dpb_tinit_2);
  // Initialize memory in order.
  nd_.os.back() = 0;
  for (Int i = 0; i < nmvp; ++i) {
    nd_.t[i].init_memory(in);
    nd_.s[i].init_memory(in);
  }
  nd_.t.back().init_memory(in);
  Timer::stop(Timer::dpb_tinit_2); Timer::start(Timer::dpb_tinit_3);
  // Initialize numerical data.
# pragma omp parallel
  {
#   pragma omp for schedule(static)
    for (Int i = i0; i < ntri; ++i)
      nd_.t[i].init_numeric(T);
    const Size ilim = nmvp*nthreads_;
#   pragma omp for schedule(static, parfor_static_size)
    for (Size i = 0; i < ilim; ++i) {
      const Int si = i / nthreads_, tid = i % nthreads_;
      nd_.s[si].init_numeric(T, tid);
    }
  }
  Timer::stop(Timer::dpb_tinit_3);
#else
  Int ti = 0, si = 0;
  typename std::vector<Box>::iterator bit = bs.begin();
  if (mvp_block_nc) {
    nd_.t[ti++].init(0, 0, mvp_block_nc);
    const Box& b = *bit; ++bit;
    nd_.s[si].init(T, b.r0, b.c0, b.nr, b.nc, in);
    nd_.os[si++] = b.c0;
  }
  while (si < nmvp) {
    { const Box& b = *bit; ++bit;
      nd_.t[ti++].init(T, b.r0, b.c0, b.nr, in); }
    { const Box& b = *bit; ++bit;
      nd_.s[si].init(T, b.r0, b.c0, b.nr, b.nc, in);
      nd_.os[si++] = b.c0; }
  }
  { const Box& b = *bit; ++bit;
    nd_.t[ti++].init(T, b.r0, b.c0, b.nr, in); }
  nd_.os.back() = 0;

  assert(bit == bs.end());
  assert(ti == static_cast<Int>(nd_.t.size()));
  assert(si == static_cast<Int>(nd_.s.size()));
#endif

  { // Initialize threading data for the inverse of the on-diag tri.
    max_diag_tri_ = 0;
    Int max_nthreads = 0;
    for (size_t i = 0; i < nd_.t.size(); ++i) {
      max_diag_tri_ = std::max(max_diag_tri_, nd_.t[i].n());
      max_nthreads = std::max(max_nthreads, nd_.t[i].nthreads());
    }
    wrk_.resize(max_diag_tri_ * in.max_nrhs);
    nd_.inv_tri_done.resize(max_nthreads);
  }

  p2p_init();

#if defined OUTPUT_RECURSIVEBLOCK_PARTITION
  do {
    if (omp_get_max_threads() != 32) break;
    const Int lc0 = mvp_block_nc + c0;
    static int ctr = 0;
    Size nnz = 0;
    for (Int r = r0; r < r0 + n; ++r) {
      const Size rn = T.ir[r+1] - T.ir[r];
      nnz += rn - find_first<Int>(T.jc + T.ir[r], rn, lc0);
    }
    if (nnz > 17500000) break;
    { // Write recursive block as a matrix-market file.
      std::stringstream ss;
      ss << "rtp" << ctr << ".mm";
      FILE* fid = fopen(ss.str().c_str(), "w");
      fprintf(fid, "%%%%MatrixMarket matrix coordinate real general\n"
              "%11d %11d %11d\n", n, n, static_cast<Int>(nnz));
      Size cnt = 0;
      for (Int r = r0; r < r0 + n; ++r) {
        const Int rn = T.ir[r+1] - T.ir[r];
        for (Size j = T.ir[r] + find_first<Int>(T.jc + T.ir[r], rn, lc0);
             j < T.ir[r+1]; ++j, ++cnt) {
          fprintf(fid, "%d %d %1.15e\n", r - r0 + 1, T.jc[j] - lc0 + 1, T.d[j]);
          assert(cnt <= nnz);
        }
      }
      fclose(fid);
    }
    { // Write the recursive block decomposition.
      std::stringstream ss;
      ss << "rtp" << ctr << ".dat";
      FILE* fid = fopen(ss.str().c_str(), "w");
      for (size_t i = 1 + (mvp_block_nc ? 1 : 0); i < bs.size(); i += 2) {
        const Box& b = bs[i];
        fprintf(fid, "%d %d %d %d\n",
                b.r0 - r0 + 1, b.c0 - lc0 + 1, b.nr, b.nc);
      }
      fclose(fid);
    }
    { // Write the cropped SerialBlocks.
      std::stringstream ss;
      ss << "rtp" << ctr << "-sb.dat";
      FILE* fid = fopen(ss.str().c_str(), "w");
      for (size_t i = mvp_block_nc ? 1 : 0; i < nd_.s.size(); ++i) {
        const TMatrix& tm = nd_.s[i];
        for (Int j = 0; j < tm.nblocks(); ++j) {
          const SerialBlock& s = *tm.block(j);
          fprintf(fid, "%d %d %d %d\n",
                  s.r0() - r0 + 1, s.c0() - lc0 + 1, s.nr(), s.nc());
        }
      }
      fclose(fid);
    }
    ++ctr;
  } while (0);
#endif
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
RecursiveTri::init_numeric (const CrsMatrix& T) {
  Int n = static_cast<Int>(nd_.t.size());
  Timer::start(Timer::numrbt);
# pragma omp parallel for schedule(dynamic)
  for (Int i = 0; i < n; ++i)
    nd_.t[i].reinit_numeric(T);
  Timer::stop(Timer::numrbt); Timer::start(Timer::numrbm);
  n = static_cast<Int>(nd_.s.size());
  for (Int i = 0; i < n; ++i)
    if ( ! nd_.s[i].empty()) nd_.s[i].reinit_numeric(T);
  Timer::stop(Timer::numrbm);
}

// Form the lists of p2p dependencies. A dependency in this algorithm is of two
// types. One is the usual dependency: a variable has to be solved for before
// the next. The second imposes an ordering to assure there is no write race
// condition when computing a row's dot product in pieces.
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
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
  std::vector<Size> w(this->n_, marker);
  // Fast way to get a set of unique elements.
  std::vector<Int> w_cnted(sn * nthreads_, 0);
  Int w_cnted_symbol = 1;

  nd_.s_done.resize(nthreads_*sn, 0);
  nd_.t_idx.resize(tn);
  std::vector<std::vector<Size> > s_ids(sn * nthreads_);
  for (Int ti = 0; ti < tn; ++ti) {
    // Consider each tri or MVP block in solution order.
    if (ti > 0) { // Tri block. First tri has no dependencies.
      const Tri& t = nd_.t[ti];
      const Int r0 = t.r0(), nr = t.n();
      Int k = ti == 1 ? 0 : nd_.t_idx[ti-1];
      for (Int r = r0, rlim = r0 + nr; r < rlim; ++r) {
        assert(r < this->n_);
        const Size wr = w[r];
        if (wr != marker && w_cnted[wr] != w_cnted_symbol) {
          const Int sid = rb_p2p_ind2sid(wr), tid = rb_p2p_ind2tid(wr);
          nd_.t_ids.push_back(rb_p2p_sub2ind(sid, tid));
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
            if (tid != bi)
              s_ids[ind].push_back(wr);
            w_cnted[wr] = w_cnted_symbol;
          }
          // I now have precedence.
          w[r] = ind;
        }
        ++w_cnted_symbol;
      }
    }
  }

  nd_.s_ids.resize(nthreads_);
  nd_.s_idx.resize(nthreads_);
  for (Int i = 0; i < nthreads_; ++i) nd_.s_idx[i].push_back(0);
  for (Size i = 0, ilim = static_cast<Size>(s_ids.size()); i < ilim; ++i) {
    const Int tid = rb_p2p_ind2tid(i);
    nd_.s_idx[tid].push_back(nd_.s_idx[tid].back() +
                             static_cast<Int>(s_ids[i].size()));
    for (Size j = 0, jlim = static_cast<Size>(s_ids[i].size()); j < jlim; ++j)
      nd_.s_ids[tid].push_back(s_ids[i][j]);
  }

  compress(nd_.t_ids);
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    compress(nd_.s_ids[tid]);
    compress(nd_.s_idx[tid]); }
  compress(nd_.s_ids);
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::RecursiveTri::p2p_reset () const {
  nd_.t_barrier = -1;
  for (size_t i = 0; i < nd_.inv_tri_done.size(); ++i)
    nd_.inv_tri_done[i] = -1;
  ++nd_.done_symbol;
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
RecursiveTri::reset_max_nrhs (const Int max_nrhs) {
  wrk_.resize(max_diag_tri_ * max_nrhs);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
LevelSetTri::init_lsets (const LevelSetter& lstr,
                         const bool save_for_reprocess) {
  ls_blk_sz_ = lstr.ls_blk_sz();
  save_for_reprocess_ = save_for_reprocess;
  lsp_.resize(lstr.size() + 1);
  lsp_[0] = 0;
  for (Int i = 0; i < lstr.size(); ++i)
    lsp_[i+1] = lsp_[i] + lstr.lset(i).size();
  nlvls_ = static_cast<Int>(lsp_.size()) - 1;
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
LevelSetTri::init (const CrsMatrix& T, const Int r0, const Int c0,
                   const Int n, const InitInfo& in) {
  n_ = n;
  t_.resize(in.nthreads);
  
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    t_[tid].lsp.push_back(0); }
  
  ps_.resize(in.nthreads);
  for (Int ils = 0; ils < static_cast<Int>(lsp_.size()) - 1; ++ils) {
    const Int
      r0 = lsp_[ils],
      c0 = 0,
      nr = lsp_[ils+1] - r0,
      nc = lsp_[ils+1] + mvp_block_nc_;

    const Int nseg = in.nthreads;
    CrsSegmenter seg(T, r0, c0, nr, nc, nseg, ls_blk_sz_);
    const std::vector<Int>& p = seg.p();
    assert(p[1] >= p[0]);

#   pragma omp parallel
    { int tid = omp_get_thread_num();
      if (tid + 1 < (int) p.size()) {
        t_[tid].p.push_back(p[tid]);
        t_[tid].lsp.push_back(t_[tid].lsp.back() + p[tid+1] - p[tid]);
        for (Int j = p[tid]; j < p[tid+1]; ++j)
          ps_[tid].push_back(j);
      } else {
        t_[tid].p.push_back(-1);
        t_[tid].lsp.push_back(t_[tid].lsp.back());
      } }
  }

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
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
LevelSetTri::update_permutation (std::vector<Int>& lsis, const Partition& p) {
  const std::vector<Int> old_lsis = lsis;
  std::vector<Int> q(p.cm->n);
# pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    Thread& t = t_[tid];
    int start = 0;
    for (int i = 0; i < tid; ++i)
      start += t_[i].m->m;
    const Int nlev = static_cast<Int>(t.lsp.size()) - 1;
    for (Int ils = 0, k = 0; ils < nlev; ++ils) {
      const Int p_ils = t.p[ils];
      t.p[ils] = start + k;
      const Int n = t.lsp[ils+1] - t.lsp[ils];
      for (Int i = 0; i < n; ++i, ++k) {
        lsis[start + k] = old_lsis[p_ils + i];
        q[p_ils + i] = start + k;
      }
    }
#   pragma omp barrier
    for (Size i = 0, nnz = t.m->ir[t.m->m]; i < nnz; ++i) {
      Int& jci = t.m->jc[i];
      if (jci < mvp_block_nc_) continue;
      jci = mvp_block_nc_ + q[jci - mvp_block_nc_];
    }
  }
}

// inside || {}
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::LevelSetTri::
find_task_responsible_for_variable (std::vector<p2p_Pair>& pairs) {
  const int tid = omp_get_thread_num();
  // Find the level and tid responsible for variable i.
  Thread& t = t_[tid];
  const std::vector<Int>& lsp = t.lsp;
  const std::vector<Int>& p = t.p;
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
template<typename Int, typename Size, typename Real>
Int Impl<Int, Size, Real>::LevelSetTri::
fill_graph (const std::vector<p2p_Pair>& pairs, std::vector<Size>& g,
            std::vector<Size>& gp, std::vector<Size>& wrk) {
  // O(nnz/nthreads) time and O(nnz) space. Reduce the entry-wise graph to the
  // (lvl, tid)-wise dependency graph.
  const int tid = omp_get_thread_num();
  const Size n = static_cast<Size>(gp.size()) - 1;
  Thread& t = t_[tid];
  const std::vector<Int>& lsp = t.lsp;
  // Max value of e + 1, used as a marker.
  const Size me = nlvls_*t_.size();
  // Build the graph g. g[e] is the list of dependencies for (level, tid) e.
  const CrsMatrix* const cm = t.m;
  const Size* const ir = cm->ir;
  const Int* const jc = cm->jc;

  // Get an unused piece of the workspace.
  Size* Te = wrk.data();
  for (int i = 0; i < tid; ++i)
    if (t_[i].m) {
      const CrsMatrix* const m = t_[i].m;
      Te += m->ir[m->m];
    }

  // Count entries in g.
  for (Int ils = 1, ils_lim = static_cast<Int>(lsp.size()) - 1;
       ils < ils_lim; ++ils) {
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
# pragma omp barrier
# pragma omp master
  for (Size i = 1; i <= n; ++i) {
    if (static_cast<Int>(gp[i]) > max_gelen) max_gelen = gp[i];
    gp[i] += gp[i-1];
  }
# pragma omp barrier

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
// O(max(an, bn)).
template<typename Int, typename Size> inline void
mark_intersection (const Size* a, const Int an, const Size* b, const Int bn,
                   Size* mark /* len(mark) == an */, const Size marker) {
  Int ai = 0, bi = 0;
  while (ai < an && bi < bn) {
    if (a[ai] < b[bi])
      ++ai;
    else if (a[ai] > b[bi])
      ++bi;
    else {
      mark[ai] = marker;
      ++ai;
      ++bi;
    }
  }
}

// inside || {}
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::LevelSetTri::
prune_graph (const std::vector<Size>& gc, const std::vector<Size>& gp,
             std::vector<Size>& g, std::vector<Int>& gsz,
             std::vector<Size>& wrk, const Int max_gelen) {
  // Time complexity is a little complicated. See Park et al 2014. Space is
  // O(#threads max |g(e)|).
  const int tid = omp_get_thread_num();
  const Size n = static_cast<Size>(gsz.size());
  const Size me = nlvls_*t_.size();
  Size* mark = &wrk[tid*max_gelen];
  // I find that it's more efficient to use a parfor than to work on a thread's
  // own e's.
# pragma omp for schedule(static,1)
  for (Size e = 0; e < n; ++e) {
    const Int gcelen = static_cast<Int>(gp[e+1] - gp[e]);
    if (gcelen == 0) continue;
    // e's dependencies.
    const Size* const gce = &gc[gp[e]];
    for (Int i = 0; i < gcelen; ++i)
      mark[i] = gce[i];
    for (Int ied = 0; ied < gcelen; ++ied) { // For each of e's deps:
      const Size ed = gce[ied];
      assert(ed >= 0 && ed < n);
      if (ls_p2p_ind2lvl(ed) == 0) continue; // No parent deps to check.
      // ed's dependencies.
      const Size* const gced = &gc[gp[ed]];
      const Int gcedlen = static_cast<Int>(gp[ed+1] - gp[ed]);
      mark_intersection(gce, gcelen, gced, gcedlen, mark, me);
    }
    // Insert the pruned set of dependencies.
    Int k = 0;
    Size* const ge = &g[gp[e]];
    const Int etid = ls_p2p_ind2tid(e);
    for (Int i = 0; i < gcelen; ++i) {
      const Size ed = mark[i];
      if (ed == me) continue;
      const Int edtid = ls_p2p_ind2tid(ed);
      if (edtid != etid) ge[k++] = ed;
    }
    gsz[e] = k;
  }
}

// inside || {}
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::LevelSetTri::
fill_dependencies (const std::vector<Size>& g, const std::vector<Size>& gp,
                   const std::vector<Int>& gsz) {
  const int tid = omp_get_thread_num();
  Thread& t = t_[tid];
  const std::vector<Int>& lsp = t.lsp;
  // Allocate inside this thread up front. I could probably do this even more
  // efficiently, but fill_dependencies is negligible compared with fill_graph
  // and prune_graph.
  t.p2p_depends_p.reserve(nlvls_ + 1);
  Int sz = 0;
  for (Int ils = 1, ils_lim = static_cast<Int>(lsp.size()) - 1;
       ils < ils_lim; ++ils)
    sz += gsz[ls_p2p_sub2ind(ils, tid)];
  t.p2p_depends.reserve(sz);
  // Make the final lists of dependencies.
  t.p2p_depends_p.push_back(0);
  for (Int ils = 1, ils_lim = static_cast<Int>(lsp.size()) - 1;
       ils < ils_lim; ++ils) {
    const Size e = ls_p2p_sub2ind(ils, tid);
    const Size* const ed = &g[gp[e]];
    const Int edsz = gsz[e];
    // Sort by increasing level number. Heuristic to speed up the
    // synchronization step in p2p_solve. Idea is to do p2p_done_ checking on
    // dependencies higher up in the tree first, as those are likely done
    // sooner.
    std::vector<p2p_SortEntry> es(edsz); //todo Remove memory alloc.
    for (Int ied = 0; ied < edsz; ++ied) {
      es[ied].lvl = ls_p2p_ind2lvl(ed[ied]);
      es[ied].tid = ls_p2p_ind2tid(ed[ied]);
    }
    std::sort(es.begin(), es.end());
    // Insert.
    t.p2p_depends_p.push_back(t.p2p_depends_p[ils-1]);
    for (Int ied = 0; ied < edsz; ++ied) {
      t.p2p_depends.push_back(ls_p2p_sub2ind(es[ied].lvl, es[ied].tid));
      ++t.p2p_depends_p.back();
    }
  }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::LevelSetTri::p2p_init () {
  if (t_[0].lsp.empty()) return;

  Timer::start(Timer::lsp2p_1);
  if (p2p_done_.empty()) {
    p2p_done_value_ = 0;
    p2p_done_.resize(t_[0].lsp.size()*t_.size(), 1);
  }
  if (t_[0].lsp.size() <= 1) return;

  Size nnz = 0;
  for (size_t i = 0, n = t_.size(); i < n; ++i)
    if (t_[i].m) {
      const CrsMatrix* const cm = t_[i].m;
      nnz += cm->ir[cm->m];
    }
  const Size n = t_.size() * t_[0].lsp.size();

  // g is a graph. nnz is an upper bound on the memory needed.
  //   g(e) is the set of dependencies (incoming edges) for node e. gp is the
  // pointer to g(e). So g(e) is written g[gp[e]]. gsz is the size of g(e).
  //   A node e represents a pair (level, thread id) that is a task.
  std::vector<Int> gsz;
  std::vector<Size> g, gc, gp, wrk;
  // Thread and level responsible for a variable.
  std::vector<p2p_Pair> pairs;
  try {
    g.resize(nnz);
    gc.resize(nnz);
    gp.resize(n+1, 0);
    gsz.resize(n);
    wrk.resize(nnz);
    pairs.resize(n_);
  } catch (...) {
    std::stringstream ss;
    ss << "p2p_init failed to resize: n = " << n << " nnz = " << nnz;
    throw hts::Exception(ss.str());
  }
  Timer::stop(Timer::lsp2p_1);

  Int max_gelen;
# pragma omp parallel
  {
    find_task_responsible_for_variable(pairs);
#   pragma omp barrier
#ifdef TIME
#   pragma omp master
    Timer::start(Timer::lsp2p_3);
#endif
    const Int max_gelen_t = fill_graph(pairs, g, gp, wrk);
#   pragma omp barrier
#   pragma omp master
    {
      Timer::stop(Timer::lsp2p_3);
      pairs.clear();
      // Tell all threads; only master's max_gelen_t is valid.
      max_gelen = max_gelen_t;
      // In the unlikely case that max_gelen * #threads > nnz, allocate more
      // workspace.
      const Size space = max_gelen * t_.size();
      if (space > nnz)
        try { wrk.resize(space); }
        catch (...) {
          std::stringstream ss;
          ss << "p2p_init failed to resize wrk: space = " << space;
          throw hts::Exception(ss.str());
        }
    }
    // Keep the original graph.
#   pragma omp for
    for (Size i = 0; i < nnz; ++i)
      gc[i] = g[i];
#ifdef TIME
#   pragma omp master
    Timer::start(Timer::lsp2p_6);
#endif
    prune_graph(gc, gp, g, gsz, wrk, max_gelen);
#   pragma omp barrier
#ifdef TIME
#   pragma omp master
    Timer::stop(Timer::lsp2p_6);
#endif
    fill_dependencies(g, gp, gsz);
  }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::LevelSetTri::p2p_reset () const {
  p2p_done_value_ = ! p2p_done_.empty() ? (p2p_done_[0] + 1) % 2 : 0;
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
set_diag_reciprocal (CrsMatrix& T) {
  for (Int i = 0; i < T.m; ++i) {
    const Size j = T.ir[i+1]-1;
    T.d[j] = 1/T.d[j];
  }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
LevelSetTri::init_numeric (const CrsMatrix& T) {
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    del(t_[tid].m);
    t_[tid].m = get_matrix_p(T, ps_[tid]);
    set_diag_reciprocal(*t_[tid].m); }
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
LevelSetTri::reinit_numeric (const CrsMatrix& T) {
  bool nthreads_ok;
  std::string msg;
# pragma omp parallel
  do {
    nthreads_ok = check_nthreads(t_.size(), omp_get_num_threads(), msg);
    if ( ! nthreads_ok) break;
    const int tid = omp_get_thread_num();
    Thread& t = t_[tid];
    assert(t.m);
    const std::vector<Int>& p = ps_[tid];
    for (size_t i = 0; i < p.size(); ++i) {
      const Int r = p[i], nc = static_cast<Int>(T.ir[r+1] - T.ir[r]);
      assert(nc == static_cast<Int>(t.m->ir[i+1] - t.m->ir[i]));
      memcpy(t.m->d + t.m->ir[i], T.d + T.ir[r], nc*sizeof(*t.m->d));
    }
    set_diag_reciprocal(*t_[tid].m);
  } while (0);
  if ( ! nthreads_ok) throw hts::Exception(msg);
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::Permuter::clear () {
  if (q_ && q_ != p_) { deln(q_); }
  deln(p_);
  deln(scale_);
  deln(px_);
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::Permuter::init (
  const Int n, const bool is_lo, const std::vector<Int>& lsis,
  const std::vector<Int>& dpis, const Int nthreads, const Int max_nrhs,
  const Int* p, const Int* q, const Real* scale)
{
  clear();
  n_ = n;
  max_nrhs_ = max_nrhs;
  is_lo_ = is_lo;
  p_ = q_ = 0; px_ = scale_ = 0;
  try {
    p_ = allocn<Int>(n_);
    px_ = allocn<Real>(n_*max_nrhs);
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

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
Permuter::reinit_numeric (const Real* scale) {
  if (scale)
    memcpy(scale_, scale, n_*sizeof(*scale_));
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
Permuter::reset_max_nrhs (const Int max_nrhs) {
  if (max_nrhs_ == max_nrhs) return;
  max_nrhs_ = max_nrhs;
  deln(px_);
  px_ = Allocnator<Real>(n_*max_nrhs_, "Permuter::reset_max_nrhs").release();
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
Permuter::check_nrhs (const Int nrhs) const {
  if (nrhs > max_nrhs_)
    throw hts::Exception("nrhs is > max_nrhs.");
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::TriSolver::clear () {}

// len(irt) == n-1 because of how the transpose algorithm offsets irt. Partition
// by nonzeros. At exit, start[i] is the starting column for the i'th thread.
template<typename Int, typename Size, typename Real> static void
partition_irt (const Int n, const Size* const irt, const Size nnz,
               const Int nparts, std::vector<Int>& start) {
  const Int nseg = std::min<Int>(nparts, n);
  Int i0 = 1, j0 = 1;
  start.resize(nseg);
  
  for (Int i = i0; i < nseg; ++i) {
    const Real d = ((Real) i / nseg)*nnz;
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

template<typename Int, typename Size, typename Real>
static void transpose (
  const Int m, const Size* const ir, const Int* const jc, const Real* const d,
  Size* const irt, Int* const jct, Real* const dt,
  std::vector<Size>* transpose_perm = 0)
{
  const Size nnz = ir[m];
  if (transpose_perm) transpose_perm->resize(nnz);
  std::vector<Int> start;
# pragma omp parallel
  {
#   pragma omp for schedule(static)
    for (Int i = 0; i <= m; ++i)
      irt[i] = 0;
#   pragma omp single
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
      partition_irt<Int, Size, Real>(m, irt + 2, nnz, omp_get_num_threads(),
                                     start);
    }
    const int tid = omp_get_thread_num();
    if (tid < static_cast<int>(start.size())) {
      const Int
        start_tid = start[tid],
        stop = tid + 1 == static_cast<int>(start.size()) ? m : start[tid+1];
      for (Int r = 0; r < m; ++r) {
        const Size irr = ir[r], jlim = ir[r+1];
        if (stop <= jc[irr] || start_tid > jc[jlim-1]) continue;
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

template<typename Int, typename Size, typename Real>
typename Impl<Int, Size, Real>::ConstCrsMatrix* Impl<Int, Size, Real>::
transpose (const ConstCrsMatrix& T, std::vector<Size>* transpose_perm) {
  const Int n = T.m;
  SparseData sd(n, T.ir[n], "transpose");
  htsimpl::transpose(n, T.ir, T.jc, T.d, sd.ir, sd.jc, sd.d, transpose_perm);
  ConstCrsMatrix* ccm = 0;
  try {
    ccm = new ConstCrsMatrix(n, n, sd.ir, sd.jc, sd.d,
                             opposite<Int, Size, Real>(T.dir), true);
  } catch (...) { throw hts::Exception("transpose failed to alloc."); }
  sd.release();
  return ccm;
}

// inside || {}
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
compose_transpose (const std::vector<Size>& transpose_perm, Partition& p) {
  std::vector<Size> A_idxs_copy(p.A_idxs);
# pragma omp for schedule(static)
  for (size_t i = 0; i < p.A_idxs.size(); ++i)
    p.A_idxs[i] = transpose_perm[A_idxs_copy[i]];
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TriSolver::init (const ConstCrsMatrix* T, Int nthreads, const Int max_nrhs,
                 const bool save_for_reprocess, const Int* p, const Int* q,
                 const Real* scale, const Options& o) {
  Timer::init();
  Timer::start(Timer::total_pre); Timer::start(Timer::setup);
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

    Timer::stop(Timer::setup); Timer::start(Timer::transpose);
    const bool make_transpose = T->dir == ConstCrsMatrix::transpose;
    std::vector<Size> transpose_perm;
    if (make_transpose) {
      T = transpose(*T, save_for_reprocess ? &transpose_perm : 0);
      delete_T = true;
    }
    Timer::stop(Timer::transpose);

    // Determine shape.
    Timer::start(Timer::tolower);
    Shape shape = determine_shape(*T);
    if ( ! shape.is_triangular) throw hts::NotTriangularException();
    if ( ! shape.has_full_diag) throw hts::NotFullDiagonal();
    is_lo_ = shape.is_lower;
    Timer::stop(Timer::tolower);

    // Find level sets.
    Timer::start(Timer::lsetfind);
    LevelSetter lstr;
    const Int lstr_threshold = o.lset_min_size * 
      (o.lset_min_size_scale_with_nthreads ? nthreads_ : 1);
    lstr.init(*T, lstr_threshold, is_lo_, o);
    Timer::stop(Timer::lsetfind);

    Timer::start(Timer::perm);
    if ( ! is_lo_) {
      const ConstCrsMatrix* Tp = delete_T ? T : 0;
      T = permute_to_other_tri(*T);
      if (delete_T) del(Tp);
      delete_T = true;
      lstr.reverse_variable_order(n_);
    }

    // Separate into three blocks: level set, scatter, and data parallel:
    //     [(1)        0
    //      (scatter) (2)],
    // where if is_lo_, (1) = level sets and (2) = data parallel; if ! is_lo_,
    // then the opposite.
    std::vector<Int> lsis, dpis;
    get_idxs(n_, lstr, lsis, dpis);
    if (o.printlvl > 0)
      std::cout << "n " << n_ << " |lsis| " << lsis.size() << " |dpis| "
                << dpis.size() << "\n";
    Timer::stop(Timer::perm);
    if (is_lo_) {
      // 1. Level-scheduling block.
      Timer::start(Timer::lsetinit);
      { PermVec lsis_pv(T->m, lsis);
        get_matrix_pp_with_covers_all(*T, lsis_pv, p_[0]);
        sort(p_[0]); }
      lst_.init_lsets(lstr, save_for_reprocess);
      lst_.init(*p_[0].cm, 0, 0, p_[0].cm->m, in);
      lst_.update_permutation(lsis, p_[0]);
      Timer::stop(Timer::lsetinit); Timer::start(Timer::lsetinitp2p);
      lst_.p2p_init();
      Timer::stop(Timer::lsetinitp2p);
      Timer::start(Timer::dpb); Timer::start(Timer::dpb_getmatrix);
      {
        PermVec dpis_pv(T->m, dpis), lsis_pv(T->m, lsis);
        get_matrix_p_qp_with_covers_all(*T, dpis_pv, lsis_pv, p_[1]);
        Timer::stop(Timer::dpb_getmatrix); Timer::start(Timer::dpb_sort);
        sort(p_[1]);
        Timer::stop(Timer::dpb_sort);
      }
      if (p_[1].cm->m > 0) {
        Timer::start(Timer::dpb_tinit);
        // 2. No MVP block. It's in the data-parallel block.
        // 3. Data-parallel block (+ MVP block).
        const Int mvp_nc = p_[1].cm->n - p_[1].cm->m;
        t_.init(*p_[1].cm, 0, 0, p_[1].cm->m, in, mvp_nc);
        Timer::stop(Timer::dpb_tinit);
      }
      Timer::stop(Timer::dpb);
    } else {
      Timer::start(Timer::dpb); Timer::start(Timer::dpb_getmatrix);
      PermVec dpis_pv(T->m, dpis);
      get_matrix_pp_with_covers_all(*T, dpis_pv, p_[1]);
      Timer::stop(Timer::dpb_getmatrix); Timer::start(Timer::dpb_sort);
      sort(p_[1]);
      Timer::stop(Timer::dpb_sort);
      if (p_[1].cm->m > 0) {
        Timer::start(Timer::dpb_tinit);
        // 3. Data-parallel block.
        t_.init(*p_[1].cm, 0, 0, p_[1].cm->m, in);
        // 2. No MVP block. It's in the level scheduling block.
        Timer::stop(Timer::dpb_tinit);
      }
      Timer::stop(Timer::dpb);
      // 1. Level-scheduling block (+ MVP block).
      Timer::start(Timer::lsetinit);
      { PermVec lsis_pv(T->m, lsis);
        get_matrix_p_qp_with_covers_all(*T, lsis_pv, dpis_pv, p_[0]);
        sort(p_[0]); }
      lst_.init_lsets(lstr, save_for_reprocess);
      lst_.set_mvp_block_nc(p_[1].cm->m);
      lst_.init(*p_[0].cm, 0, 0, p_[0].cm->m, in);
      lst_.update_permutation(lsis, p_[0]);
      Timer::stop(Timer::lsetinit); Timer::start(Timer::lsetinitp2p);
      lst_.p2p_init();
      Timer::stop(Timer::lsetinitp2p);
    }
    Timer::start(Timer::setup);
    if (delete_T) del(T);
    Timer::stop(Timer::setup); Timer::start(Timer::for_reprocess);
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
    } else {
      for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].clear();
    }
    Timer::stop(Timer::for_reprocess); Timer::start(Timer::perm);
    perm_.init(n_, is_lo_, lsis, dpis, nthreads_, max_nrhs, p, q, scale);
    Timer::stop(Timer::perm);
  } catch (...) {
    if (delete_T) del(T);
    throw;
    restore_num_threads(nthreads_state);
  }
  Timer::start(Timer::setup);
  restore_num_threads(nthreads_state);
  Timer::stop(Timer::setup); Timer::stop(Timer::total_pre);
#ifdef TIME
  Timer::print_setup();
#endif
}

// Reinitialize numbers, but keep the same structures.
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TriSolver::reinit_numeric (const ConstCrsMatrix* T, const Real* r) {
  Timer::init();
  Timer::start(Timer::total_re); Timer::start(Timer::numthr);
  NumThreads nthreads_state;
  set_num_threads(nthreads_, nthreads_state);
  Timer::stop(Timer::numthr); Timer::start(Timer::numpart);
  for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].alloc_d();
  repartition_into_2_blocks(p_, *T);
  Timer::stop(Timer::numpart); Timer::start(Timer::numls);
  lst_.reinit_numeric(*p_[0].cm);
  Timer::stop(Timer::numls);
  if (p_[1].cm->m > 0) {
    //todo Tighten up some of these init_numeric impls. Might want to do
    // reinit_numeric like for lst.
    t_.init_numeric(*p_[1].cm);
  }
  for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].clear_d();
  Timer::start(Timer::numperm);
  if (r) perm_.reinit_numeric(r);
  Timer::stop(Timer::numperm); Timer::start(Timer::numthr);
  restore_num_threads(nthreads_state);
  Timer::stop(Timer::numthr); Timer::stop(Timer::total_re);
#ifdef TIMENUM
  Timer::print_numeric();
#endif
}

template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::
TriSolver::reset_max_nrhs (const Int max_nrhs) {
  perm_.reset_max_nrhs(max_nrhs);
  t_.reset_max_nrhs(max_nrhs);
}

//> Solve code.

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
OnDiagTri::solve (const Real* b, const Int ldb, Real* x, const Int ldx,
                  const Int nrhs) const {
  if (d_) {
    return t_.empty() ? solve_dense(b, ldb, x, ldx, nrhs) :
      solve_dense_inv(b, ldb, x, ldx, nrhs);
  }
  if (m_) return solve_spars(b, ldb, x, ldx, nrhs);
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
OnDiagTri::solve_dense (const Real* b, const Int ldb, Real* x, const Int ldx,
                        const Int nrhs) const {
  for (Int irhs = 0; ; ) {
    for (Int j = 0, k = 0; j < this->n_; ++j) {
      Real a = b[j];
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
template<typename Int, typename Size, typename Real>
void Impl<Int, Size, Real>::OnDiagTri::
solve_dense_inv (const Real* b, const Int ldb, Real* x, const Int ldx,
                 const Int nrhs) const {
  const int tid = omp_get_thread_num();
  if (tid >= nthreads()) return;
  const Thread& t = t_[tid];
  for (Int irhs = 0; ; ) {
    for (Int r = t.r0, rlim = t.r0 + t.nr, k = 0; r < rlim; ++r) {
      Real a = 0;
      for (Int c = 0; c < r+1; ++c, ++k)
        a += b[c]*t.d[k];
      x[r] = a;
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldb;
  }
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
OnDiagTri::solve_spars (const Real* b, const Int ldb, Real* x, const Int ldx,
                        const Int nrhs) const {
  const CrsMatrix& T = *m_;
  const Int m = T.m;
  const Size* const ir = T.ir;
  const Int* const jc = T.jc;
  const Real* const d = T.d;
  for (int irhs = 0; ; ) {
    for (int r = 0; r < m; ++r) {
      const Size
        rp_rp1 = ir[r+1],
        jlim = rp_rp1 - 1;
      Real a = b[r];
      for (Size j = ir[r]; j < jlim; ++j)
        a -= x[jc[j]] * d[j];
      x[r] = a * d[rp_rp1 - 1];
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldb;
  }
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
SerialBlock::n1Axpy (const Real* x, const Int ldx, const Int nrhs,
                     Real* y, const Int ldy) const {
  if ( ! d_) return;
  if (ir_) n1Axpy_spars(x + coff_, ldx, nrhs, y + roff_, ldy);
  else n1Axpy_dense(x + coff_, ldx, nrhs, y + roff_, ldy);
}

template<typename Int, typename Size, typename Real>
inline void SerialBlock_n1Axpy_spars (
  const Int nr_, const Int nc_, const Size* const ir_, const Int* const jc_,
  const Real* const d_, const Real* x, const Int ldx, const Int nrhs,
  Real* y, const Int ldy)
{
  for (Int k = 0; ; ) {
    Size iri = ir_[0];
    for (Int i = 0; i < nr_; ++i) {
      const Size irip1 = ir_[i+1];
      const Int N = static_cast<Int>(irip1 - iri);
      if (N == 0) continue;
      Real a = 0; {
        const Real* const d = d_ + iri;
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

#ifdef HAVE_SHYLUHTS_MKL
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

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::SerialBlock::
n1Axpy_spars (const Real* x, const Int ldx, const Int nrhs,
              Real* y, const Int ldy) const {
  assert(ir_);
  if (ir_[nr_] == 0) return;
  SerialBlock_n1Axpy_spars(nr_, nc_, ir_, jc_, d_, x, ldx, nrhs, y, ldy);
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::SerialBlock::
n1Axpy_dense (const Real* x, const Int ldx, const Int nrhs,
              Real* y, const Int ldy) const {
  assert(d_);
#if defined(HAVE_SHYLUHTS_BLAS) || defined(HAVE_SHYLUHTS_MKL)
  gemm<Real>('t', 'n', nr_, nrhs, nc_, -1, d_, nc_, x, ldx, 1, y, ldy);
#else
  for (Int g = 0; ; ) {
    for (Int i = 0, k = 0; i < nr_; ++i) {
      Real a = 0;
      for (Int j = 0; j < nc_; ++j, ++k) a += d_[k]*x[j];
      y[i] -= a;
    }
    if (++g == nrhs) break;
    x += ldx;
    y += ldy;
  }
#endif
}

// inside || {}
template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
TMatrix::n1Axpy (const Real* x, const Int ldx, const Int nrhs, Real* y,
                 const Int ldy, const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return;
  if ((size_t) id >= bs_.size()) return;
  bs_[id].n1Axpy(x + coff_, ldx, nrhs, y + ros_[id], ldy);
}

// inside || {}
template<typename Int, typename Size, typename Real>
inline Real* Impl<Int, Size, Real>::
Permuter::from_outside (const Real* x, const Int nrhs, Int ldx) const {
  if ( ! ldx) ldx = n_;
  const int tid = omp_get_thread_num();
  const Int i0 = part_[tid], i1 = part_[tid+1];
  Real* ppx = px_;
  const Real* px = x;
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
template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::Permuter::
to_outside (Real* x, const Int nrhs, const Real a, const Real b,
            Int ldx) const {
  if ( ! ldx) ldx = n_;
  const int tid = omp_get_thread_num();
  const Int i0 = part_[tid], i1 = part_[tid+1];
  Real* ppx = px_;
  Real* px = x;
  if (b != 1) {
    if (a) {
      for (int k = 0; ; ) {
        for (Int i = i0; i < i1; ++i) {
          Real* const pxqi = px + q_[i];
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
    if (a) {
      for (int k = 0; ; ) {
        for (Int i = i0; i < i1; ++i) {
          Real* const pxqi = px + q_[i];
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
template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
LevelSetTri::solve (const Real* b, Real* x, const Int ldx,
                    const Int nrhs) const {
  if (nlvls_ == 0) return;
  const int tid = omp_get_thread_num();
  const Thread& t = t_[tid];
  const std::vector<Int>& lsp = t.lsp;
  const CrsMatrix* const T = t.m;
  const Size* const ir = T->ir;
  const Int* const jc = T->jc;
  const Real* const d = T->d;
  const std::vector<Int>& p = t.p;
  const Int lsp_size_m2 = nlvls_ - 1;
  const std::vector<Int>& p2p_depends_p = t.p2p_depends_p;
  const std::vector<Size>& p2p_depends = t.p2p_depends;
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
      // Math. 'volatile' protects the synchronization variable p2p_done_, but
      // technically nothing protects x from being inconsistent. In practice, it
      // probably always is. Nonetheless, flush.
#     pragma omp flush
      for (Int r = mvp_block_nc_ + p[ils];
           i < lsp_ilsp1;
           ++i, ++r) {
        const Size jlim = ir[i+1] - 1;
        Real a = b[r];
        for ( ; j < jlim; ++j)
          a -= x[jc[j]] * d[j];
        x[r] = a * d[j++];
      }
#     pragma omp flush
      if (ils == lsp_size_m2) break;
      // This thread and level is done.
      p2p_done_[ls_p2p_sub2ind(ils, tid)] = p2p_done_value;
    }
    x += ldx;
    b += ldx;
    // Increment the done indicator.
    ++p2p_done_value;
#   pragma omp barrier
  }
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
rbwait (volatile p2p_Done* const s_done, const Size* s_ids,
        const Int* const s_idx, const Int i, const p2p_Done done_symbol) {
  const Int si = s_idx[i], si1 = s_idx[i+1];
  if (si == si1) return;
  const Size* id = s_ids + si;
  const Size* const idn = s_ids + si1;
  while (id != idn) {
    volatile p2p_Done* const d = s_done + *id;
    while (*d != done_symbol) ;
    ++id;
  }
  // Make sure x is updated.
# pragma omp flush
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::RecursiveTri::
ondiag_solve (const OnDiagTri& t, Real* x, const Int ldx, const Int nrhs,
              const int tid, const Int step, volatile Int* const t_barrier,
              volatile Int* const inv_tri_done) const {
  const Int nthreads = t.nthreads();
  if (nthreads == 1) {
    t.solve(x, ldx, x, ldx, nrhs);
    *t_barrier = step;
  } else {
    // Solve T wrk_ = x.
    t.solve(x, ldx, wrk_.data(), t.n(), nrhs);
    { // Wait for the block row MVPs to finish.
      const Int done = (step << 1);
      inv_tri_done[tid] = done;
#   pragma omp flush
      for (Int i = 0; i < nthreads; ++i)
        while (inv_tri_done[i] < done) ;
    }
    // Copy wrk_ to x.
    const Int row_start = t.block_row_start(tid), nr = t.block_nr(tid);
    for (Int irhs = 0; irhs < nrhs; ++irhs)
      memcpy(x + irhs*ldx + row_start,
             wrk_.data() + irhs*t.n() + row_start,
             nr*sizeof(Real));
    { // Wait for the memcpy's to finish.
      const Int done = (step << 1) + 1;
      inv_tri_done[tid] = done;
#     pragma omp flush
      //todo Not every thread necessarily needs this on-diag tri's solution, but
      // our dep graph doesn't encode that yet.
      for (Int i = 0; i < nthreads; ++i)
        while (inv_tri_done[i] < done) ;
    }
    if (tid == 0)
      *t_barrier = step;
#   pragma omp flush
  }
}

// inside || {}
template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
RecursiveTri::solve (const Real* b, Real* x, const Int ldx,
                     const Int nrhs) const {
  if (nd_.t.empty()) return;
  assert(x == b);
  const int tid = omp_get_thread_num();
  Int os = 0;
  const Int sn = static_cast<Int>(nd_.s.size());
  Real* x_osi, * x_os;
  volatile Int* const t_barrier = &nd_.t_barrier;
  volatile Int* const inv_tri_done = nd_.inv_tri_done.data();
  volatile p2p_Done* const s_done = nd_.s_done.data();
  const Size* const t_ids = nd_.t_ids.data();
  const Int* const t_idx = nd_.t_idx.data();
  const Size* const s_ids = nd_.s_ids[tid].empty() ? 0 : nd_.s_ids[tid].data();
  const Int* const s_idx = nd_.s_idx[tid].empty() ? 0 : nd_.s_idx[tid].data();
  const p2p_Done done_symbol = nd_.done_symbol;
  { const OnDiagTri& t = nd_.t[0];
    if (tid < t.nthreads())
      ondiag_solve(t, x, ldx, nrhs, tid, 0, t_barrier, inv_tri_done); }
  if ( ! nd_.os.empty()) {
    os += nd_.t[0].n();
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
      s_done[rb_p2p_sub2ind(i, tid)] = done_symbol;
#     pragma omp flush
    }
    if (tid < t.nthreads()) {
      rbwait(s_done, t_ids, t_idx, i, done_symbol);
      ondiag_solve(t, x_os, ldx, nrhs, tid, i+1, t_barrier, inv_tri_done);
#     pragma omp flush
    }
    os += t.n();
    x_osi = x + nd_.os[i+1];
    x_os = x + os;
  }
}

template<typename Int, typename Size, typename Real>
inline void Impl<Int, Size, Real>::
TriSolver::solve (const Real* b, const Int nrhs, Real* x, const Real alpha,
                  const Real beta, const Int ldb, const Int ldx) const {
  t_.p2p_reset();
  lst_.p2p_reset();
  perm_.check_nrhs(nrhs);
  NumThreads nthreads_save;
  set_num_threads(nthreads_, nthreads_save);
  bool nthreads_ok;
  std::string msg;
# pragma omp parallel
  do {
    nthreads_ok = check_nthreads(nthreads_, omp_get_num_threads(), msg);
    if ( ! nthreads_ok) break;
    Real* px = perm_.from_outside(b, nrhs, ldb);
    if (is_lo_) {
#     pragma omp barrier
      lst_.solve(px, px, n_, nrhs);
      // No barrier needed here because lst_.solve does it.
      if (t_.n()) t_.solve(px, px, n_, nrhs);
#     pragma omp barrier
      perm_.to_outside(x, nrhs, alpha, beta, ldx);
    } else {
      if (t_.n()) {
#       pragma omp barrier
        t_.solve(px, px, n_, nrhs);
      }
#     pragma omp barrier
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
template<typename Int, typename Size, typename Real>
void trisolve_serial (
  const typename HTS<Int, Size, Real>::CrsMatrix& T, Real* ix, const Int nrhs,
  bool is_lo, const Int* p, const Int* q, const Real* scale, Real* w)
{
  assert(( ! p && ! q) || w);

  for (int irhs = 0; irhs < nrhs; ++irhs) {
    // x = x./scale optionally.
    if (scale)
      for (Int r = 0; r < T.m; ++r) ix[r] /= scale[r];

    // x = x(p) optionally.
    Real* x;
    if (p) {
      x = w;
      for (Int r = 0; r < T.m; ++r) x[r] = ix[p[r]];
    } else
      x = ix;

    // Solve.
    if (is_lo) {
      for (Int r = 0; r < T.m; ++r) {
        const Size rp_rp1 = T.ir[r+1];
        for (Size j = T.ir[r]; j < rp_rp1 - 1; ++j)
          x[r] -= x[T.jc[j]] * T.d[j];
        x[r] /= T.d[rp_rp1 - 1];
      }
    } else {
      for (Int r = T.m - 1; r >= 0; --r) {
        const Size rp_r = T.ir[r];
        for (Size j = rp_r + 1; j < T.ir[r+1]; ++j)
          x[r] -= x[T.jc[j]] * T.d[j];
        x[r] /= T.d[rp_r];
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

template<typename Int, typename Size, typename Real>
typename HTS<Int, Size, Real>::CrsMatrix* HTS<Int, Size, Real>::
make_CrsMatrix (const Int nrow, const Size* rowptr, const Int* col,
                const Real* val, const bool make_transpose) {
  static const bool Int_is_signed = static_cast<Int>(-1) < 1;
  if ( ! Int_is_signed)
    throw hts::Exception(
      "HTS<Int, Size, Scalar> cannot be used if Int is unsigned.");
  typedef typename htsimpl::Impl<Int, Size, Real>::ConstCrsMatrix CCM;
  CrsMatrix* m;
  try {
    m = new CrsMatrix(nrow, rowptr, col, val,
                      make_transpose ? CCM::transpose : CCM::forward);
  } catch (...) { throw hts::Exception("make_CrsMatrix failed to alloc."); }
  return m;
}

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
delete_CrsMatrix (CrsMatrix* T) { htsimpl::del(T); }

template<typename Int, typename Size, typename Real>
HTS<Int, Size, Real>::
PreprocessArgs::PreprocessArgs ()
  : T(0), max_nrhs(1), nthreads(-1), save_for_reprocess(false), p(0), q(0),
    scale_rhs(0), scale_rhs_by_division(true), scale_solution(0),
    scale_solution_by_division(true), options(0)
{}

template<typename Int, typename Size, typename Real>
typename HTS<Int, Size, Real>::Impl* HTS<Int, Size, Real>::
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
    htsimpl::Impl<Int, Size, Real>::set_options(*a.options, i->o);
  try {
    i->ts.init(a.T, a.nthreads, a.max_nrhs, a.save_for_reprocess, a.p, a.q,
               a.scale_rhs, i->o);
  } catch (hts::Exception& e) {
    htsimpl::del(i);
    throw;
  }
  return i;
}

template<typename Int, typename Size, typename Real>
typename HTS<Int, Size, Real>::Impl* HTS<Int, Size, Real>::
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

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
reprocess_numeric (Impl* impl, const CrsMatrix* T, const Real* r) {
  impl->ts.reinit_numeric(T, r);
}

template<typename Int, typename Size, typename Real>
bool HTS<Int, Size, Real>::
is_lower_tri (const Impl* impl) { return impl->ts.is_lower_tri(); }

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
delete_Impl (Impl* impl) {
  htsimpl::del(impl);
}

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
reset_max_nrhs (Impl* impl, const Int max_nrhs) {
  impl->ts.reset_max_nrhs(max_nrhs);
}

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
solve_serial (const CrsMatrix* T, const bool is_lo, Real* xb, const Int nrhs,
              const Int* p, const Int* q, const Real* r, Real* w) {
  htsimpl::trisolve_serial<Int, Size, Real>(*T, xb, nrhs, is_lo, p, q, r, w);
}

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
solve_omp (Impl* impl, Real* x, const Int nrhs, const Int ldx) {
  impl->ts.solve(x, nrhs, x, 0, 1, ldx, ldx);
}

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
solve_omp (Impl* impl, const Real* b, const Int nrhs, Real* x,
           const Int ldb, const Int ldx) {
  impl->ts.solve(b, nrhs, x, 0, 1, ldb, ldx);
}

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
solve_omp (Impl* impl, const Real* b, const Int nrhs, Real* x,
           const Real alpha, const Real beta,
           const Int ldb, const Int ldx) {
  impl->ts.solve(b, nrhs, x, alpha, beta, ldb, ldx);
}

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
print_options (const Impl* i, std::ostream& os) {
  htsimpl::Impl<Int, Size, Real>::print_options(i->o, os);
}

template<typename Int, typename Size, typename Real>
void HTS<Int, Size, Real>::
set_level_schedule_only (Options& o) {
  o.min_lset_size = 0;
}

template<typename Int, typename Size, typename Real>
HTS<Int, Size, Real>::Options::Options () {
  min_dense_density = 0.75;
  levelset_block_size = 1;
  lset_min_size_scale_with_nthreads = false;
  profile = false;
  print_level = 0;
  lset_max_bad_fraction = 0.01;
  min_lset_size = lset_min_size_scale_with_nthreads ? 1 : 10;
#ifdef __MIC__
  min_block_size = 256;
  min_parallel_rows = 64;
#else
  min_block_size = 128;
  min_parallel_rows = 32;
#endif
  pp_min_block_size = 256;
}

} // namespace Experimental

#endif
