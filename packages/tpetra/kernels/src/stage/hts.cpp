#include <omp.h>
#ifdef USE_MKL
# include <mkl.h>
#endif

#include <sys/time.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include "hts.hpp"
#include "hts_impl.hpp"

#ifdef USE_MM_MALLOC
# include <xmmintrin.h>
# include <new>
#endif

// Convenience try-catch wrapper for new.
#define ALLOC(alloc_line, catch_code, msg) do {         \
    try {                                               \
      alloc_line;                                       \
    } catch (...) {                                     \
      catch_code;                                       \
      throw Exception(msg ": failed to allocate.");     \
    }                                                   \
  } while (0)

static const int parfor_static_size = 20;

template<typename T>
inline void compress (std::vector<T>& v) { std::vector<T>(v).swap(v); }

typedef int blas_int;
#ifndef NO_BLAS
extern "C" {
  void dgemv_(
    const char*, const blas_int*, const blas_int*, const double*, const double*,
    const blas_int*, const double*, const blas_int*, const double*, double*,
    const blas_int*);
  void dgemm_(
    const char*, const char*, const blas_int*, const blas_int*, const blas_int*,
    const double*, const double*, const blas_int*, const double*,
    const blas_int*, const double*, double*, const blas_int*);
}
#endif

template<typename T> static void gemv(
  const char trans, const blas_int m, const blas_int n, const T alpha,
  const T* a, const blas_int lda, const T* x, const blas_int incx, const T beta,
  T* y, const blas_int incy);
template<typename T> static void gemm(
  const char transa, const char transb, const blas_int m, const blas_int nrhs,
  const blas_int n, const T alpha, const T* a, const blas_int lda, const T* b,
  const blas_int ldb, const T beta, T* c, const blas_int ldc);
template<typename T> static void trtri(
  const char uplo, const blas_int n, T* a, const blas_int lda);

template<> inline void
gemm<double> (
  const char transa, const char transb, blas_int m, const blas_int nrhs,
  blas_int n, const double alpha, const double* a, const blas_int lda,
  const double* b, const blas_int ldb, const double beta, double* c,
  const blas_int ldc)
{
#ifndef NO_BLAS
  dgemm_(&transa, &transb, &m, &nrhs, &n, &alpha, a, &lda, b, &ldb, &beta, c,
         &ldc);
#else
  assert(0);
#endif
}

#ifdef USE_MKL
// sparse A * dense x
inline void hts_mkl_dcsrmm (
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

// sparse T \ dense x
inline void hts_mkl_dcsrsm (
  const bool is_lower, const MKL_INT m, const double* d, const MKL_INT* ir,
  const MKL_INT* jc, const double* b, const MKL_INT ldb, double* x,
  const MKL_INT ldx, int nrhs)
{
  char transa = 'N';
  char
    uplo = is_lower ? 'L' : 'U',
    diag = 'N';
  for (int k = 0; k < nrhs; ++k)
    mkl_cspblas_dcsrtrsv(
      &uplo, &transa, &diag, const_cast<MKL_INT*>(&m),
      const_cast<double*>(d), const_cast<MKL_INT*>(ir),
      const_cast<MKL_INT*>(jc), const_cast<double*>(b + ldb*k),
      x + ldx*k);
}
#endif

namespace {
template<typename T> inline void touch (T* const p, const size_t n) {
  // 1 KB should be a safe lower bound on page size. Touch enough to touch every
  // page; I don't think there's any need to touch more memory than that. On
  // the KNC, first-touch doesn't matter.
#ifndef __MIC__
  for (size_t i = 0; i < n; i += 1024 / sizeof(T))
    p[i] = 0;
#endif
}

#ifdef USE_MM_MALLOC
#define MEM_ALIGN 64
template<typename T> inline T*
allocn (const size_t n, const bool first_touch = false) throw (std::bad_alloc) {
  T* p = (T*) _mm_malloc(n*sizeof(T), MEM_ALIGN);
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
#define MEM_ALIGN 16
template<typename T> inline T*
allocn (const size_t n, const bool first_touch = false) throw (std::bad_alloc) {
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

template<typename T> T square (const T& x) { return x*x; }
} // namespace

namespace Experimental {
namespace hts {
namespace impl {

static void set_options (const hts::Options& os, Options& od) {
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

Options::Options () { set_options(hts::Options(), *this); }

void Options::print (std::ostream& os) const {
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
#ifdef NO_BLAS
  os << " NO_BLAS";
#endif
#ifdef USE_MKL
  os << " USE_MKL";
#endif
}

void print_options (const Options& o, std::ostream& os) {
  print_compiletime_options(os);
  os << std::endl;
  o.print(os);
  os << std::endl;
}

ConstCrsMatrix::~ConstCrsMatrix () {
  if ( ! deallocate_) return;
  deln_const(ir); deln_const(jc); deln_const(d);
}

CrsMatrix::~CrsMatrix () {
  if ( ! deallocate_) return;
  deln_const(ir); deln_const(jc); deln_const(d);
}

template<typename T> static T* vec2arr (const std::vector<T>& v) {
  if (v.empty()) return NULL;
  T* a;
  ALLOC(a = allocn<T>(v.size()), ;, "vec2arr");
  memcpy(a, &v[0], sizeof(T)*v.size());
  return a;
}

inline void Partition::alloc_d () {
  assert( ! cm->d);
  ALLOC(cm->d = allocn<Real>(cm->ir[cm->m]), ;, "Partition::alloc_d");
}

inline void Partition::clear () { del(cm); }

inline void Partition::clear_d () { deln(cm->d); }

namespace {
struct SparseData {
  Int* ir, * jc;
  Real* d;
  SparseData (const Int m, const Int nnz, const bool touch = false) {
    ir = jc = 0;
    d = 0;
    try {
      ir = allocn<Int>(m+1, touch);
      if (nnz > 0) {
        jc = allocn<Int>(nnz, touch);
        d = allocn<Real>(nnz, touch);
      }
      ir[0] = 0;
      ir[m] = nnz;
    } catch (...) {
      free();
      throw Exception("SparseData failed to allocate.");
    }
  }
  void free () {
    deln(ir);
    deln(jc);
    deln(d);
  }
};

inline void
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

struct NumThreads {
  int omp, mkl;
  bool mkl_dynamic;
};

inline void set_num_threads (const int nthreads, NumThreads& save) {
  save.omp = omp_get_max_threads();
#ifdef USE_MKL
  save.mkl = mkl_get_max_threads();
  save.mkl_dynamic = mkl_get_dynamic();
#endif
  omp_set_num_threads(nthreads);
#ifdef USE_MKL
  // We never use MKL threading.
  mkl_set_dynamic(0);
  mkl_set_num_threads(1);
#endif
}

inline void restore_num_threads (const NumThreads& save) {
  // This ruins performance on especially the KNC if the number of threads is
  // actually switching. Need to think about this more. Right now, the interface
  // does not promise OMP state will remain the same.
  return;

  omp_set_num_threads(save.omp);
#ifdef USE_MKL
  mkl_set_num_threads(save.mkl);
  mkl_set_dynamic(save.mkl_dynamic);
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

void throw_if_nthreads_not_ok (const int nthreads) {
  int nr;
# pragma omp parallel
  nr = omp_get_num_threads();
  std::string msg;
  if ( ! check_nthreads(nthreads, nr, msg))
    throw Exception(msg);
}

// Return i such that ir[i] is the first value >= c. If no such i exists,
// return n.
inline Int find_first (const Int* ir, const Int n, const Int c) {
  return c == 0 ? 0 : std::lower_bound(ir, ir+n, c) - ir;
}

// Return the number of nonzeros in row r that are in [c_first, c_last). The
// corresponding indices, relative to the start of the row, are i_first:i_last.
inline Int find_first_and_last (const Int* const ir, const Int r,
                                const Int* const jc,
                                const Int c_first, const Int c_last,
                                Int& i_first, Int& i_last) {
  assert(c_last >= c_first);
  const Int
    iri = ir[r],
    irip1 = ir[r+1],
    n = irip1 - iri;
  i_first = find_first(jc + iri, n, c_first);
  if (i_first == n) {
    i_last = i_first;
    return 0;
  }
  i_last = i_first + find_first(jc + iri + i_first, n - i_first, c_last);
  // A return value of n - i_first is OK.
  return i_last - i_first;
}

// Crop the submatrix A(b) such that A(cb) has no 0 border.
Int crop_matrix (const CrsMatrix& T, const Box& b, Box& cb) {
  cb.r0 = -1;
  Int r1 = -1;
  cb.c0 = b.c0 + b.nc;
  Int c1 = b.c0;
  Int nnz = 0;
  //todo ||ize if not nested.
  for (Int r = b.r0, lrow = 0; r < b.r0 + b.nr; ++r, ++lrow) {
    Int i_first, i_last;
    const Int cnt = find_first_and_last(T.ir, r, T.jc, b.c0, b.c0 + b.nc,
                                        i_first, i_last);
    if (cnt) {
      nnz += cnt;
      const Int irr = T.ir[r];
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

typedef std::vector<Int> VI;

// Decide how many level sets to keep.
Int decide_level_set_max_index (const std::vector<Int>& N, const Int size_thr,
                                const Options& o) {
  Int N_end = (Int) N.size();
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
void alloc_lsets (const Int lsmi, const Int sns, const std::vector<Int>& level,
                  const std::vector<Int>& n, LevelSetter::LevelSets& lsets) {
  if (lsmi < 0) return;
  const Int Lm_sr = level.size();
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

namespace locrsrow {
// Serial impl for reference.
Int schedule_serial (const ConstCrsMatrix& L, const Int sns,
                     std::vector<Int>& w) {
  // Eq. 18 in Y. Saad's 1989 SIAM J Sci Stat Comput paper.
  Int max_level = -1;
  if (sns == 1) {
    w.resize(L.m, -1);
    for (Int r = 0; r < L.m; ++r) {
      Int level = -1;
      for (Int j = L.ir[r]; j < L.ir[r+1]; ++j)
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
        for (Int j = L.ir[r]; j < L.ir[r+1]; ++j) {
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

Int schedule_sns1 (const ConstCrsMatrix& L, std::vector<Int>& w,
                   const Options& o) {
  const Int
    nthreads = omp_get_max_threads(),
    blksz = nthreads*((o.pp_min_block_size + nthreads)/nthreads),
    rows_per_thread = std::max(1, blksz / nthreads);
  if (blksz > L.m)
    return schedule_serial(L, 1, w);
  std::vector<Int> frontier(blksz);
  for (Int i = 0; i < blksz; ++i) frontier[i] = L.ir[i];
  w.resize(L.m);
  Int max_level = -1;
  volatile Int done = -1;
# pragma omp parallel
  {
#   pragma omp for schedule(static, parfor_static_size)
    for (Int i = 0; i < L.m; ++i) w[i] = -1;
    const Int* const ir = L.ir;
    const Int* const jc = L.jc;
    for (Int c = 0; c < L.m; c += blksz) {
      const Int tlim = std::min(c + blksz, L.m);
      // On-diag serial triangle.
#     pragma omp single nowait
      {
        for (Int r = c; r < tlim; ++r) {
          // w[r] contains the max level seen so far.
          Int level = w[r];
          for (Int j = frontier[r - c], jlim = ir[r+1]; j < jlim; ++j)
            level = std::max(level, w[jc[j]]);
          ++level;
          w[r] = level;
          max_level = std::max(max_level, level);
        }
        done = c;
      }
      if (tlim == L.m) break;
      // Off-diag parallel block row.
      const Int rlim = std::min(tlim + blksz, L.m);
      while (done != c) ;
#     pragma omp for schedule(static, rows_per_thread)
      for (Int r = tlim; r < rlim; ++r) {
        Int level = -1;
        const Int jlim = ir[r+1];
        for (Int j = ir[r]; j < jlim; ++j) {
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

Int schedule (const ConstCrsMatrix& L, const Int sns, std::vector<Int>& w,
              const Options& o) {
  assert(L.m > 0);
  if (sns == 1) return schedule_sns1(L, w, o);
  const Int
    Lm_sr = L.m / sns,
    blksz = (o.pp_min_block_size + sns) / sns,
    bnr = sns*blksz,
    nthreads = omp_get_max_threads(),
    rows_per_thread = std::max(1, (blksz + nthreads) / nthreads);
  if (blksz > Lm_sr)
    return schedule_serial(L, sns, w);
  std::vector<Int> frontier(bnr);
  for (Int i = 0; i < bnr; ++i) frontier[i] = L.ir[i];
  w.resize(Lm_sr);
  Int max_level = -1;
# pragma omp parallel
  {
#   pragma omp for schedule(static, parfor_static_size)
    for (Int i = 0; i < Lm_sr; ++i) w[i] = -1;
    const Int* const ir = L.ir;
    const Int* const jc = L.jc;
    for (Int sc = 0; sc < Lm_sr; sc += blksz) {
      const Int stlim = std::min(sc + blksz, Lm_sr);
      // On-diag serial triangle.
#     pragma omp single nowait
      {
        const Int c = sns*sc;
        for (Int sr = sc, r = c; sr < stlim; ++sr) {
          Int level = w[sr];
          for (Int i = 0; i < sns; ++i, ++r)
            for (Int j = frontier[r - c], jlim = ir[r+1]; j < jlim; ++j) {
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
        srlim = std::min(stlim + blksz, Lm_sr),
        tlim = sns*stlim;
#     pragma omp for schedule(static, rows_per_thread)
      for (Int sr = stlim; sr < srlim; ++sr) {
        Int level = -1;
        for (Int i = 0, r = sns*sr; i < sns; ++i, ++r) {
          const Int jlim = ir[r+1];
          for (Int j = ir[r]; j < jlim; ++j) {
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
} // namespace locrsrow

void find_row_level_sets_Lcrs (const ConstCrsMatrix& L, const Int sns,
                               Int size_thr, LevelSetter::LevelSets& lsets,
                               const Options& o) {
  assert(L.m % sns == 0);

  std::vector<Int> w;
  const Int max_level = locrsrow::schedule_serial(L, sns, w);
 
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
namespace upcrscol {
//todo Need to ||ize, but this equivalent to ||izing a CSC L trisolve, which is
// harder than CSR L trisolve. Not sure yet how to do this well.

Int schedule_serial (const ConstCrsMatrix& U, const Int sns,
                     std::vector<Int>& w) {
  Int max_level = -1;
  if (sns == 1) {
    w.resize(U.m, -1);
    for (Int r = 0; r < U.m; ++r) {
      ++w[r];
      const Int level = w[r];
      max_level = std::max(level, max_level);
      for (Int j = U.ir[r] + 1; j < U.ir[r+1]; ++j) {
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
        for (Int j = U.ir[r] + 1; j < U.ir[r+1]; ++j) {
          const Int sc = U.jc[j] / sns;
          w[sc] = std::max(w[sc], level);
        }
    }
  }
  return max_level;
}
} // namespace upcrscol

void find_col_level_sets_Ucrs (const ConstCrsMatrix& U, const Int sns,
                               Int size_thr, LevelSetter::LevelSets& lsets,
                               const Options& o) {
  assert(U.m % sns == 0);

  std::vector<Int> w;
  const Int max_level = upcrscol::schedule_serial(U, sns, w);
 
  // Count level set sizes.
  std::vector<Int> n(max_level+1);
  for (size_t i = 0; i < w.size(); ++i)
    ++n[w[i]];
  // Cutoff.
  const Int lsmi = decide_level_set_max_index(n, size_thr, o);
  // Fill lsets.
  alloc_lsets(lsmi, sns, w, n, lsets);
}

inline void find_level_sets (
  const ConstCrsMatrix& T, const Int sns, const Int size_thr, const bool is_lo,
  LevelSetter::LevelSets& lsets, const Options& o)
{
  if (is_lo)
    find_row_level_sets_Lcrs(T, sns, size_thr, lsets, o);
  else
    find_col_level_sets_Ucrs(T, sns, size_thr, lsets, o);
}
} // namespace

void LevelSetter::init (const ConstCrsMatrix& T, const Int size_thr,
                        const bool is_lo, const Options& o) {
  lsets_.clear();
  is_lo_ = is_lo;
  // Guard against an invalid setting.
  ls_blk_sz_ = T.m % o.ls_blk_sz == 0 ? o.ls_blk_sz : 1;
  find_level_sets(T, ls_blk_sz_, size_thr, is_lo_, lsets_, o);
}

const VI& LevelSetter::lset (const size_t i) const {
  return is_lo_ ? lsets_[i] : lsets_[lsets_.size() - i - 1];
}

void LevelSetter::reverse_variable_order (Int n) {
  --n;
  for (size_t i = 0; i < lsets_.size(); ++i) {
    VI& ls = lsets_[i];
    for (size_t j = 0; j < ls.size(); ++j)
      ls[j] = n - ls[j];
    std::reverse(ls.begin(), ls.end());
  }
}

namespace {
CrsMatrix* get_matrix_p (const CrsMatrix& A, const std::vector<Int>& p) {
  const Int n = static_cast<Int>(p.size());
  Int nnz = 0;
  for (size_t i = 0; i < p.size(); ++i) {
    const Int r = p[i];
    nnz += A.ir[r+1] - A.ir[r];
  }

  SparseData sd(n, nnz, true);
  for (size_t i = 0; i < p.size(); ++i) {
    const Int r = p[i], nc = A.ir[r+1] - A.ir[r];
    sd.ir[i+1] = sd.ir[i] + nc;
    memcpy(sd.jc + sd.ir[i], A.jc + A.ir[r], nc*sizeof(*sd.jc));
    memcpy(sd.d + sd.ir[i], A.d + A.ir[r], nc*sizeof(*sd.d));
  }

  CrsMatrix* cm;
  ALLOC(cm = new CrsMatrix(n, A.n, sd.ir, sd.jc, sd.d, true),
        sd.free(),
        "get_matrix_p");
  return cm;
}

ConstCrsMatrix* permute_to_other_tri (const ConstCrsMatrix& U) {
  const Int n = U.m, nnz = U.ir[n];
  SparseData sd(n, nnz);
# pragma omp parallel for schedule(static)
  for (Int k = 1; k <= n; ++k)
    sd.ir[k] = nnz - U.ir[n-k];
# pragma omp parallel for schedule(static)
  for (Int k = 0; k < nnz; ++k) {
    const Int i = nnz - k - 1;
    sd.jc[k] = n - U.jc[i] - 1;
    sd.d[k] = U.d[i];
  }
  ConstCrsMatrix* ccm;
  ALLOC(ccm = new ConstCrsMatrix(n, n, sd.ir, sd.jc, sd.d, true),
        sd.free(),
        "permute_to_other_tri");
  return ccm;
}

void get_data_parallel_idxs (const Int n, const LevelSetter& lstr,
                             std::vector<Int>& dpis) {
  std::vector<char> dpisb(n, 1);
  Int cnt = n;
  for (Int i = 0; i < lstr.size(); ++i) {
    const std::vector<Int>& lset = lstr.lset(i);
    for (size_t j = 0; j < lset.size(); ++j) {
      dpisb[lset[j]] = 0;
      --cnt;
    }
  }
  dpis.resize(cnt);
  for (size_t i = 0, dk = 0; i < dpisb.size(); ++i)
    if (dpisb[i]) dpis[dk++] = i;
}

// Partition 1:n into lsis, the set of level scheduled rows, and dpis, the set
// of data-|| rows.
void get_idxs (const Int n, const LevelSetter& lstr, std::vector<Int>& lsis,
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

struct Shape {
  const bool is_lower, is_triangular;
  const bool has_full_diag; // Valid only if is_triangular.
  Shape (const bool is_lower, const bool is_triangular,
         const bool has_full_diag = false)
    : is_lower(is_lower), is_triangular(is_triangular),
      has_full_diag(has_full_diag) {}
};

Shape determine_shape (const ConstCrsMatrix& A) {
  bool tri_determined = false, is_tri = true, is_lower = true,
    has_full_diag = true;
# pragma omp parallel for schedule(static, parfor_static_size)
  for (Int r = 0; r < A.m; ++r) {
    bool diag_fnd = false;
    for (Int j = A.ir[r]; j < A.ir[r+1]; ++j) {
      const Int c = A.jc[j];
      if (c != r) {
        if ( ! tri_determined) {
#         pragma omp critical
          { is_lower = c < r;
            tri_determined = true; }
        }
        if ((is_lower && c > r) || ( ! is_lower && c < r)) {
#         pragma omp critical
          is_tri = false;
        }
      } else diag_fnd = true;
    }
    if ( ! diag_fnd) {
#     pragma omp critical
      has_full_diag = false;
    }
  }
  // If ! tri_determined, then T must be a diag matrix. Can treat as lower,
  // which is is_lower's default value.
  return Shape(is_lower, is_tri, has_full_diag);
}

class PermVec {
  const std::vector<Int>& p_;
  std::vector<Int> pi_;
public:
  // p is a set of indices into an n-vector.
  PermVec(const Int n, const std::vector<Int>& p)
    : p_(p)
  {
    pi_.resize(n, -1);
    for (size_t i = 0; i < p_.size(); ++i) pi_[p_[i]] = i;
  }
  size_t size () const { return p_.size(); }
  // Index into p.
  Int get (const Int i) const { return p_[i]; }
  // Does p contain k?
  bool has (const Int k) const { return pi_[k] >= 0; }
  // Return i for p[i] == k.
  Int to_block (const Int k) const { return pi_[k]; }
};

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
void get_matrix_pp_with_covers_all (const ConstCrsMatrix& A, const PermVec& pv,
                                    Partition& p) {
  // Count nnz.
  Int nnz = 0;
# pragma omp parallel for reduction(+:nnz)
  for (size_t i = 0; i < pv.size(); ++i) {
    const Int k = pv.get(i);
    nnz += A.ir[k+1] - A.ir[k];
  }
  SparseData sd(pv.size(), nnz);
  for (size_t ipv = 0; ipv < pv.size(); ++ipv) {
    const Int
      i = pv.get(ipv),
      nc = A.ir[i+1] - A.ir[i];
    sd.ir[ipv+1] = sd.ir[ipv] + nc;
  }
  assert(sd.ir[pv.size()] == nnz);
  p.A_idxs.resize(nnz);
  const Int ipv_lim = pv.size();
# pragma omp parallel for schedule(static, parfor_static_size)
  for (Int ipv = 0; ipv < ipv_lim; ++ipv) {
    const Int
      i = pv.get(ipv),
      nc = A.ir[i+1] - A.ir[i],
      Aj0 = A.ir[i],
      Bj0 = sd.ir[ipv];
    for (Int j = 0; j < nc; ++j) {
      const Int Aj = Aj0 + j, Bj = Bj0 + j;
      sd.jc[Bj] = pv.to_block(A.jc[Aj]);
      sd.d[Bj] = A.d[Aj];
      p.A_idxs[Bj] = Aj;
    }
  }
  ALLOC(p.cm = new CrsMatrix(pv.size(), pv.size(), sd.ir, sd.jc, sd.d, true),
        sd.free(),
        "get_matrix_pp_with_covers_all");
}

// Extract B = A(p,s), s = [q p], given that A(p,~s) is empty.
void get_matrix_p_qp_with_covers_all (
  const ConstCrsMatrix& A, const PermVec& pv, const PermVec& qv, Partition& p)
{
  Int nnz = 0;
# pragma omp parallel for reduction(+:nnz)
  for (size_t i = 0; i < pv.size(); ++i) {
    const Int r = pv.get(i);
    nnz += A.ir[r+1] - A.ir[r];
  }

  SparseData sd(pv.size(), nnz);
  for (size_t ipv = 0, lim = pv.size(); ipv < lim; ++ipv) {
    const Int
      r = pv.get(ipv),
      nc = A.ir[r+1] - A.ir[r];
    sd.ir[ipv+1] = sd.ir[ipv] + nc;
  }
  assert(sd.ir[pv.size()] == nnz);

  p.A_idxs.resize(nnz);
  const Int ipv_lim = (Int) pv.size();
# pragma omp parallel for schedule(static, parfor_static_size)
  for (Int ipv = 0; ipv < ipv_lim; ++ipv) {
    const Int
      i = pv.get(ipv),
      nc = A.ir[i+1] - A.ir[i],
      Aj0 = A.ir[i],
      Bj0 = sd.ir[ipv];
    for (Int j = 0; j < nc; ++j) {
      const Int
        Aj = Aj0 + j, Bj = Bj0 + j,
        Ac = A.jc[Aj];
      sd.jc[Bj] = qv.has(Ac) ? qv.to_block(Ac) : qv.size() + pv.to_block(Ac);
      sd.d[Bj] = A.d[Aj];
      p.A_idxs[Bj] = Aj;
    }
  }

  ALLOC(p.cm = new CrsMatrix(pv.size(), A.n, sd.ir, sd.jc, sd.d, true),
        sd.free(),
        "get_matrix_p_qp_with_covers_all");
}

template<typename T> struct SortEntry {
  Int i, A_idx; T d;
  SortEntry () {}
  bool operator< (const SortEntry& se) const { return i < se.i; }
};

void sort (Partition& p) {
  const int nthreads = omp_get_max_threads();
  CrsMatrix& A = *p.cm;

  //todo ||ize
  int max_nc = 0;
  for (Int r = 0; r < A.m; ++r)
    max_nc = std::max(max_nc, A.ir[r+1] - A.ir[r]);
  std::vector< SortEntry<Real> > sess(nthreads * max_nc);

# pragma omp parallel for schedule(static, 1)
  for (Int r = 0; r < A.m; ++r) {
    const int tid = omp_get_thread_num();
    const Int irr = A.ir[r], irrp1 = A.ir[r+1], nc = irrp1 - irr;
    SortEntry<Real>* ses = &sess[tid*max_nc];
    for (Int j = 0; j < nc; ++j) {
      const Int Aj = irr + j;
      ses[j].i = A.jc[Aj];
      ses[j].d = A.d[Aj];
      ses[j].A_idx = p.A_idxs[Aj];
    }
    std::sort(ses, ses + nc);
    for (Int j = 0; j < nc; ++j) {
      const Int Aj = irr + j;
      A.jc[Aj] = ses[j].i;
      A.d[Aj] = ses[j].d;
      p.A_idxs[Aj] = ses[j].A_idx;
    }
  }
}

void partition_into_2_blocks (
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

inline void copy_partition (const ConstCrsMatrix& A, Partition& p,
                            const bool is_lo) {
  const Int ilim = (Int) p.A_idxs.size();
  if (is_lo) {
#   pragma omp parallel for schedule(static)
    for (Int i = 0; i < ilim; ++i)
      p.cm->d[i] = A.d[p.A_idxs[i]];
  } else {
    const Int nnz = A.ir[A.m];
#   pragma omp parallel for schedule(static)
    for (Int i = 0; i < ilim; ++i)
      p.cm->d[i] = A.d[nnz - p.A_idxs[i] - 1];
  }
}

inline void
repartition_into_2_blocks (Partition* const p, const ConstCrsMatrix& A,
                           const bool is_lo) {
  for (int i = 0; i < 2; ++i)
    if (p[i].cm) copy_partition(A, p[i], is_lo);
}
} // namespace

inline void CrsSegmenter::
count_nnz_by_row_loop (const Int i, std::vector<Int>& rcnt) {
  Int i_first, i_last;
  rcnt[i] = find_first_and_last(A_.ir, r0_ + i, A_.jc, c0_, c0_ + nc_,
                                i_first, i_last);
}

void CrsSegmenter::count_nnz_by_row (std::vector<Int>& rcnt) {
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

void CrsSegmenter::init_nnz (const std::vector<Int>& rcnt) {
  const Int nseg = p_.size() - 1;
  nnz_.resize(nthreads_, 0);
  for (Int i = 0; i < nseg; ++i)
    for (Int j = p_[i] - p_[0]; j < p_[i+1] - p_[0]; ++j)
      nnz_[i] += rcnt[j];
}

// p partitions rows. It may not be a perfect partitioning in the sense of
// optimal given that each thread must handle a mutually exclusive set of
// rows. Attempt to remove spikes in #rows.
inline void
smooth_spikes (const std::vector<Int>& rcnt, std::vector<Int>& p,
               std::vector<Int>& nnz, const bool ignore_0) {
  const Int n = (Int) p.size() - 1;
  if (n <= 1) return;
  Int first = ignore_0 ? 1 : 0;
  for (int it = 0; ; ++it) {
    // Find a spike.
    Int s_val = nnz[first], s_i = 0;
    for (Int i = first + 1; i < n; ++i)
      if (nnz[i] > s_val) {
        s_i = i;
        s_val = nnz[i];
      }

    // See if it decreases the spikiness to move a row from one set to another.
    Int rnnz, d, idx, s;
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
        idx_r = p[s_i+1] - 1 - p[0],
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

void CrsSegmenter::segment () {
  assert(nr_ % ls_blk_sz_ == 0);
  const Int nseg = std::min(nthreads_, nr_ / ls_blk_sz_);
  if (nseg == 0) {
    assert(nr_ == 0);
    p_.resize(1, r0_);
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

  p_.resize(nseg + 1, r0_);
  if (nseg > 0) {
    if (cs_rcnt.back() == 0) {
      for (Int i = 1; i < nseg; ++i)
        p_[i] = r0_ + i;
    } else {
      Int i0 = 1, j0 = 1, nparts = nseg, start = 0;
      if (tid_empty) {
        // Adjust so that the first thread gets nothing and the nonzeros are
        // balanced among the rest.
        p_[1] = r0_;
        ++i0;
        --nparts;
        ++start;
      }
      for (Int i = i0; i < nseg; ++i) {
        const Real d = ((Real) (i - start) / nparts)*cs_rcnt.back();
        Int j = (std::upper_bound(cs_rcnt.begin() + j0, cs_rcnt.end(), d)
                 - cs_rcnt.begin());
        if (d - cs_rcnt[j-1] > cs_rcnt[j] - d) ++j;
        j = std::min(j, (Int) cs_rcnt.size() - nseg + i);
        if (j < j0) {
          // If not all the threads will have work, let the earlier threads have
          // the work.
          j = std::min(j0, (Int) cs_rcnt.size());
          assert(r0_ + j*ls_blk_sz_ == p_[i-1] + 1);
        }
        p_[i] = r0_ + j*ls_blk_sz_;
        j0 = j + 1;
      }
    }
  }
  p_[nseg] = r0_ + nr_;

  init_nnz(rcnt);
  if (ls_blk_sz_ == 1) smooth_spikes(rcnt, p_, nnz_, tid_empty);
}

void TMatrix::clear () {
  bs_.clear();
  ros_.clear();
  is_empty_ = true;
}

void
TMatrix::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
               const InitInfo& in, const Int block_0_nnz_os,
               const int tid_offset) {
  init_metadata(A, r0, c0, nr, nc, in, block_0_nnz_os, tid_offset);
  init_memory(in);
# pragma omp parallel
  init_numeric(A, omp_get_thread_num());
}

void
TMatrix::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
               const InitInfo& in, const CrsSegmenter& seg) {
  init_metadata(A, r0, c0, nr, nc, in, seg);
  init_memory(in);
# pragma omp parallel
  init_numeric(A, omp_get_thread_num());
}

void TMatrix::init_metadata (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
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
  const Int nseg = std::max(1, std::min((Int) in.nthreads, nt));
  CrsSegmenter seg(A, r0, c0, nr, nc, nseg, 1, block_0_nnz_os);

  if (nnz)
    init_metadata_with_seg(A, r0, c0, nr, nc, roff, in, seg);
  else
    is_parallel_ = false;
}

void TMatrix::init_metadata (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                             const InitInfo& in, const CrsSegmenter& seg) {
  clear();
  nr_ = nr;
  init_metadata_with_seg(A, r0, c0, nr, nc, 0, in, seg);
}

void TMatrix::
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
    const Int n = (Int) p.size() - 1;
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

void TMatrix::init_memory (const InitInfo& in) {
# pragma omp parallel
  { const int tid = omp_get_thread_num(), id = tid - tid_os_;
    if (id >= 0 && id < (int) bs_.size())
      bs_[id].init_memory(in);
  }
}

// possibly inside || {}
void TMatrix::init_numeric (const CrsMatrix& A, const int tid) {
  const int id = tid - tid_os_;
  if (id >= 0 && id < (int) bs_.size())
    bs_[id].init_numeric(A);
}

inline void TMatrix::reinit_numeric (const CrsMatrix& A) {
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

inline Int TMatrix::block_r0 (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return bs_[id].r0();
}

inline Int TMatrix::block_nr (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return bs_[id].nr();
}

inline const SerialBlock* TMatrix::block (const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return &bs_[id];
}
inline SerialBlock* TMatrix::block (const int tid) {
  const int id = tid - tid_os_;
  if (id < 0) return 0;
  return &bs_[id];
}

std::ostream& operator<< (std::ostream& os, const TMatrix& m) {
  if (m.bs_.empty()) {
    os << "  TMatrix empty\n";
    return os;
  }
  os << "  TMatrix (" << std::setw(5) << m.bs_[0].r0_ << ", "
     << std::setw(5) << m.bs_[0].c0_ << ") "
     << std::setw(5) << m.nr_ << ":\n";
  Int nnz = 0;
  for (size_t i = 0; i < m.bs_.size(); ++i)
    nnz += m.bs_[i].nnz();
  for (size_t i = 0; i < m.bs_.size(); ++i)
    os << "  " << std::setw(3) << m.tid_os_ + i
       << std::setw(5) << m.bs_[i].nr() << " "
       << std::setw(10) << m.bs_[i].nnz() << " "
       << (double) m.bs_[i].nnz() / nnz << "\n";
  return os;
}

inline void SerialBlock::clear () {
  if (deallocate_) {
    deln(ir_);
    deln(jc_);
    deln(d_);
  }
  ir_ = jc_ = 0;
  roff_ = coff_ = 0;
  d_ = 0;
}

void SerialBlock::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                        const InitInfo& in) {
  init_metadata(A, r0, c0, nr, nc, in);
  init_memory(in);
  init_numeric(A);
}

void SerialBlock::init_metadata (const CrsMatrix& A, Int r0, Int c0, Int nr,
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

void SerialBlock::init_memory (const InitInfo& in) {
  if (nr_ == 0 || nc_ == 0) return;
  if (is_dense_) {
    // First touch is measurably required.
    ALLOC(d_ = allocn<Real>(nr_*nc_, true), ;, "DenseSerialBlock::init");
  } else {
    try {
      ir_ = allocn<Int>(nr_+1, true);
      if (nnz_ == 0) {
        ir_[nr_] = 0;
        return;
      }
      jc_ = allocn<Int>(nnz_, true);
      d_ = allocn<Real>(nnz_, true);
    } catch (...) {
      clear();
      throw Exception("SparseSerialBlock::init: failed to allocate.");
    }
  }
}

inline void SerialBlock::init_numeric (const CrsMatrix& A) {
  if (is_dense_) memset(d_, 0, nr_*nc_*sizeof(*d_));
  reinit_numeric(A);
}

inline void SerialBlock::reinit_numeric (const CrsMatrix& A) {
  if (ir_) reinit_numeric_spars(A);
  else reinit_numeric_dense(A);
}

void SerialBlock::reinit_numeric_dense (const CrsMatrix& A) {
  const Int ilim = r0_ + nr_;
  for (Int i = r0_; i < ilim; ++i) {
    const Int
      lrow = i - r0_,
      iri = A.ir[i],
      irip1 = A.ir[i+1];
    for (Int j = iri + find_first(A.jc + iri, irip1 - iri, c0_);
         j < irip1; ++j) {
      const Int lcol = A.jc[j] - c0_;
      if (lcol >= nc_) break;
      const Int k = nc_*lrow + lcol;
      d_[k] = A.d[j];
    }
  }
}

inline void SerialBlock::reinit_numeric_spars (const CrsMatrix& A) {
  Int k = 0, pk = 0, nzero_seq = 0, prev_entry = 0;
  ir_[0] = 0;
  const Int ilim = r0_ + nr_;
  for (Int i = r0_; i < ilim; ++i) {
    const Int iri = A.ir[i], irip1 = A.ir[i+1];
    for (Int g = iri + find_first(A.jc + iri, irip1 - iri, c0_);
         g < irip1; ++g) {
      const Int jcg = A.jc[g];
      if (jcg >= c0_ + nc_) break;
      jc_[k] = jcg - c0_;
      d_[k] = A.d[g];
      ++k;
    }
    ir_[i - r0_ + 1] = k;
#ifndef USE_MKL
    // Delta to jump over 0 rows. If there is more than one zero row
    // sequentially, fill the first zero row's ir entry with a negative number
    // whose negation gives the number of rows to jump forward.
    if (k == pk)
      ++nzero_seq;
    else {
      if (nzero_seq > 1)
        ir_[prev_entry+1] = 1 - nzero_seq;
      nzero_seq = 0;
      pk = k;
      prev_entry = i - r0_ + 1;
    }
#endif
  }
  assert(k == ir_[nr_]);
}

inline Int ntri (const int n) { return (n*(n + 1))/2; }

inline Int
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

OnDiagTri::OnDiagTri ()
  : c0_(0), nnz_(0), d_(0), m_(0)
#ifdef USE_MKL
  , b_(0)
#endif
  , dense_(true)
{}

void OnDiagTri::clear () {
  deln(d_);
  del(m_);
#ifdef USE_MKL
  deln(b_);
#endif
  t_.clear();
}

void OnDiagTri::init (const CrsMatrix& T, const Int r0, const Int c0,
                      const Int n, const InitInfo& in) {
  init_metadata(T, r0, c0, n, in);
  init_memory(in);
  init_numeric(T);
}

void OnDiagTri::init (const Int r0, const Int c0, const Int n) {
  clear();
  n_ = n; r0_ = r0; c0_ = c0;
}

inline bool is_dense_tri (const CrsMatrix& T, const Int r0, const Int c0,
                          const Int n, const InitInfo& in, Int& nnz) {
  nnz = count_nnz_lotri(T, r0, c0, n);
  return nnz >= in.min_dense_density*ntri(n);
}

void OnDiagTri::init_metadata (const CrsMatrix& T, const Int r0, const Int c0,
                               const Int n, const InitInfo& in) {
  clear();
  n_ = n; r0_ = r0; c0_ = c0;
  if ( ! n_) {
    nnz_ = 0;
    return;
  }
  dense_ = is_dense_tri(T, r0, c0, n, in, nnz_);
  if (dense_) inv_init_metadata(in);
}

void OnDiagTri::init_memory (const InitInfo& in) {
  if ( ! nnz_) return;
  if (dense_) {
    ALLOC(d_ = allocn<Real>(ntri(n_), true), ;, "OnDiagTri::init dense");
    if ( ! t_.empty()) inv_init_memory();
  } else {
    SparseData sd(n_, nnz_);
    ALLOC(m_ = new CrsMatrix(n_, n_, sd.ir, sd.jc, sd.d, true),
          sd.free(); clear(),
          "OnDiagTri::init_memory");
#ifdef USE_MKL
    ALLOC(b_ = allocn<Real>(in.max_nrhs * n_, true), clear(),
          "OnDiagTri::init sparse");
#endif
  }
}

void OnDiagTri::init_numeric (const CrsMatrix& T) {
  if (dense_) {
    reinit_numeric(T);
  } else {
    Int* const ir = m_->ir;
    Int* const jc = m_->jc;
    Real* const d = m_->d;
    for (Int grow = r0_; grow < r0_ + n_; ++grow) {
      const Int
        lrow = grow - r0_,
        irg = T.ir[grow],
        irgp1 = T.ir[grow+1];
      ir[lrow+1] = ir[lrow];
      for (Int k = irg + find_first(T.jc + irg, irgp1 - irg, c0_);
           k < irgp1; ++k) {
        const Int lcol = T.jc[k] - c0_;
        if (lcol >= n_) break;
        Int& i = ir[lrow+1];
        jc[i] = lcol;
        d[i] =
#ifndef USE_MKL
          lrow == lcol ? 1/T.d[k] :
#endif
          T.d[k];
        ++i;
      }
    }
    assert(ir[n_] == nnz_);
  }
}

//todo At the cost of doubling the jc-related memory, I could store indices into
// T in init_numeric and then use them here in reinit_numeric so I wouldn't have
// to run find_first. At least record the find_first results.
void OnDiagTri::reinit_numeric (const CrsMatrix& T) {
  if (d_) {
    memset(d_, 0, ntri(n_)*sizeof(*d_));
    Int nnz = 0;
    for (Int grow = r0_; grow < r0_ + n_; ++grow) {
      const Int
        lrow = grow - r0_,
        irg = T.ir[grow],
        irgp1 = T.ir[grow+1];
      for (Int k = irg + find_first(T.jc + irg, irgp1 - irg, c0_);
           k < irgp1; ++k) {
        const Int lcol = T.jc[k] - c0_;
        if (lcol >= n_) break;
        const Int di = (lrow*(lrow + 1))/2 + lcol;
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
    for (Int grow = r0_, i = 0; grow < r0_ + n_; ++grow) {
      const Int
        lrow = grow - r0_,
        irg = T.ir[grow],
        irgp1 = T.ir[grow+1];
      for (Int k = irg + find_first(T.jc + irg, irgp1 - irg, c0_);
           k < irgp1; ++k) {
        const Int lcol = T.jc[k] - c0_;
        if (lcol >= n_) break;
        d[i] =
#ifndef USE_MKL
          lrow == lcol ? 1/T.d[k] :
#endif
          T.d[k];
        ++i;
      }
    }
  }
}

inline int OnDiagTri::nthreads () const {
  return std::max(1, static_cast<int>(t_.size()));
}

OnDiagTri::Thread::~Thread () { deln(d); }

void OnDiagTri::inv_init_metadata (const InitInfo& in) {
  const Int
    nnz = ntri(n_),
    nt = std::min((int) in.nthreads, 1 + nnz / square(in.min_parallel_rows));
  if (nt <= 1) return;

  t_.resize(nt);
  const Real nnz_per_thread = (Real) nnz / nt;
  Int r0 = 0;
  for (Int tid = 0; tid < nt; ++tid) {
    t_[tid].r0 = r0;
    const Int n_max = n_ - r0;
    if (tid+1 == nt)
      t_[tid].nr = n_max;
    else {
      // Solve for n in
      //   ((r0 + n) (r0 + n + 1))/2 - (r0 (r0 + 1))/2 = nnz_per_thread.
      const Int b = 1 + 2*r0;
      t_[tid].nr = std::min(
        (Real) n_max, round(0.5*(std::sqrt(b*b + 8*nnz_per_thread) - b)));
    }
    r0 += t_[tid].nr;
  }
}

void OnDiagTri::inv_init_memory () {
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    if (tid < nthreads()) {
      Thread& t = t_[tid];
      const Int nnz = ntri(t.r0 + t.nr) - ntri(t.r0);
      ALLOC(t.d = allocn<Real>(nnz, true), ;, "OnDiagTri::inv_init_memory");
    } }
}

// T is in row-major compressed dense tri format. The diag of T is the
// reciprocal.
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

void OnDiagTri::inv_reinit_numeric (const CrsMatrix& T) {
  { std::vector<Real> w(n_);
    invert(d_, n_, &w[0]); }
  for (int tid = 0; tid < nthreads(); ++tid) {
    Thread& t = t_[tid];
    const Int nr0 = ntri(t.r0);
    memcpy(t.d, d_ + nr0, (ntri(t.r0 + t.nr) - nr0)*sizeof(*t.d));
  }
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
namespace rt {
inline Int split (const Int n, const Int nthreads) {
  const Int s = n/2;
  if (n <= nthreads) return s;
  const Int smod = s % nthreads;
  if ( ! smod) return s;
  return s + (nthreads - smod);
}

void build_recursive_tri_r (const CrsMatrix& T, const Int r, const Int c,
                            const Int n, const InitInfo& in,
                            std::list<Box>& b) {
  if (n <= in.min_blksz) {
    Int nnz;
    if (n <= in.min_parallel_rows || is_dense_tri(T, r, c, n, in, nnz)) {
      b.push_back(Box(r, c, n, n));
      return;
    }
  }
  const Int n1 = split(n, in.nthreads), n2 = n - n1, r1 = r + n1;
  build_recursive_tri_r(T, r, c, n1, in, b);
  b.push_back(Box(r1, c, n2, n1));
  build_recursive_tri_r(T, r1, c + n1, n2, in, b);
}

void
build_recursive_tri (const CrsMatrix& T, const Int r, const Int c, const Int n,
                     const Int mvp_block_nc, const InitInfo& in,
                     std::vector<Box>& bv) {
  std::list<Box> bl;
  if (mvp_block_nc) bl.push_back(Box(r, c, n, mvp_block_nc));
  build_recursive_tri_r(T, r, c + mvp_block_nc, n, in, bl);
  bv.resize(bl.size());
  Int i = 0;
  for (std::list<Box>::const_iterator it = bl.begin(); it != bl.end(); ++it)
    bv[i++] = *it;
}
} // namespace rt

void
RecursiveTri::init (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                    const InitInfo& in, const Int mvp_block_nc) {
  clear();
  n_ = n; r0_ = r0; nthreads_ = in.nthreads;

  std::vector<Box> bs;
  rt::build_recursive_tri(T, r0, c0, n, mvp_block_nc, in, bs);
  const Int
    nblock = bs.size(),
    ntri = ((mvp_block_nc ? 2 : 1) + nblock) / 2,
    nmvp = ntri - 1;
  nd_.t.resize(ntri);
  nd_.s.resize(nmvp);
  nd_.os.resize(ntri);

  Int ti = 0, si = 0;
  std::vector<Box>::iterator bit = bs.begin();
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
  assert(ti == (Int) nd_.t.size());
  assert(si == (Int) nd_.s.size());

  { // Initialize threading data for the inverse of the on-diag tri.
    Int max_diag_tri = 0, max_nthreads = 0;
    for (size_t i = 0; i < nd_.t.size(); ++i) {
      max_diag_tri = std::max(max_diag_tri, nd_.t[i].n());
      max_nthreads = std::max(max_nthreads, nd_.t[i].nthreads());
    }
    wrk_.resize(max_diag_tri * in.max_nrhs);
    nd_.inv_tri_done.resize(max_nthreads);
  }

  p2p_init();
}

void RecursiveTri::clear () {
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

void RecursiveTri::init_numeric (const CrsMatrix& T) {
  Int n = (Int) nd_.t.size();
# pragma omp parallel for schedule(dynamic)
  for (Int i = 0; i < n; ++i)
    nd_.t[i].reinit_numeric(T);
  n = (Int) nd_.s.size();
  for (Int i = 0; i < n; ++i)
    if ( ! nd_.s[i].empty()) nd_.s[i].reinit_numeric(T);
}

void LevelSetTri::init_lsets (const LevelSetter& lstr,
                              const bool save_for_reprocess) {
  ls_blk_sz_ = lstr.ls_blk_sz();
  save_for_reprocess_ = save_for_reprocess;
  lsp_.resize(lstr.size() + 1);
  lsp_[0] = 0;
  for (Int i = 0; i < lstr.size(); ++i)
    lsp_[i+1] = lsp_[i] + lstr.lset(i).size();
}

void LevelSetTri::init (const CrsMatrix& T, const Int r0, const Int c0,
                        const Int n, const InitInfo& in) {
  n_ = n;
  t_.resize(in.nthreads);
  
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    t_[tid].lsp.push_back(0); }
  
  ps_.resize(in.nthreads);
  for (Int ils = 0; ils < (Int) lsp_.size() - 1; ++ils) {
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
void LevelSetTri::update_permutation (std::vector<Int>& lsis,
                                      const Partition& p) {
  const std::vector<Int> old_lsis = lsis;
  std::vector<Int> q(p.cm->n);
# pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    Thread& t = t_[tid];
    int start = 0;
    for (int i = 0; i < tid; ++i)
      start += t_[i].m->m;
    const Int nlev = (Int) t.lsp.size() - 1;
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
    for (Int i = 0, nnz = t.m->ir[t.m->m]; i < nnz; ++i) {
      Int& jci = t.m->jc[i];
      if (jci < mvp_block_nc_) continue;
      jci = mvp_block_nc_ + q[jci - mvp_block_nc_];
    }
  }
}

#define ls_p2p_sub2ind(lvl, tid) ((nlvls)*(tid) + (lvl))
#define ls_p2p_ind2lvl(e) ((e) % (nlvls))
#define ls_p2p_ind2tid(e) ((e) / (nlvls))

// inside || {}
void LevelSetTri::
find_task_responsible_for_variable (std::vector<p2p::Pair>& pairs) {
  const int tid = omp_get_thread_num();
  // Find the level and tid responsible for variable i.
  Thread& t = t_[tid];
  const std::vector<Int>& lsp = t.lsp;
  const std::vector<Int>& p = t.p;
  for (Int ils = 0, ils_lim = (Int) lsp.size() - 1; ils < ils_lim; ++ils) {
    const Int r0 = p[ils];
    for (Int j = 0, jlim = lsp[ils+1] - lsp[ils]; j < jlim; ++j) {
      const Int i = r0 + j;
      pairs[i].tid = tid;
      pairs[i].lvl = ils;
    }
  }
}

// inside || {}
Int LevelSetTri::
fill_graph (const std::vector<p2p::Pair>& pairs, std::vector<Int>& g,
            std::vector<Int>& gp, std::vector<Int>& wrk) {
  // O(nnz/nthreads) time and O(nnz) space. Reduce the entry-wise graph to the
  // (lvl, tid)-wise dependency graph.
  const int tid = omp_get_thread_num();
  const Int n = (Int) gp.size() - 1;
  Thread& t = t_[tid];
  const std::vector<Int>& lsp = t.lsp;
  const Int nlvls = (Int) lsp.size() - 1;
  // Build the graph g. g[e] is the list of dependencies for (level, tid) e.
  const CrsMatrix* const cm = t.m;
  const Int* const ir = cm->ir;
  const Int* const jc = cm->jc;

  // Get an unused piece of the workspace.
  Int* Te = &wrk[0];
  for (int i = 0; i < tid; ++i) {
    const CrsMatrix* const m = t_[i].m;
    if (m) Te += m->ir[m->m];
  }

  // Count entries in g.
  for (Int ils = 1, ils_lim = (Int) lsp.size() - 1; ils < ils_lim; ++ils) {
    const Int e = ls_p2p_sub2ind(ils, tid);
    Int* const Te0 = Te + ir[lsp[ils]];
    const Int Ten = ir[lsp[ils+1]] - ir[lsp[ils]];
    // Record all e = (lvl, tid) dependencies in a CRS matrix structurally
    // identical to cm.
    for (Int r = lsp[ils], rlim = lsp[ils+1]; r < rlim; ++r)
      for (Int j = ir[r], jlim = ir[r+1]; j < jlim; ++j) {
        const Int c = jc[j] - mvp_block_nc_;
        if (c < 0) {
          Te[j] = -1;
          continue;
        }
        const Int r_tid = pairs[c].tid, r_lvl = pairs[c].lvl;
        const Int ed = ls_p2p_sub2ind(r_lvl, r_tid);
        Te[j] = ed;
      }
    // Sort the values so that it's fast to ...
    std::sort(Te0, Te0 + Ten);
    // ... count the number of unique entries.
    Int k = 0, prev = -1;
    for (Int i = 0; i < Ten; ++i)
      if (Te0[i] != prev) {
        prev = Te0[i];
        const Int r_lvl = ls_p2p_ind2lvl(prev);
        if (r_lvl != ils)
          ++k;
      }
    // Now we know how much space to allocate.
    gp[e+1] = k;
  }

  // Cumsum to get row pointers.
  Int max_gelen = 0;
# pragma omp barrier
# pragma omp master
  for (Int i = 1; i <= n; ++i) {
    if (gp[i] > max_gelen) max_gelen = gp[i];
    gp[i] += gp[i-1];
  }
# pragma omp barrier

  // Fill g. Can reuse the data in Te. Everything in this loop is identical to
  // the previous except that we can now record the unique values.
  for (Int ils = 1, ils_lim = (Int) lsp.size() - 1; ils < ils_lim; ++ils) {
    const Int e = ls_p2p_sub2ind(ils, tid);
    Int* const Te0 = Te + ir[lsp[ils]];
    const Int Ten = ir[lsp[ils+1]] - ir[lsp[ils]];
    Int* const ge_g = &g[gp[e]];
    Int k = 0, prev = -1;
    for (Int i = 0; i < Ten; ++i)
      if (Te0[i] != prev) {
        prev = Te0[i];
        const Int r_lvl = ls_p2p_ind2lvl(prev);
        if (r_lvl != ils)
          ge_g[k++] = prev;
      }
    assert(k == gp[e+1] - gp[e]);
    // Sort for the prune_graph phase.
    std::sort(ge_g, ge_g + k);
  }

  return max_gelen;
}

//pre a is sorted.
//pre b is sorted.
//pre len(mark) == an.
// O(max(an, bn)).
inline void
mark_intersection (const Int* a, const Int an, const Int* b, const Int bn,
                   Int* mark /* len(mark) == an */) {
  Int ai = 0, bi = 0;
  while (ai < an && bi < bn) {
    if (a[ai] < b[bi])
      ++ai;
    else if (a[ai] > b[bi])
      ++bi;
    else {
      mark[ai] = -1;
      ++ai;
      ++bi;
    }
  }
}

// inside || {}
void LevelSetTri::
prune_graph (const std::vector<Int>& gc, const std::vector<Int>& gp,
             std::vector<Int>& g, std::vector<Int>& gsz,
             std::vector<Int>& wrk, const Int max_gelen) {
  // Time complexity is a little complicated. See Park et al 2014. Space is
  // O(#threads max |g(e)|).
  const int tid = omp_get_thread_num();
  const Int nlvls = (Int) t_[0].lsp.size() - 1;
  const Int n = (Int) gsz.size();
  Int* mark = &wrk[tid*max_gelen];
  // I find that it's more efficient to use a parfor than to work on a thread's
  // own e's.
# pragma omp for schedule(static,1)
  for (Int e = 0; e < n; ++e) {
    const Int gcelen = gp[e+1] - gp[e];
    if (gcelen == 0) continue;
    // e's dependencies.
    const Int* const gce = &gc[gp[e]];
    for (Int i = 0, ilim = gcelen; i < ilim; ++i)
      mark[i] = gce[i];
    for (Int ied = 0; ied < gcelen; ++ied) { // For each of e's deps:
      const Int ed = gce[ied], edlvl = ls_p2p_ind2lvl(ied);
      assert(ed >= 0 && ed < n);
      if (edlvl == 0) continue; // No parent deps to check.
      // ed's dependencies.
      const Int* const gced = &gc[gp[ed]];
      const Int gcedlen = gp[ed+1] - gp[ed];
      mark_intersection(gce, gcelen, gced, gcedlen, mark);
    }
    // Insert the pruned set of dependencies.
    Int k = 0;
    Int* const ge = &g[gp[e]];
    const Int etid = ls_p2p_ind2tid(e);
    for (Int i = 0, ilim = gcelen; i < ilim; ++i) {
      const Int ed = mark[i];
      if (ed == -1) continue;
      const Int edtid = ls_p2p_ind2tid(ed);
      if (edtid != etid) ge[k++] = ed;
    }
    gsz[e] = k;
  }
}

// inside || {}
void LevelSetTri::
fill_dependencies (const std::vector<Int>& g, const std::vector<Int>& gp,
                   const std::vector<Int>& gsz) {
  const int tid = omp_get_thread_num();
  Thread& t = t_[tid];
  const std::vector<Int>& lsp = t.lsp;
  const Int nlvls = (Int) lsp.size() - 1;
  // Allocate inside this thread up front. I could probably do this even more
  // efficiently, but fill_dependencies is negligible compared with fill_graph
  // and prune_graph.
  t.p2p_depends_p.reserve(nlvls + 1);
  Int sz = 0;
  for (Int ils = 1, ils_lim = (Int) lsp.size() - 1; ils < ils_lim; ++ils)
    sz += gsz[ls_p2p_sub2ind(ils, tid)];
  t.p2p_depends.reserve(sz);
  // Make the final lists of dependencies.
  t.p2p_depends_p.push_back(0);
  for (Int ils = 1, ils_lim = (Int) lsp.size() - 1; ils < ils_lim; ++ils) {
    const Int e = ls_p2p_sub2ind(ils, tid);
    const Int* const ed = &g[gp[e]];
    const Int edsz = gsz[e];
    // Sort by increasing level number. Heuristic to speed up the
    // synchronization step in p2p_solve. Idea is to do p2p_done_ checking on
    // dependencies higher up in the tree first, as those are likely done
    // sooner.
    std::vector<p2p::SortEntry> es(edsz); //todo Remove memory alloc.
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

void LevelSetTri::p2p_init () {
  if (t_[0].lsp.empty()) return;

  if (p2p_done_.empty()) {
    p2p_done_value_ = 1;
    p2p_done_.resize(t_[0].lsp.size()*t_.size(), 0);
  }
  if (t_[0].lsp.size() <= 1) return;

  Int nnz = 0;
  for (size_t i = 0, n = t_.size(); i < n; ++i) {
    const CrsMatrix* const cm = t_[i].m;
    if (cm) nnz += cm->ir[cm->m];
  }
  const Int n = t_.size() * t_[0].lsp.size();

  // g is a graph. nnz is an upper bound on the memory needed.
  //   g(e) is the set of dependencies (incoming edges) for node e. gp is the
  // pointer to g(e). So g(e) is written g[gp[e]]. gsz is the size of g(e).
  //   A node e represents a pair (level, thread id) that is a task.
  std::vector<Int> g, gp, gsz, gc;
  std::vector<Int> wrk;
  // Thread and level responsible for a variable.
  std::vector<p2p::Pair> pairs;
  ALLOC(g.resize(nnz); gp.resize(n+1, 0); gsz.resize(n); wrk.resize(nnz);
        pairs.resize(n_), ;, "p2p_init");

  Int max_gelen;
# pragma omp parallel
  {
    find_task_responsible_for_variable(pairs);
#   pragma omp barrier
    const Int max_gelen_t = fill_graph(pairs, g, gp, wrk);
#   pragma omp barrier
#   pragma omp master
    { // Keep the original graph.
      ALLOC(gc = g, ; , "p2p_init (gc = g)");
      pairs.clear();
      // Tell all threads; only master's max_gelen_t is valid.
      max_gelen = max_gelen_t;
      // In the unlikely case that max_gelen * #threads > nnz, allocate more
      // workspace.
      const Int space = max_gelen * t_.size();
      if (space > nnz)
        ALLOC(wrk.resize(space), ;, "p2p_init (space)");
    }
#   pragma omp barrier
    prune_graph(gc, gp, g, gsz, wrk, max_gelen);
#   pragma omp barrier
    fill_dependencies(g, gp, gsz);
  }
}

inline void set_diag_reciprocal (CrsMatrix& T) {
  for (Int i = 0; i < T.m; ++i) {
    const Int j = T.ir[i+1]-1;
    T.d[j] = 1/T.d[j];
  }
}

void LevelSetTri::init_numeric (const CrsMatrix& T) {
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    del(t_[tid].m);
    t_[tid].m = get_matrix_p(T, ps_[tid]);
    set_diag_reciprocal(*t_[tid].m); }
}

void LevelSetTri::reinit_numeric (const CrsMatrix& T) {
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
      const Int r = p[i], nc = T.ir[r+1] - T.ir[r];
      assert(nc == t.m->ir[i+1] - t.m->ir[i]);
      memcpy(t.m->d + t.m->ir[i], T.d + T.ir[r], nc*sizeof(*t.m->d));
    }
    set_diag_reciprocal(*t_[tid].m);
  } while (0);
  if ( ! nthreads_ok) throw Exception(msg);
}

inline void Permuter::clear () {
  if (q_ && q_ != p_) { deln(q_); }
  deln(p_);
  deln(scale_);
  deln(px_);
}

void Permuter::init (
  const Int n, const bool is_lo, const std::vector<Int>& lsis,
  const std::vector<Int>& dpis, const Int nthreads, const Int max_nrhs,
  const Int* p, const Int* q, const Real* scale)
{
  clear();
  n_ = n;
  is_lo_ = is_lo;
  p_ = q_ = 0; px_ = scale_ = 0;
  try {
    p_ = allocn<Int>(n_);
    px_ = allocn<Real>(n_*max_nrhs);
    q_ = p || q ? allocn<Int>(n_) : p_;
    if (scale) {
      scale_ = allocn<Real>(n_);
      memcpy(scale_, scale, n_*sizeof(*scale_));
    }
  } catch (...) {
    clear();
    throw Exception("Permuter::init failed to allocate.");
  }

  partition_n_uniformly(n_, nthreads, part_);

  // Load p_ and q_ with internal and possibly user permutations.
  if (p || q) {
    const Int dpis_sz = (Int) dpis.size();
    // Incorporate user's permutations.
    for (Int i = 0, ilim = (Int) lsis.size(), k = is_lo_ ? 0 : dpis_sz;
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
      memcpy(p_ + (is_lo_ ? 0 : dpis.size()), &lsis[0],
             lsis.size()*sizeof(*p_));
    if ( ! dpis.empty())
      memcpy(p_ + (is_lo_ ? n - dpis.size() : 0), &dpis[0],
             dpis.size()*sizeof(*p_));
    if ( ! is_lo_)
      for (Int i = 0; i < n_; ++i) p_[i] = n_ - p_[i] - 1;
  }
}

void TriSolver::clear () {}

// Analyze level sets in terms of size and dependencies on previous level.
void analyze (const ConstCrsMatrix& T, const LevelSetter& lsr) {
  std::vector<Int> rls(T.m, -1); // Level set responsible for row r.
  for (Int ils = 0; ils < lsr.size(); ++ils) {
    const std::vector<Int>& ls = lsr.lset(ils);
    for (size_t i = 0; i < ls.size(); ++i)
      rls[ls[i]] = ils;
  }

  std::vector<Int> cnted(T.m, -1);
  for (Int ils = 0; ils < lsr.size(); ++ils) {
    Int cnt = 0; // Number of unique dependencies in parent level set.
    const std::vector<Int>& ls = lsr.lset(ils);
    for (size_t i = 0; i < ls.size(); ++i) {
      const Int r = ls[i];
      for (Int j = T.ir[r]; j < T.ir[r+1]; ++j) {
        const Int c = T.jc[j];
        if (rls[c] == ils-1 && cnted[c] != ils) {
          ++cnt;
          cnted[c] = ils;
        }
      }
    }
    printf("%5d: %5d %5d\n", ils, (Int) ls.size(), cnt);
  }
}

void TriSolver::init (const ConstCrsMatrix* T, Int nthreads, const Int max_nrhs,
                      const bool save_for_reprocess, const Int* p, const Int* q,
                      const Real* scale, const Options& o) {
  NumThreads nthreads_state;
  set_num_threads(nthreads, nthreads_state);
  throw_if_nthreads_not_ok(nthreads);    
  bool delete_T = false;
  try {
    // Basic size parameters.
    n_ = T->m;
    nthreads_ = nthreads;
    InitInfo in;
    in.nthreads = nthreads_;
    in.min_blksz = o.min_block_size;
    in.min_dense_density = o.min_dense_density;
    in.max_nrhs = max_nrhs;
    in.min_parallel_rows = o.min_parallel_rows;

    // Determine shape.
    Shape shape = determine_shape(*T);
    if ( ! shape.is_triangular) throw NotTriangularException();
    if ( ! shape.has_full_diag) throw NotFullDiagonal();
    is_lo_ = shape.is_lower;

    // Find level sets.
    LevelSetter lstr;
    const Int lstr_threshold = o.lset_min_size * 
      (o.lset_min_size_scale_with_nthreads ? nthreads_ : 1);
    lstr.init(*T, lstr_threshold, is_lo_, o);

    if ( ! is_lo_) {
      T = permute_to_other_tri(*T);
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
    if (is_lo_) {
      // 1. Level-scheduling block.
      { PermVec lsis_pv(T->m, lsis);
        get_matrix_pp_with_covers_all(*T, lsis_pv, p_[0]);
        sort(p_[0]); }
      lst_.init_lsets(lstr, save_for_reprocess);
      lst_.init(*p_[0].cm, 0, 0, p_[0].cm->m, in);
      lst_.update_permutation(lsis, p_[0]);
      lst_.p2p_init();
      { PermVec dpis_pv(T->m, dpis), lsis_pv(T->m, lsis);
        get_matrix_p_qp_with_covers_all(*T, dpis_pv, lsis_pv, p_[1]);
        sort(p_[1]); }
      if (p_[1].cm->m > 0) {
        // 2. No MVP block. It's in the data-parallel block.
        // 3. Data-parallel block (+ MVP block).
        const Int mvp_nc = p_[1].cm->n - p_[1].cm->m;
        t_.init(*p_[1].cm, 0, 0, p_[1].cm->m, in, mvp_nc);
      }
    } else {
      PermVec dpis_pv(T->m, dpis);
      get_matrix_pp_with_covers_all(*T, dpis_pv, p_[1]);
      sort(p_[1]);
      if (p_[1].cm->m > 0) {
        // 3. Data-parallel block.
        t_.init(*p_[1].cm, 0, 0, p_[1].cm->m, in);
        // 2. No MVP block. It's in the level scheduling block.
      }
      // 1. Level-scheduling block (+ MVP block).
      { PermVec lsis_pv(T->m, lsis);
        get_matrix_p_qp_with_covers_all(*T, lsis_pv, dpis_pv, p_[0]);
        sort(p_[0]); }
      lst_.init_lsets(lstr, save_for_reprocess);
      lst_.set_mvp_block_nc(p_[1].cm->m);
      lst_.init(*p_[0].cm, 0, 0, p_[0].cm->m, in);
      lst_.update_permutation(lsis, p_[0]);
      lst_.p2p_init();
    }
    if (delete_T) del(T);
    if (save_for_reprocess) {
      // For 32-bit Int, 64-bit Real, save 2/3 of memory during solves at little
      // extra init_numeric cost.
      for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].clear_d();
    } else {
      for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].clear();
    }
    perm_.init(n_, is_lo_, lsis, dpis, nthreads_, max_nrhs, p, q, scale);
  } catch (...) {
    if (delete_T) del(T);
    throw;
    restore_num_threads(nthreads_state);
  }
  restore_num_threads(nthreads_state);
}

// Reinitialize numbers, but keep the same structures.
void TriSolver::init_numeric (const ConstCrsMatrix* T) {
  NumThreads nthreads_state;
  set_num_threads(nthreads_, nthreads_state);
  for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].alloc_d();
  repartition_into_2_blocks(p_, *T, is_lo_);
  lst_.reinit_numeric(*p_[0].cm);
  if (p_[1].cm->m > 0) {
    //todo Tighten up some of these init_numeric impls. Might want to do
    // reinit_numeric like for lst.
    t_.init_numeric(*p_[1].cm);
  }
  for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].clear_d();
  restore_num_threads(nthreads_state);
}

//> Solve code.

inline void OnDiagTri::solve (const Real* b, const Int ldb, Real* x,
                              const Int ldx, const Int nrhs) const {
  if (d_) {
    return t_.empty() ? solve_dense(b, ldb, x, ldx, nrhs) :
      solve_dense_inv(b, ldb, x, ldx, nrhs);
  }
  if (m_) return solve_spars(b, ldb, x, ldx, nrhs);
}

inline void OnDiagTri::solve_dense (const Real* b, const Int ldb, Real* x,
                                    const Int ldx, const Int nrhs) const {
  for (Int irhs = 0; ; ) {
    for (Int j = 0, k = 0; j < n_; ++j) {
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
void OnDiagTri::
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

inline void OnDiagTri::solve_spars (const Real* b, const Int ldb, Real* x,
                                    const Int ldx, const Int nrhs) const {
#ifdef USE_MKL
  memcpy(b_, b, n_*sizeof(*b_));
  for (int k = 1; k < nrhs; ++k)
    memcpy(b_ + n_*k, b + ldb*k, n_*sizeof(*b_));
  hts_mkl_dcsrsm(true, m_->m, m_->d, m_->ir, m_->jc, b_, n_, x, ldx, nrhs);
#else
  const CrsMatrix& T = *m_;
  const Int m = T.m;
  const Int* const ir = T.ir;
  const Int* const jc = T.jc;
  const Real* const d = T.d;
  for (int irhs = 0; ; ) {
    for (int r = 0; r < m; ++r) {
      const Int
        rp_rp1 = ir[r+1],
        jlim = rp_rp1 - 1;
      Real a = b[r];
      for (Int j = ir[r]; j < jlim; ++j)
        a -= x[jc[j]] * d[j];
      x[r] = a * d[rp_rp1 - 1];
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldb;
  }
#endif
}

inline void
SerialBlock::n1Axpy (const Real* RESTRICT x, const Int ldx, const Int nrhs,
                     Real* RESTRICT y, const Int ldy) const {
  if ( ! d_) return;
  if (ir_) n1Axpy_spars(x + coff_, ldx, nrhs, y + roff_, ldy);
  else n1Axpy_dense(x + coff_, ldx, nrhs, y + roff_, ldy);
}

inline void SerialBlock::
n1Axpy_spars (const Real* RESTRICT x, const Int ldx, const Int nrhs,
              Real* RESTRICT y, const Int ldy) const {
  assert(ir_);
  if (ir_[nr_] == 0) return;
#ifdef USE_MKL
  hts_mkl_dcsrmm(false, nr_, nc_, d_, ir_, jc_, x, ldx, y, ldy, nrhs);
#else
  for (Int k = 0; ; ) {
    Int iri = ir_[0];
    for (Int i = 0; i < nr_; ++i) {
      const Int irip1 = ir_[i+1];
      const Int N = irip1 - iri;
      if (N <= 0) {
        if (N < 0) i -= irip1;
        continue;
      }
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
#endif
}

inline void SerialBlock::
n1Axpy_dense (const Real* RESTRICT x, const Int ldx, const Int nrhs,
              Real* RESTRICT y, const Int ldy) const {
  assert(d_);
#ifndef NO_BLAS
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
inline void
TMatrix::n1Axpy (const Real* x, const Int ldx, const Int nrhs, Real* y,
                 const Int ldy, const int tid) const {
  const int id = tid - tid_os_;
  if (id < 0) return;
  if ((size_t) id >= bs_.size()) return;
  bs_[id].n1Axpy(x + coff_, ldx, nrhs, y + ros_[id], ldy);
}

// inside || {}
inline Real* Permuter::from_outside (const Real* x, const Int nrhs) const {
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
    px += n_;
  }
  return px_;
}

// inside || {}
inline void Permuter::
to_outside (Real* x, const Int nrhs, const Real a, const Real b) const {
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
        px += n_;
      }
    } else {
      for (int k = 0; ; ) {
        for (Int i = i0; i < i1; ++i)
          px[q_[i]] = b*ppx[i];
        if (++k == nrhs) break;
        ppx += n_;
        px += n_;
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
        px += n_;
      }
    } else {
      for (int k = 0; ; ) {
        for (Int i = i0; i < i1; ++i)
          px[q_[i]] = ppx[i];
        if (++k == nrhs) break;
        ppx += n_;
        px += n_;
      }
    }
  }
}

// inside || {}
// I've tested using a separate CRS for each level set and thread. The
// motivation is to make jc and d closer in memory. On the MIC, this is faster;
// on the CPU, it's slower.
inline void LevelSetTri::solve (const Real* b, Real* ix, const Int ldx,
                                const Int nrhs) const {
  const int tid = omp_get_thread_num();
  const CrsMatrix* const T = t_[tid].m;
  Real* x = ix;
  const Int* const ir = T->ir;
  const Int* const jc = T->jc;
  const Real* const d = T->d;
  const std::vector<Int>& p = t_[tid].p;
  const std::vector<Int>& lsp = t_[tid].lsp;
  const Int lsp_size_m1 = (Int) lsp.size() - 1;
  for (Int irhs = 0; ; ) {
    for (Int
           ils = 0,      // level set index
           i = lsp[ils], // this thread's current row
           j = ir[i];    // the usual j index into jc and d
         ils < lsp_size_m1;
         ++ils) {
      const Int lsp_ilsp1 = lsp[ils+1];
      for (Int r = mvp_block_nc_ + p[ils];
           i < lsp_ilsp1;
           ++i, ++r) {
        const Int jlim = ir[i+1] - 1;
        Real a = b[r];
        for ( ; j < jlim; ++j)
          a -= x[jc[j]] * d[j];
        x[r] = a * d[j++];
      }
#     pragma omp barrier
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldx;
  }
}

// inside || {}
inline void
LevelSetTri::p2p_solve (const Real* b, Real* x, const Int ldx,
                        const Int nrhs) const {
  const int tid = omp_get_thread_num();
  const Thread& t = t_[tid];
  const std::vector<Int>& lsp = t.lsp;
  const Int nlvls = (Int) lsp.size() - 1;
  if (nlvls == 0) return;
  const CrsMatrix* const T = t.m;
  const Int* const ir = T->ir;
  const Int* const jc = T->jc;
  const Real* const d = T->d;
  const std::vector<Int>& p = t.p;
  const Int lsp_size_m2 = (Int) lsp.size() - 2;
  const std::vector<Int>& p2p_depends_p = t.p2p_depends_p;
  const std::vector<Int>& p2p_depends = t.p2p_depends;
  p2p::Done p2p_done_value = p2p_done_value_;
  for (Int irhs = 0; irhs < nrhs; ++irhs) {
    for (Int
           ils = 0,      // level set index
           i = lsp[ils], // this thread's current row
           j = ir[i];    // the usual j index into jc and d
         ;
         ++ils) {
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
          volatile p2p::Done* const done = &p2p_done_[p2p_depends[di]];
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
        const Int jlim = ir[i+1] - 1;
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
    // Increment the done indicator. Overflow and wrap is part of its behavior.
    ++p2p_done_value;
#   pragma omp barrier
  }
# pragma omp master
  p2p_done_value_ += nrhs;
}

#define rb_p2p_sub2ind(sid, tid) (sn*(tid) + (sid))
#define rb_p2p_ind2tid(e) ((e) / (sn))
#define rb_p2p_ind2sid(e) ((e) % (sn))

// Form the lists of p2p dependencies. A dependency in this algorithm is of two
// types. One is the usual dependency: a variable has to be solved for before
// the next. The second imposes an ordering to assure there is no write race
// condition when computing a row's dot product in pieces.
void RecursiveTri::p2p_init () {
  // For L, in the case where the LS block has rows > 0 and so there is an MVP
  // block, we rely on the fact that n_ doesn't count the rows associated with
  // the space-holding triangle.
  const Int tn = (Int) nd_.t.size(), sn = (Int) nd_.s.size();
  assert(sn+1 == tn);

  // w[r] is the current MVP block that has precedence for row r.
  std::vector<Int> w(n_, -1);
  // Fast way to get a set of unique elements.
  std::vector<Int> w_cnted(sn * nthreads_, 0);
  Int w_cnted_symbol = 1;

  nd_.s_done.resize(nthreads_*(tn - 1), 0);
  nd_.t_idx.resize(tn);
  std::vector< std::vector<Int> > s_ids(sn * nthreads_);
  for (Int ti = 0; ti < tn; ++ti) {
    // Consider each tri or MVP block in solution order.
    if (ti > 0) { // Tri block. First tri has no dependencies.
      const Tri& t = nd_.t[ti];
      const Int r0 = t.r0(), nr = t.n();
      Int k = ti == 1 ? 0 : nd_.t_idx[ti-1];
      for (Int r = r0, rlim = r0 + nr; r < rlim; ++r) {
        assert(r < n_);
        const Int wr = w[r];
        if (wr >= 0 && w_cnted[wr] != w_cnted_symbol) {
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
        const Int ind = rb_p2p_sub2ind(ti, bi);
        for (Int r = r0, rlim = r0 + nr; r < rlim; ++r) {
          const Int wr = w[r];
          // If row depends on an MVP block, and that block has not yet been
          // recorded in my dependency list, record it.
          if (wr >= 0 && w_cnted[wr] != w_cnted_symbol) {
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
  for (size_t i = 0, ilim = s_ids.size(); i < ilim; ++i) {
    const Int tid = rb_p2p_ind2tid(i);
    nd_.s_idx[tid].push_back(nd_.s_idx[tid].back() + s_ids[i].size());
    for (size_t j = 0, jlim = s_ids[i].size(); j < jlim; ++j)
      nd_.s_ids[tid].push_back(s_ids[i][j]);
  }

  compress(nd_.t_ids);
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    compress(nd_.s_ids[tid]);
    compress(nd_.s_idx[tid]); }
  compress(nd_.s_ids);
}

inline void RecursiveTri::p2p_reset () const {
  nd_.t_barrier = -1;
  for (size_t i = 0; i < nd_.inv_tri_done.size(); ++i)
    nd_.inv_tri_done[i] = -1;
  ++nd_.done_symbol;
}

inline void
rbwait (volatile p2p::Done* const s_done, const Int* s_ids,
        const Int* const s_idx, const Int i, const p2p::Done done_symbol) {
  const Int si = s_idx[i], si1 = s_idx[i+1];
  if (si == si1) return;
  const Int* id = s_ids + si;
  const Int* const idn = s_ids + si1;
  while (id != idn) {
    volatile p2p::Done* const d = s_done + *id;
    while (*d != done_symbol) ;
    ++id;
  }
  // Make sure x is updated.
# pragma omp flush
}

inline void RecursiveTri::
ondiag_solve (const OnDiagTri& t, Real* x, const Int ldx, const Int nrhs,
              const int tid, const Int step, volatile Int* const t_barrier,
              volatile Int* const inv_tri_done) const {
  if (t.nthreads() == 1) {
    t.solve(x, ldx, x, ldx, nrhs);
    *t_barrier = step;
  } else {
    t.solve(x, ldx, &wrk_[0], t.n(), nrhs);
    inv_tri_done[tid] = step;
    for (int i = 0; i < t.nthreads(); ++i)
      while (inv_tri_done[i] < step) ;
    if (tid == 0) {
      //todo Could ||ize.
      for (Int irhs = 0; irhs < nrhs; ++irhs)
        memcpy(x + irhs*ldx, &wrk_[0] + irhs*t.n(), t.n()*sizeof(Real));
      *t_barrier = step;
    }
  }
}

// inside || {}
inline void RecursiveTri::solve (const Real* b, Real* x, const Int ldx,
                                 const Int nrhs) const {
  if (nd_.t.empty()) return;
  assert(x == b);
  const int tid = omp_get_thread_num();
  Int os = 0;
  const Int sn = static_cast<Int>(nd_.s.size());
  Real* x_osi, * x_os;
  volatile Int* const t_barrier = &nd_.t_barrier;
  volatile Int* const inv_tri_done = &nd_.inv_tri_done[0];
  volatile p2p::Done* const s_done = &nd_.s_done[0];
  const Int* const t_ids = &nd_.t_ids[0];
  const Int* const t_idx = &nd_.t_idx[0];
  const Int* const s_ids = nd_.s_ids[tid].empty() ? 0 : &nd_.s_ids[tid][0];
  const Int* const s_idx = nd_.s_idx[tid].empty() ? 0 : &nd_.s_idx[tid][0];
  const p2p::Done done_symbol = nd_.done_symbol;
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

inline void
TriSolver::solve (const Real* b, const Int nrhs, Real* x, const Real alpha,
                  const Real beta) const {
  t_.p2p_reset();
  NumThreads nthreads_save;
  set_num_threads(nthreads_, nthreads_save);
  bool nthreads_ok;
  std::string msg;
# pragma omp parallel
  do {
    nthreads_ok = check_nthreads(nthreads_, omp_get_num_threads(), msg);
    if ( ! nthreads_ok) break;
    Real* px = perm_.from_outside(b, nrhs);
    if (is_lo_) {
#     pragma omp barrier
      lst_.p2p_solve(px, px, n_, nrhs);
      // No barrier needed here because lst_.solve does it.
      if (t_.n()) t_.solve(px, px, n_, nrhs);
#     pragma omp barrier
      perm_.to_outside(x, nrhs, alpha, beta);
    } else {
      if (t_.n()) {
#       pragma omp barrier
        t_.solve(px, px, n_, nrhs);
      }
#     pragma omp barrier
      lst_.p2p_solve(px, px, n_, nrhs);
      // No barrier needed because of lst_.solve.
      perm_.to_outside(x, nrhs, alpha, beta);
    }
  } while (0);
  restore_num_threads(nthreads_save);
  if ( ! nthreads_ok) throw Exception(msg);
}

//< Solve code.
} // namespace impl

namespace {
// Solve T x = x with size(x, 2) = nrhs. If is_lo, T is lower tri, else upper.
void trisolve_serial (const CrsMatrix& T, Real* ix, const Int nrhs, bool is_lo,
                      const Int* p, const Int* q, const Real* scale, Real* w) {
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
        const Int rp_rp1 = T.ir[r+1];
        for (Int j = T.ir[r]; j < rp_rp1 - 1; ++j)
          x[r] -= x[T.jc[j]] * T.d[j];
        x[r] /= T.d[rp_rp1 - 1];
      }
    } else {
      for (Int r = T.m - 1; r >= 0; --r) {
        const Int rp_r = T.ir[r];
        for (Int j = rp_r + 1; j < T.ir[r+1]; ++j)
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
} // namespace

CrsMatrix* make_CrsMatrix (const Int nrow, const Int* rowptr, const Int* col,
                           const Real* val) {
  CrsMatrix* m;
  ALLOC(m = new CrsMatrix(nrow, rowptr, col, val), ;, "make_CrsMatrix");
  return m;
}

void delete_CrsMatrix (CrsMatrix* T) { del(T); }

Impl* preprocess (const CrsMatrix* T, const Int max_nrhs, const Int nthreads,
                  const bool save_for_reprocess, const Int* p, const Int* q,
                  const Real* r, const Options* options)
{
  Impl* impl;
  ALLOC(impl = new Impl(), ;, "preprocess");
  if (options) impl::set_options(*options, impl->o);
  try {
    impl->ts.init(T, nthreads, max_nrhs, save_for_reprocess, p, q, r, impl->o);
  } catch (Exception& e) {
    del(impl);
    throw;
  }
  return impl;
}

void reprocess_numeric (Impl* impl, const CrsMatrix* T) {
  impl->ts.init_numeric(T);
}

bool is_lower_tri (const Impl* impl) { return impl->ts.is_lower_tri(); }

void delete_Impl (Impl* impl) {
  del(impl);
}

void solve_serial (const CrsMatrix* T, const bool is_lo, Real* xb,
                   const Int nrhs, const Int* p, const Int* q, const Real* r,
                   Real* w) {
  trisolve_serial(*T, xb, nrhs, is_lo, p, q, r, w);
}

void solve_omp (Impl* impl, Real* x, const Int nrhs) {
  impl->ts.solve(x, nrhs, x, 0, 1);
}

void solve_omp (Impl* impl, const Real* b, const Int nrhs, Real* x) {
  impl->ts.solve(b, nrhs, x, 0, 1);
}

void solve_omp (Impl* impl, const Real* b, const Int nrhs, Real* x,
                const Real alpha, const Real beta) {
  impl->ts.solve(b, nrhs, x, alpha, beta);
}

void print_options (const Impl* impl, std::ostream& os) {
  impl::print_options(impl->o, os);
}

void set_level_schedule_only (Options& o) {
  o.min_lset_size = 0;
}

Options::Options () {
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
} // namespace hts
} // namespace Experimental
