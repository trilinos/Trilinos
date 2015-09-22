#include <omp.h>
#include <sys/time.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <Teuchos_BLAS_wrappers.hpp>
#include "shylu_gts.hpp"
#include "gts_impl.hpp"
#ifdef HAVE_SHYLUGTS_MKL
# include <mkl.h>
#endif

// Debug output control.
#if 0
#define pr(m) do {                              \
    std::stringstream ss;                       \
    ss << m << std::endl;                       \
    std::cout << ss.str();                      \
  } while (0)
//#define TIME
#define prc(m) pr(#m << " | " << (m))
#else
#define pr(m)
#define prc(m)
#endif
#define pre(m) do {                             \
    std::stringstream ss;                       \
    ss << m << std::endl;                       \
    std::cerr << ss.str();                      \
  } while (0)

// Convenience try-catch wrapper for new.
#define ALLOC(alloc_line, catch_code, msg) do {         \
    try {                                               \
      alloc_line;                                       \
    } catch (...) {                                     \
      catch_code;                                       \
      throw Exception(msg ": failed to allocate.");     \
    }                                                   \
  } while (0)

namespace details {
template<typename T>
static void prarr (const std::string& name, const T* const v, const size_t n) {
  std::cerr << name << ": ";
  for (size_t i = 0; i < n; ++i) std::cerr << " " << v[i];
  std::cerr << "\n";
}

template<typename T>
static void prvec (const std::string& name, const std::vector<T>& v) {
  prarr(name, &v[0], v.size());
}

template<typename T>
inline void compress (std::vector<T>& v) { std::vector<T>(v).swap(v); }

#ifdef HAVE_SHYLUGTS_MKL
inline void gts_mkl_dcsrmm (
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
inline void gts_mkl_dcsrsm (
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
class Timer {
public:
  enum Op { setup = 0, tolower,
            lsetcntpat, lsetfind,
            perm,
            partition1, partition2, partitionsort,
            lsetinit, 
            mvpblockinit, dpblockinit,
            NSETUPTIMERS,
            slvpf, slvls, slvmvp, slvdp, slvpt,
            NTIMERS };
  static void init () {
#ifdef TIME
    for (int i = 0; i < NTIMERS; ++i) et_[i] = 0;
#endif
  }
  static void start (const Op op) {
#ifdef TIME
    gettimeofday(&t_start_[op], 0);
#endif
  }
  static void stop (const Op op) {
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
  static void print_setup () {
#ifdef TIME
    double tot = 0;
    for (int i = 0; i < NSETUPTIMERS; ++i) tot += et_[i];
    tpr(setup); tpr(tolower);
    tpr(lsetcntpat); tpr(lsetfind);
    tpr(perm);
    tpr(partition1); tpr(partition2); tpr(partitionsort);
    tpr(lsetinit);
    tpr(mvpblockinit); tpr(dpblockinit);
    printf("%20s %10.3e %10.1f\n", "total", tot, 100.0);
#endif
  }
  static void print_solve () {
#ifdef TIME
    double tot = 0;
    for (int i = NSETUPTIMERS; i < NTIMERS; ++i) tot += et_[i];
    tpr(slvpf); tpr(slvls); tpr(slvmvp); tpr(slvdp); tpr(slvpt);
    printf("%20s %10.3e %10.1f\n", "total", tot, 100.0);
#endif
  }
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

namespace gts {
namespace impl {
static void sync_options (Options& o) {
  o.serial_block_size = o.min_block_size;
  o.serial_block_size2 = o.serial_block_size*o.serial_block_size;
}

static void set_options (const gts::Options& os, Options& od) {
  od.min_block_size = os.min_block_size;
  od.min_dense_density = os.min_dense_density;
  od.ls_blk_sz = os.levelset_block_size;
  od.max_level_set_bottlenecks = os.max_level_set_bottlenecks;
  od.min_lset_size = os.min_lset_size;
  od.min_lset_size_scale_with_nthreads = os.min_lset_size_scale_with_nthreads;
  od.profile = os.profile;
  od.printlvl = os.print_level;
  sync_options(od);
}

Options::Options () { set_options(gts::Options(), *this); }

void Options::print (std::ostream& os) const {
  os << " min_block_size " << min_block_size
     << " serial_block_size " << serial_block_size
     << " max_level_set_bottlenecks " << max_level_set_bottlenecks
     << " min_dense_density " << min_dense_density
     << " ls_blk_sz " << ls_blk_sz
     << " min_lset_size " << min_lset_size
     << " min_lset_size_scale_with_nthreads "
     << min_lset_size_scale_with_nthreads
     << " profile " << profile;
}

static void print_compiletime_options(std::ostream& os) {
#ifdef HAVE_SHYLUGTS_MKL
  os << " HAVE_SHYLUGTS_MKL";
#endif
#ifdef USE_P2P
  os << " USE_P2P";
#endif
}

const Options* options (const gts::Options* o=0) {
  static Options opts;
  if (o) set_options(*o, opts);
  return &opts;
}

void print_options (std::ostream& os) {
  print_compiletime_options(os);
  os << std::endl;
  options()->print(os);
  os << std::endl;
}

void Profiler::clear () {
  lsets_.clear(); lsets_nnz_.clear(); blocks_.clear(); tris_.clear();
}

void Profiler::set_block (const Int r0, const Int c0, const Int nr,
                          const Int nc, const std::vector<Int>& nnz_rows) {
  Block b;
  b.r0 = r0; b.c0 = c0; b.nr = nr; b.nc = nc; b.nnz_rows = nnz_rows;
  blocks_.push_back(b);
}

void Profiler::set_tri (const Int r0, const Int c0, const Int n, const Int nnz,
                        const bool is_sparse) {
  Tri t;
  t.r0 = r0; t.c0 = c0; t.n = n; t.nnz = nnz; t.is_sparse = is_sparse;
  tris_.push_back(t);
}

bool Profiler::write (const std::string& filename) const {
  FILE* fid = fopen(filename.c_str(), "w");
  if ( ! fid) return false;

  { std::stringstream ss;
    print_compiletime_options(ss);
    fprintf(fid, "p> compiletime: %s\n", ss.str().c_str()); }
  { std::stringstream ss;
    options()->print(ss);
    fprintf(fid, "p> runtime: %s\n", ss.str().c_str()); }
  fprintf(fid, "p> nthreads: %d\n", nthreads_);

  fprintf(fid, "s> blocks\n");
  fprintf(fid, "# r0 c0 nr nc, then nnz/(block row)\n");
  for (size_t i = 0; i < blocks_.size(); ++i) {
    const Block& b = blocks_[i];
    fprintf(fid, "p> block: %d %d %d %d\n", b.r0, b.c0, b.nr, b.nc);
    fprintf(fid, "p> nnz:");
    for (size_t j = 0; j < b.nnz_rows.size(); ++j)
      fprintf(fid, " %d", b.nnz_rows[j]);
    fprintf(fid, "\n");
  }

  fprintf(fid, "s> triangles\n");
  fprintf(fid, "# r0 c0 n nnz is_sparse\n");
  for (size_t i = 0; i < tris_.size(); ++i) {
    const Tri& t = tris_[i];
    fprintf(fid, "p> tri: %d %d %d %d %d\n", t.r0, t.c0, t.n, t.nnz,
            t.is_sparse);
  }

  fprintf(fid, "s> level sets\n");
  fprintf(fid, "p> used: %d\n", lsets_actual_size_);
  fprintf(fid, "p> nnz:");
  for (size_t i = 0; i < lsets_nnz_.size(); ++i)
    fprintf(fid, " %d", lsets_nnz_[i]);
  fprintf(fid, "\n");
  fprintf(fid, "p> sizes:\n");
  for (size_t i = 0; i < lsets_.size(); ++i)
    fprintf(fid, "%ld %ld\n", i, lsets_[i].size());

  fclose(fid);
  return true;
}

Profiler* profiler () {
  static Profiler p;
  if (options()->profile) return &p; else return 0;
}

ConstCrsMatrix::~ConstCrsMatrix () {
  if ( ! deallocate_) return;
  if (ir) delete[] ir;
  if (jc) delete[] jc;
  if (d) delete[] d;
}

CrsMatrix::~CrsMatrix () {
  if ( ! deallocate_) return;
  if (ir) delete[] ir;
  if (jc) delete[] jc;
  if (d) delete[] d;
}

ConstCrsMatrix* move_to_ConstCrsMatrix (CrsMatrix& A) {
  ConstCrsMatrix* Ac;
  ALLOC(Ac = new ConstCrsMatrix(A.m, A.n, A.ir, A.jc, A.d, A.get_deallocate()),
        ;, "ConstCrsMatrix");
  A.set_deallocate(false);
  return Ac;
}

template<typename T> static T* vec2arr (const std::vector<T>& v) {
  if (v.empty()) return NULL;
  T* a;
  ALLOC(a = new T[v.size()], ;, "vec2arr");
  memcpy(a, &v[0], sizeof(T)*v.size());
  return a;
}

CrsMatrix* CrsMatrixFiller::make_CrsMatrix () {
  CrsMatrix* cm = 0;
  Int* ira = 0, *jca = 0;
  Real* da = 0;
  try {
    ira = vec2arr(ir);
    jca = vec2arr(jc);
    da = vec2arr(d);
    cm = new CrsMatrix(m, n, ira, jca, da);
  } catch (...) {
    if (ira) delete[] ira;
    if (jca) delete[] jca;
    if (da) delete[] da;
    if (cm) delete cm;
    throw Exception("make_ConstCrsMatrix: failed to allocate.");
  }
  return cm;
}

void Partition::alloc_d () {
  assert( ! cm->d);
  ALLOC(cm->d = new Real[cm->ir[cm->m]], ;, "Partition::alloc_d");
}

namespace {
struct SparseData {
  Int* ir, * jc;
  Real* d;
  SparseData (const Int m, const Int nnz, const bool touch = false) {
    ir = jc = 0;
    d = 0;
    try {
      ir = new Int[m+1];
      if (nnz > 0) {
        jc = new Int[nnz];
        d = new Real[nnz];
      }
      ir[0] = 0;
      ir[m] = nnz;
    } catch (...) {
      free();
      throw Exception("SparseData failed to allocate.");
    }
    if (touch) {
      for (Int i = 0; i <= m; ++i) ir[i] = 0;
      for (Int i = 0; i < nnz; ++i) { jc[i] = 0; d[i] = 0; }
    }
  }
  void free () {
    if (ir) delete[] ir;
    if (jc) delete[] jc;
    if (d) delete[] d;    
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

inline void set_num_threads (const int nthreads) {
  omp_set_num_threads(nthreads);
#ifdef HAVE_SHYLUGTS_MKL
  // We never use MKL threading.
  mkl_set_num_threads(1);
#endif  
}

// Return i such that ir[i] is the first value >= c. If no such i exists,
// return n.
inline Int find_first (const Int* ir, const Int n, const Int c) {
  return c == 0 ? 0 : std::lower_bound(ir, ir+n, c) - ir;
}

// Crop the submatrix A(b) such that A(cb) has no 0 border.
Int crop_matrix (const CrsMatrix& T, const Box& b, Box& cb) {
  cb.r0 = -1;
  Int r1 = -1;
  cb.c0 = b.c0 + b.nc;
  Int c1 = b.c0;
  Int nnz = 0;
  for (Int r = b.r0, lrow = 0; r < b.r0 + b.nr; ++r, ++lrow) {
    const Int irr = T.ir[r], irrp1 = T.ir[r+1];
    Int cnt = 0;
    for (Int k = irr + find_first(T.jc + irr, irrp1 - irr, b.c0);
         k < irrp1; ++k){
      const Int c = T.jc[k];
      assert(c >= b.c0);
      if (c < b.c0 + b.nc) {
        cb.c0 = std::min(cb.c0, c);
        c1 = std::max(c1, c);
        ++cnt;
        ++nnz;
      } else break;
    }
    if (cnt) {
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

void find_row_level_sets_Lcrs (const ConstCrsMatrix& L, const Int sns,
                               Int size_thr, LevelSetter::LevelSets& lsets) {
  assert(L.m % sns == 0);
  Timer::start(Timer::lsetfind);

  // Find the level sets. Andrey Prokopenko suggested this algorithm. It means
  // we don't have to form the CCS matrix for L. Thanks Andrey! Erik found the
  // algorithm as eq. 18 in Y. Saad's 1989 SIAM J Sci Stat Comput paper.
  //   Implement it for a blocked matrix, where the block size is sns >= 1.
  const Int Lm_sr = L.m / sns;
  std::vector<Int> w(Lm_sr, -1);
  Int max_level = -1;
  for (Int sr = 0; sr < Lm_sr; ++sr) {
    Int level = -1;
    for (Int i = 0; i < sns; ++i) {
      const Int r = sns*sr + i;
      for (Int j = L.ir[r]; j < L.ir[r+1]; ++j) {
        const Int sc = L.jc[j] / sns;
        level = std::max(level, w[sc]);
      }
    }
    ++level;
    w[sr] = level;
    max_level = std::max(max_level, level);
  }
 
  // Count level set sizes.
  std::vector<Int> n(max_level+1);
  for (size_t i = 0; i < w.size(); ++i) ++n[w[i]];

  // Decide how many level sets to keep.
  const Int max_lset_bn = options()->max_level_set_bottlenecks;
  Int bn = 0, streak = 0;
  Int lsmi = -1; // level sets max index
  do {
    ++lsmi;
    if (n[lsmi] < size_thr) {
      ++bn;
      ++streak;
      if (bn == max_lset_bn || n[lsmi] == 0) break;
    } else streak = 0;
    if (lsmi == max_level) break;
  } while (true);
  lsmi -= streak;

  // Allocate lsets.
  lsets.resize(lsmi+1);
  for (Int i = 0; i <= lsmi; ++i)
    lsets[i].reserve(sns * n[i]);
  // Load.
  for (Int i = 0; i < Lm_sr; ++i) {
    const Int level = w[i];
    if (level <= lsmi) {
      for (Int j = 0; j < sns; ++j) lsets[level].push_back(i * sns + j);
    }
  }
  Timer::stop(Timer::lsetfind);
}

class Stack {
  std::vector<Int> d_;
  Int i_;
public:
  Stack (const Int capacity) : d_(capacity), i_(0) {}
  void push_back (const Int e) { d_[i_++] = e; }
  void clear () { i_ = 0; }
  Int size () const { return i_; }
  const Int* data () const { return &d_[0]; }
};

class ColCntPattern {
protected:
  const ConstCrsMatrix& T_;
  VI col_nnz_, touched_;
  const Int sns_;
public:
  ColCntPattern (const ConstCrsMatrix& T, const Int sns) : T_(T), sns_(sns) {
    if (sns_ == 1) {
      col_nnz_.resize(T_.n, 0);
      for (Int k = 0; k < T_.ir[T_.m]; ++k)
        ++col_nnz_[T_.jc[k]];
    } else {
      // col_nnz_ is for super columns.
      const Int n = T_.n / sns_;
      col_nnz_.resize(n, 0);
      touched_.resize(n, -1); // Has this supernode been counted yet?
      for (Int r = 0; r < T_.m; ++r) {
        const Int sr = r / sns_;
        for (Int j = T_.ir[r]; j < T_.ir[r+1]; ++j) {
          const Int sc = T_.jc[j] / sns_;
          if (touched_[sc] != sr) {
            // If this supernode has not yet been touched, count it.
            ++col_nnz_[sc];
            // Now touch it so we don't count it again. touched_ is constructed
            // and used in such a way that we don't have to do O(n) work at any
            // point to maintain it other than in the initial resize.
            touched_[sc] = sr;
          }
        }
      }
      // Reset for use in nix.
      touched_.assign(touched_.size(), n);
    }
  }
  Int n () const { return col_nnz_.size(); }
  Int cnt (const Int idx) const { return col_nnz_[idx]; }
  // Input and output are supernode indices.
  void nix (const Int r, Stack& nnz1idx) {
    if (sns_ == 1) {
      for (Int j = T_.ir[r]; j < T_.ir[r+1]; ++j) {
        const Int c = T_.jc[j];
        if (--col_nnz_[c] == 1) nnz1idx.push_back(c);
      }
    } else {
      const Int sr = r;
      for (Int r = sr*sns_; r < (sr + 1)*sns_; ++r)
        for (Int j = T_.ir[r]; j < T_.ir[r+1]; ++j) {
          const Int c = T_.jc[j], sc = c / sns_;
          if (touched_[sc] != sr) {
            if (--col_nnz_[sc] == 1) nnz1idx.push_back(sc);
            touched_[sc] = sr;
          }
        }
    }
  }
};

void find_col_level_sets_Lcrs (const ConstCrsMatrix& L, const Int sns,
                               Int size_thr, LevelSetter::LevelSets& lsets) {
  size_thr = std::max(1, size_thr);
  assert(L.m % sns == 0);
  Timer::start(Timer::lsetcntpat);
  ColCntPattern S(L, sns);
  Timer::stop(Timer::lsetcntpat);

  Timer::start(Timer::lsetfind);
  const Int max_lset_bn = options()->max_level_set_bottlenecks;
  Int bn = 0, streak = 0;
  Stack nnz1idxs(S.n());
  bool first = true;
  for (;;) {
    lsets.push_back(std::vector<Int>());
    VI& idxs = lsets.back();
    if (first) {
      // Find root (super)nodes.
      for (Int i = 0; i < S.n(); ++i)
        if (S.cnt(i) == 1) idxs.push_back(i);
      first = false;
    } else {
      idxs.resize(nnz1idxs.size());
      memcpy(&idxs[0], nnz1idxs.data(), idxs.size()*sizeof(Int));
      nnz1idxs.clear();
    }
    if ((Int) idxs.size() < size_thr) { // <<--- need to handle size_thr == 0
      ++bn;
      ++streak;
      if (bn == max_lset_bn || idxs.size() == 0) break;
    } else streak = 0;
    for (size_t i = 0; i < idxs.size(); ++i) S.nix(idxs[i], nnz1idxs);
    if (sns > 1) {
      // Supernode -> nodes.
      VI sidxs(idxs);
      idxs.clear();
      idxs.reserve(sns * sidxs.size());
      for (size_t i = 0; i < sidxs.size(); ++i)
        for (Int k = 0; k < sns; ++k) idxs.push_back(sidxs[i] * sns + k);
    }
  }
  for (Int j = 0; j < streak; ++j) lsets.pop_back();
  Timer::stop(Timer::lsetfind);
}

inline void find_level_sets (
  const ConstCrsMatrix& T, const Int sns, const Int size_thr, const bool is_lo,
  LevelSetter::LevelSets& lsets)
{
  if (is_lo)
    find_row_level_sets_Lcrs(T, sns, size_thr, lsets);
  else
    find_col_level_sets_Lcrs(T, sns, size_thr, lsets);
}
} // namespace

void LevelSetter::init (const ConstCrsMatrix& T, const Int size_thr,
                        const bool is_lo, const Int ls_blk_sz) {
  lsets_.clear();
  is_lo_ = is_lo;
  // Guard against an invalid setting.
  ls_blk_sz_ = T.m % ls_blk_sz == 0 ? ls_blk_sz : 1;
  find_level_sets(T, ls_blk_sz_, size_thr, is_lo_, lsets_);
  if (profiler()) {
    profiler()->set_lsets_size(lsets_.size());
    find_level_sets(T, ls_blk_sz_, 1, is_lo_, profiler()->get_lsets());
  }
}

const VI& LevelSetter::lset (const size_t i) const {
  return is_lo_ ? lsets_[i] : lsets_[lsets_.size() - i - 1];
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
  for (Int r = 0; r < A.m; ++r) {
    bool diag_fnd = false;
    for (Int j = A.ir[r]; j < A.ir[r+1]; ++j) {
      const Int c = A.jc[j];
      if (c != r) {
        if ( ! tri_determined) {
          is_lower = c < r;
          tri_determined = true; 
        }
        if ((is_lower && c > r) || ( ! is_lower && c < r)) {
          is_tri = false;
          return Shape(is_lower, is_tri);
        }
      } else diag_fnd = true;
    }
    if ( ! diag_fnd) has_full_diag = false;
  }
  // If ! tri_determined, then T must be a diag matrix. Can treat as lower,
  // which is is_lower's default value.
  return Shape(is_lower, is_tri, has_full_diag);
}

int are_same_ (const ConstCrsMatrix& A, const ConstCrsMatrix& B) {
#define same(expr, ret) if ( ! (expr)) return ret;
  same(A.m == B.m, 1);
  same(A.n == B.n, 2);
  same(A.ir[A.m] == B.ir[B.m], 3);
  const Int nnz = A.ir[A.m];
  for (Int k = 0; k <= A.m; ++k)
    same(B.ir[k] == A.ir[k], 4);
  for (Int k = 0; k < nnz; ++k) {
    same(B.jc[k] == A.jc[k], 5);
    same(B.d[k] == A.d[k], 6);
  }
  return 0;
#undef same
}

bool are_same (const ConstCrsMatrix& A, const ConstCrsMatrix& B) {
  const int ret = are_same_(A, B);
  if (ret) pr("are_same failed with " << ret);
  return ret == 0;
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
  Timer::start(Timer::partition1);
  // Count nnz.
  Int nnz = 0;
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
#pragma omp parallel for schedule(static,1)
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
  Timer::stop(Timer::partition1);
}

// Extract B = A(p,s), s = [q p], given that A(p,~s) is empty.
void get_matrix_p_qp_with_covers_all (
  const ConstCrsMatrix& A, const PermVec& pv, const PermVec& qv, Partition& p)
{
  Timer::start(Timer::partition2);
  Int nnz = 0;
  for (size_t i = 0, lim = pv.size(); i < lim; ++i) {
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
#pragma omp parallel for schedule(static,1)
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
  Timer::stop(Timer::partition2);
}

template<typename T> struct SortEntry {
  Int i, A_idx; T d;
  SortEntry () {}
  bool operator< (const SortEntry& se) const { return i < se.i; }
};

void sort (Partition& p) {
  Timer::start(Timer::partitionsort);
  const int nthreads = omp_get_max_threads();
  typedef std::vector< SortEntry<Real> > SortableVector;
  std::vector<SortableVector> sess(nthreads);
  CrsMatrix& A = *p.cm;
# pragma omp parallel for schedule(static,1)
  for (Int r = 0; r < A.m; ++r) {
    const int tid = omp_get_thread_num();
    const Int irr = A.ir[r], irrp1 = A.ir[r+1], nc = irrp1 - irr;
    SortableVector& ses = sess[tid];
    ses.resize(nc);
    for (Int j = 0; j < nc; ++j) {
      const Int Aj = irr + j;
      ses[j].i = A.jc[Aj];
      ses[j].d = A.d[Aj];
      ses[j].A_idx = p.A_idxs[Aj];
    }
    std::sort(ses.begin(), ses.end());
    for (Int j = 0; j < nc; ++j) {
      const Int Aj = irr + j;
      A.jc[Aj] = ses[j].i;
      A.d[Aj] = ses[j].d;
      p.A_idxs[Aj] = ses[j].A_idx;
    }
  }
  Timer::stop(Timer::partitionsort);
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
#pragma omp parallel for schedule(static)
    for (Int i = 0; i < ilim; ++i)
      p.cm->d[i] = A.d[p.A_idxs[i]];
  } else {
    const Int nnz = A.ir[A.m];
#pragma omp parallel for schedule(static)
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

void CrsSegmenter::count_nnz_by_row (std::vector<Int>& rcnt) {
  rcnt.resize(nr_);
# pragma omp parallel for schedule(guided)
  for (Int i = 0; i < nr_; ++i) {
    rcnt[i] = 0;
    const Int
      iri = A_.ir[r0_+i],
      irip1 = A_.ir[r0_+i+1];
    for (Int j = iri + find_first(A_.jc + iri, irip1 - iri, c0_);
         j < irip1; ++j) {
      if (A_.jc[j] >= c0_ + nc_) break;
      ++rcnt[i];
    }
  }
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
  if ( ! tid_empty)
    for (size_t i = 0; i < cs_rcnt.size(); ++i)
      cs_rcnt[i] += block_0_nnz_os_;

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

inline void TMatrix::clear () { bs_.clear(); }

inline bool TMatrix::empty () const {
  for (size_t i = 0; i < bs_.size(); ++i)
    if (bs_[i].nnz() != 0) return false;
  return true;
}

inline void
TMatrix::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
               const InitInfo& in, const Int block_0_nnz_os) {
  clear();
  nr_ = nr;
  
  Box b;
  const Int nnz = crop_matrix(A, Box(r0, c0, nr, nc), b);
  const Int roff = b.r0 - r0;
  coff_ = b.c0 - c0;
  r0 = b.r0; c0 = b.c0; nr = b.nr; nc = b.nc;

  // p2p sync'ing favors fewer threads. So don't run a bunch of threads on small
  // blocks.
  const Int nt = nnz / options()->serial_block_size2 + 1;
  const Int nseg = std::max(1, std::min((Int) in.nthreads, nt));
  CrsSegmenter seg(A, r0, c0, nr, nc, nseg, 1, block_0_nnz_os);

  // Serial if block is too small.
  is_parallel_ = nseg > 1;
  if (nnz > 0) {
    if ( ! is_parallel_) {
      bs_.push_back(SerialBlock());
      bs_.back().init(A, r0, c0, nr, nc, in);
      ros_.push_back(roff);
    } else {
      const std::vector<Int>& p = seg.p();
      const Int n = (Int) p.size() - 1;
      ros_.resize(n, 0);
      bs_.resize(n, SerialBlock());
#     pragma omp parallel
      { const int tid = omp_get_thread_num();
        if (tid >= 0 && tid < (int) bs_.size()) {
          ros_[tid] = p[tid] - r0 + roff;
          const Int nri = p[tid+1] - p[tid];
          bs_[tid].init(A, p[tid], c0, nri, nc, in);
        } }
    }
  }
}

inline void TMatrix::init_numeric (const CrsMatrix& A) {
  if ( ! is_parallel_) {
    if (bs_.empty()) return;
    bs_[0].init_numeric(A);
  } else {
#   pragma omp parallel
    { const int tid = omp_get_thread_num();
      if (tid < (int) bs_.size())
        bs_[tid].init_numeric(A);
    }
  }
}

inline Int TMatrix::block_r0 (const int tid) const { return bs_[tid].r0(); }

inline Int TMatrix::block_nr (const int tid) const { return bs_[tid].nr(); }

inline void SerialBlock::clear () {
  if (ir_) delete[] ir_;
  if (jc_) delete[] jc_;
  if (d_) delete[] d_;
  ir_ = jc_ = 0;
  roff_ = coff_ = 0;
  d_ = 0;
}

inline void
SerialBlock::init (const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                   const InitInfo& in) {
  clear();

  Box b;
  const Int nnz = crop_matrix(A, Box(r0, c0, nr, nc), b);
  roff_ = b.r0 - r0;
  r0 = b.r0; nr = b.nr;
  coff_ = b.c0 - c0;
  c0 = b.c0; nc = b.nc;

  r0_ = r0; c0_ = c0; nr_ = nr; nc_ = nc;
  if (nr_ == 0 || nc_ == 0) return;
  assert(nnz > 0);
  assert(nnz <= (1.0*nr_)*nc_);
  const bool is_dense = nnz >= (in.min_dense_density*nr_)*nc_;
  if (is_dense) {
    ALLOC(d_ = new Real[nr_*nc_], ;, "DenseSerialBlock::init");
    memset(d_, 0, nr_*nc_*sizeof(*d_));
  } else {
    try {
      ir_ = new Int[nr_+1];
      if (nnz == 0) {
        ir_[nr_] = 0;
        return;
      }
      jc_ = new Int[nnz]; d_ = new Real[nnz];
    } catch (...) {
      clear();
      throw Exception("SparseSerialBlock::init: failed to allocate.");
    }
  }
  init_numeric(A);
}

inline void SerialBlock::init_numeric (const CrsMatrix& A) {
  if (ir_) init_numeric_sparse(A);
  else init_numeric_dense(A);
}

inline void SerialBlock::init_numeric_sparse (const CrsMatrix& A) {
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
#ifndef HAVE_SHYLUGTS_MKL
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

void SerialBlock::init_numeric_dense (const CrsMatrix& A) {
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
      d_[nc_*lrow + lcol] = A.d[j];
    }
  }
}

inline Int SerialBlock::nnz () const {
  return d_ ? (ir_ ? ir_[nr_] : nr_*nc_) : 0;
}

inline Int ntri (const int n) { return (n*(n + 1))/2; }

void DenseTri::init (const CrsMatrix& T, const Int r0, const Int c0,
                     const Int n, const InitInfo& in) {
  clear();
  n_ = n; r0_ = r0; c0_ = c0;
  if ( ! n_) return;

  const Int N = ntri(n_);
  ALLOC(d_ = new Real[N], ;, "DenseTri::init");
  memset(d_, 0, N*sizeof(*d_));

  init_numeric(T);
  if (profiler()) profiler()->set_tri(r0_, c0_, n_, N, false);
}

void DenseTri::init_numeric (const CrsMatrix& T) {
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
      // Compressed dense triangle.
      d_[(lrow*(lrow + 1))/2 + lcol] = T.d[k];
      ++nnz;
    }
  }
}

SerialSparseTri::SerialSparseTri ()
  : c0_(0), m_(NULL)
#ifdef HAVE_SHYLUGTS_MKL
  , b_(NULL)
#endif
{}

void SerialSparseTri::clear () {
  if (m_) { delete m_; m_ = 0; }
#ifdef HAVE_SHYLUGTS_MKL
  if (b_) { delete[] b_; b_ = 0; }
#endif
}

void SerialSparseTri::init (const CrsMatrix& T, const Int r0, const Int c0,
                            const Int n, const InitInfo& in) {
  clear();
  n_ = n; r0_ = r0; c0_ = c0;
  assert(n_ >= 1);

  CrsMatrixFiller f;
  f.m = f.n = n;
  f.ir.resize(n + 1);
  f.ir[0] = 0;
  f.ir.back() = 0;
  
  for (Int grow = r0; grow < r0 + n; ++grow) {
    const Int
      lrow = grow - r0,
      irg = T.ir[grow],
      irgp1 = T.ir[grow+1];
    f.ir[lrow+1] = f.ir[lrow];
    for (Int k = irg + find_first(T.jc + irg, irgp1 - irg, c0);
         k < irgp1; ++k) {
      const Int lcol = T.jc[k] - c0;
      if (lcol >= n) break;
      f.jc.push_back(lcol);
      f.d.push_back(T.d[k]);
      ++f.ir[lrow+1];
    }
  }

  m_ = f.make_CrsMatrix();

#ifdef HAVE_SHYLUGTS_MKL
  ALLOC(b_ = new Real[in.max_nrhs * n_], clear(), "SerialSparseTri::init");
#endif

  if (profiler()) profiler()->set_tri(r0, c0, n, m_->ir[m_->m], true);
}

void SerialSparseTri::init_numeric (const CrsMatrix& T) {
  Int* const jc = m_->jc;
  Real* const d = m_->d;
  for (Int grow = r0_, i = 0; grow < r0_ + n_; ++grow) {
    const Int
      irg = T.ir[grow],
      irgp1 = T.ir[grow+1];
    for (Int k = irg + find_first(T.jc + irg, irgp1 - irg, c0_);
         k < irgp1; ++k) {
      const Int lcol = T.jc[k] - c0_;
      if (lcol >= n_) break;
      jc[i] = lcol;
      d[i] = T.d[k];
      ++i;
    }
  }
}

void RecursiveTri::clear () {
  if (t1_) delete t1_;
  s_.clear();
  if (t2_) delete t2_;
  t1_ = t2_ = 0;
}

inline Int split (const Int n, const Int nthreads) {
  const Int s = n/2;
  if (n <= nthreads) return s;
  const Int smod = s % nthreads;
  if ( ! smod) return s;
  return s + (nthreads - smod);
}

inline Int
count_nnz_lotri (const CrsMatrix& T, const Int r0, const Int c0,
                 const Int n) {
  Int cnt = 0;
  for (Int r = r0, lrow = 0; r < r0 + n; ++r, ++lrow) {
    const Int irr = T.ir[r], irrp1 = T.ir[r+1];
    for (Int k = irr + find_first(T.jc + irr, irrp1 - irr, c0);
         k < irrp1; ++k)
      if (T.jc[k] <= c0 + lrow) ++cnt;
  }
  return cnt;
}

inline SerialTri*
make_serial_Tri (const CrsMatrix& T, const Int r0, const Int c0,
                 const Int n, const InitInfo& in) {
  SerialTri* t;
  const int nnz = count_nnz_lotri(T, r0, c0, n);
  if (nnz >= in.min_dense_density*ntri(n))
    ALLOC(t = new DenseTri(), ;, "make_serial_Tri");
  else
    ALLOC(t = new SerialSparseTri(), ;, "make_serial_Tri");
  t->init(T, r0, c0, n, in);
  return t;
}

inline Tri*
make_Tri (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
          const InitInfo& in) {
  if (n == 0) return NULL;
  Tri* t = 0;
  // Use 2 rather than 1 min_blksz here for the following reason. If 1 is used,
  // then the tri is subdivided. But then each block is serial. So might as well
  // make the whole thing one serial block. The one tradeoff is sparse/dense
  // fineness.
  if (n > 2*in.min_blksz) {
    ALLOC(t = new RecursiveTri(), ;, "RecursiveTri::make_Tri");
    t->init(T, r0, c0, n, in);
  } else t = make_serial_Tri(T, r0, c0, n, in);
  return t;
}

void RecursiveTri::
init_toplevel (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
               const InitInfo& in, const Int mvp_block_nc) {
  nthreads_ = in.nthreads;
  if (mvp_block_nc) {
    // Special case. The t1 triangle is empty. This handles the MVP block in the
    // gts 3-block decomposition. Input sizes are inconsistent with obvious
    // usage b/c T is not a square matrix; rather, it looks like this: [B T],
    // where B is the MVP block and T is the true top-level recursive tri.
    clear();
    n_ = n;
    t1_ = new SpaceHolderTri();
    t1_->init(T, 0, 0, mvp_block_nc, in);
    s_.init(T, r0, c0, n, mvp_block_nc, in);
    t2_ = make_Tri(T, r0, c0 + mvp_block_nc, n, in);
  } else
    init(T, r0, c0, n, in);
  init_nd(nd_, 0);
  // Placeholder one past what is actually used for calculations. This lets us
  // avoid an 'if' on whether to read.
  nd_.os.push_back(0);
  init_p2p();
}

void RecursiveTri::init (const CrsMatrix& T, const Int r0, const Int c0,
                         const Int n, const InitInfo& in) {
  clear();
  n_ = n;
  const Int n1 = split(n, in.nthreads), n2 = n - n1;
  assert((n <= in.nthreads || n1 % in.nthreads == 0) && n1 >= 0 && n1 <= n);
  t1_ = make_Tri(T, r0, c0, n1, in);
  s_.init(T, r0 + n1, c0, n2, n1, in);
  t2_ = make_Tri(T, r0 + n1, c0 + n1, n2, in);
}

void RecursiveTri::init_nd (NonrecursiveData& nd, Int os) {
  assert(t1_ || t2_);
  if ( ! t1_ || ! t2_) {
    //todo This situation arises when LS claims all but one row. Instead of
    // special casing this here, I should special case adding that one row into
    // the LS block.
    nd.t.push_back(static_cast<SerialTri*>(t1_ ? t1_ : t2_));
    return;
  }
  { RecursiveTri* rt = dynamic_cast<RecursiveTri*>(t1_);
    if (rt) rt->init_nd(nd, os);
    else nd.t.push_back(static_cast<SerialTri*>(t1_)); }
  nd.s.push_back(s_.empty() ? 0 : &s_);
  nd.os.push_back(os);
  if (t2_) {
    RecursiveTri* rt = dynamic_cast<RecursiveTri*>(t2_);
    if (rt) rt->init_nd(nd, os + t1_->n());
    else nd.t.push_back(static_cast<SerialTri*>(t2_));
  }
}

void RecursiveTri::init_numeric (const CrsMatrix& T) {
  int n = (Int) nd_.t.size();
#pragma omp parallel for schedule(dynamic)
  for (Int i = 0; i < n; ++i)
    nd_.t[i]->init_numeric(T);
  n = nd_.s.size();
  for (Int i = 0; i < n; ++i)
    if (nd_.s[i]) nd_.s[i]->init_numeric(T);
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

  if (profiler()) {
    std::vector<Int> nnz(t_.size(), 0);
    for (size_t i = 0; i < t_.size(); ++i) {
      const CrsMatrix* const A = t_[i].m;
      if (A) nnz[i] = A->ir[A->m];
    }
    profiler()->set_lsets_nnz(nnz);
  }
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
// particular, t.m->d and t.m->ir are unchanged. A cool thing in my
// implementation is that the diagonal element remains the last element in the
// row, which favors the streaming of data through memory.
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

#ifdef USE_P2P
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
#endif // USE_P2P

void LevelSetTri::init_numeric (const CrsMatrix& T) {
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    if (t_[tid].m) delete t_[tid].m;
    t_[tid].m = get_matrix_p(T, ps_[tid]); }
}

void LevelSetTri::reinit_numeric (const CrsMatrix& T) {
# pragma omp parallel
  { const int tid = omp_get_thread_num();
    Thread& t = t_[tid];
    assert(t.m);
    const std::vector<Int>& p = ps_[tid];
    for (size_t i = 0; i < p.size(); ++i) {
      const Int r = p[i], nc = T.ir[r+1] - T.ir[r];
      assert(nc == t.m->ir[i+1] - t.m->ir[i]);
      memcpy(t.m->d + t.m->ir[i], T.d + T.ir[r], nc*sizeof(*t.m->d));
    } }
}

inline void Permuter::clear () {
  if (q_ && q_ != p_) { delete[] q_; q_ = 0; }
  if (p_) { delete[] p_; p_ = 0; }
  if (scale_) { delete[] scale_; scale_ = 0; }
  if (px_) { delete[] px_; px_ = 0; }
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
    p_ = new Int[n_];
    px_ = new Real[n_*max_nrhs];
    q_ = p || q ? new Int[n_] : p_;
    if (scale) {
      scale_ = new Real[n_];
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

// This is the top-implementation-level preprocessing method. Every step of this
// method needs to be redone, as development of this prototype has been solely
// focused on a fast solve phase, with the (fairly safe) assumption that the
// preprocessing phase can be made reasonably fast, perhaps 10-20 times the time
// of a single solve. However, the preprocessing phase is the more complicated
// of the two phases, so this will take some work. The timers separate
// preprocessing into logical segments so we can see how this work is going.
void TriSolver::init (const ConstCrsMatrix* T, const Int nthreads,
                      const Int max_nrhs, const Int min_blksz,
                      const bool save_for_reprocess, const Int* p, const Int* q,
                      const Real* scale) {
  bool delete_T = false;
  try {
    Timer::init();
    Timer::start(Timer::setup);
    set_num_threads(nthreads);

    // Basic size parameters.
    n_ = T->m;
    nthreads_ = nthreads;
    InitInfo in;
    in.nthreads = nthreads_;
    in.min_blksz = min_blksz;
    in.min_dense_density = options()->min_dense_density;
    in.max_nrhs = max_nrhs;

    // Determine shape.
    Shape shape = determine_shape(*T);
    if ( ! shape.is_triangular) throw NotTriangularException();
    if ( ! shape.has_full_diag) throw NotFullDiagonal();
    is_lo_ = shape.is_lower;
    Timer::stop(Timer::setup); Timer::start(Timer::tolower);
    if ( ! is_lo_) {
      T = permute_to_other_tri(*T);
      delete_T = true;
    }
    Timer::stop(Timer::tolower);

    // Find level sets.
    LevelSetter lstr;
    const Int lstr_threshold = options()->min_lset_size * 
      (options()->min_lset_size_scale_with_nthreads ? nthreads_ : 1);
    lstr.init(*T, lstr_threshold, is_lo_, options()->ls_blk_sz);

    // Separate into three blocks: level set, scatter, and data parallel:
    //     [(1)        0
    //      (scatter) (2)],
    // where if is_lo_, (1) = level sets and (2) = data parallel; if ! is_lo_,
    // then the opposite.
    std::vector<Int> lsis, dpis;
    get_idxs(n_, lstr, lsis, dpis);
    Timer::stop(Timer::perm);
    if (options()->printlvl > 0)
      pr("n " << n_ << " |lsis| " << lsis.size() << " |dpis| " << dpis.size());
    if (is_lo_) {
      // 1. Level-scheduling block.
      { PermVec lsis_pv(T->m, lsis);
        get_matrix_pp_with_covers_all(*T, lsis_pv, p_[0]);
        sort(p_[0]); }
      Timer::start(Timer::lsetinit);
      lst_.init_lsets(lstr, save_for_reprocess);
      lst_.init(*p_[0].cm, 0, 0, p_[0].cm->m, in);
      lst_.update_permutation(lsis, p_[0]);
#ifdef USE_P2P
      lst_.p2p_init();
#endif
      Timer::stop(Timer::lsetinit);
      { PermVec dpis_pv(T->m, dpis), lsis_pv(T->m, lsis);
        get_matrix_p_qp_with_covers_all(*T, dpis_pv, lsis_pv, p_[1]);
        sort(p_[1]); }
      if (p_[1].cm->m > 0) {
        // 2. No MVP block. It's in the data-parallel block.
        // 3. Data-parallel block (+ MVP block).
        Timer::start(Timer::dpblockinit);
        const Int mvp_nc = p_[1].cm->n - p_[1].cm->m;
        t_.init_toplevel(*p_[1].cm, 0, 0, p_[1].cm->m, in, mvp_nc);
        Timer::stop(Timer::dpblockinit);
      }
    } else {
      PermVec dpis_pv(T->m, dpis);
      get_matrix_pp_with_covers_all(*T, dpis_pv, p_[1]);
      sort(p_[1]);
      if (p_[1].cm->m > 0) {
        // 3. Data-parallel block.
        Timer::start(Timer::dpblockinit);
        t_.init_toplevel(*p_[1].cm, 0, 0, p_[1].cm->m, in);
        Timer::stop(Timer::dpblockinit);
        // 2. No MVP block. It's in the level scheduling block.
      }
      // 1. Level-scheduling block (+ MVP block).
      { PermVec lsis_pv(T->m, lsis);
        get_matrix_p_qp_with_covers_all(*T, lsis_pv, dpis_pv, p_[0]);
        sort(p_[0]); }
      Timer::start(Timer::lsetinit);
      lst_.init_lsets(lstr, save_for_reprocess);
      lst_.set_mvp_block_nc(p_[1].cm->m);
      lst_.init(*p_[0].cm, 0, 0, p_[0].cm->m, in);
      lst_.update_permutation(lsis, p_[0]);
#ifdef USE_P2P
      lst_.p2p_init();
#endif
      Timer::stop(Timer::lsetinit);
    }
    if (delete_T) delete T;
    if (save_for_reprocess) {
      // For 32-bit Int, 64-bit Real, save 2/3 of memory during solves at little
      // extra init_numeric cost.
      for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].clear_d();
    } else {
      for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].clear();
    }
    Timer::start(Timer::perm);
    perm_.init(n_, is_lo_, lsis, dpis, nthreads_, max_nrhs, p, q, scale);
    Timer::stop(Timer::perm);
    if (options()->printlvl > 0) Timer::print_setup();
  } catch (...) {
    if (delete_T) delete T;
    throw;
  }
}

// Reinitialize numbers, but keep the same structures.
void TriSolver::init_numeric (const ConstCrsMatrix* T) {
  for (Int i = 0; i < 2; ++i) if (p_[i].cm) p_[i].alloc_d();
  Timer::init(); Timer::start(Timer::partition1);
  repartition_into_2_blocks(p_, *T, is_lo_);
  Timer::stop(Timer::partition1); Timer::start(Timer::lsetinit);
  lst_.reinit_numeric(*p_[0].cm);
  Timer::stop(Timer::lsetinit);
  if (p_[1].cm->m > 0) {
    Timer::start(Timer::dpblockinit);
    //todo Tighten up some of these init_numeric impls. Might want to do
    // reinit_numeric like for lst.
    t_.init_numeric(*p_[1].cm);
    Timer::stop(Timer::dpblockinit);
  }
  if (options()->printlvl > 0) Timer::print_setup();
}

//> Solve code.

inline void DenseTri::solve (const Real* b, Real* x, const Int ldx,
                             const Int nrhs) const {
  for (Int irhs = 0; ; ) {
    for (Int j = 0, k = 0; j < n_; ++j) {
      Real a = b[j];
      const Int ilim = j;
      for (Int i = 0; i < ilim; ++i, ++k)
        a -= x[i]*d_[k];
      x[j] = a / d_[k++];
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldx;
  }
}

inline void
SerialSparseTri::solve (const Real* b, Real* x, const Int ldx, const Int nrhs)
  const {
  // No OpenMP: block is too small, so serialize.
#ifdef HAVE_SHYLUGTS_MKL
  memcpy(b_, x, n_*sizeof(*b_));
  for (int k = 1; k < nrhs; ++k)
    memcpy(b_ + n_*k, x + ldx*k, n_*sizeof(*b_));
  gts_mkl_dcsrsm(true, m_->m, m_->d, m_->ir, m_->jc, b_, m_->m, x, ldx, nrhs);
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
      x[r] = a / d[rp_rp1 - 1];
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldx;
  }
#endif
}

inline void
SerialBlock::n1Axpy (const Real* RESTRICT x, const Int ldx, const Int nrhs,
                     Real* RESTRICT y, const Int ldy) const {
  if ( ! d_) return;
  if (ir_) n1Axpy_sparse(x + coff_, ldx, nrhs, y + roff_, ldy);
  else n1Axpy_dense(x + coff_, ldx, nrhs, y + roff_, ldy);
}

inline void SerialBlock::
n1Axpy_sparse (const Real* RESTRICT x, const Int ldx, const Int nrhs,
               Real* RESTRICT y, const Int ldy) const {
  assert(ir_);
  if (ir_[nr_] == 0) return;
#ifdef HAVE_SHYLUGTS_MKL
  gts_mkl_dcsrmm(false, nr_, nc_, d_, ir_, jc_, x, ldx, y, ldy, nrhs);
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

#if defined(HAVE_SHYLUGTS_BLAS) && ! defined(HAVE_SHYLUGTS_MKL)
typedef Int blas_int;
template<typename T> static void gemm(
  const char transa, const char transb, const blas_int m, const blas_int nrhs,
  const blas_int n, const T alpha, const T* a, const blas_int lda, const T* b,
  const blas_int ldb, const T beta, T* c, const blas_int ldc);

template<> inline void gemm<double> (
  const char transa, const char transb, blas_int m, const blas_int nrhs,
  blas_int n, const double alpha, const double* a, const blas_int lda,
  const double* b, const blas_int ldb, const double beta, double* c,
  const blas_int ldc)
{
  // Calling Teuchos::BLAS<Int, Real>::GEMM is measurably slower, so use
  // Teuchos_BLAS_wrappers directly. I'm copying the following macro from
  // Teuchos_BLAS.cpp. Perhaps one should be exposed in Teuchos_BLAS.hpp or
  // Teuchos_BLAS_wrappers.hpp.
# ifdef INTEL_CXML
  unsigned int cone = 1;
# define CHAR_MACRO(char_var) &char_var, cone
# else
# define CHAR_MACRO(char_var) &char_var
# endif
  DGEMM_F77(CHAR_MACRO(transa), CHAR_MACRO(transb), &m, &nrhs, &n, &alpha, a,
            &lda, b, &ldb, &beta, c, &ldc);
}
#endif // HAVE_SHYLUGTS_BLAS

inline void SerialBlock::
n1Axpy_dense (const Real* RESTRICT x, const Int ldx, const Int nrhs,
              Real* RESTRICT y, const Int ldy) const {
  assert(d_);
#if defined(HAVE_SHYLUGTS_MKL)
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nr_, nrhs, nc_, -1, d_,
              nc_, x, ldx, 1, y, ldy);
#elif defined(HAVE_SHYLUGTS_BLAS)
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
  if ((size_t) tid >= bs_.size()) return;
  bs_[tid].n1Axpy(x + coff_, ldx, nrhs, y + ros_[tid], ldy);
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
        x[r] = a / d[j++];
      }
#     pragma omp barrier
    }
    if (++irhs == nrhs) break;
    x += ldx;
    b += ldx;
  }
}

#ifdef USE_P2P
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
        x[r] = a / d[j++];
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
#endif

#define rb_p2p_sub2ind(sid, tid) (sn*(tid) + (sid))
#define rb_p2p_ind2tid(e) ((e) / (sn))
#define rb_p2p_ind2sid(e) ((e) % (sn))

// Form the lists of p2p dependencies. A dependency in this algorithm is of two
// types. One is the usual dependency: a variable has to be solved for before
// the next. The second imposes an ordering to assure there is no write race
// condition when computing a row's dot product in pieces.
void RecursiveTri::init_p2p () {
  // For L, we rely on the fact that t_ doesn't count SpaceHolderTri's #rows.
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
      const Tri* const t = nd_.t[ti];
      const Int r0 = t->r0(), nr = t->n();
      Int k = ti == 1 ? 0 : nd_.t_idx[ti-1];
      for (Int r = r0, rlim = r0 + nr; r < rlim; ++r) {
        assert(r < n_);
        const Int wr = w[r];
        if (wr >= 0 && w_cnted[wr] != w_cnted_symbol) {
          const Int sid = rb_p2p_ind2sid(wr), tid = rb_p2p_ind2tid(wr);
          if (tid != 0) {
            nd_.t_ids.push_back(rb_p2p_sub2ind(sid, tid));
            ++k;
          }
          w_cnted[wr] = w_cnted_symbol;
        }
      }
      nd_.t_idx[ti] = k;
      ++w_cnted_symbol;
    }

    if (ti+1 < tn) { // MVP block.
      const TMatrix* const s = nd_.s[ti];
      if ( ! s) continue;
      for (Int bi = 0, bn = s->nblocks(); bi < bn; ++bi) {
        const Int r0 = s->block_r0(bi), nr = s->block_nr(bi);
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
#pragma omp parallel
  { const int tid = omp_get_thread_num();
    compress(nd_.s_ids[tid]);
    compress(nd_.s_idx[tid]); }
  compress(nd_.s_ids);
}

inline void RecursiveTri::p2p_reset () const {
  nd_.t_barrier = -1;
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
  volatile p2p::Done* const s_done = &nd_.s_done[0];
  const Int* const t_ids = &nd_.t_ids[0];
  const Int* const t_idx = &nd_.t_idx[0];
  const Int* const s_ids = nd_.s_ids[tid].empty() ? 0 : &nd_.s_ids[tid][0];
  const Int* const s_idx = nd_.s_idx[tid].empty() ? 0 : &nd_.s_idx[tid][0];
  const p2p::Done done_symbol = nd_.done_symbol;
# pragma omp master
  { nd_.t[0]->solve(x, x, ldx, nrhs);
    *t_barrier = 0; }
  if ( ! nd_.os.empty()) {
    os += nd_.t[0]->n();
    x_osi = x + nd_.os[0];
    x_os = x + os;
  } else return;
  for (Int i = 0; i < sn; ++i) {
    const TMatrix* const s = nd_.s[i];
    if (s) {
      if (s->parallel()) {
        if (tid) // Master doesn't need to wait on the tri.
          while (*t_barrier < i) ;
        rbwait(s_done, s_ids, s_idx, i, done_symbol);
        s->n1Axpy(x_osi, ldx, nrhs, x_os, ldx, tid);
#       pragma omp flush
        s_done[rb_p2p_sub2ind(i, tid)] = done_symbol;
#       pragma omp master
        {
          Real* const xos = x + os;
          rbwait(s_done, t_ids, t_idx, i, done_symbol);
          nd_.t[i+1]->solve(xos, xos, ldx, nrhs);
          *t_barrier = i+1;
#         pragma omp flush
        }
      } else {
#       pragma omp master
        {
          rbwait(s_done, s_ids, s_idx, i, done_symbol);
          s->n1Axpy(x_osi, ldx, nrhs, x_os, ldx, 0);
#         pragma omp flush
          s_done[rb_p2p_sub2ind(i, 0)] = done_symbol;
          Real* const xos = x + os;
          rbwait(s_done, t_ids, t_idx, i, done_symbol);
          nd_.t[i+1]->solve(xos, xos, ldx, nrhs);
          *t_barrier = i+1;
#         pragma omp flush
        }
      }
    } else {
#     pragma omp master
      {
        Real* const xos = x + os;
        rbwait(s_done, t_ids, t_idx, i, done_symbol);
        nd_.t[i+1]->solve(xos, xos, ldx, nrhs);
        *t_barrier = i+1;
#       pragma omp flush
      }
    }
    os += nd_.t[i+1]->n();
    x_osi = x + nd_.os[i+1];
    x_os = x + os;
  }
}

inline void
TriSolver::solve (const Real* b, const Int nrhs, Real* x, const Real alpha,
                  const Real beta) const {
  set_num_threads(nthreads_);
  t_.p2p_reset();
# pragma omp parallel
  {
    Real* px = perm_.from_outside(b, nrhs);
    if (is_lo_) {
#     pragma omp barrier
#ifdef USE_P2P
      lst_.p2p_solve(px, px, n_, nrhs);
#else
      lst_.solve(px, px, n_, nrhs);
#endif
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
#ifdef USE_P2P
      lst_.p2p_solve(px, px, n_, nrhs);
#else
      lst_.solve(px, px, n_, nrhs);
#endif
      // No barrier needed because of lst_.solve.
      perm_.to_outside(x, nrhs, alpha, beta);
    }
  }
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

void delete_CrsMatrix (CrsMatrix* T) { delete T; }

Impl* preprocess (const CrsMatrix* T, const Int max_nrhs, const Int nthreads,
                  const bool save_for_reprocess, const Int* p, const Int* q,
                  const Real* r)
{
  Impl* impl;
  ALLOC(impl = new Impl(), ;, "preprocess");
  try {
    impl->ts.init(T, nthreads, max_nrhs, impl::options()->serial_block_size,
                  save_for_reprocess, p, q, r);
  } catch (Exception& e) {
    delete impl;
    throw;
  }
  if (impl::profiler()) impl::profiler()->set_nthreads(nthreads);
  return impl;
}

void reprocess_numeric (Impl* impl, const CrsMatrix* T) {
  impl->ts.init_numeric(T);
}

bool is_lower_tri (const Impl* impl) { return impl->ts.is_lower_tri(); }

void delete_Impl (Impl* impl) {
  delete impl;
  if (impl::profiler()) impl::profiler()->clear();
  if (impl::options()->printlvl > 0) Timer::print_solve();
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

void print_options (std::ostream& os) {
  impl::print_options(os);
}

Options::Options ()
  : min_block_size(32), min_dense_density(0.75), levelset_block_size(1),
    max_level_set_bottlenecks(1<<30), min_lset_size_scale_with_nthreads(false),
    profile(false), print_level(0)
{
  min_lset_size = min_lset_size_scale_with_nthreads ?
#ifdef USE_P2P
    1 : 10
#else
    5 : 40
#endif
    ;
}

void set_options (const Options& o) { impl::options(&o); }

void write_profile (const std::string& filename) {
  if ( ! impl::profiler()) throw Exception("Not profiled.");
  if ( ! impl::profiler()->write(filename))
    throw Exception("Could not write " + filename);
}

void write_levelsets (const std::string& filename) {
  if ( ! impl::profiler()) throw Exception("Not profiled.");
  const impl::LevelSetter::LevelSets& lsets = impl::profiler()->get_lsets();
  FILE* fid = fopen(filename.c_str(), "w");
  if ( ! fid) throw Exception("Could not write " + filename);
  for (Int lsi = 0; lsi < impl::profiler()->get_lsets_n_used(); ++lsi)
    for (size_t i = 0; i < lsets[lsi].size(); ++i)
      fprintf(fid, "%d %d\n", static_cast<int>(lsi), lsets[lsi][i] + 1);
  fclose(fid);
}

void write_matrixmarket (const CrsMatrix* T, const std::string& filename) {
  FILE* fid = fopen(filename.c_str(), "w");
  if ( ! fid) throw Exception("Could not write " + filename);
  fprintf(fid, "%%%%MatrixMarket matrix coordinate real general\n"
          "%11d %11d %11d\n", T->m, T->n, T->ir[T->m]);
  for (Int r = 0; r < T->m; ++r)
    for (Int j = T->ir[r]; j < T->ir[r+1]; ++j)
      fprintf(fid, "%d %d %1.15e\n", r+1, T->jc[j] + 1, T->d[j]);
  fclose(fid);
}

void write_binary (const CrsMatrix* T, const std::string& filename) {
  FILE* fid = fopen(filename.c_str(), "w");
  if ( ! fid) throw Exception("Could not write " + filename);
  fwrite(&T->m, sizeof(Int), 1, fid);
  fwrite(&T->n, sizeof(Int), 1, fid);
  fwrite(T->ir, sizeof(Int), T->m + 1, fid);
  fwrite(T->jc, sizeof(Int), T->ir[T->m], fid);
  fwrite(T->d, sizeof(Real), T->ir[T->m], fid);
  fclose(fid);
}
} // namespace gts
} // namespace details
