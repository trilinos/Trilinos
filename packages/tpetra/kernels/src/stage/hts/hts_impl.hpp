#ifndef HTS_IMPL_HPP
#define HTS_IMPL_HPP

#include <string.h>
#include <assert.h>
#include "hts.hpp"

#ifdef __MIC__
# define USE_MM_MALLOC
#endif

#ifdef USE_MKL
# undef NO_BLAS
#endif

namespace Experimental {
namespace hts {
namespace impl {

#define RESTRICT

struct Box {
  Int r0, c0, nr, nc;
  Box () : r0(0), c0(0), nr(0), nc(0) {}
  Box (const Int r0, const Int c0, const Int nr, const Int nc)
    : r0(r0), c0(c0), nr(nr), nc(nc) {}
};

struct ConstCrsMatrix {
  enum Direction { forward = 0, transpose };

  const Int m, n;
  const Size* const ir; // row pointer
  const Int* const jc;  // col
  const Real* const d;
  Direction dir;

  ConstCrsMatrix (const Int nrow, const Int ncol, const Size* ir, const Int* jc,
                  const Real* d, Direction dir, bool deallocate)
    : m(nrow), n(ncol), ir(ir), jc(jc), d(d), dir(dir),
      deallocate_(deallocate)
  {}
  ~ConstCrsMatrix();
private:
  const bool deallocate_;
};

struct CrsMatrix {
  Int m, n;
  Size* const ir; // rowptr
  Int* const jc; // col
  Real* d;

  CrsMatrix (const Int nrow, const Int ncol, Size* ir, Int* jc, Real* d)
    : m(nrow), n(ncol), ir(ir), jc(jc), d(d)
  {}
  ~CrsMatrix();
};

struct Options {
  Int min_block_size, ls_blk_sz, lset_min_size, printlvl, pp_min_block_size,
    min_parallel_rows;
  Real min_dense_density, lset_max_bad_fraction;
  bool lset_min_size_scale_with_nthreads, profile;
  Options();
  void print(std::ostream& os) const;
};

struct Partition {
  CrsMatrix* cm;
  std::vector<Size> A_idxs;
  Partition () : cm(0) {}
  ~Partition () { clear(); }
  void clear();
  void clear_d();
  void alloc_d();
};

// Finds level sets.
class LevelSetter {
public:
  typedef std::vector<Int> VI;
  typedef std::vector<VI> LevelSets;

  void init(const ConstCrsMatrix& T, const Int size_thr, const bool is_lo,
            const Options& o);
  Int size () const { return lsets_.size(); }
  Int ls_blk_sz () const { return ls_blk_sz_; }
  const VI& lset(const size_t i) const;
  void reverse_variable_order(Int n);

private:
  LevelSets lsets_;
  bool is_lo_;
  Int ls_blk_sz_;
};

class Segmenter {
protected:
  std::vector<Size> nnz_;
  std::vector<Int> p_;
public:
  virtual ~Segmenter () {}
  Size nnz (const Int idx) const { return nnz_[idx]; }
  const std::vector<Int>& p () const { return p_; }
};

// Segment a CRS matrix into blocks for threaded MVP.
class CrsSegmenter : public Segmenter {
  const CrsMatrix& A_;
  const Int r0_, c0_, nr_, nc_, nthreads_, ls_blk_sz_, block_0_nnz_os_;
public:
  CrsSegmenter (const CrsMatrix& A, const Int r0, const Int c0, const Int nr,
                const Int nc, const Int nthreads, const Int ls_blk_sz = 1,
                const Int block_0_nnz_os = 0)
    : A_(A), r0_(r0), c0_(c0), nr_(nr), nc_(nc), nthreads_(nthreads),
      ls_blk_sz_(ls_blk_sz), block_0_nnz_os_(block_0_nnz_os)
  { segment(); }
  // nnz in segment idx.
  Size nnz (const Int idx) const { return nnz_[idx]; }
  const std::vector<Int>& p () const { return p_; }
private:
  void count_nnz_by_row(std::vector<Int>& rcnt);
  void count_nnz_by_row_loop(const Int i, std::vector<Int>& rcnt);
  void init_nnz(const std::vector<Int>& rcnt);
  void segment();
};

struct InitInfo {
  Int nthreads, min_blksz, max_nrhs, min_parallel_rows;
  Real min_dense_density;
};

namespace p2p {
typedef int Done;

struct Pair { Int lvl, tid; };

struct SortEntry {
  Int lvl, tid;
  bool operator< (const p2p::SortEntry& se) const {
    if (lvl < se.lvl) return true;
    if (lvl > se.lvl) return false;
    return tid < se.tid;
  }
};
} // namespace p2p

class TMatrix;

// Represent both sparse and dense blocks.
class SerialBlock {
  friend std::ostream& operator<<(std::ostream& os, const TMatrix& m);

  Int r0_, c0_, nr_, nc_, roff_, coff_;
  Size nnz_;
  Real* d_;
  Size* ir_;
  Int* jc_;
  bool deallocate_, is_dense_;

  void reinit_numeric_spars(const CrsMatrix& A);
  void reinit_numeric_dense(const CrsMatrix& A);
  void n1Axpy_spars(const Real* RESTRICT x, const Int ldx, const Int nrhs,
                    Real* RESTRICT y, const Int ldy) const;
  void n1Axpy_dense(const Real* RESTRICT x, const Int ldx, const Int nrhs,
                    Real* RESTRICT y, const Int ldy) const;

public:
  SerialBlock ()
    : r0_(0), c0_(0), nr_(0), nc_(0), roff_(0), coff_(0), nnz_(0), d_(0),
      ir_(0), jc_(0), deallocate_(true), is_dense_(true)
  {}
  ~SerialBlock () { clear(); }
  void clear();
  Int r0 () const { return r0_; }
  // *Cropped* nr.
  Int nr () const { return nr_; }
  void init(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
            const InitInfo& in);
  void init_metadata(const CrsMatrix& T, const Int r0, const Int c0,
                     const Int nr, const Int nc, const InitInfo& in);
  void init_memory(const InitInfo& in);
  void init_numeric(const CrsMatrix& T);
  void reinit_numeric(const CrsMatrix& A);
  bool inited () const { return d_; }
  void n1Axpy(const Real* RESTRICT x, const Int ldx, const Int nrhs,
              Real* RESTRICT y, const Int ldy) const;
  Size nnz () const { return nnz_; }
  bool is_sparse () const { return ir_; }

  // Experimental.
  Int Int_size() const;
  Int Real_size() const;
};

// Threaded block matrix.
class TMatrix {
  friend std::ostream& operator<<(std::ostream& os, const TMatrix& m);

  Int nr_, coff_, tid_os_;
  std::vector<Int> ros_; // Block row offsets.
  std::vector<SerialBlock> bs_;
  bool is_parallel_, is_empty_;

  void init_metadata_with_seg(const CrsMatrix& A, Int r0, Int c0, Int nr,
                              Int nc, const Int roff, const InitInfo& in,
                              const CrsSegmenter& seg);

public:
  TMatrix () : nr_(0), coff_(0), tid_os_(0), is_parallel_(false) {}
  ~TMatrix () { clear(); }
  bool empty () const { return is_empty_; }
  bool parallel () const { return is_parallel_; }
  Int nr () const { return nr_; }
  void init(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
            const InitInfo& in, const Int block_0_nnz_os = 0,
            const int tid_offset = 0);
  void init(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
            const InitInfo& in, const CrsSegmenter& seg);
  void clear();
  void init_metadata(const CrsMatrix& T, const Int r0, const Int c0,
                     const Int nr, const Int nc, const InitInfo& in,
                     // In calculating loads, assume thread 0 is already
                     // handling this many other nonzeros.
                     const Int block_0_nnz_os = 0,
                     // Offset tid by this number. It doesn't make sense to have
                     // block_0_nnz_os > 0 if this argument is > 0.
                     const int tid_offset = 0);
  // Version of the above in which a segmentation is imposed.
  void init_metadata(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                     const InitInfo& in, const CrsSegmenter& seg);
  void init_memory(const InitInfo& in);
  void init_numeric(const CrsMatrix& T, const int tid);
  void reinit_numeric(const CrsMatrix& A);
  // inside || {}
  void n1Axpy(const Real* x, const Int ldx, const Int nrhs, Real* y,
              const Int ldy, const int tid) const;
  // Raw block-level numbers; 0 padding is excluded.
  Int nblocks () const { return tid_os_ + bs_.size(); }
  // Valid only for tid < nblocks().
  Int block_r0(const int tid) const;
  // *Cropped* nr.
  Int block_nr(const int tid) const;
  const SerialBlock* block(const int tid) const;
  SerialBlock* block(const int tid);
};

// A triangle.
class Tri {
public:
  Tri () : n_(0), r0_(0) {}
  virtual ~Tri () {}
  Int n () const { return n_; }
  Int r0 () const { return r0_; }
protected:
  Int n_, r0_;
  static const bool is_lo_ = true;
};

class OnDiagTri : public Tri {
  Int c0_;
  Size nnz_;
  Real* d_; // Packed by column in dense tri format.
  CrsMatrix* m_;
#ifdef USE_MKL
  Real* b_; // x and b must be separate for MKL.
#endif
  bool dense_;

  void solve_spars(const Real* b, const Int ldb, Real* x, const Int ldx,
                   const Int nrhs) const;
  void solve_dense(const Real* b, const Int ldb, Real* x, const Int ldx,
                   const Int nrhs) const;
 
public:
  OnDiagTri();
  ~OnDiagTri () { clear(); }
  void clear();
  void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
            const InitInfo& in);
  // Initialization that takes up space as a block but has no data.
  void init(const Int r0, const Int c0, const Int n);
  Int nthreads() const;
  void init_metadata(const CrsMatrix& T, const Int r0, const Int c0,
                     const Int n, const InitInfo& in);
  void init_memory(const InitInfo& in);
  void init_numeric(const CrsMatrix& T);
  void reinit_numeric(const CrsMatrix& T);
  void solve(const Real* b, const Int ldb, Real* x, const Int ldx,
             const Int nrhs) const;
  Size nnz () const { return nnz_; }
  bool inited () const { return d_ || m_; }

private: // For inverse of on-diag triangle.
  struct Thread {
    Int r0, nr;
    Real* d;
    Thread() : r0(0), nr(0), d(0) {}
    ~Thread();
  };
  std::vector<Thread> t_;
  mutable int done_ctr_;
  void inv_init_metadata(const InitInfo& in);
  void inv_init_memory();
  void inv_reinit_numeric(const CrsMatrix& T);
  void solve_dense_inv(const Real* b, const Int ldb, Real* x, const Int ldx,
                       const Int nrhs) const;
};

// Recursive representation. A lower tri looks like this:
//     [Tri        ]
//     [TMatrix Tri].
class RecursiveTri : public Tri {
  Int nthreads_;

  struct NonrecursiveData {
    std::vector<OnDiagTri> t;
    std::vector<TMatrix> s;
    std::vector<Int> os;

    // s k waits until t_barrier >= k.
    mutable Int t_barrier;
    mutable p2p::Done done_symbol;
    mutable std::vector<p2p::Done> s_done;
    // idx is pointers into ids.
    std::vector<Int> t_ids, t_idx;
    std::vector< std::vector<Int> > s_ids, s_idx;
    // For the inverse on-diag tri.
    mutable std::vector<Int> inv_tri_done;
  } nd_;

  // For the inverse on-diag tri.
  mutable std::vector<Real> wrk_;

  void p2p_init();
  void ondiag_solve(const OnDiagTri& t, Real* x, const Int ldx, const Int nrhs,
                    const int tid, const Int stp, volatile Int* const t_barrier,
                    volatile Int* const inv_tri_done) const;

public:
  RecursiveTri () {}
  ~RecursiveTri () { clear(); }
  void clear();
  // Client should call this:
  void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
            const InitInfo& in, const Int mvp_block_nc = 0);
  void init_numeric(const CrsMatrix& T);
  void p2p_reset() const;
  void solve(const Real* b, Real* x, const Int ldx, const Int nrhs) const;
};

class LevelSetTri {
  Int n_, ls_blk_sz_;
  Int mvp_block_nc_; // Optional scatter block incorporated into this block.
  std::vector< std::vector<Int> > ps_; // Rows for a thread. For init_numeric.
  std::vector<Int> lsp_; // Level set i is lsp_[i] : lsp_[i+1]-1.
  // Use point-to-point synchronization and the tasking procedure described in
  // Park et al, Sparsifying Synchronizations for High-Performance Shared-Memory
  // Sparse Triangular Solve, ISC 2014.
  struct Thread {
    CrsMatrix* m;
    std::vector<Int> p;   // Rows this thread owns.
    std::vector<Int> lsp; // Pointer to start of level set.
    std::vector<Int> p2p_depends;
    std::vector<Int> p2p_depends_p;
    Thread () : m(NULL) {}
    ~Thread () { if (m) delete m; }
  };
  std::vector<Thread> t_;
  bool save_for_reprocess_;

public:
  LevelSetTri () : mvp_block_nc_(0) {}
  // Call this before init().
  void init_lsets(const LevelSetter& lstr, const bool save_for_reprocess);
  // If > 0, the scatter block is attached.
  void set_mvp_block_nc (const Int n) { mvp_block_nc_ = n; }
  void update_permutation(std::vector<Int>& lsis, const Partition& p);
  void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
            const InitInfo& in);
  void init_numeric(const CrsMatrix& T);
  void reinit_numeric(const CrsMatrix& T);
  void solve(const Real* b, Real* x, const Int ldx, const Int nrhs) const;

private:
  mutable std::vector<p2p::Done> p2p_done_;
  mutable p2p::Done p2p_done_value_;
  void find_task_responsible_for_variable(std::vector<p2p::Pair>& pairs);
  Int fill_graph(const std::vector<p2p::Pair>& pairs, std::vector<Int>& g,
                 std::vector<Size>& gp, std::vector<Int>& wrk);
  void prune_graph(const std::vector<Int>& gc, const std::vector<Size>& gp,
                   std::vector<Int>& g, std::vector<Int>& gsz,
                   std::vector<Int>& wrk, const Int max_gelen);
  void fill_dependencies(const std::vector<Int>& g, const std::vector<Size>& gp,
                         const std::vector<Int>& gsz);
public:
  void p2p_init();
  void p2p_solve(const Real* b, Real* x, const Int ldx, const Int nrhs) const;
};

// Change between internal and external ordering.
class Permuter {
  Int n_;
  Int* p_, * q_;
  Real* scale_;
  bool is_lo_;
  mutable Real* px_;
  std::vector<Int> part_;

public:
  Permuter () : n_(0), p_(NULL), q_(NULL), scale_(NULL), px_(NULL) {}
  ~Permuter () { clear(); }
  void clear();
  void init(const Int n, const bool is_lo, const std::vector<Int>& lsis,
            const std::vector<Int>& dpis, const Int nthreads,
            const Int max_nrhs, const Int* p, const Int* q, const Real* r);
  void reinit_numeric(const Real* r);
  Real* from_outside(const Real* x, const Int nrhs) const;
  void to_outside(Real* x, const Int nrhs, const Real a, const Real b) const;
};

// Top-level solver.
class TriSolver {
  Int n_;
  bool is_lo_;
  Int nthreads_;
  Partition p_[2];

  Permuter perm_;
  LevelSetTri lst_;
  RecursiveTri t_;

  void clear();

public:
  TriSolver () : n_(0), is_lo_(true), nthreads_(1) {}
  ~TriSolver () { clear(); }
  void init(const ConstCrsMatrix* T, Int nthreads, const Int max_nrhs,
            const bool save_for_reprocess, const Int* p, const Int* q,
            const Real* r, const Options& o);
  void reinit_numeric(const ConstCrsMatrix* T, const Real* r);
  bool is_lower_tri () const { return is_lo_; }
  // x and b can be the same pointers.
  void solve(const Real* b, const Int nrhs, Real* x, const Real alpha,
             const Real beta) const;
};

} // namespace impl

struct CrsMatrix : public impl::ConstCrsMatrix {
  CrsMatrix (const Int nrow, const Size* ir, const Int* jc, const Real* d,
             const impl::ConstCrsMatrix::Direction dir)
    : impl::ConstCrsMatrix(nrow, nrow, ir, jc, d, dir, false)
  {}
};

struct Impl {
  impl::TriSolver ts;
  impl::Options o;
};
} // namespace hts
} // namespace Experimental

#endif
