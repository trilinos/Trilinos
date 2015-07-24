#ifndef GTS_TRISOLVE_IMPL_HPP
#define GTS_TRISOLVE_IMPL_HPP

#include <string.h>
#include <assert.h>
#include "shylu_gts.hpp"
#include "ShyLUGTS_config.h"

// Experimental P2P(-sparse) impl following Park et al 2014.
#define USE_P2P

namespace details {
namespace gts {
namespace impl {

//#define RESTRICT __restrict__
#define RESTRICT

struct NotTriangularException : public std::exception {
  virtual const char* what () const throw () {
    return "Not a triangular matrix.";
  }
};

void print_options(std::ostream& os);

struct Box {
  Int r0, c0, nr, nc;
  Box () : r0(0), c0(0), nr(0), nc(0) {}
  Box (const Int r0, const Int c0, const Int nr, const Int nc)
    : r0(r0), c0(c0), nr(nr), nc(nc) {}
};

struct ConstCrsMatrix {
  const Int m, n;
  const Int* const ir; // rowptr
  const Int* const jc; // col
  const Real* const d;

  ConstCrsMatrix (const Int nrow, const Int ncol, const Int* ir, const Int* jc,
                  const Real* d, bool deallocate = true)
    : m(nrow), n(ncol), ir(ir), jc(jc), d(d), deallocate_(deallocate)
  {}
  ~ConstCrsMatrix();
private:
  const bool deallocate_;
};

struct CrsMatrix {
  Int m, n;
  Int* const ir; // rowptr
  Int* const jc; // col
  Real* d;

  CrsMatrix (const Int nrow, const Int ncol, Int* ir, Int* jc, Real* d,
             bool deallocate = true)
    : m(nrow), n(ncol), ir(ir), jc(jc), d(d), deallocate_(deallocate)
  {}
  ~CrsMatrix();
  bool get_deallocate () const { return deallocate_; }
  void set_deallocate (const bool val) { deallocate_ = val; }
private:
  bool deallocate_;
};

struct CrsMatrixFiller {
  Int m, n;
  std::vector<Int> ir, jc;
  std::vector<Real> d;

  CrsMatrix* make_CrsMatrix ();
};

struct Partition {
  CrsMatrix* cm;
  std::vector<Int> A_idxs;
  Partition () : cm(0) {}
  ~Partition () { clear(); }
  void clear () { if (cm) { delete cm; cm = 0; } }
  void clear_d () { delete[] cm->d; cm->d = 0; }
  void alloc_d();
};

// Finds level sets.
class LevelSetter {
public:
  typedef std::vector<Int> VI;
  typedef std::vector<VI> LevelSets;

  void init(const ConstCrsMatrix& T, const Int size_thr, const bool is_lo,
            const Int ls_blk_sz);
  Int size () const { return lsets_.size(); }
  Int ls_blk_sz () const { return ls_blk_sz_; }
  const VI& lset(const size_t i) const;

private:
  LevelSets lsets_;
  bool is_lo_;
  Int ls_blk_sz_;
};

class Segmenter {
protected:
  std::vector<Int> nnz_, p_;
public:
  virtual ~Segmenter () {}
  Int nnz (const Int idx) const { return nnz_[idx]; }
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
  Int nnz (const Int idx) const { return nnz_[idx]; }
  const std::vector<Int>& p () const { return p_; }
private:
  void count_nnz_by_row(std::vector<Int>& rcnt);
  void init_nnz(const std::vector<Int>& rcnt);
  void segment();
};

struct InitInfo {
  unsigned short nthreads, min_blksz, max_nrhs;
  double min_dense_density;
};

// Represent both sparse and dense blocks. I'm choosing not to use runtime
// polymorphism. I'm also choosing not to use compile-time polymorphism because
// I want to avoid one level of pointer indirection. This has a particularly
// measurable effect on the symbolic phase's performance.
class SerialBlock {
private:
  Int r0_, c0_, nr_, nc_, roff_, coff_;
  Real* d_;
  Int* ir_, * jc_;

  void clear();
  void init_numeric_sparse(const CrsMatrix& A);
  void init_numeric_dense(const CrsMatrix& A);
  void n1Axpy_sparse(const Real* RESTRICT x, const Int ldx, const Int nrhs,
                     Real* RESTRICT y, const Int ldy) const;
  void n1Axpy_dense(const Real* RESTRICT x, const Int ldx, const Int nrhs,
                    Real* RESTRICT y, const Int ldy) const;

public:
  SerialBlock ()
    : r0_(0), c0_(0), nr_(0), nc_(0), roff_(0), coff_(0), d_(0), ir_(0), jc_(0)
  {}
  ~SerialBlock () { clear(); }
  Int r0 () const { return r0_; }
  // *Cropped* nr.
  Int nr () const { return nr_; }
  void init(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
            const InitInfo& in);
  void init_numeric(const CrsMatrix& A);
  bool inited () const { return d_; }
  void n1Axpy(const Real* RESTRICT x, const Int ldx, const Int nrhs,
              Real* RESTRICT y, const Int ldy) const;
  Int nnz() const;
};

// Threaded block matrix.
class TMatrix {
  Int nr_, coff_;
  std::vector<Int> ros_; // Block row offsets.
  std::vector<SerialBlock> bs_;
  bool is_parallel_;

public:
  TMatrix () : nr_(0), coff_(0), is_parallel_(false) {}
  ~TMatrix () { clear(); }
  bool empty() const;
  bool parallel () const { return is_parallel_; }
  Int nr () const { return nr_; }
  void init(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
            const InitInfo& in, const Int block_0_nnz_os = 0);
  void clear();
  void init_numeric(const CrsMatrix& A);
  // inside || {}
  void n1Axpy(const Real* x, const Int ldx, const Int nrhs, Real* y,
              const Int ldy, const int tid) const;
  // Raw block-level numbers; 0 padding is excluded.
  Int nblocks () const { return bs_.size(); }
  // Valid only for tid < nblocks().
  Int block_r0(const int tid) const;
  // *Cropped* nr.
  Int block_nr(const int tid) const;
};

// A triangle.
class Tri {
public:
  Tri () : n_(0), r0_(0) {}
  virtual ~Tri () {}
  Int n () const { return n_; }
  Int r0 () const { return r0_; }
  virtual void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                    const InitInfo& in) = 0;
  virtual void init_numeric(const CrsMatrix& T) = 0;
  virtual void solve(const Real* b, Real* x, const Int ldx, const Int nrhs)
    const = 0;
  virtual Int nnz() const { assert(0); return -1; }
protected:
  Int n_, r0_;
  static const bool is_lo_ = true;
};

class SpaceHolderTri : public Tri {
  virtual void init (const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                     const InitInfo& in) { n_ = n; r0_ = r0; }
  virtual void init_numeric (const CrsMatrix& T) {}
  virtual void solve (const Real* b, Real* x, const Int ldx, const Int nrhs)
    const {}
};

// Dense representation.
class DenseTri : public Tri {
  Int c0_;
  Real* d_; // Packed by column in dense tri format.
  void clear () { if (d_) delete[] d_; }
public:
  DenseTri () : c0_(0), d_(NULL) {}
  ~DenseTri () { clear(); }
  virtual void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                    const InitInfo& in);
  virtual void init_numeric(const CrsMatrix& T);
  virtual void solve(const Real* b, Real* x, const Int ldx, const Int nrhs)
    const;
  virtual Int nnz () const { return d_ ? ((n_ + 1)*n_) / 2 : 0; }
};

// Serial sparse representation.
class SerialSparseTri : public Tri {
  Int c0_;
  CrsMatrix* m_;
#ifdef HAVE_SHYLUGTS_MKL
  Real* b_; // x and b must be separate for MKL.
#endif
  void clear();
public:
  SerialSparseTri();
  ~SerialSparseTri () { clear(); }
  virtual void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                    const InitInfo& in);
  virtual void init_numeric(const CrsMatrix& T);
  virtual void solve(const Real* b, Real* x, const Int ldx, const Int nrhs)
    const;
  virtual Int nnz () const { return m_ ? m_->ir[m_->m] : 0; }
};

typedef Tri SerialTri;

namespace p2p {
typedef int Done;
}

// Recursive representation. A lower tri looks like this:
//     [Tri        ]
//     [TMatrix Tri].
// Siva R. suggested this to me. It's based on an idea possibly by Jonsson and
// Kagstrom, but perhaps it goes back further in time.
class RecursiveTri : public Tri {
  // Represent the blocks in a nonrecursive format.
  struct NonrecursiveData {
    std::vector<SerialTri*> t;
    std::vector<TMatrix*> s;
    std::vector<Int> os;

    // s k waits until t_barrier >= k.
    mutable Int t_barrier;
    mutable p2p::Done done_symbol;
    mutable std::vector<p2p::Done> s_done;
    // idx is pointers into ids.
    std::vector<Int> t_ids, t_idx;
    std::vector< std::vector<Int> > s_ids, s_idx;
  };

  Tri* t1_;
  TMatrix s_;
  Tri* t2_;
  Int nthreads_;
  RecursiveTri::NonrecursiveData nd_;

  void clear();
  void init_nd(NonrecursiveData& nd, Int os);
  void init_p2p();

  virtual void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                    const InitInfo& in);

public:
  RecursiveTri () : t1_(NULL), t2_(NULL) {}
  ~RecursiveTri () { clear(); }
  // Client should call this:
  void init_toplevel(const CrsMatrix& T, const Int r0, const Int c0,
                     const Int n, const InitInfo& in,
                     const Int mvp_block_nc = 0);
  virtual void init_numeric(const CrsMatrix& T);
  void p2p_reset() const;
  virtual void solve(const Real* b, Real* x, const Int ldx, const Int nrhs)
    const;
};

namespace p2p {
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

class LevelSetTri : public Tri {
  Int ls_blk_sz_;
  Int mvp_block_nc_; // Optional scatter block incorporated into this block.
  std::vector< std::vector<Int> > ps_; // Rows for a thread. For init_numeric.
  std::vector<Int> lsp_; // Level set i is lsp_[i] : lsp_[i+1]-1.
  struct Thread {
    CrsMatrix* m;
    std::vector<Int> p;   // Rows this thread owns.
    std::vector<Int> lsp; // Pointer to start of level set.
#ifdef USE_P2P
    std::vector<Int> p2p_depends;
    std::vector<Int> p2p_depends_p;
#endif
    Thread () : m(NULL) {}
    ~Thread () { if (m) { delete m; } }
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
  virtual void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
                    const InitInfo& in);
  virtual void init_numeric(const CrsMatrix& T);
  void reinit_numeric(const CrsMatrix& T);
  virtual void solve(const Real* b, Real* x, const Int ldx, const Int nrhs)
    const;

#ifdef USE_P2P
private:
  mutable std::vector<p2p::Done> p2p_done_;
  mutable p2p::Done p2p_done_value_;
  void find_task_responsible_for_variable(std::vector<p2p::Pair>& pairs);
  Int fill_graph(const std::vector<p2p::Pair>& pairs, std::vector<Int>& g,
                 std::vector<Int>& gp, std::vector<Int>& wrk);
  void prune_graph(const std::vector<Int>& gc, const std::vector<Int>& gp,
                   std::vector<Int>& g, std::vector<Int>& gsz,
                   std::vector<Int>& wrk, const Int max_gelen);
  void fill_dependencies(const std::vector<Int>& g, const std::vector<Int>& gp,
                         const std::vector<Int>& gsz);
public:
  void p2p_init();
  void p2p_solve(const Real* b, Real* x, const Int ldx, const Int nrhs) const;
#endif
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
  Real* from_outside(const Real* x, const Int nrhs) const;
  void to_outside(Real* x, const Int nrhs, const Real a, const Real b) const;
};

// Top-level solver.
class TriSolver {
  Int n_;
  bool is_lo_;
  unsigned short nthreads_;
  Partition p_[2];

  Permuter perm_;
  LevelSetTri lst_;
  RecursiveTri t_;

  void clear();

public:
  TriSolver () : n_(0), is_lo_(true), nthreads_(1) {}
  ~TriSolver () { clear(); }
  void init(const ConstCrsMatrix* T, const Int nthreads, const Int max_nrhs,
            const Int min_blksz, const bool save_for_reprocess, const Int* p,
            const Int* q, const Real* r);
  void init_numeric(const ConstCrsMatrix* T);
  bool is_lower_tri () const { return is_lo_; }
  // x and b can be the same pointers.
  void solve(const Real* b, const Int nrhs, Real* x, const Real alpha,
             const Real beta) const;
};

struct Options {
  Int min_block_size, serial_block_size, serial_block_size2,
    ls_blk_sz, max_level_set_bottlenecks, min_lset_size, printlvl;
  double min_dense_density;
  bool min_lset_size_scale_with_nthreads, profile;
  Options();
  void print(std::ostream& os) const;
};

class Profiler {
  Int nthreads_;
  Int lsets_actual_size_;
  LevelSetter::LevelSets lsets_;
  std::vector<Int> lsets_nnz_;

  struct Block {
    Int r0, c0, nr, nc;
    std::vector<Int> nnz_rows;
  };
  std::vector<Block> blocks_;

  struct Tri {
    Int r0, c0, n, nnz;
    bool is_sparse;
  };
  std::vector<Tri> tris_;

public:
  void clear();
  void set_nthreads (const int nt) { nthreads_ = nt; }
  void set_lsets_size (const int sz) { lsets_actual_size_ = sz; }
  void set_lsets_nnz (const std::vector<Int>& nnz) { lsets_nnz_ = nnz; }
  LevelSetter::LevelSets& get_lsets () { return lsets_; }
  void set_block(const Int r0, const Int c0, const Int nr, const Int nc,
                 const std::vector<Int>& nnz_rows);
  void set_tri(const Int r0, const Int c0, const Int n, const Int nnz,
               bool is_sparse);
  bool write(const std::string& filename) const;
  const LevelSetter::LevelSets& get_lsets () const { return lsets_; }
  Int get_lsets_n_used () const { return lsets_actual_size_; }
};

} // namespace impl

struct CrsMatrix : impl::ConstCrsMatrix {
  CrsMatrix (const Int nrow, const Int* ir, const Int* jc, const Real* d)
    : impl::ConstCrsMatrix(nrow, nrow, ir, jc, d, false)
  {}
};

struct Impl {
  impl::TriSolver ts;
};
} // namespace gts
} // namespace details

#endif
