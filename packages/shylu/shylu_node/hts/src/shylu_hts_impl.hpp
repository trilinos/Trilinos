// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef INCLUDE_SHYLU_HTS_IMPL_HPP
#define INCLUDE_SHYLU_HTS_IMPL_HPP

#include <cassert>
#include <cstring>
#include <limits>
#include "shylu_hts.hpp"

namespace Experimental {
namespace htsimpl {

struct Direction {
  enum Enum { forward = 0, transpose };
  static Enum opposite (const Enum dir)
  { return static_cast<Enum>((static_cast<int>(dir) + 1) % 2); }
};

// RAII array that does not touch all entries on construction like
// std::vector. It does not support realloc-like operations.
template <typename T> class Array {
  T* p_;
  std::size_t n_, cap_;
public:
  Array () { init(); }
  Array(std::size_t n);
  Array(std::size_t n, const T& init);
  ~Array () { clear(); }
  // Initialize the object with the assumption that all variables are uninit'ed
  // prior to calling.
  void init();
  void clear();
  // optclear means optionally clear. The function has the semantics of
  // clearing, but it may not actually release the memory.
  void optclear_and_resize(std::size_t n);
  // _ft indicates first touch.
  void optclear_and_resize_ft(std::size_t n);
  void optclear_and_resize(std::size_t n, const T& i);
  void optclear_and_reserve(std::size_t n);
  void optclear_and_reserve_ft(std::size_t n);
  T& operator[] (std::size_t i) { return p_[i]; }
  const T& operator[] (std::size_t i) const { return p_[i]; }
  T& back () { return p_[n_-1]; }
  const T& back () const { return p_[n_-1]; }
  std::size_t size () const { return n_; }
  bool empty () const { return size() == 0; }
  T* data () const { return p_; }
  // This does not realloc; reserve must provide the necessary memory. It does
  // not throw, either. It asserts.
  void unsafe_push_back(const T& e);
  T* begin () { return p_; }
  T* end () { return p_ + n_; }
};

template<typename Int, typename Size, typename Sclr> struct Impl {
  // User API type.
  typedef HTS<Int, Size, Sclr> ihts;
  typedef typename hts::ScalarTraits<Sclr>::Real Real;

  struct Box {
    Int r0, c0, nr, nc;
    Box () : r0(0), c0(0), nr(0), nc(0) {}
    Box (const Int ir0, const Int ic0, const Int inr, const Int inc)
      : r0(ir0), c0(ic0), nr(inr), nc(inc) {}
  };

  struct ConstCrsMatrix {
    // Callers's specification.
    const Int m, n;
    const Size* const ir; // row pointer
    const Int* const jc;  // col
    const Sclr* const d;
    const Direction::Enum dir;
    const bool conj;
    mutable typename HTS<Int, Size, Sclr>::Deallocator* deallocator;
    // Either inferred by determine_shape from the user's matrix, or these data
    // are propagated.
    mutable bool unitdiag, is_lo;

    ConstCrsMatrix (const Int inrow, const Int incol, const Size* iir,
                    const Int* ijc, const Sclr* id, Direction::Enum idir,
                    const bool iconj, bool ideallocate,
                    const bool unit_diag, const bool is_lower)
      : m(inrow), n(incol), ir(iir), jc(ijc), d(id), dir(idir), conj(iconj),
        deallocator(0), unitdiag(unit_diag), is_lo(is_lower),
        deallocate_(ideallocate)
    {}
    ~ConstCrsMatrix();

    void deallocate() const;

  private:
    const bool deallocate_;
  };

  struct CrsMatrix {
    Int m, n;
    Size* const ir; // rowptr
    Int* const jc; // col
    Sclr* d;

    CrsMatrix (const Int nrow, const Int ncol, Size* iir, Int* ijc, Sclr* id)
      : m(nrow), n(ncol), ir(iir), jc(ijc), d(id)
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
    Size* A_idxs;
    Size nnz;

    Partition () : cm(0), A_idxs(0) {}
    ~Partition () { clear(); }
    void clear();

    void alloc_d();
    void clear_d();

    void alloc_A_idxs(const Size nnz);
    // Handle implicit unit diag. A_idxs[i] can't point to anything if the i'th
    // element is a diagonal entry and the matrix has an implicit unit diag.
    void A_invalidate (const Size& i)
    { if (A_idxs) A_idxs[i] = std::numeric_limits<Size>::max(); }
    bool A_valid (const Size& i) const
    { return A_idxs[i] != std::numeric_limits<Size>::max(); }
  };

  // Finds level sets.
  class LevelSetter {
  public:
    typedef Array<Array<Int> > LevelSets;

    void init(const ConstCrsMatrix& T, const Int size_thr, const bool is_lo,
              const Options& o);
    Int size () const { return lsets_.size(); }
    Int ls_blk_sz () const { return ls_blk_sz_; }
    const Array<Int>& lset(const size_t i) const;
    void reverse_variable_order(Int n);

  private:
    LevelSets lsets_;
    bool is_lo_;
    Int ls_blk_sz_;
  };

  class Segmenter {
  protected:
    Array<Size> nnz_;
    Array<Int> p_;
  public:
    virtual ~Segmenter () {}
    Size nnz (const Int idx) const { return nnz_[idx]; }
    const Array<Int>& get_p () const { return p_; }
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
    Size nnz (const Int idx) const { return this->nnz_[idx]; }
  private:
    void count_nnz_by_row(Array<Int>& rcnt);
    void count_nnz_by_row_loop(const Int i, Array<Int>& rcnt);
    void init_nnz(const Array<Int>& rcnt);
    void segment();
  };

  struct InitInfo {
    Int nthreads, min_blksz, max_nrhs, min_parallel_rows;
    Real min_dense_density;
  };

  typedef int p2p_Done;

  struct p2p_Pair {
    Int lvl, tid;
    p2p_Pair& operator= (const Int& v) {
      lvl = tid = v;
      return *this;
    }
  };

  struct p2p_SortEntry {
    Int lvl, tid;
    bool operator< (const p2p_SortEntry& se) const {
      if (lvl < se.lvl) return true;
      if (lvl > se.lvl) return false;
      return tid < se.tid;
    }
  };

  class TMatrix;

  // Represent both sparse and dense blocks.
  class SerialBlock {
    Int r0_, c0_, nr_, nc_, roff_, coff_;
    Size nnz_;
    Sclr* d_;
    Size* ir_;
    Int* jc_;
    bool deallocate_, is_dense_;

    void reinit_numeric_spars(const CrsMatrix& A);
    void reinit_numeric_dense(const CrsMatrix& A);

  public:
    SerialBlock ()
      : r0_(0), c0_(0), nr_(0), nc_(0), roff_(0), coff_(0), nnz_(0), d_(0),
        ir_(0), jc_(0), deallocate_(true), is_dense_(true)
    {}
    ~SerialBlock () { clear(); }
    void clear();
    // Cropped r0.
    Int get_r0 () const { return r0_; }
    // Cropped nr.
    Int get_nr () const { return nr_; }
    void init(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
              const InitInfo& in);
    void init_metadata(const CrsMatrix& T, Int r0, Int c0, Int nr, Int nc,
                       const InitInfo& in);
    void init_memory(const InitInfo& in);
    void init_numeric(const CrsMatrix& T);
    void reinit_numeric(const CrsMatrix& A);
    bool inited () const { return d_; }
    void n1Axpy(const Sclr* x, const Int ldx, const Int nrhs,
                Sclr* y, const Int ldy) const;
    Size get_nnz () const { return nnz_; }
    bool is_sparse () const { return ir_; }

    // For analysis; not algorithms.
    // Cropped c0 and nc.
    Int get_c0 () const { return c0_; }
    Int get_nc () const { return nc_; }
  };

  // Threaded block matrix.
  class TMatrix {
    Int nr_, coff_, tid_os_;
    Array<Int> ros_; // Block row offsets.
    Array<SerialBlock> bs_;
    bool is_parallel_, is_empty_;

    void init_metadata_with_seg(const CrsMatrix& A, Int r0, Int c0, Int nr,
                                Int nc, const Int roff, const InitInfo& in,
                                const CrsSegmenter& seg);

  public:
    TMatrix () : nr_(0), coff_(0), tid_os_(0), is_parallel_(false) {}
    ~TMatrix () { clear(); }
    bool empty () const { return is_empty_; }
    bool parallel () const { return is_parallel_; }
    Int get_nr () const { return nr_; }
    void init(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
              const InitInfo& in, const Int block_0_nnz_os = 0,
              const int tid_offset = 0);
    void init(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
              const InitInfo& in, const CrsSegmenter& seg);
    void clear();
    void init_metadata(const CrsMatrix& T, Int r0, Int c0, Int nr, Int nc,
                       const InitInfo& in,
                       // In calculating loads, assume thread 0 is already
                       // handling this many other nonzeros.
                       const Int block_0_nnz_os = 0,
                       // Offset tid by this number. It doesn't make sense to
                       // have block_0_nnz_os > 0 if this argument is > 0.
                       const int tid_offset = 0);
    // Version of the above in which a segmentation is imposed.
    void init_metadata(const CrsMatrix& A, Int r0, Int c0, Int nr, Int nc,
                       const InitInfo& in, const CrsSegmenter& seg);
    void init_memory(const InitInfo& in);
    void init_numeric(const CrsMatrix& T, const int tid);
    void reinit_numeric(const CrsMatrix& A, const int tid);
    // inside || {}
    void n1Axpy(const Sclr* x, const Int ldx, const Int nrhs, Sclr* y,
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
    Int get_n () const { return n_; }
    Int get_r0 () const { return r0_; }
  protected:
    Int n_, r0_;
    static const bool is_lo_ = true;
  };

  class OnDiagTri : public Tri {
    Int c0_;
    Size nnz_;
    Sclr* d_; // Packed by column in dense tri format.
    CrsMatrix* m_;
    bool dense_;

    void solve_spars(const Sclr* b, const Int ldb, Sclr* x, const Int ldx,
                     const Int nrhs) const;
    void solve_dense(const Sclr* b, const Int ldb, Sclr* x, const Int ldx,
                     const Int nrhs) const;
 
  public:
    OnDiagTri();
    ~OnDiagTri () { clear(); }
    void clear();
    void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
              const InitInfo& in);
    // Initialization that takes up space as a block but has no data.
    void init(const Int r0, const Int c0, const Int n);
    // > 1 only if (C) the solve is implemented as an MVP with inv(T).
    Int nthreads() const;
    // Start of tid's block row. > 0 only if C.
    Int block_row_start(const int tid) const;
    // #rows in tid's block row. < n() only if C.
    Int block_nr(const int tid) const;
    void init_metadata(const CrsMatrix& T, const Int r0, const Int c0,
                       const Int n, const InitInfo& in);
    void init_memory(const InitInfo& in);
    void init_numeric(const CrsMatrix& T, const bool invert = false);
    void reinit_numeric(const CrsMatrix& T, const bool invert = false);
    void solve(const Sclr* b, const Int ldb, Sclr* x, const Int ldx,
               const Int nrhs) const;
    Size get_nnz () const { return nnz_; }
    bool inited () const { return d_ || m_; }

  private: // For inverse of on-diag triangle.
    struct Thread {
      Int r0, nr;
      Sclr* d;
      Thread() : r0(0), nr(0), d(0) {}
      ~Thread();
    };
    Array<Thread> t_;
    mutable int done_ctr_;
    void inv_init_metadata(const InitInfo& in);
    void inv_init_memory();
    void inv_reinit_numeric();
    void solve_dense_inv(const Sclr* b, const Int ldb, Sclr* x, const Int ldx,
                         const Int nrhs) const;
  public:
    // Support for external computation of inverses.
    bool is_dense_inverted () const { return ! t_.empty(); }
    Sclr* dense_tri () const { return d_; }
    void inv_copy();
  };

  // Recursive representation. A lower tri looks like this:
  //     [Tri        ]
  //     [TMatrix Tri].
  class RecursiveTri : public Tri {
    Int nthreads_;

    struct NonrecursiveData {
      Array<OnDiagTri> t;
      Array<TMatrix> s;
      Array<Int> os;

      // s k waits until t_barrier >= k.
      mutable Int t_barrier;
      mutable p2p_Done done_symbol;
      mutable Array<p2p_Done> s_done;
      // idx is pointers into ids.
      Array<Size> t_ids;
      Array<Int> t_idx;
      Array< Array<Size> > s_ids;
      Array< Array<Int> > s_idx;
      // For the inverse on-diag tri.
      mutable Array<Int> inv_tri_done;
    } nd_;

    Size rb_p2p_sub2ind (const Int sid, const int tid) const
    { return nd_.s.size()*tid + sid; }
    Int rb_p2p_ind2tid (const Size e) const { return e / nd_.s.size(); }
    Int rb_p2p_ind2sid (const Size e) const { return e % nd_.s.size(); }

    // For the inverse on-diag tri.
    Int max_diag_tri_;
    bool invert_separately_;
    mutable Array<Sclr> wrk_;

    void p2p_init();
    void init_invert_ondiag_tris_separately();
    void invert_ondiag_tris();
    void ondiag_solve(const OnDiagTri& t, Sclr* x, const Int ldx, const Int nrhs,
                      const int tid, const Int stp, volatile Int* const t_barrier,
                      volatile Int* const inv_tri_done) const;

  public:
    RecursiveTri () {}
    ~RecursiveTri () { clear(); }
    void clear();
    // Client should call this:
    void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
              const InitInfo& in, const Int mvp_block_nc = 0);
    void reinit_numeric(const CrsMatrix& T);
    void p2p_reset() const;
    void solve(const Sclr* b, Sclr* x, const Int ldx, const Int nrhs) const;
    void reset_max_nrhs(const Int max_nrhs);
  };

  class LevelSetTri {
    Int n_, ls_blk_sz_;
    Int mvp_block_nc_; // Optional scatter block incorporated into this block.
    Array<Array<Int> > ps_; // Rows for a thread. For init_numeric.
    Array<Int> lsp_; // Level set i is lsp_[i] : lsp_[i+1]-1.
    // Use point-to-point synchronization and the tasking procedure described in
    // Park et al, Sparsifying Synchronizations for High-Performance Shared-Memory
    // Sparse Triangular Solve, ISC 2014.
    struct Thread {
      CrsMatrix* m;
      Array<Int> p;   // Rows this thread owns.
      Array<Int> lsp; // Pointer to start of level set.
      Array<Size> p2p_depends;
      Array<Int> p2p_depends_p;
      Thread () : m(0) {}
      ~Thread () { if (m) delete m; }
    };
    Array<Thread> t_;
    Int nlvls_;
    bool save_for_reprocess_;

    Size ls_p2p_sub2ind (const Int lvl, const int tid) const
    { return nlvls_*tid + lvl; }
    Int ls_p2p_ind2lvl (const Size e) const { return e % nlvls_; }
    Int ls_p2p_ind2tid (const Size e) const { return e / nlvls_; }

  public:
    LevelSetTri () : mvp_block_nc_(0) {}
    // Call this before init().
    void init_lsets(const LevelSetter& lstr, const bool save_for_reprocess);
    // If > 0, the scatter block is attached.
    void set_mvp_block_nc (const Int n) { mvp_block_nc_ = n; }
    void update_permutation(Array<Int>& lsis, const Partition& p);
    void init(const CrsMatrix& T, const Int r0, const Int c0, const Int n,
              const InitInfo& in);
    void init_numeric(const CrsMatrix& T);
    void reinit_numeric(const CrsMatrix& T);
    void solve(const Sclr* b, Sclr* x, const Int ldx, const Int nrhs) const;

  private:
    mutable Array<p2p_Done> p2p_done_;
    mutable p2p_Done p2p_done_value_;
    void find_task_responsible_for_variable(p2p_Pair* const pairs);
    Int fill_graph(const p2p_Pair* const pairs, Size* const g, Size* const gp,
                   Size* const wrk);
    void prune_graph(const Size* const gc, const Size* const gp, Size* const g,
                     Int* const gsz, Size* const wrk, const Int max_gelen);
    void fill_dependencies(const Size* const g, const Size* const gp,
                           const Int* const gsz);
  public:
    void p2p_init();
    void p2p_reset() const;
  };

  // Change between internal and external ordering.
  class Permuter {
    Int n_, max_nrhs_;
    Int* p_, * q_;
    Real* scale_;
    bool is_lo_;
    mutable Sclr* px_;
    Array<Int> part_;

  public:
    Permuter ()
      : n_(0), max_nrhs_(0), p_(0), q_(0), scale_(0), px_(0) {}
    ~Permuter () { clear(); }
    void clear();
    void init(const Int n, const bool is_lo, const Array<Int>& lsis,
              const Array<Int>& dpis, const Int nthreads,
              const Int max_nrhs, const Int* p, const Int* q, const Real* r);
    void reinit_numeric(const Real* r);
    void reset_max_nrhs(const Int max_nrhs);
    void check_nrhs(const Int nrhs) const;
    Sclr* from_outside(const Sclr* x, const Int nrhs, Int ldx=0) const;
    void to_outside(Sclr* x, const Int nrhs, const Sclr a, const Sclr b,
                    Int ldx=0) const;
  };

  // Top-level solver.
  class TriSolver {
    Int n_;
    bool is_lo_, unitdiag_; // Properties to be discovered of the user's T.
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
    void reset_max_nrhs(const Int max_nrhs);
    bool is_lower_tri () const { return is_lo_; }
    bool has_implicit_unit_diag () const { return unitdiag_; }
    // x and b can be the same pointers.
    void solve(const Sclr* b, const Int nrhs, Sclr* x, const Sclr alpha,
               const Sclr beta, const Int ldb, const Int ldx) const;
  };

  struct SparseData {
    Size* ir;
    Int* jc;
    Sclr* d;
    bool dealloc_;
    SparseData () : ir(0), jc(0), d(0), dealloc_(false) {}
    SparseData (const Int m, const Size nnz, const char* fail_msg,
                const bool touch = false)
      : ir(0), jc(0), d(0), dealloc_(false)
    { init(m, nnz, fail_msg, touch); }
    ~SparseData () { if (dealloc_) free(); }
    void init(const Int m, const Size nnz, const char* fail_msg,
              const bool touch = false);
    // Tell SparseData not to deallocate on destruction.
    void release () { dealloc_ = false; }
    void free();
  };

  struct Shape {
    const bool is_lower, is_triangular;
    const bool has_full_diag; // Valid only if is_triangular.
    const bool has_no_diag;   // Ditto.
    Shape (const bool iis_lower, const bool iis_triangular,
           const bool ihas_full_diag, const bool ihas_no_diag)
      : is_lower(iis_lower), is_triangular(iis_triangular),
        has_full_diag(ihas_full_diag), has_no_diag(ihas_no_diag) {}
  };

  class PermVec {
    const Array<Int>& p_;
    Array<Int> pi_;
  public:
    // p is a set of indices into an n-vector.
    PermVec(const Int n, const Array<Int>& p)
      : p_(p)
    {
      pi_.optclear_and_resize(n, -1);
      for (std::size_t i = 0; i < p_.size(); ++i) pi_[p_[i]] = i;
    }
    std::size_t size () const { return p_.size(); }
    // Index into p.
    Int get (const Int i) const { return p_[i]; }
    // Does p contain k?
    bool has (const Int k) const { return pi_[k] >= 0; }
    // Return i for p[i] == k.
    Int to_block (const Int k) const { return pi_[k]; }
  };

  class DenseTrisInverter {
  public:
    struct Tri {
      Sclr* d;
      Int idx, n;
      Tri () : d(0), idx(0), n(0) {}
      Tri (const Int iidx, Sclr* const id, const Int in) : d(id), idx(iidx), n(in) {}
    };
    DenseTrisInverter(Array<Tri>& tris);
    void compute();
  private:
    Array<Tri>& tris_;
    Int nthreads_;
    Array<Sclr> w_;
    void invert(Sclr* const T, const Int n, Sclr* const w);
    void copy(Sclr* const T, const Int n, Sclr* const w);
  };

  static void set_options(const typename ihts::Options& os, Options& od);
  static void print_options(const Options& o, std::ostream& os);
  static void partition_n_uniformly(const Int n, const Int nparts,
                                    Array<Int>& p);
  static void throw_if_nthreads_not_ok(const int nthreads);
  static Int find_first_and_last(
    const Size* const ir, const Int r, const Int* const jc,
    const Int c_first, const Int c_last, Int& i_first, Int& i_last);
  static Size crop_matrix(const CrsMatrix& T, const Box& b, Box& cb);
  static Int decide_level_set_max_index(
    const Array<Int>& N, const Int size_thr, const Options& o);
  static void alloc_lsets(
    const Int lsmi, const Int sns, const Array<Int>& level,
    const Array<Int>& n, typename LevelSetter::LevelSets& lsets);
  static Int locrsrow_schedule_serial (
    const ConstCrsMatrix& L, const Int sns, Array<Int>& w);
  static Int locrsrow_schedule_sns1(
    const ConstCrsMatrix& L, Array<Int>& w, const Options& o);
  static Int locrsrow_schedule(
    const ConstCrsMatrix& L, const Int sns, Array<Int>& w,
    const Options& o);
  static void find_row_level_sets_Lcrs(
    const ConstCrsMatrix& L, const Int sns, Int size_thr,
    typename LevelSetter::LevelSets& lsets, const Options& o);
  static Int upcrscol_schedule_serial(
    const ConstCrsMatrix& U, const Int sns, Array<Int>& w);
  static void find_col_level_sets_Ucrs(
    const ConstCrsMatrix& U, const Int sns, Int size_thr,
    typename LevelSetter::LevelSets& lsets, const Options& o);
  static inline void find_level_sets(
    const ConstCrsMatrix& T, const Int sns, const Int size_thr,
    const bool is_lo, typename LevelSetter::LevelSets& lsets, const Options& o);
  static CrsMatrix* get_matrix_p(const CrsMatrix& A, const Array<Int>& p,
                                 const bool set_diag_reciprocal = false);
  static ConstCrsMatrix* permute_to_other_tri(const ConstCrsMatrix& U);
  static void get_idxs(const Int n, const LevelSetter& lstr,
                       Array<Int>& lsis, Array<Int>& dpis);
  static Shape determine_shape(const ConstCrsMatrix& A);
  static void get_matrix_common_with_covers_all(
    const ConstCrsMatrix& A, const PermVec& pv, const PermVec& qv, Partition& p,
    const bool get_A_idxs, const bool pp);
  static void get_matrix_pp_with_covers_all(
    const ConstCrsMatrix& A, const PermVec& pv, Partition& p,
    const bool get_A_idxs);
  static void get_matrix_p_qp_with_covers_all(
    const ConstCrsMatrix& A, const PermVec& pv, const PermVec& qv, Partition& p,
    const bool get_A_idxs);
  static void sort(Partition& p, const Array<Int>& start);
  static void reverse_A_idxs(const Size nnz, Partition& p);
  static void copy_partition(const ConstCrsMatrix& A, Partition& p);
  static void repartition_into_2_blocks(
    Partition* const p, const ConstCrsMatrix& A);
  static Int count_nnz_lotri(const CrsMatrix& T, const Int r0, const Int c0,
                             const Int n);
  static bool is_dense_tri(const CrsMatrix& T, const Int r0, const Int c0,
                           const Int n, const InitInfo& in, Size* nnz = 0);
  static void find_split_rows(
    const CrsMatrix& T, const Int r0, const Int c0, const Int n,
    const InitInfo& in, std::vector<Int>& split_rows);
  static void build_recursive_tri_r(
    const CrsMatrix& T, const Int r, const Int c, const Int n,
    const InitInfo& in, std::vector<Int>& split_rows,
    std::list<Box>& b);
  static void build_recursive_tri(
    const CrsMatrix& T, const Int r, const Int c, const Int n,
    const Int mvp_block_nc, const InitInfo& in, Array<Box>& bv);
  static ConstCrsMatrix* transpose(
    const ConstCrsMatrix& T, Array<Size>* transpose_perm);
  static void compose_transpose(const Array<Size>& transp_perm, Partition& p);
  static void rbwait(
    volatile p2p_Done* const s_done, const Size* s_ids, const Int* const s_idx,
    const Int i, const p2p_Done done_symbol);
}; // Impl<Int, Size, Sclr>

} // namespace htsimpl

template<typename Int, typename Size, typename Sclr>
struct HTS<Int, Size, Sclr>::CrsMatrix
  : public htsimpl::Impl<Int, Size, Sclr>::ConstCrsMatrix
{
  typedef typename htsimpl::Impl<Int, Size, Sclr>::ConstCrsMatrix CCM;
  CrsMatrix (const Int inrow, const Size* iir, const Int* ijc, const Sclr* id,
             const typename htsimpl::Direction::Enum idir, const bool iconj)
    : CCM(inrow, inrow, iir, ijc, id, idir, iconj,
          false /* don't dealloc */,
          false /* unknown whether implicit unit diag */,
          false /* unknown whether lower tri*/)
  {}
};

template<typename Int, typename Size, typename Sclr>
struct HTS<Int, Size, Sclr>::Impl {
  typename htsimpl::Impl<Int, Size, Sclr>::TriSolver ts;
  typename htsimpl::Impl<Int, Size, Sclr>::Options o;
};

} // namespace Experimental

#endif
