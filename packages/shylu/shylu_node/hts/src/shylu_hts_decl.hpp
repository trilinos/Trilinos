// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef INCLUDE_SHYLU_HTS_DECL_HPP
#define INCLUDE_SHYLU_HTS_DECL_HPP

#include <exception>
#include <string>

namespace Experimental {

namespace hts {
/*! \brief Wraps std::exception.
 */
class Exception : public std::exception {
  std::string msg_;
public:
  Exception (const std::string& msg) : msg_(msg) {}
  virtual ~Exception () throw () {}
  virtual const char* what () const throw () {
    return msg_.c_str();
  }
};

/*! \brief The matrix is neither lower nor upper triangular.
 */
struct NotTriangularException : public Exception {
  NotTriangularException () : Exception("Not a triangular matrix.") {}
};
/*! \brief The matrix lacks a full diagonal.
 */
struct NotFullDiagonalException : public Exception {
  NotFullDiagonalException () : Exception("Lacks a full diagonal.") {}
};

template <typename T> struct ScalarTraits { typedef typename T::value_type Real; };
template <> struct ScalarTraits<double> { typedef double Real; };
template <> struct ScalarTraits<float> { typedef float Real; };
} // namespace hts

template<typename IntT=int, typename SizeT=IntT, typename SclrT=double>
struct HTS {
  /*! \struct HTS
   *  \brief OpenMP-parallelized sparse triangular solver based on the symbolic
   *         analysis of the nonzero structure of T.
   *
   *   Solve problems of the form
   *     \code
   *     c <- r ? b : b./r
   *     T x = c
   *     P' T x = c
   *     T Q' x = c
   *     P' T Q' x = c
   *     \endcode
   * where T is a lower or upper triangular matrix in compressed row storage (CRS)
   * format, and P and Q are permutations, b is in dense format and can be
   * multiple r.h.s, and x is in dense format.
   *   In terms of Matlab's sparse lu, this solver can be used as follows. In
   * Matlab, probably the most efficient way to handle a sparse LU factorization
   * is like this:
   *     \code
   *     % Factorize.
   *     [L U p q R] = lu(A, 'vector');
   *     r = full(diag(R));
   *     % Solve.
   *     b = b./r;
   *     b = b(p);
   *     x = U \ (L \ b);
   *     x(q) = x;
   *     \endcode
   *   For this factorization, use HTS as follows.
   *     \code
   *     typedef HTS<int, int, double> ihts;
   *     Impl* Limpl = ihts::preprocess(L, nrhs, nthreads, p,    NULL, r);
   *     Impl* Uimpl = ihts::preprocess(U, nrhs, nthreads, NULL, q);
   *     ihts::solve_omp(Limpl, xb, nrhs); // xb is the rhs b on input.
   *     ihts::solve_omp(Uimpl, xb, nrhs); // xb is the solution x on output.
   *     \endcode
   *
   * Int is the index type. Size is the array pointer type in the CRS/CCS data
   * structure; it is appropriate for quantities that have magnitude
   * O(nnz). Sclr is the scalar type, including
   * complex. hts::ScalarTraits<Sclr>::Real is the real type. Int must be
   * signed. Size may be unsigned.
   *
   * Algorithm sketch.
   *   HTS automatically combines two basic algorithms, level scheduling and
   * recursive blocking, to address the spectrum of extremely sparse matrices
   * (level scheduling) to matrices just sparse enough to be efficiently stored in
   * a sparse data structure (recursive blocking).
   *   I write the following for a lower triangular matrix T = L, but an upper one
   * U is almost the same.
   *   The matrix L is permuted to have three blocks:
   *     \code
   *     [ (level sets)  0              ]
   *     [ (MVP)        (recursive tri) ].
   *     \endcode
   * The MVP block does a matrix-vector product. The sparse data in that block are
   * stored in CRS format, and so parallelization uses no extra memory and has
   * just one synchronization point, the end of the MVP. The recursive tri looks
   * like this:
   *     \code
   *     (recursive tri) =
   *           (on-diagonal small tri)
   *       |
   *         [ (recursive tri)  0                ]
   *         [ (MVP)             (recursive tri) ].
   *     \endcode
   *
   * For more, see
   *   A. M. Bradley, A hybrid multithreaded direct sparse triangular solver,
   *   Proc. SIAM CSC, 2016, doi:10.1137/1.9781611974690.ch2.
   */

  typedef IntT Int;
  typedef SizeT Size;
  typedef SclrT Sclr;
  typedef typename hts::ScalarTraits<Sclr>::Real Real;
  typedef hts::Exception Exception;
  typedef hts::NotTriangularException NotTriangularException;
  typedef hts::NotFullDiagonalException NotFullDiagonalException;

  //! \brief Opaque CRS type to interact with the solver.
  struct CrsMatrix;
  //! \brief Opaque handle containing symbolic analysis data.
  struct Impl;

  /*! \brief Construct a shallow wrapper to the user's C[R|S]S matrix.
   *
   * The data passed in must persist at least as long as the CrsMatrix.
   * 
   * This interface is written in terms of a CrsMatrix; use make_transpose in
   * this function and is_lo in the solve functions to solve the various
   * combinations of upper/lower triangles stored in CRS/CCS format.
   *
   * An implicit unit diagonal is permitted. If there are no diagonal elements,
   * a unit diagonal is assumed. If only some diagonal elements are missing,
   * preprocess will throw the exception NotFullDiagonalException.
   *
   * If the matrix is not triangular, preprocess will throw
   * NotTriangularException.
   *
   * \param n [in] Dimension of the matrix.
   *
   * \param rowptr [in] Base-0 index pointers into col. Array has length nrow+1,
   * where rowptr[nrow] is the number of nonzeros.
   *
   * \param col [in] Ordered base-0 index pointers to columns.
   *
   * \param val [in] Matrix values.
   *
   * \param make_transpose [in] Make the transpose of this matrix. Set this to
   * true, for example, to make a CRS matrix from a CCS one.
   *
   * \param make_transpose [in] Take the complex conjugate of the entries of
   * this matrix. Set make_transpose = true, make_conj = true to get the
   * conjugate transpose of this matrix.
   */
  static CrsMatrix* make_CrsMatrix(
    const Int n, const Size* rowptr, const Int* col, const Sclr* val,
    const bool make_transpose = false, const bool make_conj = false);

  /*! \brief Delete the CrsMatrix. Does not destroy the user's data that T wraps.
   *
   * \param T [in] The matrix.
   */
  static void delete_CrsMatrix(CrsMatrix* T);

  struct Deallocator {
    int counter;
    Deallocator () : counter(0) {}
    virtual ~Deallocator () {}
    virtual void free_CrsMatrix_data() = 0;
  };
  static void register_Deallocator(CrsMatrix* mtx, Deallocator* d);

  struct Options;

  //! Print compile-time and run-time options.
  static void print_options(const Impl* impl, std::ostream& os);

  /*! \brief Use level scheduling only.
   *
   * For very sparse problems, it can be efficient to turn off the symbolic
   * analysis's block finder and instead instruct it to level schedule only. A
   * typical use case is Gauss-Seidel smoothing of a volume discretization.
   */
  static void set_level_schedule_only(Options& o);

  struct PreprocessArgs {
    const CrsMatrix* T;
    Int max_nrhs;
    Int nthreads;
    // - Optional
    bool save_for_reprocess;
    const Int* p;
    const Int* q;
    const Real* scale_rhs;
    bool scale_rhs_by_division;
    const Real* scale_solution;
    bool scale_solution_by_division;
    const Options* options;

    PreprocessArgs();
  };

  static Impl* preprocess(const PreprocessArgs& args);

  // Preprocess T. After the Impl is constructed, T (and the user's data) can be
  // deleted if desired.
  static Impl*
  preprocess(
    const CrsMatrix* T,
    // Max number of r.h.s. to be used with this Impl.
    const Int max_nrhs,
    // Number of OpenMP threads this Impl will use.
    const Int nthreads,
    // Optionally save nonzero structural data for fast injection of new numbers
    // with the same nonzero structure. On problems large relative to system
    // memory, this can slow down solves.
    const bool save_for_reprocess = false,
    // Optional row permutation P. On input, b <- b(p).
    const Int* p = 0,
    // Optional column permutation Q. On output, x(q) <- x.
    const Int* q = 0,
    // Optional diagonal scaling R. On input, b <- b./r.
    const Real* r = 0,
    // Options for profiling and performance tuning.
    const Options* options = 0);

  // Use T to replace the numbers in impl. The symbolic structure of T must be
  // the same as in the call to preprocess.
  static void reprocess_numeric(
    Impl* impl, const CrsMatrix* T, const Real* r = 0);

  //! \brief Ask whether T is lower or upper triangular.
  static bool is_lower_tri(const Impl* impl);

  //! \brief Ask whether T has an implicit unit diagonal.
  static bool has_implicit_unit_diag(const Impl* impl);

  //! \brief Delete the Impl after the solves are completed.
  static void delete_Impl(Impl* impl);

  static void reset_max_nrhs(Impl* impl, const Int max_nrhs);

  // Preprocessor-based solve.
  static void solve_omp(
    Impl* impl,
    // On input, the r.h.s. b; on output, the solution x in T x = b.
    Sclr* xb,
    // Number of r.h.s.
    const Int nrhs,
    const Int ldxb=0);
  // x = T \ b. b and x can be the same.
  static void solve_omp(Impl* impl, const Sclr* b, const Int nrhs, Sclr* x,
                        const Int ldb=0, const Int ldx=0);
  // x = alpha x + beta (T \ b). b and x can be the same.
  static void solve_omp(Impl* impl, const Sclr* b, const Int nrhs, Sclr* x,
                        const Sclr alpha, const Sclr beta,
                        const Int ldb=0, const Int ldx=0);

  struct Options {
    // - Obvious user parameters.
    // Tile rows and cols into nxn blocks. For scalar problems, n=1 (default) is
    // sensible. For PDEs with multiple equations per node (e.g. 3 in 3D
    // elasticity), n = ndof/node is useful.
    Int levelset_block_size;
    // Minimum level set size. Any smaller is considered a bottleneck. Set to 0
    // to make hts solve the whole problem using level scheduling.
    Int min_lset_size;

    // - Nonobvious performance tuning parameters.
    // Recurse until blocks are this size.
    Int min_block_size;
    // Guide to how to break up blocks for parallel execution.
    Int min_parallel_rows;
    // In the symbolic analysis phase, serial block sizes need to be larger
    // because memory access patterns are terrible.
    Int pp_min_block_size;
    // Where sparse and dense representations are decided dynamically, choose
    // dense if nonzero density is >= min_dense_density.
    Real min_dense_density;
    // Max fraction of rows in the LS block that are in too-small level sets.
    Real lset_max_bad_fraction;
    // Scale min_lset_size with nthreads. If true, effective minimum size is
    //     nthreads*lset_min_size
    // rather than simply min_lset_size.
    bool lset_min_size_scale_with_nthreads;

    // - Debug/output parameters.
    // Turn on profiling. Increases computation time.
    bool profile;
    // 0 for no printing.
    int print_level;

    // Sets default values.
    Options();
  };

  // For testing.

  // Straightforward serial solve implementation for checks.
  static void solve_serial(
    const CrsMatrix* T,
    // Is lower triangular? If not, then it's upper triangular.
    const bool is_lo,
    // Has an implicit unit diagonal?
    const bool implicit_unit_diag,
    // On input, the r.h.s. b; on output, the solution x in T x = b.
    Sclr* xb,
    // Number of r.h.s.
    const Int nrhs,
    // Optional row permutation P. On input, b <- b(p).
    const Int* p = 0,
    // Optional column permutation Q. On output, x(q) <- x.
    const Int* q = 0,
    // Optional diagonal scaling R. On input, b <- b./r.
    const Real* r = 0,
    // If p or q are provided, then also provide size(T,1) workspace.
    Sclr* w = 0);
};

} // namespace Experimental

#endif
