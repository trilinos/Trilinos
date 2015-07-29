#ifndef SHYLU_GTS_TRISOLVE_HPP
#define SHYLU_GTS_TRISOLVE_HPP

#include <exception>
#include <string>

namespace details {
namespace gts {

/* API for a prototype OpenMP-parallelized triangular solver based on
 * preprocessing of the matrix T.
 *   Solve problems of the form
 *     c <- r ? b : b./r.
 *     T x = c
 *     P' T x = c
 *     T Q' x = c
 *     P' T Q' x = c
 * where T is a lower or upper triangular matrix in compressed row storage (CRS)
 * format, and P and Q are permutations, b is in dense format and can be
 * multiple r.h.s, and x is in dense format.
 *   In terms of Matlab's sparse lu, this solver can be used as follows. In
 * Matlab, probably the most efficient way to handle a sparse LU factorization
 * is like this:
 *     % Factorize.
 *     [L U p q R] = lu(A, 'vector');
 *     r = full(diag(R));
 *     % Solve.
 *     b = b./r;
 *     b = b(p);
 *     x = U \ (L \ b);
 *     x(q) = x;
 *   A simple serial solver (solve_serial) is available for sanity checking. For
 * this factorization, it can be used like this:
 *     solve_serial(L, true,  xb, nrhs, p,    NULL, r,    work);
 *     solve_serial(U, false, xb, nrhs, NULL, q,    NULL, work);
 *   An OpenMP-parallelized solver (preprocess + solve_omp) is the main
 * functionality. For this factorization, use this solver as follows.
 *     Impl* Limpl = preprocess(L, nrhs, nthreads, p,    NULL, r);
 *     Impl* Uimpl = preprocess(U, nrhs, nthreads, NULL, q);
 *     solve_omp(Limpl, xb, nrhs); // xb is the rhs b on input.
 *     solve_omp(Uimpl, xb, nrhs); // xb is the solution x on output.
 *
 * Algorithm sketch.
 *   I write the following for a lower triangular matrix T = L, but an upper one
 * U is almost the same.
 *   Although this solver is meant to be indifferent to the source of T, the
 * level set extraction algorithm uses a heuristic that supposes T comes from a
 * factorization. In a factorization that has fill, L typically fills with
 * increasing row number, and U typically fills with increasing column
 * number. Therefore, it makes sense to start a breadth-first search on rows of
 * L at element L(1,1) and columns of U at element U(1,1). In the case of L, the
 * level sets levs{i} so obtained give a sequence of solves on blocks
 * L(levs{i}, levs{i}). In the case of U, the level sets are in reverse order.
 *   The matrix L is permuted to have three blocks:
 *     [ (level sets)  0              ]
 *     [ (MVP)        (recursive tri) ].
 * The MVP block does a matrix-vector product. The sparse data in that block are
 * stored in CRS format, and so parallelization uses no extra memory and has
 * just one synchronization point, the end of the MVP. The recursive tri looks
 * like this:
 *     (recursive tri) =
 *           (dense tri)
 *       |
 *         [ (recursive tri)  0                ]
 *         [ (MVP)             (recursive tri) ].
 * A dense triangle is small and therefore solved serially.
 */

// Integer type of the row pointer and column index arrays.
typedef int Int;
// Real type of the matrix entries.
typedef double Real;

// Opaque types to interact with the solver.
class CrsMatrix;
class Impl;

class Exception : public std::exception {
  std::string msg_;
public:
  Exception (const std::string& msg) : msg_(msg) {}
  virtual ~Exception () throw () {}
  virtual const char* what () const throw () {
    return msg_.c_str();
  }
};

// In the future, it would be nice to handle psychologically triangular
// matrices. For now, the matrix must be lower or upper triangular.
struct NotTriangularException : public Exception {
  NotTriangularException () : Exception("Not a triangular matrix.") {}
};
// The matrix must have a full diagonal.
struct NotFullDiagonal : public Exception {
  NotFullDiagonal () : Exception("Lacks a full diagonal.") {}
};

// Interface between the user's and the solver's data. It's just a shallow
// wrapper to the user's data. The data passed in must persist at least as long
// as the CrsMatrix.
CrsMatrix*
make_CrsMatrix(
  // Number of rows in the matrix.
  const Int nrow,
  // Base-0 index pointers into col. Array has length nrow+1, where rowptr[nrow]
  // is the number of nonzeros.
  const Int* rowptr,
  // Base-0 index pointers to columns.
  const Int* col,
  // Values.
  const Real* val);

// Delete the CrsMatrix. This does not destroy the user's data that T wraps.
void delete_CrsMatrix(CrsMatrix* T);

// Straightforward serial solve implementation for checks.
void solve_serial(
  const CrsMatrix* T,
  // Is lower triangular? If not, then it's upper triangular.
  const bool is_lo,
  // On input, the r.h.s. b; on output, the solution x in T x = b.
  Real* xb,
  // Number of r.h.s.
  const Int nrhs,
  // Optional row permutation P. On input, b <- b(p).
  const Int* p = 0,
  // Optional column permutation Q. On output, x(q) <- x.
  const Int* q = 0,
  // Optional diagonal scaling R. On input, b <- b./r.
  const Real* r = 0,
  // If p or q are provided, then also provide size(T,1) workspace.
  Real* w = 0);

// Preprocess T. After the Impl is constructed, T (and the user's data) can be
// deleted if desired.
Impl*
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
  const Real* r = 0);

// Use T to replace the numbers in impl. The symbolic structure of T must be the
// same as in the call to preprocess.
void reprocess_numeric(Impl* impl, const CrsMatrix* T);

// Ask whether T is lower or upper triangular.
bool is_lower_tri(const Impl* impl);

// Delete the Impl after the solves are completed.
void delete_Impl(Impl* impl);

// Preprocessor-based solve.
void solve_omp(
  Impl* impl,
  // On input, the r.h.s. b; on output, the solution x in T x = b.
  Real* xb,
  // Number of r.h.s.
  const Int nrhs);
// x = T \ b.
void solve_omp(Impl* impl, const Real* b, const Int nrhs, Real* x);
// x = alpha x + beta (T \ b).
void solve_omp(Impl* impl, const Real* b, const Int nrhs, Real* x,
               const Real alpha, const Real beta);

// Print compile-time options defined in gts_trisolve_impl.hpp.
void print_options(std::ostream& os);

struct Options {
  // Recurse until blocks are this size.
  Int min_block_size;
  // Where sparse and dense representations are decided dynamically, choose
  // dense if nonzero density is >= min_dense_density.
  double min_dense_density;
  // Tile rows and cols into nxn blocks. For scalar problems, n=1 (default) is
  // sensible. For PDEs with multiple equations per node (e.g. 3 in 3D
  // elasticity), n = ndof/node is useful.
  Int levelset_block_size;
  // Maximum number of bottlenecks in the level-set solver before switching to
  // the data-parallel solver.
  Int max_level_set_bottlenecks;
  // Minimum level set size. Any smaller is considered a bottleneck.
  Int min_lset_size;
  // Scale min_lset_size with nthreads. If true, effective minimum size is
  //     nthreads*min_lset_size
  // rather than simply min_lset_size.
  bool min_lset_size_scale_with_nthreads;
  // Turn on profiling. Increases computation time.
  bool profile;
  // 0 for no printing.
  int print_level;
  // Sets default values.
  Options();
};

// Optionally set options; not calling this is equivalent to
// set_options(Options()).
void set_options(const Options& o);

// Profiling functions. Write a file containing profile data. If not
// Options::profile, then throw exception.
void write_profile(const std::string& filename);
void write_levelsets(const std::string& filename);
void write_matrixmarket(const CrsMatrix* T, const std::string& filename);
void write_binary(const CrsMatrix* T, const std::string& filename);

} // namespace gts
} // namespace details

#endif
