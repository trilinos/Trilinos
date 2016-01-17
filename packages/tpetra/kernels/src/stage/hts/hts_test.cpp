// Unit test for HTS based on randomly generated matrices and no external data.

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>

#include <omp.h>

#include "hts_test_util.hpp"

namespace Experimental {
namespace htsimpl {
template<typename Int, typename Size, typename Real> class Tester {
  typedef HTS<Int, Size, Real> ihts;
  typedef util<Int, Size, Real> ut;
  typedef typename ut::TestOptions TestOptions;
  typedef typename ut::Data Data;

  static int are_same (const Data& A, const Data& B) {
#define same(expr, ret) if ( ! (expr)) return ret;
    same(A.m == B.m, 1);
    same(A.ir[A.m] == B.ir[B.m], 3);
    const Size nnz = A.ir[A.m];
    for (Int k = 0; k <= A.m; ++k)
      same(B.ir[k] == A.ir[k], 4);
    for (Size k = 0; k < nnz; ++k) {
      same(B.jc[k] == A.jc[k], 5);
      same(B.v[k] == A.v[k], 6);
    }
    return 0;
#undef same
  }

  static int test_transpose (const int verbose, const Int n) {
    Data ao, at, a;
  
    TestOptions to;
    to.n = n;
    to.matrix_type = TestOptions::sparse;
    ut::gen_tri_matrix(to, ao);
  
    ut::transpose(ao, at);
    ut::transpose(at, a);

    const int nerr = are_same(ao, a);
    if (nerr && verbose)
      std::cout << "failed: test_transpose\n";
    return nerr;
  }

  // Run one valid triangle. Solve for x in
  //     P' T Q' x = R \ b,
  // where P = diag(p) and similarly for Q and R, and T is a triangle.
  static int test (const TestOptions& to, const bool print_options=false,
                   const bool exception_expected=false) {
    int nerr = 0;

    const Int max_nrhs = 3;

    // Generate matrix data.
    Data d;
    {
      ut::gen_tri_matrix(to, d);     // triangle
      ut::gen_rand_perm(d.m, d.p);   // row permutation vector
      ut::gen_rand_perm(d.m, d.q);   // col permutation vector
      ut::gen_rand_vector(d.m, d.r); // row scaling
    }

    std::vector<Real> b(d.m*max_nrhs), xt(b.size()), x(b.size());
    {
      // True x.
      ut::gen_rand_vector(xt.size(), xt);
      // Generate the rhs b.
      for (Int irhs = 0; irhs < max_nrhs; ++irhs) {
        const Real* const xtp =xt.data() + irhs*d.m;
        Real* const bp =b.data() + irhs*d.m;
        std::vector<Real> y(d.m);
        for (Int i = 0; i < d.m; ++i)
          x[i] = xtp[d.q[i]];
        ut::mvp(d, to.transpose,x.data(),y.data());
        for (Int i = 0; i < d.m; ++i)
          bp[d.p[i]] = y[i];
        for (Int i = 0; i < d.m; ++i)
          bp[i] *= d.r[i];
      }
    }

    typename ihts::CrsMatrix* T;
    typename ihts::Options opts;
    {
      T = ihts::make_CrsMatrix(d.m, d.ir.data(), d.jc.data(), d.v.data(),
                               to.transpose);
      if (to.solve_type == TestOptions::ls_only)
        ihts::set_level_schedule_only(opts);
      else if (to.solve_type == TestOptions::rb_only)
        opts.min_lset_size = d.m + 1;
      // To really test things well, choose very small block sizes. (This is bad
      // for performance.) These parameters are not meant to be set by the user
      // ordinarily, but they are exposed in case an expert is tuning performance
      // or, as here, for testing.
      opts.min_block_size = 6;
      opts.min_parallel_rows = 2;
      opts.pp_min_block_size = 12;
      if (to.matrix_type == TestOptions::block_sparse)
        opts.levelset_block_size = to.block_size;
    }

    {
      typename ihts::Impl* impl;
      try {
        impl = ihts::preprocess(T, max_nrhs - 1 /* For testing; see below. */,
                                to.nthreads, to.reprocess, d.p.data(), d.q.data(),
                                d.r.data(), &opts);
      } catch (...) {
        if ( ! exception_expected) {
          std::cerr << "Unexpected exception on ";
          to.print(std::cerr);
          std::cerr << "\n";
          ut::write_matrixmarket(d, "unexpected_exception.mm");
        }
        ihts::delete_CrsMatrix(T);
        throw;
      }
      if (print_options)
        ihts::print_options(impl, std::cout);
      if (to.reprocess) {
        // This isn't necessary since we aren't changing the numbers, but pretend
        // we are to test the numerical phase. Do it 3 times to test idempotency.
        for (int rep = 0; rep < 3; ++rep)
          ihts::reprocess_numeric(impl, T, d.r.data());
      }
      // Exercise reset_max_nrhs.
      ihts::reset_max_nrhs(impl, max_nrhs);
      if (ihts::is_lower_tri(impl) &&
          ((to.upper && ! to.transpose) || ( ! to.upper && to.transpose)) &&
          d.m > 1 && d.ir.back() > static_cast<Size>(d.m) /* not diag */)
        ++nerr;
      for (int rep = 0; rep < 3; ++rep)
        ihts::solve_omp(impl, b.data(), to.nrhs, x.data());
      ihts::delete_Impl(impl);
      const double rd = ut::reldif(xt.data(), x.data(), to.nrhs*d.m);
      if (rd >= std::numeric_limits<Real>::epsilon()*1e6) {
        ++nerr;
        if (to.verbose) std::cout << "rd: " << rd << "\n";
      }
    }
    ihts::delete_CrsMatrix(T);

    if (to.verbose) {
      const bool print = to.verbose == 2 || (to.verbose == 1 && nerr);
      if (print) {
        std::cout << (nerr ? "fail" : "pass") << "ed: ";
        to.print(std::cout);
        std::cout << "\n";
      }
    }

    return nerr;
  }

  // Run one invalid triangle and make sure an exception is thrown.
  static int test_for_exception (const TestOptions& to) {
    bool correct_exception_thrown = false;
    try {
      test(to, false, true);
    } catch (const hts::NotFullDiagonal&) {
      correct_exception_thrown = to.matrix_type == TestOptions::missing_diag;
    } catch (const hts::NotTriangularException&) {
      correct_exception_thrown = to.matrix_type == TestOptions::not_tri;
    }
    const int nerr = correct_exception_thrown ? 0 : 1;
    if (nerr && to.verbose) {
      std::cout << "test_for_exception failed for ";
      to.print(std::cout);
      std::cout << "\n";      
    }
    return nerr;
  }

public:
  // Run all the tests.
  static int test (const int verbose) {
    const bool on_mic =
#ifdef __MIC__
      true
#else
      false
#endif
      ;

    int nerr = 0;

    // Test our own transpose to make sure it's OK for subsequent use.
    nerr += test_transpose(verbose, 277);

    const int ns[] = {1, 2, 11, 300};
    const int max_nthreads = std::min(on_mic ? 240 : 32, omp_get_max_threads());
    const int nthreads_step = on_mic ? 11 : 3;

    TestOptions to;
    to.block_size = 3; // only for block_sparse
    to.verbose = verbose;
    bool print_options = to.verbose;
    for (int nthreads = 1; ; ) {
      to.nthreads = nthreads;
      if (to.verbose >= 1) {
        std::cout << " " << nthreads;
        std::cout.flush();
      }
      for (int ti = 0; ti < 2; ++ti) {
        to.transpose = ti;
        for (size_t ni = 0; ni < sizeof(ns)/sizeof(*ns); ++ni) {
          to.n = ns[ni];
          for (size_t si = 0; si < TestOptions::n_solve_types; ++si) {
            to.solve_type = static_cast<typename ut::TestOptions::SolveType>(si);
            for (int ui = 0; ui < 2; ++ui) {
              to.upper = ui;
              for (int ri = 0; ri < 2; ++ri) {
                for (int nrhs = 1; nrhs <= 3; nrhs += 2) {
                  to.nrhs = nrhs;
                  to.reprocess = ri;
                  to.matrix_type = TestOptions::diag;
                  nerr += test(to, print_options);
                  print_options = false;
                  to.matrix_type = TestOptions::dense;  nerr += test(to);
                  to.matrix_type = TestOptions::sparse; nerr += test(to);
                  to.matrix_type = TestOptions::block_sparse; nerr += test(to);
                  if (to.n > 2) {
                    to.matrix_type = TestOptions::not_tri; nerr += test_for_exception(to);
                    to.matrix_type = TestOptions::missing_diag; nerr += test_for_exception(to);
                  }
                }
              }
            }
          }
        }
      }
      if (nthreads == max_nthreads) break;
      nthreads = std::min(max_nthreads, nthreads + nthreads_step);
    }

    return nerr;
  }
};
} // namespace htsimpl
} // namespace Experimental

int main (int argc, char** argv) {
  int verbose = 0;
  if (argc > 1) {
    const std::string s(argv[1]);
    verbose = s == "-v" ? 1 : s == "-V" ? 2 : 0;
  }
  if (verbose >= 1) std::cout << "<int, int, double>\n";
  int nerr = Experimental::htsimpl::Tester<int, int, double>::test(verbose);
  if (verbose >= 1) std::cout << "\n<short, size_t, float>\n";
  nerr += Experimental::htsimpl::Tester<short, size_t, float>::test(verbose);
  if (verbose >= 1) std::cout << "\n";
  std::cout << "HTS Test: " << (nerr ? "Failed" : "Passed") << "\n";
  return 0;
}
