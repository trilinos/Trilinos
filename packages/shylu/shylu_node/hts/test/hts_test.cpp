// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Unit test for HTS based on randomly generated matrices and no external data.

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <limits>

#include <complex>

#include "hts_test_util.hpp"

namespace Experimental {
namespace htsimpl {
template <typename T> bool is_complex () { return false; }
template <> bool is_complex<std::complex<double> >() { return true; }
template <> bool is_complex<std::complex<float> >() { return true; }

template<typename Int, typename Size, typename Sclr> class Tester {
  typedef HTS<Int, Size, Sclr> ihts;
  typedef typename ihts::Real Real;
  typedef util<Int, Size, Sclr> ut;
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
    const Real tol = std::numeric_limits<Real>::epsilon()*1e6;

    // Generate matrix data.
    Data d;
    Size nnz;
    {
      ut::gen_tri_matrix(to, d);     // triangle
      nnz = d.ir.back();
      ut::gen_rand_perm(d.m, d.p);   // row permutation vector
      ut::gen_rand_perm(d.m, d.q);   // col permutation vector
      ut::gen_rand_vector(d.m, d.r); // row scaling
    }

    const Int ldb = d.m + 3, ldx = d.m + 4;
    std::vector<Sclr> b(ldb*max_nrhs), xt(d.m*max_nrhs), x(ldx*max_nrhs);
    {
      // True x.
      ut::gen_rand_vector(xt.size(), xt);
      // Generate the rhs b.
      std::vector<Sclr> y(d.m);
      for (Int irhs = 0; irhs < max_nrhs; ++irhs) {
        const Sclr* const xtp = xt.data() + irhs*d.m;
        Sclr* const bp = b.data() + irhs*ldb;
        for (Int i = 0; i < d.m; ++i)
          x[i] = xtp[d.q[i]];
        ut::mvp(d, to.transpose, to.conjugate, x.data(), y.data());
        if (to.has_unit_diag())
          for (Int i = 0; i < d.m; ++i)
            y[i] += x[i];
        for (Int i = 0; i < d.m; ++i)
          bp[d.p[i]] = y[i];
        for (Int i = 0; i < d.m; ++i)
          bp[i] *= d.r[i];
      }
    }
    std::vector<Sclr> bo(b);

    typename ihts::CrsMatrix* T;
    typename ihts::Options opts;
    {
      T = ihts::make_CrsMatrix(d.m, d.ir.data(), d.jc.data(), d.v.data(),
                               to.transpose, to.conjugate);
      if ( ! to.reprocess && to.nthreads > 1)
        ihts::register_Deallocator(T, &d);
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
        std::vector<Sclr> Td;
        if (to.reprocess) {
          // Save the true values.
          Td = d.v;
          // Zero the matrix values to simulate reprocessing.
          d.v.assign(d.v.size(), 1);
        }
        impl = ihts::preprocess(T, max_nrhs - 1 /* For testing; see below. */,
                                to.nthreads, to.reprocess, d.p.data(), d.q.data(),
                                d.r.data(), &opts);
        if (to.reprocess) {
          // Restore the values.
          d.v = Td;
        }
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
        // Run 2 times to test idempotency.
        for (int rep = 0; rep < 2; ++rep)
          ihts::reprocess_numeric(impl, T, d.r.data());
      }
      // Exercise reset_max_nrhs.
      ihts::reset_max_nrhs(impl, max_nrhs);
      if (ihts::is_lower_tri(impl) &&
          ((to.upper && ! to.transpose) || ( ! to.upper && to.transpose)) &&
          d.m > 1 && nnz > static_cast<Size>(d.m) /* not diag */)
        ++nerr;
      for (int slv = 0; slv <= 2; ++slv) {
        // Check each solve interface.
        switch (slv) {
        case 0:
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), ldb, ldx);
          break;
        case 1:
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), 0.0,  1.0,
                          ldb, ldx);
          ihts::solve_omp(impl, b.data(), to.nrhs, x.data(), 2.0, -1.0,
                          ldb, ldx);
          break;
        case 2:
          ihts::solve_omp(impl, b.data(), to.nrhs, ldb);
          break;
        }
        double rd = 0;
        for (Int i = 0; i < to.nrhs; ++i)
          rd = std::max(
            rd, ut::reldif(xt.data() + i*d.m,
                           (slv == 2 ? b.data() + i*ldb : x.data() + i*ldx),
                           d.m));
        if (slv == 2) b = bo;
        if (rd >= tol) {
          ++nerr;
          if (to.verbose) std::cout << "rd " << slv << ": " << rd << "\n";
        }
      }
      ihts::delete_Impl(impl);
    }

    if (to.nthreads == 1) {
      std::vector<Sclr> xb(d.m*to.nrhs), w(d.m);
      for (Int irhs = 0; irhs < to.nrhs; ++irhs) {
        const Sclr* const bp = b.data() + irhs*ldb;
        Sclr* const xbp = xb.data() + irhs*d.m;
        for (Int i = 0; i < d.m; ++i)
          xbp[i] = bp[i];
      }
      ihts::solve_serial(T, ! to.upper, to.has_unit_diag(), xb.data(), to.nrhs,
                         d.p.data(), d.q.data(), d.r.data(), w.data());
      const double rd = ut::reldif(xt.data(), xb.data(), d.m*to.nrhs);
      if (rd >= tol) {
        ++nerr;
        if (to.verbose) std::cout << "serial rd: " << rd << "\n";
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
    } catch (const hts::NotFullDiagonalException&) {
      correct_exception_thrown = to.matrix_type == TestOptions::missing_some_diag;
    } catch (const hts::NotTriangularException&) {
      correct_exception_thrown = ut::is_not_tri(to.matrix_type);
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
    int nerr = 0;

    // Test that we throw on an unsigned Int.
    {
      bool caught = false;
      try {
        HTS<unsigned int, int, double>::make_CrsMatrix(1, 0, 0, 0);
      } catch (const hts::Exception& e) {
        caught = true;
      }
      if ( ! caught) ++nerr;
    }

    // Test our own transpose to make sure it's OK for subsequent use.
    nerr += test_transpose(verbose, 277);

    const int ns[] = {1, 2, 3, 21, 300};
    const int max_nthreads = omp_get_max_threads();
    const int nthreads_step = max_nthreads > 40 ? 11 : 3;

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
      for (int ti = 0; ti < 4; ++ti) {
        to.transpose = ti % 2;
        if ( ! is_complex<Sclr>() && ti >= 2) break;
        if (ti >= 2) to.conjugate = true;
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
                  to.matrix_type = TestOptions::dense; nerr += test(to);
                  to.matrix_type = TestOptions::sparse; nerr += test(to);
                  to.matrix_type = TestOptions::block_sparse; nerr += test(to);
                  to.matrix_type = TestOptions::implicit_unit_diag; nerr += test(to);
                  to.matrix_type = TestOptions::block_sparse_implicit_unit_diag; nerr += test(to);
                  if (to.n > 2) {
                    to.matrix_type = TestOptions::not_tri;
                    nerr += test_for_exception(to);
                    to.matrix_type = TestOptions::missing_some_diag;
                    nerr += test_for_exception(to);
                    to.matrix_type = TestOptions::not_tri_almost_diag;
                    nerr += test_for_exception(to);
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
#ifdef HAVE_SHYLU_NODEHTS_COMPLEX
  if (verbose >= 1) std::cout << "\n<int, size_t, std::complex<double> >\n";
  nerr += Experimental::htsimpl::Tester<int, size_t, std::complex<double> >::test(verbose);
#endif
  if (verbose >= 1) std::cout << "\n<short, size_t, float>\n";
  nerr += Experimental::htsimpl::Tester<short, size_t, float>::test(verbose);
  if (verbose >= 1) std::cout << "\n";
  std::cout << "HTS Test: " << (nerr ? "Failed" : "Passed") << "\n";
  return 0;
}
