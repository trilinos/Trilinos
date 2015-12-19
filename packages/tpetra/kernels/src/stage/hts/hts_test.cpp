// Unit test for HTS based on randomly generated matrices and no external data.

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <omp.h>

#include "hts_test_util.hpp"

namespace Experimental {
using namespace hts::util;

int are_same (const Data& A, const Data& B) {
#define same(expr, ret) if ( ! (expr)) return ret;
  same(A.m == B.m, 1);
  same(A.ir[A.m] == B.ir[B.m], 3);
  const hts::Size nnz = A.ir[A.m];
  for (hts::Int k = 0; k <= A.m; ++k)
    same(B.ir[k] == A.ir[k], 4);
  for (hts::Size k = 0; k < nnz; ++k) {
    same(B.jc[k] == A.jc[k], 5);
    same(B.v[k] == A.v[k], 6);
  }
  return 0;
#undef same
}

int test_transpose (const bool verbose, const hts::Int n) {
  Data ao, at, a;
  
  TestOptions to;
  to.n = n;
  to.matrix_type = TestOptions::sparse;
  gen_tri_matrix(to, ao);
  
  transpose(ao, at);
  transpose(at, a);

  const int nerr = are_same(ao, a);
  if (nerr && verbose)
    std::cout << "failed: test_transpose\n";
  return nerr;
}

// Run one valid triangle. Solve for x in
//     P' T Q' x = R \ b,
// where P = diag(p) and similarly for Q and R, and T is a triangle.
int test (const TestOptions& to, const bool print_options=false) {
  int nerr = 0;

  const hts::Int max_nrhs = 3;

  // Generate matrix data.
  Data d;
  {
    gen_tri_matrix(to, d);     // triangle
    gen_rand_perm(d.m, d.p);   // row permutation vector
    gen_rand_perm(d.m, d.q);   // col permutation vector
    gen_rand_vector(d.m, d.r); // row scaling
  }

  std::vector<hts::Real> b(d.m*max_nrhs), xt(b.size()), x(b.size());
  {
    // True x.
    gen_rand_vector(xt.size(), xt);
    // Generate the rhs b.
    for (hts::Int irhs = 0; irhs < max_nrhs; ++irhs) {
      const hts::Real* const xtp =xt.data() + irhs*d.m;
      hts::Real* const bp =b.data() + irhs*d.m;
      std::vector<hts::Real> y(d.m);
      for (hts::Int i = 0; i < d.m; ++i)
        x[i] = xtp[d.q[i]];
      mvp(d, to.transpose,x.data(),y.data());
      for (hts::Int i = 0; i < d.m; ++i)
        bp[d.p[i]] = y[i];
      for (hts::Int i = 0; i < d.m; ++i)
        bp[i] *= d.r[i];
    }
  }

  hts::CrsMatrix* T;
  hts::Options opts;
  {
    T = hts::make_CrsMatrix(d.m, d.ir.data(), d.jc.data(), d.v.data(),
                            to.transpose);
    if (to.solve_type == TestOptions::ls_only)
      hts::set_level_schedule_only(opts);
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
    hts::Impl* impl;
    try {
      impl = hts::preprocess(T, max_nrhs, to.nthreads, to.reprocess,
                             d.p.data(), d.q.data(), d.r.data(), &opts);
    } catch (...) {
      hts::delete_CrsMatrix(T);
      throw;
    }
    if (print_options)
      hts::print_options(impl, std::cout);
    if (to.reprocess) {
      // This isn't necessary since we aren't changing the numbers, but pretend
      // we are to test the numerical phase. Do it 3 times to test idempotency.
      for (int rep = 0; rep < 3; ++rep)
        hts::reprocess_numeric(impl, T, d.r.data());
    }
    if (hts::is_lower_tri(impl) &&
        ((to.upper && ! to.transpose) || ( ! to.upper && to.transpose)) &&
        d.m > 1 && d.ir.back() > static_cast<hts::Size>(d.m) /* not diag */)
      ++nerr;
    hts::solve_omp(impl, b.data(), to.nrhs, x.data());
    hts::delete_Impl(impl);
    const double rd = reldif(xt.data(), x.data(), to.nrhs*d.m);
    if (rd >= 1e-10) {
      ++nerr;
      if (to.verbose) std::cout << "rd: " << rd << "\n";
    }
  }
  hts::delete_CrsMatrix(T);

  if (to.verbose)
    std::cout << (nerr ? "fail" : "pass") << "ed: " << to << "\n";

  return nerr;
}

// Run one invalid triangle and make sure an exception is thrown.
int test_for_exception (const TestOptions& to) {
  bool correct_exception_thrown = false;
  try {
    test(to);
  } catch (const hts::NotFullDiagonal&) {
    correct_exception_thrown = to.matrix_type == TestOptions::missing_diag;
  } catch (const hts::NotTriangularException&) {
    correct_exception_thrown = to.matrix_type == TestOptions::not_tri;
  }
  const int nerr = correct_exception_thrown ? 0 : 1;
  if (nerr && to.verbose)
    std::cout << "test_for_exception failed for " << to << "\n";
  return nerr;
}

// Run all the tests.
int test (const bool verbose) {
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
  const int max_nthreads = std::min(on_mic ? 228 : 19, omp_get_max_threads());
  const int nthreads_step = on_mic ? 11 : 1;

  TestOptions to;
  to.block_size = 3; // only for block_sparse
  to.verbose = verbose;
  bool print_options = to.verbose;
  for (int nthreads = 1; nthreads <= max_nthreads; nthreads += nthreads_step) {
    to.nthreads = nthreads;
    for (int ti = 0; ti < 2; ++ti) {
      to.transpose = ti;
      for (size_t ni = 0; ni < sizeof(ns)/sizeof(*ns); ++ni) {
        to.n = ns[ni];
        for (size_t si = 0; si < TestOptions::n_solve_types; ++si) {
          to.solve_type = static_cast<TestOptions::SolveType>(si);
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
  }

  return nerr;
}
} // namespace Experimental

int main (int argc, char** argv) {
  const bool verbose = argc > 1 && std::string(argv[1]) == "-v";
  const int nerr = Experimental::test(verbose);
  std::cout << "HTS Test: " << (nerr ? "Failed" : "Passed") << "\n";
  return 0;
}
