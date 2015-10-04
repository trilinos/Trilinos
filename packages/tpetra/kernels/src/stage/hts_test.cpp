// Unit test for HTS based on randomly generated matrices and no external data.

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <omp.h>

#include "hts.hpp"

namespace Experimental {
struct TestOptions {
  enum MatrixType { diag, dense, sparse, missing_diag, not_tri };
  enum SolveType { hybrid = 0, ls_only, rb_only, n_solve_types };

  int n, nthreads;
  bool verbose, reprocess, upper;
  MatrixType matrix_type;
  SolveType solve_type;

  TestOptions ()
    : n(0), nthreads(1), verbose(false), reprocess(false), upper(false),
      matrix_type(sparse), solve_type(hybrid)
  {}
};

std::ostream& operator<< (std::ostream& os, const TestOptions& to) {
  os << "[n " << to.n << " reprocess " << to.reprocess << " upper " << to.upper
     << " matrix_type " << to.matrix_type << " solve_type " << to.solve_type
     << " nthreads " << to.nthreads << "]";
  return os;
}

struct Data {
  hts::Int m;
  std::vector<hts::Int> ir, jc, p, q;
  std::vector<hts::Real> v, r;
  void clear () { ir.clear(); jc.clear(); v.clear(); }
};

// Return a uniform random [0, 1).
double urand () {
#if 0 // Not all compilers have this, it seems.
  static std::default_random_engine generator;
  static std::uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
#else
  return rand() / ((double) RAND_MAX + 1.0);
#endif
}

// Return a random diagonal entry.
inline double gen_diag () {
  double diag = urand() - 0.5;
  diag += (diag > 0 ? 1 : -1);
  return diag;
}

// Make a uniform [-0.5, 0.5) random vector.
void gen_rand_vector (const size_t n, std::vector<hts::Real>& x) {
  x.resize(n);
  for (size_t i = 0; i < n; ++i)
    x[i] = urand() - 0.5;
}

// Make a random permutation vector.
void gen_rand_perm (const size_t n, std::vector<hts::Int>& p) {
  p.resize(n);
  for (size_t i = 0; i < n; ++i)
    p[i] = i;
  for (size_t i = 0; i < n; ++i) {
    const int j = urand()*n, k = urand()*n;
    std::swap(p[j], p[k]);
  }
}

// Make a random diagonal matrix.
void gen_diag_matrix (const TestOptions& to, Data& d) {
  for (int r = 0; r < to.n; ++r) {
    d.jc.push_back(r);
    d.v.push_back(gen_diag());
    d.ir[r+1] = d.ir[r] + 1;
  }
}

// y = a*x.
void mvp (const Data& a, const hts::Real* const x, hts::Real* const y) {
  using hts::Int;
  using hts::Real;
  const Int* const ir = &a.ir[0];
  const Int* const jc = &a.jc[0];
  const Real* const d = &a.v[0];
  for (Int i = 0; i < a.m; ++i) {
    Real acc = 0;
    const Int iri = ir[i], irip1 = ir[i+1];
    for (Int j = iri; j < irip1; ++j)
      acc += d[j]*x[jc[j]];
    y[i] = acc;
  }
}

template<typename T> inline double square (const T& x) { return x*x; }

// norm(a - b, 2)/norm(a, 2).
double reldif (const std::vector<hts::Real>& a,
               const std::vector<hts::Real>& b) {
  assert(a.size() == b.size());
  double num = 0, den = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    num += square(a[i] - b[i]);
    den += square(a[i]);
  }
  return std::sqrt(num/den);
}

// Remove or add an entry to make the matrix non-triangular. This will test the
// shape-detection code.
void test_shape_checking (const TestOptions& to, Data& d, const int r) {
  if (to.matrix_type == TestOptions::missing_diag && r == to.n/2) {
    // Remove the previous diagonal entry.
    --d.ir[r+1];
    d.jc.pop_back();
    d.v.pop_back();
  } else if (to.matrix_type == TestOptions::not_tri && r == to.n/2) {
    // Add an entry on the wrong tri side.
    d.jc.push_back(to.upper ? std::max(0, r-1) : std::min(to.n-1, r+1));
    d.v.push_back(0.5);
    ++d.ir[r+1];
  }
}

// Make a random sparse matrix with probability of an off-diag entry set to
// p. If p is not specified or is negative, set it to min(1, 20/n).
void gen_tri_sparse_matrix (const TestOptions& to, Data& d,
                            double p = -1 /* density */) {
  if (p < 0) p = std::min(1.0, 20.0/to.n);
  std::vector<double> row_val;
  std::vector<hts::Int> row_col;
  for (int r = 0; r < to.n; ++r) {
    double diag;
    {
      // Generate the off-diag row and a diag. Scale the off-diag row to keep
      // the condition number reasonable.
      row_val.clear();
      row_col.clear();
      const int ntrial = to.upper ? to.n - r - 1 : r;
      const int base = to.upper ? r + 1 : 0;
      for (int i = 0; i < ntrial; ++i)
        if (urand() < p) {
          row_val.push_back(urand() - 0.5);
          row_col.push_back(static_cast<hts::Int>(base + i));
        }
      double sum = 0;
      for (std::size_t i = 0; i < row_val.size(); ++i)
        sum += std::abs(row_val[i]);
      diag = gen_diag();
      const double scale = 0.1*std::abs(diag)/std::max(1.0, sum);
      for (std::size_t i = 0; i < row_val.size(); ++i)
        row_val[i] *= scale;
    }
    d.ir[r+1] = d.ir[r];
    if (to.upper) {
      ++d.ir[r+1];
      d.jc.push_back(r);
      d.v.push_back(diag);
      test_shape_checking(to, d, r);
      d.ir[r+1] += static_cast<hts::Int>(row_val.size());
      d.jc.insert(d.jc.end(), row_col.begin(), row_col.end());
      d.v.insert(d.v.end(), row_val.begin(), row_val.end());
    } else {
      d.ir[r+1] += static_cast<hts::Int>(row_val.size());
      d.jc.insert(d.jc.end(), row_col.begin(), row_col.end());
      d.v.insert(d.v.end(), row_val.begin(), row_val.end());
      ++d.ir[r+1];
      d.jc.push_back(r);
      d.v.push_back(diag);
      test_shape_checking(to, d, r);
    }
  }
}

// Make a random dense matrix (but in sparse format).
void gen_tri_dense_matrix (const TestOptions& to, Data& d) {
  gen_tri_sparse_matrix(to, d, 1);
}

// Make a random triangular matrix according to the TestOptions.
void gen_tri_matrix (const TestOptions& to, Data& d) {
  d.clear();
  d.m = to.n;
  d.ir.resize(d.m+1, 0);
  switch (to.matrix_type) {
  case TestOptions::diag: gen_diag_matrix(to, d); break;
  case TestOptions::dense: gen_tri_dense_matrix(to, d); break;
  case TestOptions::sparse:
  case TestOptions::missing_diag:
  case TestOptions::not_tri:
    gen_tri_sparse_matrix(to, d); break;
  }
  assert((std::size_t) d.ir.back() == d.v.size());
  assert((std::size_t) d.ir.back() == d.jc.size());
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
    gen_tri_matrix(to, d);      // triangle
    gen_rand_perm(to.n, d.p);   // row permutation vector
    gen_rand_perm(to.n, d.q);   // col permutation vector
    gen_rand_vector(to.n, d.r); // row scaling
  }

  std::vector<hts::Real> b(to.n*max_nrhs), xt(b.size()), x(b.size());
  {
    // True x.
    gen_rand_vector(xt.size(), xt);
    // Generate the rhs b.
    for (hts::Int irhs = 0; irhs < max_nrhs; ++irhs) {
      const hts::Real* const xtp = &xt[0] + irhs*d.m;
      hts::Real* const bp = &b[0] + irhs*d.m;
      std::vector<hts::Real> y(d.m);
      for (hts::Int i = 0; i < d.m; ++i)
        x[i] = xtp[d.q[i]];
      mvp(d, &x[0], &y[0]);
      for (hts::Int i = 0; i < d.m; ++i)
        bp[d.p[i]] = y[i];
      for (hts::Int i = 0; i < d.m; ++i)
        bp[i] *= d.r[i];
    }
  }

  hts::CrsMatrix* T;
  hts::Options opts;
  {
    T = hts::make_CrsMatrix(to.n, &d.ir[0], &d.jc[0], &d.v[0]);
    if (to.solve_type == TestOptions::ls_only)
      hts::set_level_schedule_only(opts);
    else if (to.solve_type == TestOptions::rb_only)
      opts.min_lset_size = to.n + 1;
    // To really test things well, choose very small block sizes. (This is bad
    // for performance.) These parameters are not meant to be set by the user
    // ordinarily, but they are exposed in case an expert is tuning performance
    // or, as here, for testing.
    opts.min_block_size = 6;
    opts.min_parallel_rows = 2;
    opts.pp_min_block_size = 12;
  }
  
  {
    hts::Impl* impl;
    try {
      impl = hts::preprocess(T, max_nrhs, to.nthreads, to.reprocess,
                             &d.p[0], &d.q[0], &d.r[0], &opts);
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
        hts::reprocess_numeric(impl, T);
    }
    if (hts::is_lower_tri(impl) && to.upper && to.n > 1)
      ++nerr;
    hts::solve_omp(impl, &b[0], max_nrhs, &x[0]);
    hts::delete_Impl(impl);
    if (reldif(xt, x) >= 1e-13)
      ++nerr;
  }
  hts::delete_CrsMatrix(T);

  if (nerr && to.verbose)
    std::cout << "test failed for " << to << "\n";
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

  const int ns[] = {1, 2, 10, 300};
  const int max_nthreads = std::min(on_mic ? 228 : 19, omp_get_max_threads());
  const int nthreads_step = on_mic ? 11 : 1;
  
  int nerr = 0;
  TestOptions to;
  to.verbose = verbose;
  bool print_options = to.verbose;
  for (size_t ni = 0; ni < sizeof(ns)/sizeof(*ns); ++ni) {
    to.n = ns[ni];
    for (int nthreads = 1; nthreads <= max_nthreads; nthreads += nthreads_step) {
      to.nthreads = nthreads;
      for (size_t si = 0; si < TestOptions::n_solve_types; ++si) {
        to.solve_type = static_cast<TestOptions::SolveType>(si);
        for (int ri = 0; ri < 2; ++ri) {
          to.reprocess = ri;
          to.matrix_type = TestOptions::diag;   to.upper = false;
          nerr += test(to, print_options);
          print_options = false;
          to.matrix_type = TestOptions::dense;  to.upper = false; nerr += test(to);
          to.matrix_type = TestOptions::dense;  to.upper = true;  nerr += test(to);
          to.matrix_type = TestOptions::sparse; to.upper = false; nerr += test(to);
          to.matrix_type = TestOptions::sparse; to.upper = true;  nerr += test(to);
          if (to.n > 2) {
            to.matrix_type = TestOptions::not_tri; to.upper = false;
            nerr += test_for_exception(to);
            to.matrix_type = TestOptions::not_tri; to.upper = true;
            nerr += test_for_exception(to);
            to.matrix_type = TestOptions::missing_diag; to.upper = false;
            nerr += test_for_exception(to);
            to.matrix_type = TestOptions::missing_diag; to.upper = true;
            nerr += test_for_exception(to);
          }
        }
      }
    }
  }

  return nerr;
}

// If I want to write out the generated matrix. Not used in practice.
void write_matrixmarket (const Data& T, const std::string& filename) {
  FILE* fid = fopen(filename.c_str(), "w");
  if ( ! fid) return;
  fprintf(fid, "%%%%MatrixMarket matrix coordinate real general\n"
          "%11d %11d %11d\n", T.m, T.m, T.ir[T.m]);
  for (hts::Int r = 0; r < T.m; ++r)
    for (hts::Int j = T.ir[r]; j < T.ir[r+1]; ++j)
      fprintf(fid, "%d %d %1.15e\n", r+1, T.jc[j] + 1, T.v[j]);
  fclose(fid);
}
} // namespace Experimental

int main (int argc, char** argv) {
  const bool verbose = argc > 1 && std::string(argv[1]) == "-v";
  const int nerr = Experimental::test(verbose);
  std::cout << "HTS Test: " << (nerr ? "Failed" : "Passed") << "\n";
  return 0;
}
