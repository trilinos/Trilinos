#ifndef INCLUDE_HTS_TEST_UTIL_HPP
#define INCLUDE_HTS_TEST_UTIL_HPP

#include <vector>
#include <iostream>
#include <sstream>
#include <sys/time.h>

#include "hts.hpp"

namespace Experimental {
namespace hts {
namespace util {

// Simple RAII timer.
class Timer {
  std::string name_;
  timeval t1_;
  static double calc_et (const timeval& t1, const timeval& t2) {
    static const double us = 1.0e6;
    return (t2.tv_sec*us + t2.tv_usec - t1.tv_sec*us - t1.tv_usec)/us;
  }
public:
  Timer (const std::string& name) : name_(name) { gettimeofday(&t1_, 0); }
  ~Timer () {
    timeval t2;
    gettimeofday(&t2, 0);
    std::cout << "T " << name_ << "\nT   " << calc_et(t1_, t2) << "\n";
  }
};

struct TestOptions {
  enum MatrixType { diag, dense, sparse, missing_diag, not_tri, block_sparse };
  enum SolveType { hybrid = 0, ls_only, rb_only, n_solve_types };

  int n, nthreads, block_size, nrhs;
  bool verbose, reprocess, upper, transpose;
  MatrixType matrix_type;
  SolveType solve_type;
  double density;

  TestOptions ()
    : n(0), nthreads(1), block_size(1), nrhs(1), verbose(false),
      reprocess(false), upper(false), transpose(false), matrix_type(sparse),
      solve_type(hybrid), density(-1)
  {}
};

struct Data {
  Int m;
  std::vector<Size> ir;
  std::vector<Int> jc, p, q;
  std::vector<Real> v, r;
  void clear () {
    m = 0;
    ir.clear(); jc.clear(); v.clear();
    p.clear(); q.clear(); r.clear();
  }
};

std::ostream& operator<<(std::ostream& os, const TestOptions& to);

// Return a uniform random [0, 1).
double urand();

// Return a random diagonal entry.
inline double gen_diag();

// Make a uniform [-0.5, 0.5) random vector.
void gen_rand_vector(const size_t n, std::vector<Real>& x);

// Make a random permutation vector.
void gen_rand_perm(const size_t n, std::vector<Int>& p);

// Make a random diagonal matrix.
void gen_diag_matrix(const TestOptions& to, Data& d);

// y = a*x.
void mvp(const Data& a, const bool transp, const Real* const x, Real* const y);

// norm(a - b, 2)/norm(a, 2).
double reldif(const Real* const a, const Real* const b, const Size n);

// Remove or add an entry to make the matrix non-triangular. This will test the
// shape-detection code.
void test_shape_checking(const TestOptions& to, Data& d, const int r);

// Make a random triangular matrix according to the TestOptions.
void gen_tri_matrix(const TestOptions& to, Data& d);

void write_matrixmarket(const Data& T, const std::string& filename);

std::ostream& operator<<(std::ostream& os, const Data& d);

void transpose(const Data& src, Data& dst);

}}}

#endif
