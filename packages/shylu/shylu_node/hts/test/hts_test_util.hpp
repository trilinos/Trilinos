// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef INCLUDE_HTS_TEST_UTIL_HPP
#define INCLUDE_HTS_TEST_UTIL_HPP

#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>
#include <sys/time.h>

#include "shylu_hts.hpp"

namespace Experimental {
namespace htsimpl {
static double urand () {
#if 0 // Not all compilers have this, it seems.
  static std::default_random_engine generator;
  static std::uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
#else
  return rand() / ((double) RAND_MAX + 1.0);
#endif
}

template <typename T> inline T rand_scalar () { return urand() - 0.5; }
template<> inline std::complex<double> rand_scalar<std::complex<double> > ()
{ return std::complex<double>(rand_scalar<double>(), rand_scalar<double>()); }
template<> inline std::complex<float> rand_scalar<std::complex<float> > ()
{ return std::complex<float>(rand_scalar<float>(), rand_scalar<float>()); }

template <typename T> inline T real_part (const T& v) { return v; }
template <typename T> inline T real_part (const std::complex<T>& v)
{ return v.real(); }

template<typename Int, typename Size, typename Sclr> struct util {
  typedef typename hts::ScalarTraits<Sclr>::Real Real;

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
    enum MatrixType { diag, dense, sparse, block_sparse, implicit_unit_diag,
                      block_sparse_implicit_unit_diag, missing_some_diag,
                      not_tri, not_tri_almost_diag };
    enum SolveType { hybrid = 0, ls_only, rb_only, n_solve_types };

    int n, nthreads, block_size, nrhs;
    bool reprocess, upper, transpose, conjugate;
    int verbose;
    MatrixType matrix_type;
    SolveType solve_type;
    double density;

    TestOptions ()
      : n(0), nthreads(1), block_size(1), nrhs(1), reprocess(false),
        upper(false), transpose(false), conjugate(false), verbose(0),
        matrix_type(sparse), solve_type(hybrid), density(-1)
    {}

    void print (std::ostream& os) const {
      os << "[n " << n << " reprocess " << reprocess << " upper " << upper
         << " matrix_type " << matrix_type << " transpose " << transpose
         << " conjugate " << conjugate << " solve_type " << solve_type
         << " nrhs " << nrhs << " nthreads " << nthreads << "]";
    }

    bool has_unit_diag () const {
      return (matrix_type == implicit_unit_diag ||
              matrix_type == block_sparse_implicit_unit_diag);
    }
  };

  struct Data : public HTS<Int, Size, Sclr>::Deallocator {
    Int m;
    std::vector<Size> ir;
    std::vector<Int> jc, p, q;
    std::vector<Sclr> v;
    std::vector<Real> r;
    void clear () {
      m = 0;
      ir.clear(); jc.clear(); v.clear();
      p.clear(); q.clear(); r.clear();
    }
    virtual void free_CrsMatrix_data () {
      ir.clear(); jc.clear(); v.clear();
    }
  };

  // Return a random diagonal entry.
  static Sclr gen_diag () {
    Sclr diag = rand_scalar<Sclr>();
    diag += (real_part(diag) > 0 ? 1.0 : -1.0);
    return diag;
  }

  // Make a uniform [-0.5, 0.5) random vector.
  static void gen_rand_vector (const size_t n, std::vector<Real>& x) {
    x.resize(n);
    for (size_t i = 0; i < n; ++i)
      x[i] = urand() - 0.5;
  }

  template <typename T> static void
  gen_rand_vector (const size_t n, std::vector<std::complex<T> >& x) {
    x.resize(n);
    for (size_t i = 0; i < n; ++i)
      x[i] = rand_scalar<std::complex<T> >();
  }

  // Make a random permutation vector.
  static void gen_rand_perm (const size_t n, std::vector<Int>& p) {
    p.resize(n);
    for (size_t i = 0; i < n; ++i)
      p[i] = i;
    for (size_t i = 0; i < n; ++i) {
      const int j = urand()*n, k = urand()*n;
      std::swap(p[j], p[k]);
    }
  }

  // Make a random diagonal matrix.
  static void gen_diag_matrix (const TestOptions& to, Data& d) {
    for (int r = 0; r < to.n; ++r) {
      d.jc.push_back(r);
      d.v.push_back(gen_diag());
      d.ir[r+1] = d.ir[r] + 1;
    }
  }

  // y = a*x.
  static void mvp (const Data& a, const bool transp, const bool conj,
                   const Sclr* const x, Sclr* const y) {
    const Size* const ir = &a.ir[0];
    const Int* const jc = &a.jc[0];
    const Sclr* const d = &a.v[0];
    if (transp) {
      for (Int i = 0; i < a.m; ++i)
        y[i] = 0;
      for (Int i = 0; i < a.m; ++i) {
        const Sclr xi = x[i];
        const Size iri = ir[i], irip1 = ir[i+1];
        for (Size j = iri; j < irip1; ++j) {
          Sclr dj = d[j];
          if (conj) conjugate(dj);
          y[jc[j]] += dj*xi;
        }
      }
    } else {
      for (Int i = 0; i < a.m; ++i) {
        Sclr acc = 0;
        const Size iri = ir[i], irip1 = ir[i+1];
        for (Size j = iri; j < irip1; ++j) {
          Sclr dj = d[j];
          if (conj) conjugate(dj);
          acc += dj*x[jc[j]];
        }
        y[i] = acc;
      }
    }
  }

  template<typename T> static inline double square (const T& x) { return x*x; }
  template<typename T> static inline double square (const std::complex<T>& x)
  { return std::norm(x); }

  // norm(a - b, 2)/norm(a, 2).
  static double
  reldif (const Sclr* const a, const Sclr* const b, const Size n) {
    double num = 0, den = 0;
    for (Size i = 0; i < n; ++i) {
      num += square(a[i] - b[i]);
      den += square(a[i]);
    }
    return std::sqrt(num/den);
  }

  static bool is_not_tri (const typename TestOptions::MatrixType& mt) {
    return mt == TestOptions::not_tri || mt == TestOptions::not_tri_almost_diag;
  }

  // Remove or add an entry to make the matrix non-triangular. This will test
  // the shape-detection code.
  static void
  test_shape_checking (const TestOptions& to, Data& d, const int r) {
    if (to.matrix_type == TestOptions::missing_some_diag && r == to.n/2) {
      // Remove the previous diagonal entry.
      --d.ir[r+1];
      d.jc.pop_back();
      d.v.pop_back();
    } else if (is_not_tri(to.matrix_type) && r == to.n/2) {
      // Add an entry on the wrong tri side.
      d.jc.push_back(to.upper ? std::max(0, r-1) : std::min(to.n-1, r+1));
      d.v.push_back(0.5);
      ++d.ir[r+1];
    }
  }

  // Make a random sparse matrix with probability of an off-diag entry set to
  // p. If p is not specified or is negative, set it to min(1, 20/n).
  static void gen_tri_sparse_matrix (const TestOptions& to, Data& d, double p) {
    if (p < 0) p = std::min(0.4, 20.0/to.n);
    std::vector<Sclr> row_val;
    std::vector<Int> row_col;
    for (int r = 0; r < to.n; ++r) {
      Sclr diag;
      {
        // Generate the off-diag row and a diag. Scale the off-diag row to keep
        // the condition number reasonable.
        row_val.clear();
        row_col.clear();
        const int ntrial = to.upper ? to.n - r - 1 : r;
        const int base = to.upper ? r + 1 : 0;
        for (int i = 0; i < ntrial; ++i)
          if (urand() < p) {
            row_val.push_back(rand_scalar<Sclr>());
            row_col.push_back(static_cast<Int>(base + i));
          }
        // If testing shape checking, make sure there's at least one entry on
        // the correct side.
        if (is_not_tri(to.matrix_type) && to.n > 2 && r == to.n/2 &&
            row_val.empty()) {
          row_val.push_back(rand_scalar<Sclr>());
          row_col.push_back(static_cast<Int>(base + ntrial - 1));
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
        if (is_not_tri(to.matrix_type))
          test_shape_checking(to, d, r);
        ++d.ir[r+1];
        d.jc.push_back(r);
        d.v.push_back(diag);
        if (to.matrix_type == TestOptions::missing_some_diag)
          test_shape_checking(to, d, r);
        d.ir[r+1] += static_cast<Int>(row_val.size());
        d.jc.insert(d.jc.end(), row_col.begin(), row_col.end());
        d.v.insert(d.v.end(), row_val.begin(), row_val.end());
      } else {
        d.ir[r+1] += static_cast<Int>(row_val.size());
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
  static void gen_tri_dense_matrix (const TestOptions& to, Data& d) {
    gen_tri_sparse_matrix(to, d, 1);
  }

  static void gen_tri_block_matrix_from_tri_matrix (
    const Data& src, const Int blksz, const bool upper, Data& dst)
  {
    assert(blksz > 1);
    const Size nnz =
      blksz*blksz*src.ir.back() - ((blksz*(blksz-1))/2)*src.m // full blocks
      - src.ir.back(); // remove 1 from each block
    dst.m = blksz*src.m;
    dst.ir.resize(dst.m+1);
    dst.jc.resize(nnz);
    dst.v.resize(nnz);

    dst.ir[0] = 0;
    Size dk = 0;
    for (Int sr = 0, dr = 0; sr < src.m; ++sr) {
      for (Int br = 0; br < blksz; ++br, ++dr) {
        dst.ir[dr+1] = dst.ir[dr];
        const Real scale_odb = 1.0/(src.ir[sr+1] - src.ir[sr]), scale = 1.0/blksz;
        for (Size j = src.ir[sr]; j < src.ir[sr+1]; ++j) {
          const bool on_diag = ((   upper && j == src.ir[sr]) ||
                                ( ! upper && j == src.ir[sr+1] - 1));
          const Int sc = src.jc[j];
          const Sclr& val = src.v[j];
          const Int
            bcs = on_diag && upper ? br : 0,
            bce = on_diag ? (upper ? blksz : br+1) : blksz;
          for (Int bc = bcs; bc < bce; ++bc) {
            // Skip one entry in the block.
            if (on_diag && upper) {
              if (br == 0 && bc+1 == bce) continue;
            } else if (br+1 == blksz && bc == bcs) continue;

            ++dst.ir[dr+1];
            dst.jc[dk] = blksz*sc + bc;
            dst.v[dk] = val*2.0*(urand() - 0.5);
            // Scale the off-diag elements of the on-diag block.
            if (on_diag && bc != br)
              dst.v[dk] *= scale_odb;
            // Scale the off-diag elements of the tri.
            if ( ! (on_diag && bc == br))
              dst.v[dk] *= scale;
            ++dk;
          }
        }
      }
    }
    assert(dk == nnz);
  }

  static void gen_implicit_unit_diag_from_tri_matrix (
    const bool upper, const Data& src, Data& dst)
  {
    dst.m = src.m;
    dst.ir = src.ir;
    dst.jc.resize(src.jc.size() - dst.m);
    dst.v.resize(dst.jc.size());
    dst.ir[0] = 0;
    Size z = 0;
    for (Int r = 0; r < dst.m; ++r) {
      dst.ir[r+1] = src.ir[r+1] - (r + 1);
      const Sclr diag = upper ? src.v[src.ir[r]] : src.v[src.ir[r+1] - 1];
      Size
        K_beg = src.ir[r]   + (upper ? 1 : 0),
        K_end = src.ir[r+1] - (upper ? 0 : 1);
      for (Size K = K_beg; K < K_end; ++K) {
        dst.jc[z] = src.jc[K];
        dst.v[z] = src.v[K] / diag;
        ++z;
      }
    }
  }

  // Make a random triangular matrix according to the TestOptions.
  static void gen_tri_matrix (const TestOptions& to, Data& d) {
    d.clear();
    d.m = to.n;
    d.ir.resize(d.m+1, 0);
    switch (to.matrix_type) {
    case TestOptions::diag: gen_diag_matrix(to, d); break;
    case TestOptions::dense: gen_tri_dense_matrix(to, d); break;
    case TestOptions::sparse:
    case TestOptions::implicit_unit_diag:
    case TestOptions::missing_some_diag:
    case TestOptions::not_tri:
      gen_tri_sparse_matrix(to, d, to.density); break;
    case TestOptions::block_sparse:
    case TestOptions::block_sparse_implicit_unit_diag: {
      Data point_sparse(d);
      gen_tri_sparse_matrix(to, point_sparse, to.density);
      gen_tri_block_matrix_from_tri_matrix(point_sparse, to.block_size,
                                           to.upper, d);
    } break;
    case TestOptions::not_tri_almost_diag:
      gen_tri_sparse_matrix(to, d, 0); break;
    }
    if (to.has_unit_diag()) {
      Data src(d);
      d.clear();
      gen_implicit_unit_diag_from_tri_matrix(to.upper, src, d);
    }
    assert((std::size_t) d.ir.back() == d.v.size());
    assert((std::size_t) d.ir.back() == d.jc.size());
  }

  static void write_matrixmarket (const Data& T, const std::string& filename) {
    std::ofstream f(filename.c_str());
    if ( ! f.is_open()) return;
    f << "%%MatrixMarket matrix coordinate real general\n"
      << T.m << " " << T.m << " " << T.ir[T.m] << "\n";
    for (Int r = 0; r < T.m; ++r)
      for (Size j = T.ir[r]; j < T.ir[r+1]; ++j)
        f << r+1 << " " << T.jc[j] + 1 << " " << T.v[j] << "\n";
  }

  static void transpose (
    const std::vector<Size>& ir, const std::vector<Int>& jc,
    const std::vector<Sclr>& d, std::vector<Size>& irt,
    std::vector<Int>& jct, std::vector<Sclr>& dt)
  {
    const Int m = static_cast<Int>(ir.size()) - 1;
    // 1. Count the number of entries in each col.
    for (Int i = 0; i <= m; ++i) irt[i] = 0;
    const Size nnz = ir[m];
    for (Size k = 0; k < nnz; ++k) {
      // Store shifted up by 1. This is the trick that makes extra workspace
      // unnecessary.
      ++irt[jc[k]+1];
    }
    // 2. Cumsum to get the col pointers.
    Size sum = 0;
    for (Int i = 0; i < m; ++i) {
      const Size incr = irt[i+1];
      irt[i+1] = sum;
      sum += incr;
    }
    assert(sum == nnz);
    // At this point, At.ir[i+1] gives what is normally in At.ir[i].
    // 3. Fill in jc and d.
    for (Int r = 0; r < m; ++r)
      for (Size j = ir[r], jlim = ir[r+1]; j < jlim; ++j) {
        const Int c = jc[j];
        // Points to next entry to be filled.
        // 3a. Increment to next available entry.
        const Size p = irt[c+1]++;
        jct[p] = r;
        dt[p] = d[j];
      }
    assert(irt[m] == ir[m]);
    // 4. Step 3a made it such that At.ir[i] is now what we need except for
    // i = 0. At.ir[0] was already set to 0 in step 1, so we're done.
  }

  static void transpose (const Data& src, Data& dst) {
    dst.m = src.m;
    dst.ir.resize(src.m+1);
    dst.jc.resize(src.ir.back());
    dst.v.resize(src.ir.back());
    transpose(src.ir, src.jc, src.v, dst.ir, dst.jc, dst.v);
  }
};
} // namespace htsimpl
} // namespace Experimental

#endif
