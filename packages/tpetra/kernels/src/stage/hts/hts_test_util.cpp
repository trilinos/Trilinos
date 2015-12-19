#include "hts_test_util.hpp"

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

namespace Experimental {
namespace hts {
namespace util {

std::ostream& operator<< (std::ostream& os, const TestOptions& to) {
  os << "[n " << to.n << " reprocess " << to.reprocess << " upper " << to.upper
     << " matrix_type " << to.matrix_type << " transpose " << to.transpose
     << " solve_type " << to.solve_type << " nrhs " << to.nrhs
     << " nthreads " << to.nthreads << "]";
  return os;
}

double urand () {
#if 0 // Not all compilers have this, it seems.
  static std::default_random_engine generator;
  static std::uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
#else
  return rand() / ((double) RAND_MAX + 1.0);
#endif
}

inline double gen_diag () {
  double diag = urand() - 0.5;
  diag += (diag > 0 ? 1 : -1);
  return diag;
}

void gen_rand_vector (const size_t n, std::vector<hts::Real>& x) {
  x.resize(n);
  for (size_t i = 0; i < n; ++i)
    x[i] = urand() - 0.5;
}

void gen_rand_perm (const size_t n, std::vector<hts::Int>& p) {
  p.resize(n);
  for (size_t i = 0; i < n; ++i)
    p[i] = i;
  for (size_t i = 0; i < n; ++i) {
    const int j = urand()*n, k = urand()*n;
    std::swap(p[j], p[k]);
  }
}

void gen_diag_matrix (const TestOptions& to, Data& d) {
  for (int r = 0; r < to.n; ++r) {
    d.jc.push_back(r);
    d.v.push_back(gen_diag());
    d.ir[r+1] = d.ir[r] + 1;
  }
}

void mvp (const Data& a, const bool transp, const hts::Real* const x,
          hts::Real* const y) {
  using hts::Size;
  using hts::Int;
  using hts::Real;
  const Size* const ir = &a.ir[0];
  const Int* const jc = &a.jc[0];
  const Real* const d = &a.v[0];
  if (transp) {
    for (Int i = 0; i < a.m; ++i)
      y[i] = 0;
    for (Int i = 0; i < a.m; ++i) {
      const Real xi = x[i];
      const Size iri = ir[i], irip1 = ir[i+1];
      for (Size j = iri; j < irip1; ++j)
        y[jc[j]] += d[j]*xi;
    }
  } else {
    for (Int i = 0; i < a.m; ++i) {
      Real acc = 0;
      const Size iri = ir[i], irip1 = ir[i+1];
      for (Size j = iri; j < irip1; ++j)
        acc += d[j]*x[jc[j]];
      y[i] = acc;
    }
  }
}

template<typename T> inline double square (const T& x) { return x*x; }

// norm(a - b, 2)/norm(a, 2).
double reldif (const Real* const a, const Real* const b, const Size n) {
  double num = 0, den = 0;
  for (Size i = 0; i < n; ++i) {
    num += square(a[i] - b[i]);
    den += square(a[i]);
  }
  return std::sqrt(num/den);
}

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
void gen_tri_sparse_matrix (const TestOptions& to, Data& d, double p) {
  if (p < 0) p = std::min(0.4, 20.0/to.n);
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
      if (to.matrix_type == TestOptions::not_tri)
        test_shape_checking(to, d, r);
      ++d.ir[r+1];
      d.jc.push_back(r);
      d.v.push_back(diag);
      if (to.matrix_type == TestOptions::missing_diag)
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

void gen_tri_block_matrix_from_tri_matrix (
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
        const Real val = src.v[j];
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
          dst.v[dk] = val*2*(urand() - 0.5);
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
    gen_tri_sparse_matrix(to, d, to.density); break;
  case TestOptions::block_sparse:
    Data point_sparse(d);
    gen_tri_sparse_matrix(to, point_sparse, to.density);
    gen_tri_block_matrix_from_tri_matrix(point_sparse, to.block_size,
                                         to.upper, d);
    break;
  }
  assert((std::size_t) d.ir.back() == d.v.size());
  assert((std::size_t) d.ir.back() == d.jc.size());
}

void write_matrixmarket (const Data& T, const std::string& filename) {
  FILE* fid = fopen(filename.c_str(), "w");
  if ( ! fid) return;
  fprintf(fid, "%%%%MatrixMarket matrix coordinate real general\n"
          "%11d %11d %11d\n", (int) T.m, (int) T.m, (int) T.ir[T.m]);
  for (hts::Int r = 0; r < T.m; ++r)
    for (hts::Size j = T.ir[r]; j < T.ir[r+1]; ++j)
      fprintf(fid, "%d %d %1.15e\n", r+1, T.jc[j] + 1, T.v[j]);
  fclose(fid);
}

std::ostream& operator<< (std::ostream& os, const Data& d) {
  os << "m: " << d.m << " nnz: " << d.ir.back() << "\n";
  for (Int r = 0; r < d.m; ++r) {
    Int pc = -1;
    for (Size j = d.ir[r]; j < d.ir[r+1]; ++j) {
      for (Int c = pc; c < d.jc[j]-1; ++c)
        os << " ";
      pc = d.jc[j];
      os << "*";
    }
    os << "\n";
  }
  return os;
}

void transpose (
  const std::vector<Size>& ir, const std::vector<Int>& jc,
  const std::vector<Real>& d, std::vector<Size>& irt,
  std::vector<Int>& jct, std::vector<Real>& dt)
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

void transpose (const Data& src, Data& dst) {
  dst.m = src.m;
  dst.ir.resize(src.m+1);
  dst.jc.resize(src.ir.back());
  dst.v.resize(src.ir.back());
  transpose(src.ir, src.jc, src.v, dst.ir, dst.jc, dst.v);
}

}}}
