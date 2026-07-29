//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <iostream>
#include <chrono>
#include <iomanip>

#include "KokkosGraph_Merge.hpp"
#include "Kokkos_Random.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_IOUtils.hpp"

struct Result {
  bool ok;
  double us;
  size_t bytes;
  Result() : ok(false) {}
};

struct Stats {
  size_t numRows;
  size_t nnz;
  float nnzPerRow;
  float nnzPerRowStddev;
};

struct MatrixResult {
  Stats stats;
  Result spmv[2];
  Result spmvmv[2][2][2];
};

template <typename Crs>
float nnz_per_row_stddev(const Crs &crs) {
  auto rmh = Kokkos::create_mirror_view(crs.graph.row_map);
  Kokkos::deep_copy(rmh, crs.graph.row_map);

  float mean = 0;
  for (size_t i = 0; i < size_t(crs.numRows()); ++i) {
    mean += rmh(i + 1) - rmh(i);
  }
  mean /= crs.numRows();

  float var = 0;
  for (size_t i = 0; i < size_t(crs.numRows()); ++i) {
    float v = rmh(i + 1) - rmh(i) - mean;
    var += v * v;
  }
  var /= crs.numRows();
  return std::sqrt(var);
}

template <typename Crs>
Stats get_stats(const Crs &crs) {
  Stats stats;
  stats.numRows         = crs.numRows();
  stats.nnz             = crs.nnz();
  stats.nnzPerRow       = float(stats.nnz) / stats.numRows;
  stats.nnzPerRowStddev = nnz_per_row_stddev(crs);
  return stats;
}

template <bool CONJ, typename Matrix>
Result bench_spmv(const Matrix &crs, int nWarmup, int nIters) {
  using exec_space = typename Matrix::execution_space;
  typedef typename Matrix::memory_space mem_space;
  using Scalar = typename Matrix::non_const_value_type;

  typedef Kokkos::View<Scalar *, mem_space> view_t;

  using Spmv = KokkosSparse::Impl::SpmvMergeHierarchical<exec_space, Matrix, view_t, view_t>;

  Result res;

  view_t y_init("y_init", crs.numRows());
  view_t y_exp("y", crs.numRows());
  view_t y_act("y", crs.numRows());
  view_t x("x", crs.numCols());

  res.bytes = crs.nnz() * sizeof(Scalar)                                       // A values
              + crs.nnz() * sizeof(typename Matrix::ordinal_type)              // A cols
              + crs.graph.row_map.size() * sizeof(typename Matrix::size_type)  // A row-map
              + 2 * y_act.size() * sizeof(Scalar)                              // load / store y
              + x.size() * sizeof(Scalar)                                      // load x
      ;

  Kokkos::Random_XorShift64_Pool<exec_space> random_pool(12345);
  fill_random(y_init, random_pool, 0.0, 1.0);

#if 0  // SIMPLE CASE, output should be number of non-zeros in each row
  Kokkos::parallel_for("x", 
    Kokkos::RangePolicy<exec_space>(0, x.size()), 
    KOKKOS_LAMBDA(const int &i) {
      x(i) = 1.0;
    }
  );
  Scalar alpha = 1;
  Scalar beta = 0;
#else
  fill_random(x, random_pool, 0.0, 1.0);
  Scalar alpha = -1;
  Scalar beta  = 0.273;
#endif

  Kokkos::deep_copy(y_act, y_init);
  Kokkos::deep_copy(y_exp, y_init);

  const char *mode = CONJ ? KokkosSparse::Conjugate : KokkosSparse::NoTranspose;

  // test non-transpose
  KokkosSparse::spmv(mode, alpha, crs, x, beta, y_exp);
  Kokkos::fence();
  Spmv::spmv(exec_space(), mode, alpha, crs, x, beta, y_act);
  Kokkos::fence();

  {
    using AT       = KokkosKernels::ArithTraits<Scalar>;
    using mag_type = typename AT::mag_type;
    using ATM      = KokkosKernels::ArithTraits<mag_type>;

    // max value in A
    mag_type maxA = 0;
    Kokkos::parallel_reduce(
        "maxA", Kokkos::RangePolicy<exec_space>(0, crs.nnz()),
        KOKKOS_LAMBDA(const int &i, mag_type &lmax) {
          mag_type v = AT::abs(crs.values(i));
          lmax       = lmax > v ? lmax : v;
        },
        maxA);

    double eps           = AT::epsilon();
    const double max_val = AT::abs(beta * 1.0 + crs.numCols() * alpha * maxA * 1.0);

    auto h_exp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y_exp);
    auto h_act = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y_act);

#define DIE_ON_ERROR 0
#define CHECK_CORRECTNESS_FIRST 1
#if CHECK_CORRECTNESS_FIRST
    size_t err = 0;
    std::vector<size_t> errIdx;
    for (size_t i = 0; i < h_exp.size(); ++i) {
      const mag_type error = ATM::abs(h_exp(i) - h_act(i));
      if (error > eps * max_val) {
        ++err;
        errIdx.push_back(i);
      }
    }
    if (err > 0) {
      size_t errLimit = 100;
      std::cerr << "first " << errLimit << " errors...\n";
      std::cerr << "i\texp\tact" << std::endl;
      std::cerr << "-\t---\t---" << std::endl;
      for (auto i : errIdx) {
        std::cerr << i << "\t" << h_exp(i) << "\t" << h_act(i) << std::endl;
        if (0 == --errLimit) {
          break;
        }
      }
      std::cerr << __FILE__ << ":" << __LINE__ << ": ERROR: correctness failed " << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << ": threshold was " << eps * max_val << std::endl;
#if DIE_ON_ERROR
      exit(EXIT_FAILURE);
#else
      res.ok = false;
      return res;
#endif
    }

#endif
#undef DIE_ON_ERROR
#undef CHECK_CORRECTNESS_FIRST
  }

  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::duration<double> Duration;

  exec_space space;

  for (int i = 0; i < nWarmup; ++i) {
    Spmv::spmv(exec_space(), mode, alpha, crs, x, beta, y_exp);
  }
  Kokkos::fence();
  auto start = Clock::now();
  for (int i = 0; i < nIters; ++i) {
    Spmv::spmv(exec_space(), mode, alpha, crs, x, beta, y_exp);
  }
  Kokkos::fence();
  Duration elapsed(Clock::now() - start);

  res.us = elapsed.count() * 1e6 / nIters;
  res.ok = true;
  return res;
}

/*! \brief bench multivector spmv

    \tparam Layout of X multivector
    \tparam Layout of Y multivector

    \param k number of vectors
*/
template <bool CONJ, typename XLayout, typename YLayout, typename Matrix>
Result bench_spmvmv(const Matrix &crs, int nWarmup, int nIters, int k) {
  typedef typename Matrix::execution_space exec_space;
  typedef typename Matrix::memory_space mem_space;
  using Scalar = typename Matrix::non_const_value_type;

  typedef Kokkos::View<Scalar **, XLayout, mem_space> x_view_t;
  typedef Kokkos::View<Scalar **, YLayout, mem_space> y_view_t;

  using Spmv = KokkosSparse::Impl::SpmvMvMergeHierarchical<Matrix, x_view_t, y_view_t>;

  Result res;

  y_view_t y_init("y_init", crs.numRows(), k);
  y_view_t y_exp("y", crs.numRows(), k);
  y_view_t y_act("y", crs.numRows(), k);
  x_view_t x("x", crs.numCols(), k);

  res.bytes = crs.nnz() * sizeof(Scalar)                                       // A values
              + crs.nnz() * sizeof(typename Matrix::ordinal_type)              // A cols
              + crs.graph.row_map.size() * sizeof(typename Matrix::size_type)  // A row-map
              + 2 * y_act.extent(0) * y_act.extent(1) * sizeof(Scalar)         // load / store y
              + x.extent(0) * x.extent(1) * sizeof(Scalar)                     // load x
      ;

  Kokkos::Random_XorShift64_Pool<exec_space> random_pool(12345);
  fill_random(y_init, random_pool, 0.0, 1.0);

#if 0  // SIMPLE CASE, output should be number of non-zeros in each row
  Kokkos::parallel_for("x", 
    Kokkos::MDRangePolicy<exec_space, Kokkos::Rank<2>>({0, 0}, {x.extent(0), x.extent(1)}), 
    KOKKOS_LAMBDA(const int &i, const int &j) {
      x(i,j) = 1.0;
    }
  );
  Scalar alpha = 1;
  Scalar beta = 0;
#else
  fill_random(x, random_pool, 0.0, 1.0);
  Scalar alpha = -1;
  Scalar beta  = 0.273;
#endif

  Kokkos::deep_copy(y_act, y_init);
  Kokkos::deep_copy(y_exp, y_init);

  const char *mode = CONJ ? KokkosSparse::Conjugate : KokkosSparse::NoTranspose;

  // test non-transpose
  KokkosSparse::spmv(mode, alpha, crs, x, beta, y_exp);
  Kokkos::fence();
  Spmv::spmv(mode, alpha, crs, x, beta, y_act);
  Kokkos::fence();

  {
    using AT       = KokkosKernels::ArithTraits<Scalar>;
    using mag_type = typename AT::mag_type;
    using ATM      = KokkosKernels::ArithTraits<mag_type>;

    // max value in A
    mag_type maxA = 0;
    Kokkos::parallel_reduce(
        "maxA", Kokkos::RangePolicy<exec_space>(0, crs.nnz()),
        KOKKOS_LAMBDA(const int &i, mag_type &lmax) {
          mag_type v = AT::abs(crs.values(i));
          lmax       = lmax > v ? lmax : v;
        },
        maxA);

    double eps           = AT::epsilon();
    const double max_val = AT::abs(beta * 1.0 + crs.numCols() * alpha * maxA * 1.0);

    auto h_exp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y_exp);
    auto h_act = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y_act);

#define DIE_ON_ERROR 0
#define CHECK_CORRECTNESS_FIRST 1
#if CHECK_CORRECTNESS_FIRST
    size_t err = 0;
    std::vector<std::pair<size_t, size_t>> errIdx;
    for (size_t i = 0; i < h_exp.extent(0); ++i) {
      for (size_t j = 0; j < h_exp.extent(1); ++j) {
        const mag_type error = ATM::abs(h_exp(i, j) - h_act(i, j));
        if (error > eps * max_val) {
          ++err;
          errIdx.push_back(std::make_pair(i, j));
        }
      }
    }
    if (err > 0) {
      size_t errLimit = 100;
      std::cerr << "first " << errLimit << " errors...\n";
      std::cerr << "i,j\texp\tact" << std::endl;
      std::cerr << "-\t---\t---" << std::endl;
      for (auto p : errIdx) {
        const auto i = p.first;
        const auto j = p.second;
        std::cerr << i << "," << j << "\t" << h_exp(i, j) << "\t" << h_act(i, j) << std::endl;
        if (0 == --errLimit) {
          break;
        }
      }
      std::cerr << __FILE__ << ":" << __LINE__ << ": ERROR: correctness failed" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << ": threshold was " << eps * max_val << std::endl;
#if DIE_ON_ERROR
      exit(EXIT_FAILURE);
#else
      res.ok = false;
      return res;
#endif
    }
#endif
#undef DIE_ON_ERROR
#undef CHECK_CORRECTNESS_FIRST
  }

  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::duration<double> Duration;

  exec_space space;

  for (int i = 0; i < nWarmup; ++i) {
    Spmv::spmv(mode, alpha, crs, x, beta, y_exp);
  }
  Kokkos::fence();
  auto start = Clock::now();
  for (int i = 0; i < nIters; ++i) {
    Spmv::spmv(mode, alpha, crs, x, beta, y_exp);
  }
  Kokkos::fence();
  Duration elapsed(Clock::now() - start);

  res.us = elapsed.count() * 1e6 / nIters;
  res.ok = true;
  return res;
}

template <typename T, typename Device>
MatrixResult read_and_bench(const std::string &path, int nWarmup, int nIters) {
  typedef int Ordinal;
  typedef KokkosSparse::CrsMatrix<T, Ordinal, Device> Matrix;
  Matrix crs = KokkosSparse::Impl::read_kokkos_crst_matrix<Matrix>(path.c_str());

  MatrixResult mr;

  mr.stats = get_stats(crs);

  for (bool conj : {true, false}) {
    if (conj) {
      mr.spmv[conj] = bench_spmv<true>(crs, nWarmup, nIters);
    } else {
      mr.spmv[conj] = bench_spmv<false>(crs, nWarmup, nIters);
    }

    for (bool xl : {true, /*, false*/}) {
      for (bool yl : {true /*, false*/}) {
        if (conj && xl && yl) {
          mr.spmvmv[conj][xl][yl] = bench_spmvmv<true, Kokkos::LayoutLeft, Kokkos::LayoutLeft>(crs, nWarmup, nIters, 4);
        } else if (conj && xl && !yl) {
          mr.spmvmv[conj][xl][yl] =
              bench_spmvmv<true, Kokkos::LayoutLeft, Kokkos::LayoutRight>(crs, nWarmup, nIters, 4);
        } else if (conj && !xl && yl) {
          mr.spmvmv[conj][xl][yl] =
              bench_spmvmv<true, Kokkos::LayoutRight, Kokkos::LayoutLeft>(crs, nWarmup, nIters, 4);
        } else if (conj && !xl && !yl) {
          mr.spmvmv[conj][xl][yl] =
              bench_spmvmv<true, Kokkos::LayoutRight, Kokkos::LayoutRight>(crs, nWarmup, nIters, 4);
        } else if (!conj && xl && yl) {
          mr.spmvmv[conj][xl][yl] =
              bench_spmvmv<false, Kokkos::LayoutLeft, Kokkos::LayoutLeft>(crs, nWarmup, nIters, 4);
        } else if (!conj && xl && !yl) {
          mr.spmvmv[conj][xl][yl] =
              bench_spmvmv<false, Kokkos::LayoutLeft, Kokkos::LayoutRight>(crs, nWarmup, nIters, 4);
        } else if (!conj && !xl && yl) {
          mr.spmvmv[conj][xl][yl] =
              bench_spmvmv<false, Kokkos::LayoutRight, Kokkos::LayoutLeft>(crs, nWarmup, nIters, 4);
        } else if (!conj && !xl && !yl) {
          mr.spmvmv[conj][xl][yl] =
              bench_spmvmv<false, Kokkos::LayoutRight, Kokkos::LayoutRight>(crs, nWarmup, nIters, 4);
        }
      }
    }
  }
  return mr;
}

void usage(char **argv) { std::cerr << argv[0] << " matrix.mtx [... matrix.mtx]" << std::endl; }

std::string header() {
  std::string s("path,rows,nnz,nnz/row,stddev nnz/row");

  for (bool conj : {true, false}) {
    std::string pfx = conj ? "spmv[c]" : "spmv[n]";
    s += "," + pfx + " us";
    s += "," + pfx + " GFLOPS";
    s += "," + pfx + " GB/s";

    for (bool xl : {true, false}) {
      for (bool yl : {true, false}) {
        pfx = "spmvmv";
        pfx += "[";
        pfx += conj ? "c" : "n";
        pfx += "|";
        pfx += xl ? "xl" : "xr";
        pfx += "|";
        pfx += yl ? "yl" : "yr";
        pfx += "]";
        s += "," + pfx + " us";
        s += "," + pfx + " GFLOPS";
        s += "," + pfx + " GB/s";
      }
    }
  }
  return s;
}

int main(int argc, char **argv) {
  Kokkos::initialize();

  if (argc < 2) {
    usage(argv);
    Kokkos::finalize();
    return EXIT_FAILURE;
  }

#define DO_BENCH(DEVICE, TYPE)                                            \
  {                                                                       \
    for (int i = 1; i < argc; ++i) {                                      \
      MatrixResult mr = read_and_bench<TYPE, DEVICE>(argv[i], 10, 500);   \
      if (1 == i) std::cout << header() << std::endl;                     \
      std::cout << argv[i];                                               \
      std::cout << "," << mr.stats.numRows;                               \
      std::cout << "," << mr.stats.nnz;                                   \
      std::cout << "," << mr.stats.nnzPerRow;                             \
      std::cout << "," << mr.stats.nnzPerRowStddev;                       \
      for (bool conj : {true, false}) {                                   \
        const Result &r = mr.spmv[conj];                                  \
        if (r.ok) {                                                       \
          std::cout << "," << r.us;                                       \
          std::cout << "," << mr.stats.nnz * 2 / (r.us / 1e6) / 1e9;      \
          std::cout << "," << r.bytes / (r.us / 1e6) / 1e9;               \
        } else {                                                          \
          std::cout << ",,,";                                             \
        }                                                                 \
        for (bool xl : {true, false}) {                                   \
          for (bool yl : {true, false}) {                                 \
            const Result &r2 = mr.spmvmv[conj][xl][yl];                   \
            if (r2.ok) {                                                  \
              std::cout << "," << r2.us;                                  \
              std::cout << "," << mr.stats.nnz * 2 / (r2.us / 1e6) / 1e9; \
              std::cout << "," << r2.bytes / (r2.us / 1e6) / 1e9;         \
            } else {                                                      \
              std::cout << ",,,";                                         \
            }                                                             \
          }                                                               \
        }                                                                 \
      }                                                                   \
      std::cout << std::endl;                                             \
    }                                                                     \
  }

#ifdef KOKKOS_ENABLE_CUDA
  std::cout << "CUDA SpMV (double) c=" << Kokkos::Cuda().concurrency() << std::endl;
  DO_BENCH(Kokkos::Cuda, double);

  std::cout << "CUDA SpMV (float) c=" << Kokkos::Cuda().concurrency() << std::endl;
  DO_BENCH(Kokkos::Cuda, float);
#endif

#ifdef KOKKOS_ENABLE_OPENMP
  std::cout << "OpenMP SpMV (double) c=" << Kokkos::OpenMP().concurrency() << std::endl;
  DO_BENCH(Kokkos::OpenMP, double);

  std::cout << "OpenMP SpMV (float) c=" << Kokkos::OpenMP().concurrency() << std::endl;
  DO_BENCH(Kokkos::OpenMP, float);
#endif

#ifdef KOKKOS_ENABLE_HIP
  std::cout << "HIP SpMV (double) c=" << Kokkos::HIP().concurrency() << std::endl;
  DO_BENCH(Kokkos::HIP, double);

  std::cout << "HIP SpMV (float) c=" << Kokkos::HIP().concurrency() << std::endl;
  DO_BENCH(Kokkos::HIP, float);
#endif

  Kokkos::finalize();

  return 0;
}
