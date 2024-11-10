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

/*! \file KokkosSparse_spmv_bsr_benchmark.cpp

    Read a matrix market file, choose a block size and a number of multivectors,
   and compare Bsr SpMV implementations
*/

#include <typeindex>
#include <memory>

#include <Kokkos_Core.hpp>

/* Some versions of clang that hipcc is basedoff of haven't stabilized
 * std::filesystem yet */
#if defined(KOKKOS_ENABLE_HIP) && __HIPCC__
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#include <rocsparse/rocsparse.h>
#endif

#include <benchmark/benchmark.h>

#include <Benchmark_Context.hpp>

#include <Kokkos_ArithTraits.hpp>

#include "Benchmark_Utils.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_IOUtils.hpp"
#include "KokkosSparse_spmv.hpp"
#include "KokkosSparse_crs_to_bsr_impl.hpp"
#include "KokkosSparse_crs_detect_block_size.hpp"

using namespace KokkosKernelsBenchmark;

/*  Since benchmarks have to be defined before they are executed, the file IO
   for each benchmark needs to be in the execution itself, otherwise every
   matrix would have to be resident in memory before any benchmark can run.

    If multiple benchmarks need the same file, it would be read over and over
   again. This is especially painful on network file systems, so this executable
   has a global cache to store the most recently-read matrix.

    Despite that the matrix is always read with the same precision, we don't
   know the Device at this time, so we can't define the value type of the cache
   yet. Instead, we'll erase the type, and use a pointer to void. The cache will
   be keyed on a combination of the path and the requested type, so we know if
   the actual CrsMatrix behind the void pointer matches the requested type or
   not
*/
using Key = std::tuple<fs::path, std::type_index>;
using Val = std::shared_ptr<void>;  // type-erased Crs matrix (since we don't
                                    // know the template params)
static Key CACHE_KEY = {"", std::type_index(typeid(void))};
static Val CACHE_VAL = nullptr;

// This can be called before Kokkos::finalize to kill the matrix that is living
// in the cache
void drop_cache() {
  CACHE_KEY = {"", std::type_index(typeid(void))};
  CACHE_VAL = nullptr;
}

/// cache repeated reads to \c path
template <typename Crs>
Crs cached_read(const fs::path &path) {
  // check if the cached matrix is a Crs from path
  const Key key(path, std::type_index(typeid(Crs)));

  // if this is not the cached matrix, overwrite the cache
  if (CACHE_KEY != key) {
    CACHE_KEY = key;
    CACHE_VAL = std::make_shared<Crs>(KokkosSparse::Impl::read_kokkos_crst_matrix<Crs>(path.c_str()));
  }

  // the Crs type is part of the key, so we know this cast is safe
  return *std::static_pointer_cast<Crs>(CACHE_VAL);
}

/* Cache a map of path -> matrix block size so that scanning the matrix to
 * register the benchmark and then actually running the becnchmark don't both
 * need to run the matrix */
template <typename Device>
size_t detect_block_size(const fs::path &path) {
  using ReadScalar  = double;
  using ReadOrdinal = int64_t;
  using ReadOffset  = uint64_t;
  using Crs         = KokkosSparse::CrsMatrix<ReadScalar, ReadOrdinal, Device, void, ReadOffset>;

  static std::map<fs::path, size_t> cache;

  if (0 == cache.count(path)) {
    std::cerr << "read " << path << "...\n";
    const Crs crs       = cached_read<Crs>(path);
    size_t detectedSize = KokkosSparse::Impl::detect_block_size(crs);
    std::cerr << "detected block size = " << detectedSize << "\n";
    cache[path] = detectedSize;
  }
  return cache.at(path);
}

// Test that y_act is close to y_exp.
// This needs the matrix, alpha, and beta to compute the error tolerance
// properly
template <typename View, typename Matrix, typename Alpha, typename Beta>
void check_correctness(benchmark::State &state, const View &y_exp, const View &y_act, const Matrix &crs,
                       const Alpha &alpha, const Beta &beta, const DieOnError &die, const SkipOnError &skip) {
  using execution_space = typename View::execution_space;
  using scalar_type     = typename View::non_const_value_type;
  using AT              = Kokkos::ArithTraits<scalar_type>;
  using mag_type        = typename AT::mag_type;
  using ATM             = Kokkos::ArithTraits<mag_type>;

  // max value in A
  mag_type maxA = 0;
  Kokkos::parallel_reduce(
      "maxA", Kokkos::RangePolicy<execution_space>(0, crs.nnz()),
      KOKKOS_LAMBDA(const int &i, mag_type &lmax) {
        mag_type v = AT::abs(crs.values(i));
        lmax       = lmax > v ? lmax : v;
      },
      maxA);

  double eps           = AT::epsilon();
  const double max_val = AT::abs(beta * 1.0 + crs.numCols() * alpha * maxA * 1.0);

  auto h_exp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y_exp);
  auto h_act = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), y_act);

  size_t err = 0;
  std::vector<std::pair<size_t, size_t>> errIdx;
  for (size_t i = 0; i < h_exp.extent(0); ++i) {
    for (size_t k = 0; k < h_exp.extent(1); ++k) {
      const mag_type error = ATM::abs(h_exp(i, k) - h_act(i, k));
      if (error > eps * max_val) {
        ++err;
        errIdx.push_back({i, k});
      }
    }
  }
  if (err > 0) {
    size_t errLimit = 100;  // how many errors to print
    std::cerr << "first " << errLimit << " errors...\n";
    std::cerr << "i\tk\texp\tact" << std::endl;
    std::cerr << "-\t-\t---\t---" << std::endl;
    for (auto [i, k] : errIdx) {
      std::cerr << i << "\t" << k << "\t" << h_exp(i, k) << "\t" << h_act(i, k) << std::endl;
      if (0 == --errLimit) {
        break;
      }
    }
    std::cerr << __FILE__ << ":" << __LINE__ << ": ERROR: correctness failed " << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << ": threshold was " << eps * max_val << std::endl;

    if (die) {
      exit(EXIT_FAILURE);
    } else if (skip) {
      state.SkipWithError("correctness check failed");
    }
  }
}

// Wrapper to create a common interface for all SpMVs to benchmark
struct SpmvDefault {
  template <typename Alpha, typename Matrix, typename XView, typename Beta, typename YView>
  static void spmv(const char *mode, const Alpha &alpha, const Matrix &crs, const XView &x, const Beta &beta,
                   const YView &y) {
    return KokkosSparse::spmv(mode, alpha, crs, x, beta, y);
  }

  static std::string name() { return "default"; }
};

// Wrapper to create a common interface for all SpMVs to benchmark
struct SpmvNative {
  template <typename Alpha, typename Matrix, typename XView, typename Beta, typename YView>
  static void spmv(const char *mode, const Alpha &alpha, const Matrix &crs, const XView &x, const Beta &beta,
                   const YView &y) {
    KokkosSparse::SPMVHandle<typename Matrix::execution_space, Matrix, XView, YView> handle(KokkosSparse::SPMV_NATIVE);
    return KokkosSparse::spmv(&handle, mode, alpha, crs, x, beta, y);
  }

  static std::string name() { return "native"; }
};

// Wrapper to create a common interface for all SpMVs to benchmark
struct SpmvV41 {
  template <typename Alpha, typename Matrix, typename XView, typename Beta, typename YView>
  static void spmv(const char *mode, const Alpha &alpha, const Matrix &crs, const XView &x, const Beta &beta,
                   const YView &y) {
    KokkosSparse::SPMVHandle<typename Matrix::execution_space, Matrix, XView, YView> handle(KokkosSparse::SPMV_BSR_V41);
    return KokkosSparse::spmv(&handle, mode, alpha, crs, x, beta, y);
  }

  static std::string name() { return "v4.1"; }
};

template <typename Spmv, typename Bsr>
void run(benchmark::State &state, const Bsr &bsr, const size_t k) {
  using execution_space = typename Bsr::execution_space;
  using memory_space    = typename Bsr::memory_space;
  using scalar_type     = typename Bsr::non_const_value_type;
  using ordinal_type    = typename Bsr::non_const_ordinal_type;
  using size_type       = typename Bsr::non_const_size_type;

  // multivector should be layoutleft for CPU, makes
  // slices of a single vector contiguous
  using view_t = Kokkos::View<scalar_type **, Kokkos::LayoutLeft, memory_space>;

  state.counters["nnz"]        = bsr.nnz() * bsr.blockDim() * bsr.blockDim();
  state.counters["num_rows"]   = bsr.numRows() * bsr.blockDim();
  state.counters["block_size"] = bsr.blockDim();
  state.counters["num_vecs"]   = k;

  view_t y_init("y_init", bsr.numRows() * bsr.blockDim(), k);
  view_t y_exp("ye", bsr.numRows() * bsr.blockDim(), k);
  view_t y_act("ya", bsr.numRows() * bsr.blockDim(), k);
  view_t x("x", bsr.numCols() * bsr.blockDim(), k);

  Kokkos::Random_XorShift64_Pool<execution_space> random_pool(12345);
  fill_random(y_init, random_pool, 0.0, 1.0);
  fill_random(x, random_pool, 0.0, 1.0);
  scalar_type alpha{1.17};
  scalar_type beta{-0.3};

  Kokkos::deep_copy(y_act, y_init);
  Kokkos::deep_copy(y_exp, y_init);

  const char *mode = KokkosSparse::NoTranspose;

  // test the SpMV against whatever the default is
  Spmv::spmv(mode, alpha, bsr, x, beta, y_act);
  Kokkos::fence();
  KokkosSparse::spmv(mode, alpha, bsr, x, beta, y_exp);
  Kokkos::fence();

  check_correctness(state, y_exp, y_act, bsr, alpha, beta, DieOnError(false), SkipOnError(true));

  Kokkos::fence();
  for (auto _ : state) {
    Spmv::spmv(mode, alpha, bsr, x, beta, y_exp);
    Kokkos::fence();
  }

  const size_t bytesPerSpmv = bsr.nnz() * bsr.blockDim() * bsr.blockDim() * sizeof(scalar_type)  // A values
                              + bsr.nnz() * sizeof(ordinal_type)                                 // A col indices
                              + (bsr.numRows() + 1) * sizeof(size_type)                          // A row-map
                              + 2 * bsr.numRows() * bsr.blockDim() * k * sizeof(scalar_type)     // load / store y
                              + bsr.numCols() * bsr.blockDim() * k * sizeof(scalar_type)         // load x
      ;

  state.SetBytesProcessed(bytesPerSpmv * state.iterations());
}

template <typename Bsr, typename Spmv>
void read_expand_run(benchmark::State &state, const fs::path &path, const size_t blockSize, const size_t k) {
  using scalar_type  = typename Bsr::non_const_value_type;
  using ordinal_type = typename Bsr::non_const_ordinal_type;

  // read Crs into host memory
  using Crs = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, Kokkos::HostSpace>;

  const Crs crs = cached_read<Crs>(path);
  Bsr bsr;
  try {
    bsr = KokkosSparse::Impl::expand_crs_to_bsr<Bsr>(crs, blockSize);
  } catch (std::exception &e) {
    state.SkipWithError(e.what());
    return;
  }

  run<Spmv>(state, bsr, k);
}

template <typename Bsr, typename Spmv>
void read_convert_run(benchmark::State &state, const fs::path &path, const size_t blockSize, const size_t k) {
  using scalar_type  = typename Bsr::non_const_value_type;
  using ordinal_type = typename Bsr::non_const_ordinal_type;

  using Crs = KokkosSparse::CrsMatrix<scalar_type, ordinal_type, Kokkos::HostSpace>;

  const Crs crs = cached_read<Crs>(path);
  Bsr bsr;
  try {
    bsr = KokkosSparse::Impl::blocked_crs_to_bsr<Bsr>(crs, blockSize);
  } catch (std::exception &e) {
    state.SkipWithError(e.what());
    return;
  }

  run<Spmv>(state, bsr, k);
}

template <typename Ordinal, typename Scalar, typename Offset, typename Device, typename Spmv>
void register_expand_type(const fs::path &path) {
  using Bsr              = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device, void, Offset>;
  std::vector<size_t> ks = {1, 3};
  for (size_t bs : {4, 7, 10, 16}) {  // block sizes
    for (size_t k : ks) {             // multivector sizes
      std::string name = std::string("MatrixMarketExpanded") + "/" + std::string(path.stem()) + "/" +
                         Kokkos::ArithTraits<Scalar>::name() + "/" + Kokkos::ArithTraits<Ordinal>::name() + "/" +
                         Kokkos::ArithTraits<Offset>::name() + "/" + std::to_string(bs) + "/" + std::to_string(k) +
                         "/" + Spmv::name() + "/" + Device::name();
      benchmark::RegisterBenchmark(name.c_str(), read_expand_run<Bsr, Spmv>, path, bs, k)->UseRealTime();
    }
  }
}

template <typename Ordinal, typename Scalar, typename Offset, typename Device, typename Spmv>
void register_convert_type(const fs::path &path, size_t bs) {
  using Bsr              = KokkosSparse::Experimental::BsrMatrix<Scalar, Ordinal, Device, void, Offset>;
  std::vector<size_t> ks = {1, 3};

  for (size_t k : ks) {  // multivector sizes
    std::string name = std::string("MatrixMarketConvert") + "/" + std::string(path.stem()) + "/" +
                       Kokkos::ArithTraits<Scalar>::name() + "/" + Kokkos::ArithTraits<Ordinal>::name() + "/" +
                       Kokkos::ArithTraits<Offset>::name() + "/" + std::to_string(bs) + "/" + std::to_string(k) + "/" +
                       Spmv::name() + "/" + Device::name();
    benchmark::RegisterBenchmark(name.c_str(), read_convert_run<Bsr, Spmv>, path, bs, k)->UseRealTime();
  }
}

template <typename Device>
void register_converts(const fs::path &path, const size_t bs) {
  std::cerr << "benchmarks will use detected blocksize\n";
  // clang-format off
  register_convert_type<int, float, int, Device, SpmvDefault>(path, bs);
  register_convert_type<int, float, int, Device, SpmvNative>(path, bs);
  register_convert_type<int, float, int, Device, SpmvV41>(path, bs);

  register_convert_type<int, double, int, Device, SpmvDefault>(path, bs);
  register_convert_type<int, double, int, Device, SpmvNative>(path, bs);
  register_convert_type<int, double, int, Device, SpmvV41>(path, bs);

  register_convert_type<int, float, unsigned, Device, SpmvDefault>(path, bs);
  register_convert_type<int, float, unsigned, Device, SpmvNative>(path, bs);
  register_convert_type<int, float, unsigned, Device, SpmvV41>(path, bs);

  register_convert_type<int64_t, double, size_t, Device, SpmvDefault>(path, bs);
  register_convert_type<int64_t, double, size_t, Device, SpmvNative>(path, bs);
  register_convert_type<int64_t, double, size_t, Device, SpmvV41>(path, bs);

  register_convert_type<int64_t, double, int64_t, Device, SpmvDefault>(path, bs);
  register_convert_type<int64_t, double, int64_t, Device, SpmvNative>(path, bs);
  register_convert_type<int64_t, double, int64_t, Device, SpmvV41>(path, bs);

  // clang-format on
}

template <typename Device>
void register_expands(const fs::path &path) {
  std::cerr << "benchmarks will expand each non-zero into a larger block\n";
  // clang-format off
  register_expand_type<int, float, int, Device, SpmvDefault>(path);
  register_expand_type<int, float, int, Device, SpmvNative>(path);
  register_expand_type<int, float, int, Device, SpmvV41>(path);

  register_expand_type<int, double, int, Device, SpmvDefault>(path);
  register_expand_type<int, double, int, Device, SpmvNative>(path);
  register_expand_type<int, double, int, Device, SpmvV41>(path);

  register_expand_type<int, float, unsigned, Device, SpmvDefault>(path);
  register_expand_type<int, float, unsigned, Device, SpmvNative>(path);
  register_expand_type<int, float, unsigned, Device, SpmvV41>(path);

  register_expand_type<int64_t, double, uint64_t, Device, SpmvDefault>(path);
  register_expand_type<int64_t, double, uint64_t, Device, SpmvNative>(path);
  register_expand_type<int64_t, double, uint64_t, Device, SpmvV41>(path);

  register_expand_type<int64_t, double, int64_t, Device, SpmvDefault>(path);
  register_expand_type<int64_t, double, int64_t, Device, SpmvNative>(path);
  register_expand_type<int64_t, double, int64_t, Device, SpmvV41>(path);
  // clang-format on
}

template <typename Device>
void register_path(const fs::path &path) {
  size_t detectedSize;
  try {
    detectedSize = detect_block_size<Device>(path);
  } catch (const std::exception &e) {
    std::cerr << "ERROR while reading: " << e.what() << "\n"
              << "skipping!\n";
    return;
  }

  /* If a block size can be detected, just use that block size without
     expanding the matrix.
     Otherwise, expand the matrix to some arbitrary block sizes to test BSR
  */
  if (detectedSize != 1) {
    register_converts<Device>(path, detectedSize);
  } else {
    register_expands<Device>(path);
  }
}

int main(int argc, char **argv) {
  Kokkos::initialize(argc, argv);
  benchmark::Initialize(&argc, argv);
  benchmark::SetDefaultTimeUnit(benchmark::kMicrosecond);
  KokkosKernelsBenchmark::add_benchmark_context(true);

  for (int i = 1; i < argc; ++i) {
#if defined(KOKKOS_ENABLE_CUDA)
    register_path<Kokkos::Cuda>(argv[i]);
#endif
#if defined(KOKKOS_ENABLE_HIP)
    register_path<Kokkos::HIP>(argv[i]);
#endif
#if defined(KOKKOS_ENABLE_SERIAL)
    register_path<Kokkos::Serial>(argv[i]);
#endif
  }

  benchmark::RunSpecifiedBenchmarks();

  benchmark::Shutdown();
  drop_cache();
  Kokkos::finalize();
  return 0;
}
