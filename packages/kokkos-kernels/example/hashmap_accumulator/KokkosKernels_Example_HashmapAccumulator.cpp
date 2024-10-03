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

#include <stdlib.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <sys/time.h>

// Kokkos Includes
#include <Kokkos_Core.hpp>
#include <Kokkos_UniqueToken.hpp>
#include <Kokkos_Timer.hpp>

// Kokkos Kernels Includes
#include <KokkosKernels_HashmapAccumulator.hpp>
#include <KokkosKernels_Uniform_Initialized_MemoryPool.hpp>

// Command Line Parameters structure
typedef struct params {
  uint32_t use_serial  = false;
  uint32_t use_threads = false;
  uint32_t use_cuda    = false;
  uint32_t use_openmp  = false;
  bool verbose         = false;

  size_t problem_size = 20;
  size_t repeat       = 1;
} parameters_t;

namespace KokkosKernels {
namespace Experiment {

template <typename ExecutionSpace, typename uniform_memory_pool_t, typename scalar_t>
struct functorTestHashmapAccumulator {
  typedef ExecutionSpace execution_space;
  typedef typename Kokkos::View<scalar_t*> data_view_t;

  const size_t _num_entries;
  const data_view_t _data;
  uniform_memory_pool_t _memory_pool;
  const size_t _hash_size;
  const size_t _max_hash_entries;
  const parameters_t& _params;

  typedef Kokkos::Experimental::UniqueToken<execution_space, Kokkos::Experimental::UniqueTokenScope::Global>
      unique_token_t;
  unique_token_t tokens;

  functorTestHashmapAccumulator(const size_t num_entries, const data_view_t& data, uniform_memory_pool_t memory_pool,
                                const size_t hash_size, const size_t max_hash_entries, const parameters_t& params)
      : _num_entries(num_entries),
        _data(data),
        _memory_pool(memory_pool),
        _hash_size(hash_size),
        _max_hash_entries(max_hash_entries),
        _params(params),
        tokens(ExecutionSpace()) {
    if (_params.verbose) {
      std::cout << "UniqueToken.size: " << tokens.size() << std::endl;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const scalar_t idx) const {
    typedef scalar_t hash_size_type;
    typedef scalar_t hash_key_type;
    typedef scalar_t hash_value_type;

    // Alternative to team_policy thread id
    auto tid = tokens.acquire();

    // Acquire a chunk from the memory pool using a spin-loop.
    volatile scalar_t* ptr_temp = nullptr;
    while (nullptr == ptr_temp) {
      ptr_temp = (volatile scalar_t*)(_memory_pool.allocate_chunk(tid));
    }
    scalar_t* ptr_memory_pool_chunk = (scalar_t*)(ptr_temp);

    KokkosKernels::Experimental::HashmapAccumulator<hash_size_type, hash_key_type, hash_value_type> hash_map;

    // Set pointer to hash indices
    scalar_t* used_hash_indices = (scalar_t*)(ptr_temp);
    ptr_temp += _hash_size;

    // Set pointer to hash begins
    hash_map.hash_begins = (scalar_t*)(ptr_temp);
    ptr_temp += _hash_size;

    // Set pointer to hash nexts
    hash_map.hash_nexts = (scalar_t*)(ptr_temp);
    ptr_temp += _max_hash_entries;

    // Set pointer to hash keys
    hash_map.keys = (scalar_t*)(ptr_temp);
    // ptr_temp += _max_hash_entries;

    // Set pointer to hash values
    // hash_map.values = (scalar_t*)(ptr_temp);

    // Set up limits in Hashmap_Accumulator
    hash_map.hash_key_size  = _max_hash_entries;
    hash_map.max_value_size = _max_hash_entries;

    // hash function is hash_size-1 (note: hash_size must be a power of 2)
    scalar_t hash_func_pow2 = _hash_size - 1;

    // These are updated by Hashmap_Accumulator insert functions.
    scalar_t used_hash_size  = 0;
    scalar_t used_hash_count = 0;

    // Loop over stuff
    for (size_t i = 0; i < _num_entries; i++) {
      scalar_t key = _data(i);

      // Compute the hash index using & instead of % (modulus is slower).
      scalar_t hash = key & hash_func_pow2;

      int r = hash_map.sequential_insert_into_hash_TrackHashes(hash, key, &used_hash_size, hash_map.max_value_size,
                                                               &used_hash_count, used_hash_indices);

      // Check return code
      if (r) {
        // insert should return nonzero if the insert failed, but for
        // sequential_insert_into_hash_TrackHashes the 'full' case is currently
        // ignored, so r will always be 0.
      }
    }

    // TODO: Get the # of unique values inserted and return that out of the
    // functor.

    // Reset the Begins values to -1 before releasing the memory pool chunk.
    // If you don't do this the next thread that grabs this memory chunk will
    // not work properly.
    for (scalar_t i = 0; i < used_hash_count; i++) {
      scalar_t dirty_hash              = used_hash_indices[i];
      hash_map.hash_begins[dirty_hash] = -1;
    }

    // Release the memory pool chunk back to the pool
    _memory_pool.release_chunk(ptr_memory_pool_chunk);

    // Release the UniqueToken
    tokens.release(tid);

  }  // operator()

};  // functorTestHashmapAccumulator

template <typename execution_space, typename scalar_t = int>
void experiment(const parameters_t& params) {
  typedef typename KokkosKernels::Impl::UniformMemoryPool<execution_space, scalar_t> uniform_memory_pool_t;
  typedef typename Kokkos::View<scalar_t*> data_view_t;
  typedef typename data_view_t::HostMirror data_view_hostmirror_t;

  size_t num_entries = params.problem_size;

  // Set max value in the list
  size_t max_value = 100;

  // Get the concurrecny
  size_t concurrency = execution_space().concurrency();

  // Set up random number generator
  std::random_device rd;
  std::mt19937 eng(rd());
  std::uniform_int_distribution<scalar_t> distr(1, max_value);

  // Create a view of random values
  data_view_t d_data("data", num_entries);
  data_view_hostmirror_t h_data = Kokkos::create_mirror_view(d_data);

  for (size_t i = 0; i < num_entries; i++) {
    h_data(i) = distr(eng);
  }

  // Print out the array of random numbers if the list size is small.
  if (num_entries <= 50 || params.verbose) {
    std::cout << "Data: ";
    for (size_t i = 0; i < num_entries; i++) {
      std::cout << h_data(i) << " ";
    }
    std::cout << std::endl;
  }

  Kokkos::Timer timer;

  // Deep copy initialized values to device memory.
  Kokkos::deep_copy(d_data, h_data);

  // Set Hash Table Parameters
  size_t max_hash_entries = max_value;  // Max number of entries that can be
                                        // inserted (values allowed are 1..100)
  size_t hash_size_hint = max_value;    // How many hash keys are allowed. The actual hash size will
                                        // be set to the next power of 2 bigger than hash_size_hint.

  // Set the hash_size as the next power of 2 bigger than hash_size_hint.
  // - hash_size must be a power of two since we use & rather than % (which is
  // slower) for computing the hash value for HashmapAccumulator.
  size_t hash_size = 1;
  while (hash_size < hash_size_hint) {
    hash_size *= 2;
  }

  // Create Uniform Initialized Memory Pool
  KokkosKernels::Impl::PoolType pool_type = KokkosKernels::Impl::OneThread2OneChunk;

  // Determine memory chunk size for UniformMemoryPool
  size_t mem_chunk_size = hash_size;   // for hash indices
  mem_chunk_size += hash_size;         // for hash begins
  mem_chunk_size += max_hash_entries;  // for hash nexts
  mem_chunk_size += max_hash_entries;  // for hash keys
  // mem_chunk_size += max_entries;          // for hash values

  // Set a cap on # of chunks to 32.  In application something else should be
  // done here differently if we're OpenMP vs. GPU but for this example we can
  // just cap our number of chunks at 32.
  size_t mem_chunk_count = KOKKOSKERNELS_MACRO_MIN(32, concurrency);

  // KokkosKernels::Impl::UniformMemoryPool<Kokkos::DefaultExecutionSpace,
  // size_t> m_space(mem_chunk_count, mem_chunk_size, -1, pool_type);
  uniform_memory_pool_t memory_pool(mem_chunk_count, mem_chunk_size, -1, pool_type);

  functorTestHashmapAccumulator<execution_space, uniform_memory_pool_t, scalar_t> testHashmapAccumulator(
      num_entries, d_data, memory_pool, hash_size, max_hash_entries, params);

  Kokkos::parallel_for("testHashmapAccumulator", num_entries, testHashmapAccumulator);

  if (params.verbose) {
    double t = timer.seconds();
    std::cout << "Execution Time: " << std::setw(-2) << t << std::endl;
    timer.reset();
  }
}

}  // namespace Experiment
}  // namespace KokkosKernels

void print_options(std::ostream& os, const char* app_name, unsigned int indent = 0) {
  std::string spaces(indent, ' ');
  os << "Usage:" << std::endl
     << spaces << "  " << app_name << " [parameters]" << std::endl
     << std::endl
     << spaces << "Parameters:" << std::endl
     << spaces << "  Parallelism (select one of the following):" << std::endl
     << spaces << "      --serial <N>        Execute serially." << std::endl
     << spaces << "      --threads <N>       Use N posix threads." << std::endl
     << spaces << "      --openmp <N>        Use OpenMP with N threads." << std::endl
     << spaces << "      --cuda              Use CUDA" << std::endl
     << spaces << "  Optional Parameters:" << std::endl
     << spaces << "      --problem-size <N>  Problem Size (Default: 20)" << std::endl
     << spaces << "      --verbose           Verbose output" << std::endl
     << spaces << "      --help              Print out command line help." << std::endl
     << spaces << " " << std::endl;
}  // print_options

// int parse_inputs(KokkosKernels::Experiment::Parameters &params, int argc,
// char **argv)
int parse_inputs(parameters_t& params, int argc, char** argv) {
  if (argc == 1) {
    print_options(std::cout, argv[0]);
    return 1;
  }

  for (int i = 1; i < argc; ++i) {
    if (0 == Test::string_compare_no_case(argv[i], "--threads")) {
      params.use_threads = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--serial")) {
      params.use_serial = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--openmp")) {
      params.use_openmp = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--cuda")) {
      params.use_cuda = 1;
    } else if (0 == Test::string_compare_no_case(argv[i], "--repeat")) {
      params.repeat = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--problem-size")) {
      params.problem_size = atoi(argv[++i]);
    } else if (0 == Test::string_compare_no_case(argv[i], "--verbose") ||
               0 == Test::string_compare_no_case(argv[i], "-V")) {
      params.verbose = true;
    } else if (0 == Test::string_compare_no_case(argv[i], "help") || 0 == Test::string_compare_no_case(argv[i], "-h")) {
      print_options(std::cout, argv[0]);
      return 1;
    } else {
      std::cerr << "3-Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options(std::cout, argv[0]);
      return 1;
    }
  }
  if (!params.use_serial && !params.use_threads && !params.use_openmp && !params.use_cuda) {
    print_options(std::cout, argv[0]);
    return 1;
  }
  return 0;
}  // parse_inputs

int main(int argc, char* argv[]) {
  // KokkosKernels::Experiment::Parameters params;
  parameters_t params;

  // Override default repeats (default is 6)
  params.repeat = 1;

  if (parse_inputs(params, argc, argv)) {
    return 1;
  }

  const int device_id   = 0;
  const int num_threads = params.use_openmp;  // Assumption is that use_openmp variable is provided
                                              // as number of threads

  Kokkos::initialize(Kokkos::InitializationSettings().set_num_threads(num_threads).set_device_id(device_id));

  if (params.verbose) {
    Kokkos::print_configuration(std::cout);
  }

  // Work goes here.
  KokkosKernels::Experiment::experiment<Kokkos::DefaultExecutionSpace>(params);

  Kokkos::finalize();
  std::cout << "Done." << std::endl;
  return 0;
}
