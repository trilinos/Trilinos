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

#include "KokkosODE_BDF.hpp"

#include "KokkosKernels_TestUtils.hpp"
#include "KokkosKernels_perf_test_utilities.hpp"

#include <benchmark/benchmark.h>
#include "Benchmark_Context.hpp"

namespace {
// Robertson's autocatalytic chemical reaction:
// H. H. Robertson, The solution of a set of reaction rate equations,
// in J. Walsh (Ed.), Numerical Analysis: An Introduction,
// pp. 178–182, Academic Press, London (1966).
//
// Equations: y0' = -0.04*y0 + 10e4*y1*y2
//            y1' =  0.04*y0 - 10e4*y1*y2 - 3e7 * y1**2
//            y2' = 3e7 * y1**2
struct StiffChemistry {
  static constexpr int neqs = 3;

  StiffChemistry() {}

  template <class vec_type1, class vec_type2>
  KOKKOS_FUNCTION void evaluate_function(const double /*t*/, const double /*dt*/, const vec_type1& y,
                                         const vec_type2& f) const {
    f(0) = -0.04 * y(0) + 1.e4 * y(1) * y(2);
    f(1) = 0.04 * y(0) - 1.e4 * y(1) * y(2) - 3.e7 * y(1) * y(1);
    f(2) = 3.e7 * y(1) * y(1);
  }

  template <class vec_type, class mat_type>
  KOKKOS_FUNCTION void evaluate_jacobian(const double /*t*/, const double /*dt*/, const vec_type& y,
                                         const mat_type& jac) const {
    jac(0, 0) = -0.04;
    jac(0, 1) = 1.e4 * y(2);
    jac(0, 2) = 1.e4 * y(1);
    jac(1, 0) = 0.04;
    jac(1, 1) = -1.e4 * y(2) - 3.e7 * 2.0 * y(1);
    jac(1, 2) = -1.e4 * y(1);
    jac(2, 0) = 0.0;
    jac(2, 1) = 3.e7 * 2.0 * y(1);
    jac(2, 2) = 0.0;
  }
};

template <class ode_type, class mat_type, class vec_type, class scalar_type>
struct BDF_Solve_wrapper {
  const ode_type my_ode;
  const scalar_type t_start, t_end, dt, max_step;
  const vec_type y0, y_new;
  const mat_type temp, temp2;

  BDF_Solve_wrapper(const ode_type& my_ode_, const scalar_type& t_start_, const scalar_type& t_end_,
                    const scalar_type& dt_, const scalar_type& max_step_, const vec_type& y0_, const vec_type& y_new_,
                    const mat_type& temp_, const mat_type& temp2_)
      : my_ode(my_ode_),
        t_start(t_start_),
        t_end(t_end_),
        dt(dt_),
        max_step(max_step_),
        y0(y0_),
        y_new(y_new_),
        temp(temp_),
        temp2(temp2_) {}

  KOKKOS_FUNCTION void operator()(const int idx) const {
    auto subTemp  = Kokkos::subview(temp, Kokkos::ALL(), Kokkos::ALL(), idx);
    auto subTemp2 = Kokkos::subview(temp2, Kokkos::ALL(), Kokkos::ALL(), idx);
    auto subY0    = Kokkos::subview(y0, Kokkos::ALL(), idx);
    auto subYnew  = Kokkos::subview(y_new, Kokkos::ALL(), idx);

    KokkosODE::Experimental::BDFSolve(my_ode, t_start, t_end, dt, max_step, subY0, subYnew, subTemp, subTemp2);
  }
};

}  // namespace

struct bdf_input_parameters {
  int num_odes;
  int repeat;
  bool verbose;

  bdf_input_parameters(const int num_odes_, const int repeat_, const bool verbose_)
      : num_odes(num_odes_), repeat(repeat_), verbose(verbose_){};
};

template <class execution_space>
void run_ode_chem(benchmark::State& state, const bdf_input_parameters& inputs) {
  using scalar_type = double;
  using KAT         = Kokkos::ArithTraits<scalar_type>;
  using vec_type    = Kokkos::View<scalar_type**, execution_space>;
  using mat_type    = Kokkos::View<scalar_type***, execution_space>;

  StiffChemistry mySys{};

  const bool verbose = inputs.verbose;
  const int num_odes = inputs.num_odes;
  const int neqs     = mySys.neqs;

  const scalar_type t_start = KAT::zero(), t_end = 350 * KAT::one();
  scalar_type dt = KAT::zero();
  vec_type y0("initial conditions", neqs, num_odes);
  vec_type y_new("solution", neqs, num_odes);

  // Set initial conditions
  auto y0_h = Kokkos::create_mirror(y0);
  for (int sysIdx = 0; sysIdx < num_odes; ++sysIdx) {
    y0_h(0, sysIdx) = KAT::one();
    y0_h(1, sysIdx) = KAT::zero();
    y0_h(2, sysIdx) = KAT::zero();
  }

  mat_type temp("buffer1", neqs, 23 + 2 * neqs + 4, num_odes), temp2("buffer2", 6, 7, num_odes);

  if (verbose) {
    std::cout << "Number of problems solved in parallel: " << num_odes << std::endl;
  }

  Kokkos::RangePolicy<execution_space> policy(0, num_odes);

  Kokkos::Timer time;
  time.reset();
  for (auto _ : state) {
    (void)_;

    // Set initial conditions for each test iteration
    state.PauseTiming();
    dt = KAT::zero();
    Kokkos::deep_copy(y0, y0_h);
    Kokkos::deep_copy(y_new, KAT::zero());
    Kokkos::deep_copy(temp, KAT::zero());
    Kokkos::deep_copy(temp2, KAT::zero());
    BDF_Solve_wrapper bdf_wrapper(mySys, t_start, t_end, dt, (t_end - t_start) / 10, y0, y_new, temp, temp2);
    state.ResumeTiming();

    // Actually run the time integrator
    Kokkos::parallel_for(policy, bdf_wrapper);
    Kokkos::fence();
  }
  double run_time = time.seconds();
  std::cout << "Run time: " << run_time << std::endl;

  Kokkos::deep_copy(y0_h, y0);
  double error;
  for (int odeIdx = 0; odeIdx < num_odes; ++odeIdx) {
    error = 0;
    // error += Kokkos::abs(y0_h(0, odeIdx) - 0.4193639) / 0.4193639;
    // error += Kokkos::abs(y0_h(1, odeIdx) - 0.000002843646) / 0.000002843646;
    // error += Kokkos::abs(y0_h(2, odeIdx) - 0.5806333) / 0.5806333;
    error += Kokkos::abs(y0_h(0, odeIdx) - 0.462966) / 0.462966;
    error += Kokkos::abs(y0_h(1, odeIdx) - 3.42699e-06) / 3.42699e-06;
    error += Kokkos::abs(y0_h(2, odeIdx) - 0.537030) / 0.537030;
    error = error / 3;

    if (error > 1e-6) {
      std::cout << "Large error in problem " << odeIdx << ": " << error << std::endl;
    }
  }
}

void print_options() {
  std::cerr << "Options\n" << std::endl;

  std::cerr << perf_test::list_common_options();

  std::cerr << "\t[Optional] --repeat      :: how many times to repeat overall test" << std::endl;
  std::cerr << "\t[Optional] --verbose     :: enable verbose output" << std::endl;
  std::cerr << "\t[Optional] --n           :: number of ode problems to solve" << std::endl;
}  // print_options

int parse_inputs(bdf_input_parameters& params, int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (perf_test::check_arg_int(i, argc, argv, "--n", params.num_odes)) {
      ++i;
    } else if (perf_test::check_arg_int(i, argc, argv, "--repeat", params.repeat)) {
      ++i;
    } else if (perf_test::check_arg_bool(i, argc, argv, "--verbose", params.verbose)) {
    } else {
      std::cerr << "Unrecognized command line argument #" << i << ": " << argv[i] << std::endl;
      print_options();
      return 1;
    }
  }
  return 0;
}  // parse_inputs

template <class execution_space>
void run_benchmark_wrapper(benchmark::State& state, bdf_input_parameters params) {
  run_ode_chem<execution_space>(state, params);
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  {
    benchmark::Initialize(&argc, argv);
    benchmark::SetDefaultTimeUnit(benchmark::kMillisecond);
    KokkosKernelsBenchmark::add_benchmark_context(true);

    perf_test::CommonInputParams common_params;
    perf_test::parse_common_options(argc, argv, common_params);

    std::string bench_name = "KokkosODE_BDF_Stiff_Chem";
    bdf_input_parameters params(1000, 1, false);
    parse_inputs(params, argc, argv);

    if (0 < common_params.repeat) {
      benchmark::RegisterBenchmark(bench_name.c_str(), run_benchmark_wrapper<Kokkos::DefaultExecutionSpace>, params)
          ->UseRealTime()
          ->ArgNames({"n"})
          ->Args({params.num_odes})
          ->Iterations(common_params.repeat);
    } else {
      benchmark::RegisterBenchmark(bench_name.c_str(), run_benchmark_wrapper<Kokkos::DefaultExecutionSpace>, params)
          ->UseRealTime()
          ->ArgNames({"n"})
          ->Args({params.num_odes});
    }

    benchmark::RunSpecifiedBenchmarks();

    benchmark::Shutdown();
  }
  Kokkos::finalize();

  return 0;
}
