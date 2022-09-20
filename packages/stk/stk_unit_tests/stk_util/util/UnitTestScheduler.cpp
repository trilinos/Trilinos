// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <random>
#include <stk_util/environment/Scheduler.hpp>
#include <stk_util/environment/WallTime.hpp>

#include "gtest/gtest.h"

TEST(Scheduler, AdjustDt)
{
  stk::util::Scheduler scheduler;
  scheduler.set_lookahead(4);
  scheduler.add_interval(0.0, 1.0);

  const double eps = std::numeric_limits<double>::epsilon();

  {
    const double t = 0.9 - 10 * eps;
    const double dt = 0.1;
    const double dt_new = scheduler.adjust_dt(dt, t);
    EXPECT_NEAR(dt + 10 * eps, dt_new, eps);
  }

  {
    const double t = 1.0 - 10 * eps;
    const double dt = 0.1;
    const double dt_new = scheduler.adjust_dt(dt, t);
    EXPECT_DOUBLE_EQ(dt, dt_new);
  }

  {
    const double t = 0.85;
    const double dt = 0.1;
    const double dt_new = scheduler.adjust_dt(dt, t);
    EXPECT_NEAR(0.075, dt_new, eps);
  }

  {
    const double t = 0.95;
    const double dt = 0.1;
    const double dt_new = scheduler.adjust_dt(dt, t);
    EXPECT_NEAR(0.05, dt_new, eps);
  }
}

TEST(Scheduler, LogarithmicOutput)
{
  stk::util::Scheduler scheduler;
  scheduler.set_lookahead(4);
  scheduler.add_interval(0.0, 1e-5);
  scheduler.add_interval(1e-4, 1e-4);
  scheduler.add_interval(1e-3, 1e-3);
  scheduler.add_interval(1e-2, 1e-2);
  scheduler.add_interval(1e-1, 1e-1);
  scheduler.add_interval(1.0, 1.0);
  scheduler.add_interval(10.0, 10.0);
  scheduler.add_interval(100.0, 100.0);

  double t = 0.0;
  double dt = 1e-6;
  const double dt_max = 100.0;

  std::mt19937 rng;
  rng.seed(666);

  std::uniform_real_distribution<double> noise(-1.0, 1.0);

  int step = 0;
  std::vector<double> output_times;
  while (t < 3600.0) {
    const double fac = 1.1 + (0.1 - std::numeric_limits<double>::epsilon()) * noise(rng);
    dt = std::min(dt_max, fac * dt);
    dt = scheduler.adjust_dt(dt, t);

    EXPECT_TRUE(dt >= 1e-6);

    if (scheduler.is_it_time(t, step)) {
      output_times.push_back(t);
    }

    t += dt;
    ++step;
  }

  std::vector<double> expected_times = {0.0, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 1e-4, 2e-4, 3e-4,
      4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2, 2e-2, 3e-2, 4e-2,
      5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
      7.0, 8.0, 9.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
      1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100,
      3200, 3300, 3400, 3500};

  const double eps = 10 * std::numeric_limits<double>::epsilon();
  for (unsigned i = 0; i < expected_times.size(); ++i) {
    EXPECT_NEAR(expected_times[i], output_times[i], expected_times[i] * eps);
  }
}
