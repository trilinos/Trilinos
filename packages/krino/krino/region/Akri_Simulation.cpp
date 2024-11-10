// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Simulation.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_Region.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/human_bytes.hpp>

#include <cmath>

namespace krino{

Simulation::Simulation(const std::string & in_name)
  : my_name(in_name),
    my_timer("krino", sierra::Diag::sierraTimer()),
    my_is_transient(false),
    my_stop_time(0.0),
    my_step_count(0),
    my_current_time(0.0),
    my_time_step_size(0.0)
{
  stk::diag::setEnabledTimerMetricsMask(stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME);
}

Simulation::~Simulation() {}

void Simulation::add_region(Region * region)
{
  my_regions.emplace_back(region);
}

void Simulation::commit()
{
  for (auto && region : my_regions)
  {
    region->commit();
  }
}

void Simulation::execute()
{
  stk::diag::TimeBlock timer__(my_timer);

  //=====================================
  // start-up procedure
  //=====================================

  // initial conditions
  for (auto && region : my_regions)
  {
    region->initialize();
  }

  for (auto && region : my_regions)
  {
    region->process_output(false);
  }

  //=====================================
  // time integration
  //=====================================

  while ( my_current_time < my_stop_time )
  {
    static const double sqrt_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
    my_old_time = my_current_time;
    if (my_current_time+my_time_step_size > my_stop_time*(1.-sqrt_epsilon))
    {
      my_time_step_size = my_stop_time - my_current_time;
      my_current_time = my_stop_time;
    }
    else
    {
      my_current_time += my_time_step_size;
    }
    ++my_step_count;

    const double percentComplete = 100.0*my_current_time/my_stop_time;

    krinolog << "Krino Simulation " << my_name
        << ", step " << my_step_count
        << ", time " << my_current_time
        << ", time step " << my_time_step_size
        << ", " << percentComplete << "% complete"
        << stk::diag::dendl;

    for (auto && region : my_regions)
    {
      region->execute();
    }

    for (auto && region : my_regions)
    {
      region->process_output(false);
    }

    for (auto && region : my_regions)
    {
      region->mesh_bulk_data().update_field_data_states();
    }
  }

  timer__.stop();
  print_performance_info();
}

void
Simulation::print_performance_info() const
{
  sierra::Env::outputP0() << sierra::Env::section_separator() << std::endl
                  << "Timing summary running on " << sierra::Env::parallel_size() << " processor" << (sierra::Env::parallel_size() == 1 ? "" : "s") << std::endl;
  stk::diag::printTimersTable(sierra::Env::outputP0(), sierra::Diag::sierraTimer(), stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME, false, sierra::Env::parallel_comm());
  sierra::Env::outputP0() << std::endl << std::endl;

  {
    size_t now, hwm;
    stk::get_memory_usage(now, hwm);
    // min, max, sum
    size_t global_now[3] = {now,now,now};
    size_t global_hwm[3] = {hwm,hwm,hwm};

    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceSum<1>( &global_now[2] ) );
    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceMin<1>( &global_now[0] ) );
    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceMax<1>( &global_now[1] ) );

    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceSum<1>( &global_hwm[2] ) );
    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceMin<1>( &global_hwm[0] ) );
    stk::all_reduce(sierra::Env::parallel_comm(), stk::ReduceMax<1>( &global_hwm[1] ) );

    sierra::Env::outputP0() << "Memory Overview: " << std::endl;

    sierra::Env::outputP0() << "  total (over all cores) current/high-water mark = "
                            << std::setw(15) << stk::human_bytes(global_now[2])
                            << std::setw(15) << stk::human_bytes(global_hwm[2])
                            << std::endl;

    sierra::Env::outputP0() << "    min (over all cores) current/high-water mark = "
                            << std::setw(15) << stk::human_bytes(global_now[0])
                            << std::setw(15) << stk::human_bytes(global_hwm[0])
                            << std::endl;

    sierra::Env::outputP0() << "    max (over all cores) current/high-water mark = "
                            << std::setw(15) << stk::human_bytes(global_now[1])
                            << std::setw(15) << stk::human_bytes(global_hwm[1])
                            << std::endl;
  }
}

} // namespace krino
