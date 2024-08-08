// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_Simulation_h
#define Akri_Simulation_h

#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/diag/Timer.hpp>
#include <memory>
#include <string>
#include <vector>

namespace krino{

class Region;

class Simulation {
public:
  Simulation(const std::string & in_name);
  ~Simulation();

  void commit();
  void initialize() {}
  void execute();

  void add_region(Region * region);

  const std::string & get_name() const { return my_name; }
  bool is_transient() const { return my_is_transient; }
  double get_time_step() const { return my_time_step_size; }
  double get_old_time() const { return my_old_time; }
  double get_current_time() const { return my_current_time; }
  double get_stop_time() const { return my_stop_time; }
  int get_time_step_count() const { return my_step_count; }

  void set_time_step(const double in_time_step) { my_time_step_size = in_time_step; }
  void set_current_time(const double in_time) { my_current_time = in_time; }
  void set_stop_time(const double in_time) { my_is_transient = true; my_stop_time = in_time; }

  stk::diag::Timer & get_timer() const { return my_timer; }
  void print_performance_info() const;

private:
  std::string my_name;
  mutable stk::diag::Timer my_timer;
  bool my_is_transient;
  double my_stop_time;
  unsigned my_step_count;
  double my_old_time;
  double my_current_time;
  double my_time_step_size;
  std::vector<std::unique_ptr<Region>> my_regions;
};

} // namespace krino

#endif // Akri_Simulation_h

