#ifndef STK_IO_DBStepTimeInterval
#define STK_IO_DBStepTimeInterval

#include <limits>
#include <Ioss_Region.h>

namespace stk {
  namespace io {
    struct DBStepTimeInterval
    {
      DBStepTimeInterval(Ioss::Region *region, double time)
	: region(region),
	  t_analysis(time),
	  t_before(-std::numeric_limits<double>::max()),
	  t_after(std::numeric_limits<double>::max()),
	  s_before(0),
	  s_after(0),
	  exists_before(false),
	  exists_after(false)
      {
	// If no steps on database, exists_before and exists_after will be false.
	size_t step_count = region->get_property("state_count").get_int();
	for (size_t istep = 0; istep < step_count; istep++) {
	  double state_time = region->get_state_time(istep+1);
	  if (state_time <= time) {
	    if (state_time > t_before) {
	      t_before = state_time;
	      s_before = istep+1;
	      exists_before = true;
	    }
	  }

	  if (state_time >= time) {
	    if (state_time < t_after) {
	      t_after = state_time;
	      s_after = istep+1;
	      exists_after = true;
	    }
	  }
	}
      }
    
      size_t get_closest_step() const
      {
	size_t step = 0;
	if (exists_before && !exists_after) {
	  step = s_before;
	}
	else if (!exists_before && exists_after) {
	  step = s_after;
	}
	else {
	  double delta_b = t_analysis - t_before;
	  double delta_a = t_after    - t_analysis;
	  step = delta_b < delta_a ? s_before : s_after;
	}
	return step;
      }

      Ioss::Region *region; 
      double t_analysis;  // Analysis time.
      double t_before;    // Time at step s_before; valid only if exists_before == true
      double t_after;     // Time at step s_after; valid only if exists_after == true
      size_t s_before;    // Step before passed in time (1-based).
      size_t s_after;     // Step after passed in time (1-based).
      bool exists_before; // True if there is a db step with time <= passed in time
      bool exists_after;  // True if there is a db step with time >= passed in time
    };
  }
}
#endif
