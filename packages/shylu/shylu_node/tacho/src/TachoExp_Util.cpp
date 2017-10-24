#include <stdio.h>
#include <string.h>

namespace Tacho {

  double g_time_per_thread[2048] = {};
#if defined( TACHO_PROFILE_TIME_PER_THREAD )
  void resetTimePerThread() {
    memset((void*)&g_time_per_thread[0], 0, sizeof(double)*2048);
  }
  void getTimePerThread(const int nthreads, 
                        double &t_total,
                        double &t_avg,
                        double &t_min,
                        double &t_max) {
    t_total = 0; t_min = g_time_per_thread[0]; t_max = g_time_per_thread[0];
    for (int i=0;i<nthreads;++i) {
      const double t = g_time_per_thread[i];
      t_total += t;
      t_min = t_min < t ? t_min : t;
      t_max = t_max > t ? t_max : t;
    }
    t_avg = t_total / double(nthreads);
  }
#else  
  void resetTimePerThread() {
  // no-op
  }
  void getTimePerThread(const int nthreads, 
                        double &t_total,
                        double &t_avg,
                        double &t_min,
                        double &t_max) {
    t_total = 0; t_avg = 0; t_min = 0; t_max = 0;
  }
#endif  
  const char* Version() {
    return "Tacho:: Trilinos Git";
  }


}


