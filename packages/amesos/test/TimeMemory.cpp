
#include "time.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>

void  gettimes(double *user_time, double *system_time, 
	       double *wall_time ) {

#ifdef TFLOP 
  *user_time = 0;
  *system_time = 0 ; 
#else
  struct rusage self_usage ; 
  int usage_status = getrusage( RUSAGE_SELF, &self_usage ) ; 
  *user_time =  (double) self_usage.ru_utime.tv_sec 
    + (double) self_usage.ru_utime.tv_usec /1000000.0;
  *system_time = (double) self_usage.ru_stime.tv_sec 
    + (double) self_usage.ru_stime.tv_usec /1000000.0 ;
#endif

  time_t now = time( NULL ) ; 
  *wall_time = (double) now ; 

}
void  getmemusage(long *virtual_mem, long *resident_mem, 
		  long *page_faults ) {

#if (!defined(SOLARIS ) && !defined(TFLOP) ) 
  char procfile[100];
  char comm[100], state[10] ;

  int pid, ppid, pgrp, session, tty, tpgid, flags, 
    miknflt, cminflt, majflt,
    cmajflt, utime, stime, cutime, cstime, counter, priority, timeout, 
    itrealvalue, starttime, vsize, rss, rlim, startcode, endcode, 
    startstack, kstkesp, kstkeip, signal, blocked, sigignore, sigcatch, 
    wchan ;

  sprintf( procfile, "/proc/%d/stat", getpid() ) ; 
  FILE *open_fd = fopen ( procfile, "r" ) ; 
  assert( open_fd ) ; 

  fscanf( open_fd, "%d %s %c %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d ",
	  &pid, &comm, &state, & ppid, &pgrp, &session, &tty, &tpgid, &flags, 
	  &miknflt, &cminflt, &majflt, &cmajflt, &utime, &stime, &cutime, 
	  &cstime, &counter, &priority, &timeout, &itrealvalue, &starttime, 
	  &vsize, &rss, &rlim, &startcode, &endcode, &startstack, &kstkesp, 
	  &kstkeip, &signal, &blocked, &sigignore, &sigcatch, &wchan ) ; 

  *virtual_mem = vsize ; 
  *resident_mem = rss; 
  *page_faults = majflt; 
#endif
}
