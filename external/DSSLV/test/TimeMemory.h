#ifndef TIMEMEMORY
#define TIMEMEMORY

void  gettimes(double *user_start_time, double *system_start_time, 
	       double *wall_start_time ) ; 
void  getmemusage(long *virtual_mem_start, long *resident_mem_start, 
		  long *page_faults_start ); 
#endif 
