/* NOTE: This timer assumes it is being called from Fortran on a 
         machine that appends underscores.  
*/
#if defined(UseClock)

#include <time.h>
double second_(void)
{
   clock_t t1;
   static clock_t t0=0;
   static double CPS = CLOCKS_PER_SEC;
   double d;

   if (t0 == 0) t0 = clock();
   t1 = clock() - t0;
   d = t1 / CPS;
   return(d);
}

#elif defined(WALL)

#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
double second_(void)
{
   struct timeval tp;
   static long start=0, startu;
   if (!start)
   {
      gettimeofday(&tp, NULL);
      start = tp.tv_sec;
      startu = tp.tv_usec;
      return(0.0);
   }
   gettimeofday(&tp, NULL);
   return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );
}

#elif defined(UseTimes)

#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
double second_(void)
{
   struct tms ts;
   static double ClockTick=0.0;

   if (ClockTick == 0.0) ClockTick = (double) sysconf(_SC_CLK_TCK);
   times(&ts);
   return( (double) ts.tms_utime / ClockTick );
}

#else

#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
double second_(void)
{
   struct rusage ruse;
   getrusage(RUSAGE_SELF, &ruse);
   return( (double)(ruse.ru_utime.tv_sec+ruse.ru_utime.tv_usec / 1000000.0) );
}

#endif
