#include <sys/time.h>
double second_()
{
hrtime_t gethrtime(void);
hrtime_t ticks;
double seconds;
/* get current counter in integer nanoseconds */
ticks = gethrtime();
seconds = ticks / 1000000000.;
return(seconds);
}
