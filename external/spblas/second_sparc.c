#include <sys/time.h>
void second_sparc_(double *seconds)
{
hrtime_t gethrtime(void);
hrtime_t ticks;
/* get current counter in integer nanoseconds */
ticks = gethrtime();
*seconds = ticks / 1000000000.;
return;
}
