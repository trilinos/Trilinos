#include <time.h>
void second_linux(double *seconds)
{
clock_t time_now;
/* get current counter in integer seconds */
time_now = clock();
*seconds = ((double)  time_now )/ CLOCKS_PER_SEC;
return;
}
