/* timer for win 2000 pc.
   author: rbl */

#include <time.h>

double second(void) {
   clock_t start;
   double duration;

   start = clock();
  return (double)( start ) / CLOCKS_PER_SEC;
}
