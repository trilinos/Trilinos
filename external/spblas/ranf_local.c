#include <stdlib.h>

/* #define DEBUG */

/* Prototypes */

double drand48(void);

#ifdef pclinux
double ranf_local__()
#else
double ranf_local_()
#endif
{
#ifdef DEBUG
double val = 1.0;
return(val);
#else
return(drand48());
#endif
}


