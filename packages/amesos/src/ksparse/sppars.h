#ifndef SPPARS_H
#define SPPARS_H
#ifdef SHARED_MEM
#define MAX_STRIPS 4
#else
#define MAX_STRIPS 1
#endif

/* Only works with MIN_PES_SOLVE=2, but this is the only reasonable value! */
#define MIN_PES_SOLVE 2

#define OF_THRESHOLD 10000
#endif
