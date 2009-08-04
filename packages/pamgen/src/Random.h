// $Id$


// the system defined random number generator
#ifdef _MSC_VER
# define SRANDOM(i) srand(i)
# define RANDOM() rand()
#else
# define SRANDOM(i) srandom(i)
# define RANDOM() random()
#endif
