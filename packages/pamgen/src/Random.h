// $Id$


// the system defined random number generator
#ifdef _WIN32
# define SRANDOM(i) srand(i)
# define RANDOM() rand()
#else
# define SRANDOM(i) srandom(i)
# define RANDOM() random()
#endif
