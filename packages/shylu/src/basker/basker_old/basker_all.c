
#define Int int
#define BASKER(name) basker_ ## name

#include "basker.c"

#undef Int
#undef BASKER
#define Int long
#define BASKER(name) basker_ ## name ## _l

#include "basker.c"
