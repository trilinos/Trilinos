#ifndef TPETRA_DETAILS_DEBUG_CWP
#define TPETRA_DETAILS_DEBUG_CWP

#if CWP_SHUTUP == 0
#define CWP_CERR(x)
#else
#define CWP_CERR(x) (std::cerr << x)
#endif

#if CWP_SHUTUP == 0
#define CWP_PRINTF(x, ...)
#else
#define CWP_PRINTF(x, ...) printf(x, __VA_ARGS__)
#endif

#endif