#ifndef SUPLIB_C_COPY_STRING
#define SUPLIB_C_COPY_STRING
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
char *copy_string(char *dest, char const *source, size_t elements);
#ifdef __cplusplus
} /* close brackets on extern "C" declaration */
#endif
#endif
