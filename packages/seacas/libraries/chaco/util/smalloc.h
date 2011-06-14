#include <stddef.h>

/* Safe version of malloc.  Does not initialize memory .*/
extern void *smalloc(size_t n);

/* Safe version of malloc.  Does not initialize memory .*/
/* Returns instead of dying if it fails. */
extern void *smalloc_ret(size_t n);

/* Safe version of realloc */
extern void *srealloc(void *ptr, size_t n);

/* Safe version of realloc */
/* Returns instead of dying if it fails. */
extern void *srealloc_ret(void *ptr, size_t n);

/* Safe version of free. */
extern void sfree(void *ptr);

void smalloc_stats(void);




