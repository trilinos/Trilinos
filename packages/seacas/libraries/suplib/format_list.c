#include <stdint.h>
#include <stdio.h>
#if defined(_WIN32) && defined(_MSC_VER) && _MSC_VER < 1900
#define PRId64 "I64d"
#else
#include <inttypes.h>
#endif

#if defined(ADDC_)
void formlist_(int64_t *ids, int64_t *len)
#else
void formlist(int64_t *ids, int64_t *len)
#endif
{
  const char *rng_sep = "..";
  const char *seq_sep = ", ";

  printf("\t\t");
  int64_t num = 0;
  while (num < *len) {
    if (num == 0) {
      printf("%" PRId64 "", ids[num]);
    }
    else {
      printf("%s%" PRId64 "", seq_sep, ids[num]);
    }

    int64_t begin    = ids[num]; // first id in range of 1 or more ids
    int64_t previous = ids[num]; // last id in range of 1 or more ids
    // Gather a range or single value... (begin .. previous)
    while (previous == ids[num] && ++num < *len && ids[num] == previous + 1) {
      previous++;
    }

    if (begin != previous) {
      printf("%s%" PRId64 "", previous == begin + 1 ? seq_sep : rng_sep, previous);
    }
  }
  printf("\n\n");
}
