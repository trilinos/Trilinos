/*
 * Copyright(C) 1999-2021, 2024, 2025 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 */

/*
 *  NOTES ON MODULE EXMEMY:
 *
 * -2. The big assumption on which this is based is that a FTNINT
 *     which is either an int or a long int, is the same size as
 *     the void* returned by malloc. There are both compile-time
 *     and run-time assertions that attempt to catch this.
 *
 * -1. The information below is the original implementation. At that
 *     time, the only C routines linked with most executables were the
 *     supes library itself and there were no other routines calling
 *     'sbrk' and life was good. This situation has changed with the
 *     addition of the exodusII, netcdf, nemesis, and other libraries.
 *     Some of these are C-based and do 'malloc' and 'free'
 *     calls. Most man pages for 'malloc' and 'sbrk' warn of dire
 *     consequences of mixing calls to 'malloc' and 'sbrk' and we have
 *     had instances of memory corruption under the old model. To
 *     eliminate this undefined and undesirable behavior, the EXMEMY
 *     routine has been rewritten to use 'malloc' and 'free' instead
 *     of 'sbrk' calls. The PAGESIZE manipulations have also been
 *     removed since the underlying 'malloc' implementation is better
 *     at dealing with this than we are (hopefully).
 *
 *     The implementation of 'free' is somewhat of a kluge and
 *     required changes to the calling 'mem_mgr' routines.  Specifically,
 *
 *     1) A few routines were disabled -- mdget which is supposed to
 *     get a large block of memory and parcel it out in smaller chunks.
 *     This minimizes system calls, but is incompatible with 'free' since
 *     memory passed to 'mddel' may not be on a boundary obtained
 *     my malloc.
 *
 *     2) 'mdlong' is mapped to 'realloc'.  A flag is passed into the
 *     exmemy call which tells whether we are freeing memory or
 *     changing size.  The '*memret' argument is used to pass the flag
 *     in to avoid changing the interface to exmemy.
 *
 *     3) 'deferred' mode is eliminated; all memory requests are
 *     satisfied immediately.
 *
 *     4) 'mdcomp' is also eliminated.
 *
 *     The interface has not been changed; the disabled calls are
 *     just noops.
 *
 *     ---Below here is historical notes---
 *
 *  0. On a request for more memory, we ask for it in multiples of
 *     PAGES of memory of size PAGESIZE.  The reason is that on allocation,
 *     the memory manager does keep track the voids in memory.  Unfortunately,
 *     the same cannot be said for the voids created on deallocation, so...
 *
 *  1. The SUPES memory handler, EXMEMY, normally deallocates memory
 *     only if one asks for it a a location that is an exact page multiple.
 *     Further, the memory manager used to try to release memory
 *     from each memory block---this, of course, assumes that the OS
 *     allows for gaps in a users data space.  It so happens that VMS
 *     does.  I've changed the memory manager (cf. mem_mgr/mxgive) so that
 *     it now only tries to release a void if it is at the top of the
 *     program space.  I release as much as I can on any given release
 *     request.
 *
 *  2. Although it's not clear from the comments, or from the documentation
 *     that I've read, on deallocation *memrtn is the amount of memory that
 *     exmemy was *NOT* able to return to the system.  Still unclear?
 *     An example:  Suppose that I tried to return 524 bytes to the system.
 *     Ultimately, I would do this through a call to exmemy with a value
 *     for *memreq of -524 / Numsize (remember we ask to give in numeric
 *     storage units).  Assuming that the amount actually returned
 *     was 512 then *memrtn should be set to 12 divided by the NumSize
 *     for most systems---see the comments below for the Cray strategy.
 *     (memrtn is now the amount of memory released/allocated.  If it is
 *      less than zero, an error occurred)
 */

/*
 * Assumptions:
 * 1. All systems have stdlib.h and define the prototype for malloc
 * 2. All systems are ANSI C compliant (have function prototypes)
 * 3. sizeof(FTNINT) == sizeof(void*)
 */

#include <assert.h>
#include <stdlib.h>
/*
 * Define the Fortran/C interface for the system in question
 * See also itools/config/cf/fortranc.h, which defaults to float, int  --pw
 */

#include <fortranc.h>

#define CT_ASSERT(e) extern char(*ct_assert(void))[sizeof(char[1 - 2 * !(e)])]
/* If the following line causes a compile-time error, then there is a problem
 * which will cause the supes memory manager to not work correctly on this
 * platform. The error will be something similar to:
 *
 * exmemy.c(141): error: the size of an array must be greater than zero
 * CT_ASSERT(sizeof(FTNINT) == sizeof(void*));
 *
 * Contact sierra-help@sandia.gov for assistance.
 */
CT_ASSERT(sizeof(FTNINT) == sizeof(void *));

#if defined(ADDC_)
void exmemy_(FTNINT *memreq, FTNINT *locblk, FTNINT *memrtn)
#else
void exmemy(FTNINT *memreq, FTNINT *locblk, FTNINT *memrtn)
#endif /* ADDC_ */
{
  size_t numbytes;
  size_t NumSize; /* Get multiplier which allows us to */
  /* convert everything in terms of Numeric Storage Units. */
  NumSize = sizeof(FTNREAL) / sizeof(char);

  /* Normal allocation call */
  if (*memreq > 0) { /* Then we need more memory. */
    numbytes = (*memreq) * NumSize;

    size_t *block_location = (size_t *)malloc(numbytes);

    /* printf("%11x %11d\n", block_location, numbytes); */
    /* convert back to numerical storage units */
    *locblk = (FTNINT)((size_t)block_location / NumSize);

    /* See if we have lost any information in the conversion. For example, if
     * the pointer is a 64-bit quantity and a FTNINT is 32-bits, we may not be
     * able to recover the block_location...  This should have been caught in
     * the ct_assert above, but if not, we check again here...
     */
    assert(block_location == (size_t *)((size_t)(*locblk) * NumSize));

    if (block_location == NULL) {
      /*
       * Then the call to 'malloc' has failed, most likely due to
       * asking for more memory than what's available.
       */
      *memrtn = -1;
    }
    else {
      /* The last call for memory was successful! */
      *memrtn = (FTNINT)(numbytes / NumSize);
    }
  }
  else {
    /* Otherwise, if *memrtn is negative, we want to give some back. */

    /*
     * There are a few kludges used here to try to get the memory
     * system to work correctly with malloc/free instead of sbrk.
     * A flag is passed in via the memrtn variable:
     * -999 == free the memory block
     * -998 == realloc
     */

    /* Adjust '*locblk' to be the raw address.  It is passed
     * in in terms of Numeric Storage Units.
     */
    size_t *block_location = (size_t *)((size_t)(*locblk) * NumSize);

    if (*memrtn == -999 || *memreq == 0) {
      /* Handle normal 'free' */
      /* printf("FREE:  %11x %11d\n", block_location, *memreq); */
      free(block_location);
      *memrtn = -(*memreq);
    }
    else if (*memrtn == -998) {
      /* realloc call to shrink/grow memory */
      numbytes = -(*memreq) * NumSize;
      /* printf("PRE:  %11x %11d (realloc)\n", block_location, numbytes); */
      size_t *new_location = realloc(block_location, numbytes);
      /* printf("POST: %11x %11d\n", new_location, numbytes); */

      if (new_location == NULL && *memreq > 0) {
        /*
         * Then the call to 'realloc' has failed, most likely due to
         * asking for more memory than what's available.
         */
        *memrtn = -1;
      }
      else {
        *memrtn = -(*memreq);

        /* convert back to numerical storage units */
        *locblk = (FTNINT)((size_t)new_location / NumSize);
      }
      /* If '*memrtn' flag was not set to -999 or -998,
       * ignore free request */
    }
  }
}

#ifdef TEST
#include <stdio.h>

int main(int argc, char **argv)
{
  FTNINT locblk;
  FTNINT memrtn;
  size_t first = 0;
  size_t last  = 0;

  assert(sizeof(FTNREAL) == sizeof(FTNINT));

  printf("Size of FORTRAN int/real  = %lu bytes\n", sizeof(FTNINT));

  FTNINT memreq = 10000000; /* 10 Million words */
  for (int i = 0; i < 1000; i++) {
#if defined(ADDC_)
    exmemy_(&memreq, &locblk, &memrtn);
#else
    exmemy(&memreq, &locblk, &memrtn);
#endif /* ADDC_ */
    printf("%4d: %11x %11d\n", i, locblk, memrtn);

    if (memrtn < 0) {
      FTNINT allocated;
      FTNINT overhead;
      allocated = i * (memreq / 1000000);
      overhead  = (last - first - (i - 1) * memreq) * sizeof(FTNINT);
      printf("Total Allocation = %d MegaWords, %d MB\n", allocated, allocated * sizeof(FTNINT));

      /* NOTE: 'last' is pointer to beginning of last block, not end */
      printf("Overhead = %d bytes, %d bytes/allocation\n", overhead, overhead / (i - 1));

      break;
    }
    else {
      last = locblk;
    }
    if (i == 0) {
      first = locblk;
    }
  }
}
#endif /* TEST */
