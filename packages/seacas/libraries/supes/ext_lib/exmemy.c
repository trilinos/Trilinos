/*
 * Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *           
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *                         
 * * Neither the name of Sandia Corporation nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *                                                 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/*
 * $Id: exmemy.c,v 1.19 2008/12/17 22:47:19 gdsjaar Exp $
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

#define ct_assert(e) extern char (*ct_assert(void)) [sizeof(char[1 - 2*!(e)])]
/* If the following line causes a compile-time error, then there is a problem
 * which will cause the supes memory manager to not work correctly on this
 * platform. The error will be something similar to:
 *
 * exmemy.c(141): error: the size of an array must be greater than zero
 * ct_assert(sizeof(FTNINT) == sizeof(void*));
 *
 * Contact Greg Sjaardema, gdsjaar@sandia.gov for asisstance.
 */
ct_assert(sizeof(FTNINT) == sizeof(void*));

#if defined(ADDC_)
void exmemy_(FTNINT *memreq, FTNINT *locblk, FTNINT *memrtn)
#else
void exmemy( FTNINT *memreq, FTNINT *locblk, FTNINT *memrtn)
#endif /* ADDC_ */
{
  size_t numbytes;
  size_t NumSize;	/* Get multiplier which allows us to */
  /* convert everything in terms of Numeric Storage Units. */
  size_t *block_location;
  size_t *new_location;
  
  NumSize = sizeof( FTNREAL ) / sizeof( char ); 

  /* Normal allocation call */
  if(*memreq > 0) {		/* Then we need more memory. */
    numbytes = (*memreq) * NumSize;
    
    block_location = (size_t*)malloc(numbytes);
    
    /* printf("%11x %11d\n", block_location, numbytes); */
    /* convert back to numerical storage units */
    *locblk = (FTNINT)((size_t)block_location / NumSize);
    
    /* See if we have lost any information in the conversion. For example, if
     * the pointer is a 64-bit quantity and a FTNINT is 32-bits, we may not be
     * able to recover the block_location...  This should have been caught in
     * the ct_assert above, but if not, we check again here...
     */
    assert (block_location == (size_t*)((size_t)(*locblk) * NumSize));

    if (block_location == 0) {
      /*
       * Then the call to 'malloc' has failed, most likely due to
       * asking for more memory than what's available.
       */
      *memrtn = -1;
    } else {
      /* The last call for memory was successful! */
      *memrtn = (FTNINT)(numbytes / NumSize); 
    }
  } else {
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
    block_location = (size_t*)((size_t)(*locblk) * NumSize);

    if (*memrtn == -999 || *memreq == 0) {
      /* Handle normal 'free' */
      /* printf("FREE:  %11x %11d\n", block_location, *memreq); */
      free(block_location);
      *memrtn = -(*memreq);
	
    } else if (*memrtn == -998) {
      /* realloc call to shrink/grow memory */
      numbytes = -(*memreq) * NumSize;
      /* printf("PRE:  %11x %11d (realloc)\n", block_location, numbytes); */
      new_location = realloc(block_location, numbytes);
      /* printf("POST: %11x %11d\n", new_location, numbytes); */

      if (new_location == 0 && *memreq > 0) {
	/*
	 * Then the call to 'realloc' has failed, most likely due to
	 * asking for more memory than what's available.
	 */
	*memrtn = -1;
      } else {
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
#include <assert.h>
#include <stdio.h>

int main(int argc, char **argv)
{
  int i;
  FTNINT memreq;
  FTNINT locblk;
  FTNINT memrtn;
  size_t first;
  size_t last;

  assert(sizeof(FTNREAL) == sizeof(FTNINT));

  printf("Size of FORTRAN int/real  = %d bytes\n", sizeof(FTNINT));

  memreq = 10000000; /* 10 Million words */
  for (i=0; i < 1000; i++) {
#if defined(ADDC_)
    exmemy_(&memreq, &locblk, &memrtn);
#else
    exmemy(&memreq, &locblk, &memrtn);
#endif /* ADDC_ */
    printf("%4d: %11x %11d\n", i, locblk, memrtn);

    if (memrtn < 0) {
      FTNINT allocated;
      FTNINT overhead;
      allocated = i*(memreq/1000000);
      overhead  = (last-first-(i-1)*memreq)*sizeof(FTNINT);
      printf("Total Allocation = %d MegaWords, %d MB\n",
	     allocated, allocated*sizeof(FTNINT));

      /* NOTE: 'last' is pointer to beginning of last block, not end */
      printf("Overhead = %d bytes, %d bytes/allocation\n",
	     overhead, overhead/(i-1));
	     
      break;
    } else {
      last = locblk;
    }
    if (i==0) {
      first = locblk;
    }
  }
}
#endif /* TEST */
