/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include<stdio.h>
#include<stdlib.h>
#include "malloc.h"

#define ALIGNED_TO   ((unsigned)32) /* number of bytes to align to */
#define EXTRA        (ALIGNED_TO)
#define MASK         (~(EXTRA-1))

/*
** On the Intel Paragon we are CHAMELEON aligned. No problem.
*/

void *MALLOC_8(int size)
{
    return((void *) malloc(size));
}  /* end of MALLOC_8() */



/* --------------------------------------------------------------- */
/* return a void pointer that is aligned along ALIGNED_TO boundary */
/* heaven help you if you free one of these...                     */
/* --------------------------------------------------------------- */

static long ALIGN_64 = ~0x3f;

void *MALLOC_64(len)
size_t len;
{
  void *p;

  p=malloc(len+63);

  if(p==NULL) {
    return(NULL);
  }
  else{
     p = (void *)(((unsigned long)p+63) & ALIGN_64);
  }
  return p;
}


void *MALLOC_32(len)
size_t len;
{
  void *p;

  p=malloc(len+EXTRA);
  if(p==NULL) {
    return(NULL);
  }
  else{
    return((void *)((((unsigned long)p) + EXTRA) & MASK));
  }
}

