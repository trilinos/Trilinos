/*
// @HEADER
// ***********************************************************************
// 
//                Pliris: Parallel Dense Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
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

