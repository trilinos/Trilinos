/* 
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
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

/* liballoc - flex run-time memory allocation */

/* $Header: /usr/sharelan/cvs/SEACAS/prepost/aprepro/liballoc.c,v 1.6 2006/02/16 20:29:26 gdsjaar Exp $ */

#include <stdio.h>
#include <stdlib.h>

void *yy_flex_alloc(int size);
void *yy_flex_realloc(void *ptr, int size);
void yy_flex_free(void *ptr);
void *yy_flex_xmalloc(int size);

void *yy_flex_alloc( int size )
{
  return (void *) malloc( (unsigned) size );
}

void *yy_flex_realloc( void *ptr, int size )
{
  return (void *) realloc( ptr, (unsigned) size );
}

void yy_flex_free( void *ptr )
{
  free( ptr );
}

/* The following is only used by bison/alloca. */
void *yy_flex_xmalloc( int size )
{
  void *result = yy_flex_alloc( size );

  if ( ! result  )
    {
      fprintf( stderr,
	       "flex memory allocation failed in yy_flex_xmalloc()" );
      exit( 1 );
    }
  return result;
}
