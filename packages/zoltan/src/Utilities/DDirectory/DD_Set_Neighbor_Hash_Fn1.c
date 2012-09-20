/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"



#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

struct dd_nh1_struct {
 int max_gid ;
 int groupsize ;
};

static unsigned int dd_nh1 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc, struct dd_nh1_struct* hashdata, ZOLTAN_HASH_FN *fn) ;



/*************  Zoltan_DD_Set_Hash_Fn1() ***********************/
/*
**  These routines associate the first n=groupsize GIDs to proc 0, the
**  next n to proc 1, etc.  It assumes the GIDs are consecutive numbers.
**  It assumes that GIDs primarily stay near their original owner. The
**  GID length is assumed to be 1. GIDs outside of range are evenly
**  distributed among the processors via modulo(nproc).
*/


int Zoltan_DD_Set_Neighbor_Hash_Fn1 (
 Zoltan_DD_Directory *dd,          /* directory state information */
 int size)                         /* number of reserved GIDs per CPU */
   {
   char *yo = "Zoltan_DD_Set_Neighbor_Hash_Fn1";
   struct dd_nh1_struct *hashdata;

   if (dd == NULL) {
     ZOLTAN_PRINT_ERROR (0, yo, "NULL DDirectory pointer");
     return ZOLTAN_FATAL;
   }
   if (size < 1) {
     ZOLTAN_PRINT_WARN (0, yo, "Invalid input argument; size < 1");
     return ZOLTAN_WARN;
   }

   hashdata = (struct dd_nh1_struct*) ZOLTAN_MALLOC(sizeof(struct dd_nh1_struct));
   if (hashdata == NULL) {
     ZOLTAN_PRINT_ERROR (0, yo, "Memory error");
     return ZOLTAN_FATAL;
   }

   hashdata->groupsize   = size;
   dd->hash    = (DD_Hash_fn*) &dd_nh1;
   dd->hashfn  = NULL;
   dd->hashdata = hashdata;
   dd->cleanup = (DD_Cleanup_fn*) &Zoltan_DD_default_cleanup;
   hashdata->max_gid = size * dd->nproc;     /* larger GIDs out of range */

   return ZOLTAN_OK;
   }



static unsigned int dd_nh1 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc, struct dd_nh1_struct* hashdata, ZOLTAN_HASH_FN *fn)
   {
   int id = (signed) *gid;
   return (unsigned int) ((id < hashdata->max_gid) ? (id / hashdata->groupsize)
                                                   : (id % nproc));
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
