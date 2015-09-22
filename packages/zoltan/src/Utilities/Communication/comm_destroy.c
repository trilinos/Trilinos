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
#include "comm.h"
#include "zoltan_mem.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

int       Zoltan_Comm_Destroy(
ZOLTAN_COMM_OBJ **plan)		/* communication data structure pointer */
{
    if (*plan == NULL) return ZOLTAN_OK;

    /* Free fields of the communication object. */
    ZOLTAN_FREE(&((*plan)->status));
    ZOLTAN_FREE(&((*plan)->request));
    ZOLTAN_FREE(&((*plan)->sizes));
    ZOLTAN_FREE(&((*plan)->sizes_to));
    ZOLTAN_FREE(&((*plan)->sizes_from));
    ZOLTAN_FREE(&((*plan)->starts_to_ptr));
    ZOLTAN_FREE(&((*plan)->starts_from_ptr));
    ZOLTAN_FREE(&((*plan)->indices_to_ptr));
    ZOLTAN_FREE(&((*plan)->indices_from_ptr));
    ZOLTAN_FREE(&((*plan)->indices_from));
    ZOLTAN_FREE(&((*plan)->indices_to));
    ZOLTAN_FREE(&((*plan)->lengths_from));
    ZOLTAN_FREE(&((*plan)->starts_to));
    ZOLTAN_FREE(&((*plan)->starts_from));
    ZOLTAN_FREE(&((*plan)->lengths_to));
    ZOLTAN_FREE(&((*plan)->procs_from));
    ZOLTAN_FREE(&((*plan)->procs_to));

    /* Free the communication object itself */
    ZOLTAN_FREE(plan);

    return(ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
