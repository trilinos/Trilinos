/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
/**
 *   @file nssi_types.h
 *
 *   @brief Type definitions and method prototypes for
 *   the LWFS.
 *
 *   This file includes the necessary data structures
 *   required by an application that uses the LWFS.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   $Revision: 1064 $
 *   $Date: 2007-01-19 10:05:13 -0700 (Fri, 19 Jan 2007) $
 *
 */

#ifndef TRIOS_NSSI_TYPES_H_
#define TRIOS_NSSI_TYPES_H_

#include "Trios_config.h"

#include <stdint.h>

/* If sockaddr_in is not defined, create some bogus
 * definition. (We don't use it anyway, but sys/socket.h
 * references it.  This works around a Cray bug.
 */
#ifndef HAVE_TRIOS_STRUCT_SOCKADDR_IN
#ifndef HAVE_STRUCT_SOCKADDR_IN
#define HAVE_STRUCT_SOCKADDR_IN
struct sockaddr_in {
    int a;
};
#endif
#endif



typedef struct {
    uint8_t  use_buffer_queue;
    uint32_t buffer_queue_initial_size;
    uint32_t buffer_queue_max_size;
    uint8_t  buffer_queue_create_if_empty;
} nssi_config_t;



/* Treat the rpc-generated includes as system includes to avoid warnings */
#include <Trios_nssi_types_xdr.h>
#include "Trios_nssi_fprint_types.h"

#include <Trios_nnti_xdr.h>
#include "Trios_nnti_fprint_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define UUID_TO_UINT64(u) \
    ((uint64_t)( \
            (*((uint64_t *)&(u[0]))) \
            ^(*((uint64_t *)&(u[8]))) ))

#define UUID_TO_UINT32(u) \
    ((uint32_t)(\
            (*((uint32_t *)&(u[0]))) \
            ^(*((uint32_t *)&(u[4]))) \
            ^(*((uint32_t *)&(u[8]))) \
            ^(*((uint32_t *)&(u[12]))) ))


#if defined(__STDC__) || defined(__cplusplus)

#else /* K&R C */

#endif /* K&R C */

#ifdef __cplusplus
}
#endif

#endif
