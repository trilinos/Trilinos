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
/*-------------------------------------------------------------------------*/
/**  @file Trios_ring_buffer.h
 *
 *   @brief API for a circular list of NNTI_buffer_t.
 *
 *   @author Todd Kordenbrock (thkorde\@sandia.gov).
 *
 */

#ifndef _TRIOS_RING_BUFFER_H_
#define _TRIOS_RING_BUFFER_H_


#include "Trios_config.h"

#include "Trios_threads.h"

#include "nnti.h"

#include <deque>

typedef std::deque<NNTI_buffer_t *>  ring_buffer_queue_t;

typedef struct trios_ring_buffer {
    nthread_mutex_t         mutex;
    ring_buffer_queue_t     ring;
    uint32_t                ring_size;
    const NNTI_transport_t *trans_hdl;
} trios_ring_buffer_t;


#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)

    extern int trios_ring_buffer_init(
            trios_ring_buffer_t    *rb,
            const uint32_t          ring_size,
            const NNTI_transport_t *trans_hdl,
            const NNTI_buf_ops_t    op,
            const uint32_t          buffer_size);
    extern NNTI_buffer_t *trios_ring_buffer_pop(
            trios_ring_buffer_t *rb);
    extern void trios_ring_buffer_push(
            trios_ring_buffer_t *rb,
            NNTI_buffer_t       *ring_buffer);
    extern int trios_ring_buffer_fini(
            trios_ring_buffer_t *rb);

#endif


#ifdef __cplusplus
}
#endif

#endif
