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
 * This API is modeled after pthreads.  pthreads is a POSIX API.  TRIOS runs in non-POSIX
 * environments such as Cray MTA.  These systems don't have pthreads, but may have other
 * threading libs.  This is lib creates a single API for TRIOS internal usage.
 */

#include "Trios_config.h"

#include "Trios_ring_buffer.h"

#include "assert.h"

int trios_ring_buffer_init(
        trios_ring_buffer_t    *rb,
        const uint32_t          ring_size,
        const NNTI_transport_t *trans_hdl,
        const NNTI_buf_ops_t    op,
        const uint32_t          buffer_size)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    int i;
    char *b=NULL;
    NNTI_buffer_t *ring_buffer=NULL;

    nthread_mutex_init(&rb->mutex, NTHREAD_MUTEX_NORMAL);

    nthread_lock(&rb->mutex);
    rb->ring_size=ring_size;
    rb->trans_hdl=trans_hdl;
    for (i=0;i<ring_size;i++) {
        b=(char *)malloc(buffer_size);
        assert(b);
        memset(b, 0, buffer_size);
        ring_buffer=(NNTI_buffer_t *)malloc(sizeof(NNTI_buffer_t));
        assert(ring_buffer);
        nnti_rc=NNTI_register_memory(
                trans_hdl,
                b,
                buffer_size,
                1,
                op,
                NULL,
                ring_buffer);
        rb->ring.push_back(ring_buffer);
    }
    nthread_unlock(&rb->mutex);

    return((int)nnti_rc);
}

NNTI_buffer_t *trios_ring_buffer_pop(
        trios_ring_buffer_t *rb)
{
    NNTI_buffer_t *ring_buffer=NULL;

    nthread_lock(&rb->mutex);
    if (!rb->ring.empty()) {
        ring_buffer=rb->ring.front();
        rb->ring.pop_front();
    }
    nthread_unlock(&rb->mutex);

    return(ring_buffer);
}

void trios_ring_buffer_push(
        trios_ring_buffer_t *rb,
        NNTI_buffer_t       *ring_buffer)
{
    nthread_lock(&rb->mutex);
    rb->ring.push_back(ring_buffer);
    nthread_unlock(&rb->mutex);
}

int trios_ring_buffer_fini(
        trios_ring_buffer_t *rb)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    char *b=NULL;
    NNTI_buffer_t *ring_buffer=NULL;

    nthread_lock(&rb->mutex);
    while (!rb->ring.empty()) {
        ring_buffer=rb->ring.front();
        b=NNTI_BUFFER_C_POINTER(ring_buffer);
        nnti_rc=NNTI_unregister_memory(ring_buffer);
        rb->ring.pop_front();
        free(ring_buffer);
        free(b);
    }
    nthread_unlock(&rb->mutex);

    return((int)nnti_rc);
}
