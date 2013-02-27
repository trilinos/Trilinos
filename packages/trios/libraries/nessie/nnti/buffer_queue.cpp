/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
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
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */

#include "Trios_config.h"

#include "Trios_logger.h"
#include "buffer_queue.h"

#include "assert.h"
#include "string.h"

log_level bq_debug_level = LOG_UNDEFINED;


static NNTI_result_t create_buffer(
        const NNTI_transport_t *trans_hdl,
        const NNTI_buf_ops_t    op,
        const uint32_t          buffer_size,
        NNTI_buffer_t         **buffer)
{
    NNTI_result_t nnti_rc;
    char *b;

    log_debug(bq_debug_level, "enter");

    b=(char *)malloc(buffer_size);
    assert(b);
    memset(b, 0, buffer_size);
    *buffer=(NNTI_buffer_t *)malloc(sizeof(NNTI_buffer_t));
    assert(*buffer);
    nnti_rc=NNTI_register_memory(
            trans_hdl,
            b,
            buffer_size,
            1,
            op,
            NULL,
            *buffer);
    if (nnti_rc != NNTI_OK) {
        log_error(bq_debug_level, "failed registering queue buffer: %d", nnti_rc);
    }

    log_debug(bq_debug_level, "exit");

    return(nnti_rc);
}

static NNTI_result_t destroy_buffer(NNTI_buffer_t **buffer)
{
    NNTI_result_t nnti_rc;
    char *b;

    log_debug(bq_debug_level, "enter");

    b=NNTI_BUFFER_C_POINTER(*buffer);
    nnti_rc=NNTI_unregister_memory(*buffer);
    if (nnti_rc != NNTI_OK) {
        log_error(bq_debug_level, "failed unregistering queue buffer: %d", nnti_rc);
    }
    free(*buffer);
    *buffer=NULL;
    free(b);

    log_debug(bq_debug_level, "exit");

    return(nnti_rc);
}

int trios_buffer_queue_init(
        trios_buffer_queue_t *bq,
        uint32_t              initial_size,
        uint32_t              max_size,
        uint8_t               create_on_fly,
        NNTI_transport_t     *trans_hdl,
        NNTI_buf_ops_t        op,
        uint32_t              buffer_size)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    NNTI_buffer_t *buffer=NULL;

    log_debug(bq_debug_level, "enter");

    nthread_lock_init(&bq->mutex);

    if (nthread_lock(&bq->mutex)) log_warn(bq_debug_level, "failed to get thread lock");

    bq->current_size=0;
    bq->initial_size=initial_size;
    bq->max_size=max_size;
    bq->create_on_fly=create_on_fly;
    bq->trans_hdl=trans_hdl;
    bq->op=op;
    bq->buffer_size=buffer_size;

    for (uint32_t i=0;i<bq->initial_size;i++) {
        log_debug(bq_debug_level, "creating queue buffer");
        nnti_rc=create_buffer(
                bq->trans_hdl,
                bq->op,
                bq->buffer_size,
                &buffer);
        if (nnti_rc==NNTI_OK) {
            log_debug(bq_debug_level, "pushing queue buffer");
            bq->queue.push_back(buffer);
            bq->current_size++;
        } else {
            log_error(bq_debug_level, "failed creating queue buffer: %d", nnti_rc);
            break;
        }
    }
    nthread_unlock(&bq->mutex);

    log_debug(bq_debug_level, "exit");

    return((int)nnti_rc);
}

NNTI_buffer_t *trios_buffer_queue_pop(
        trios_buffer_queue_t *bq)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    NNTI_buffer_t *buffer=NULL;

    log_debug(bq_debug_level, "enter");

    if (nthread_lock(&bq->mutex)) log_warn(bq_debug_level, "failed to get thread lock");
    if (!bq->queue.empty()) {
        log_debug(bq_debug_level, "getting buffer from queue");
        buffer=bq->queue.front();
        bq->queue.pop_front();
    } else {
        if (bq->current_size < bq->max_size) {
            nnti_rc=create_buffer(
                    bq->trans_hdl,
                    bq->op,
                    bq->buffer_size,
                    &buffer);
            if (nnti_rc==NNTI_OK) {
                log_debug(bq_debug_level, "expanding buffer queue");
                bq->current_size++;
            } else {
                log_error(bq_debug_level, "failed creating queue buffer to expand the queue: %d", nnti_rc);
            }
        } else if ((bq->current_size == bq->max_size) && (bq->create_on_fly==TRUE)) {
            nnti_rc=create_buffer(
                    bq->trans_hdl,
                    bq->op,
                    bq->buffer_size,
                    &buffer);
            if (nnti_rc!=NNTI_OK) {
                log_error(bq_debug_level, "failed creating on the fly queue buffer: %d", nnti_rc);
                buffer=NULL;
            } else {
                log_debug(bq_debug_level, "creating on the fly queue buffer");
            }
        }
    }
    nthread_unlock(&bq->mutex);

    log_debug(bq_debug_level, "exit");

    return(buffer);
}

void trios_buffer_queue_push(
        trios_buffer_queue_t *bq,
        NNTI_buffer_t       *buffer)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    log_debug(bq_debug_level, "enter");

    if (nthread_lock(&bq->mutex)) log_warn(bq_debug_level, "failed to get lock");
    if (bq->queue.size() < bq->max_size) {
        log_debug(bq_debug_level, "returning buffer to queue");
        bq->queue.push_front(buffer);
    } else {
        nnti_rc=destroy_buffer(&buffer);
        if (nnti_rc!=NNTI_OK) {
            log_error(bq_debug_level, "failed destroying on the fly queue buffer: %d", nnti_rc);
        } else {
            log_debug(bq_debug_level, "destroying on the fly queue buffer");
        }
    }
    nthread_unlock(&bq->mutex);

    log_debug(bq_debug_level, "exit");
}

int trios_buffer_queue_fini(
        trios_buffer_queue_t *bq)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    NNTI_buffer_t *buffer=NULL;

    log_debug(bq_debug_level, "enter");

    if (nthread_lock(&bq->mutex)) log_warn(bq_debug_level, "failed to get lock");
    if (bq->queue.size() != bq->current_size) {
        log_warn(bq_debug_level, "buffer queue (%p) has missing entries (bq->queue.size(%llu) != bq->current_size(%llu))",
                bq, (uint64_t)bq->queue.size(), (uint64_t)bq->current_size);
    }
    if (bq->queue.size() < bq->initial_size) {
        log_warn(bq_debug_level, "buffer queue (%p) has missing entries (bq->queue.size(%llu) < bq->initial_size(%llu))",
                bq, (uint64_t)bq->queue.size(), (uint64_t)bq->initial_size);
    }
    if (bq->queue.size() > bq->max_size) {
        log_warn(bq_debug_level, "buffer queue (%p) has extra entries (bq->queue.size(%llu) > bq->max_size(%llu))",
                bq, (uint64_t)bq->queue.size(), (uint64_t)bq->max_size);
    }
    while (!bq->queue.empty()) {
        buffer=bq->queue.front();
        bq->queue.pop_front();
        nnti_rc=destroy_buffer(&buffer);
        if (nnti_rc!=NNTI_OK) {
            log_error(bq_debug_level, "failed destroying queue buffer: %d", nnti_rc);
            break;
        }
    }
    nthread_unlock(&bq->mutex);

    nthread_lock_fini(&bq->mutex);

    log_debug(bq_debug_level, "exit");

    return((int)nnti_rc);
}
