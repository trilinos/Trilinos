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
/*-------------------------------------------------------------------------*/
/**  @file buffer_queue.h
 *
 *   @brief API for a circular list of NNTI_buffer_t.
 *
 *   @author Todd Kordenbrock (thkorde\@sandia.gov).
 *
 */

#ifndef _TRIOS_BUFFER_QUEUE_H_
#define _TRIOS_BUFFER_QUEUE_H_


#include "Trios_config.h"

#include "Trios_threads.h"

#include "Trios_nnti.h"

#include <deque>

typedef std::deque<NNTI_buffer_t *>  buffer_queue_t;

typedef struct trios_buffer_queue {
    nthread_lock_t    mutex;
    buffer_queue_t    queue;
    uint32_t          current_size;
    uint32_t          initial_size;
    uint32_t          max_size;
    uint8_t           create_on_fly;
    NNTI_transport_t *trans_hdl;
    NNTI_buf_ops_t    op;
    uint32_t          buffer_size;
} trios_buffer_queue_t;


#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)

    extern int trios_buffer_queue_init(
            trios_buffer_queue_t *bq,
            uint32_t              initial_size,
            uint32_t              max_size,
            uint8_t               create_on_fly,
            NNTI_transport_t     *trans_hdl,
            NNTI_buf_ops_t        op,
            uint32_t              buffer_size);
    extern NNTI_buffer_t *trios_buffer_queue_pop(
            trios_buffer_queue_t *bq);
    extern void trios_buffer_queue_push(
            trios_buffer_queue_t *bq,
            NNTI_buffer_t       *buffer);
    extern int trios_buffer_queue_fini(
            trios_buffer_queue_t *bq);

#endif


#ifdef __cplusplus
}
#endif

#endif
