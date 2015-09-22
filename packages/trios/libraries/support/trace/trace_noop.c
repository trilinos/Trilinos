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
/**  @file trace_noop.c
 *
 *   @brief Noop implementation of trace API.
 *
 *   @author Ron Oldfield (raoldfi@sandia.gov)
 *   @version $Revision: 406 $
 *   @date $Date: 2005-10-07 15:08:29 -0600 (Fri, 07 Oct 2005) $
 *
 */

#include "Trios_trace.h"

/**
 * A noop function that always returns 0.  This function
 * replaces all trace functions when TRACING_ENABLED is
 * not set at compile time.
 */
int trace_noop(int a, ...)
{
    return 0;
}

int trace_init(const char *trace_fname, const int trace_ftype)
{
    return 0;
}

int trace_reset(
        const char *trace_fname,
        const int trace_ftype,
        const char *enable)
{
    return 0;
}

int trace_set_buffer_size(const unsigned long int bufsize)
{
    return 0;
}

int trace_fini(void)
{
    return 0;
}

int trace_register_group(const char *name, int *gid)
{
    return 0;
}

int trace_register_subgroup(const char *name, int parentid, int *gid)
{
    return 0;
}

int trace_get_gid(const char *name)
{
    return 0;
}

const char *trace_get_gname(const int gid)
{
    return 0;
}

int tracing_enabled(const int gid)
{
    return 0;
}

int trace_enable(const char *name)
{
    return 0;
}

int trace_enable_gid(const int gid)
{
    return 0;
}

int trace_disable(const int gid)
{
    return 0;
}

int trace_enable_all()
{
    return 0;
}

int trace_disable_all()
{
    return 0;
}

int trace_event(
        const int gid,
        const int eventID,
        const int pid,
        const char *data)
{
    return 0;
}

int trace_start_interval(
        const int gid,
        const int pid)
{
    return 0;
}

int trace_end_interval(
        const int gid,
        const int event_id,
        const int pid,
        const char *data)
{
    return 0;
}

int trace_start_tput_interval(const int gid, const int pid)
{
    return 0;
}

int trace_end_tput_interval(
        const int gid,
        const int event_id,
        const int pid,
        const long num_processed,
        const char *data)
{
    return 0;
}

int trace_get_count(
        const int gid,
        const int cid)
{
    return 0;
}

int trace_inc_count(
        const int gid,
        const int cid,
        const int pid,
        const char *data)
{
    return 0;
}


int trace_dec_count(
        const int gid,
        const int cid,
        const int pid,
        const char *data)
{
    return 0;
}

int trace_set_count(
        const int gid,
        const int cid,
        const int pid,
        const char *data,
        const int newval)
{
    return 0;
}

int trace_reset_count(
        const int gid,
        const int cid,
        const int pid,
        const char *data)
{
    return 0;
}

int trace_put_all_counts(
        const int gid,
        const int pid,
        const char *data)
{
    return 0;
}
