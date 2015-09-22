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
/**  @file Trios_trace.h
 *
 *   @brief A simple tracing API.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *   $Revision: 406 $.
 *   $Date: 2005-10-07 15:08:29 -0600 (Fri, 07 Oct 2005) $.
 *
 */

#ifndef _TRACE_H_
#define _TRACE_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

    enum trace_ftype {
        TRACE_UNDEFINED_FTYPE = -1,
        TRACE_SDDF_BINARY,
        TRACE_SDDF_ASCII
    };

#define TRACE_DEFAULT_FTYPE TRACE_SDDF_BINARY;

#if defined(__STDC__) || defined(__cplusplus)

    /* initialization and wrapup functions */
    extern int trace_init(const char *trace_fname, const int trace_ftype);
    extern int trace_reset(
            const char *trace_fname,
            const int trace_ftype,
            const char *enable);
    extern int trace_set_buffer_size(const unsigned long int bufsize);
    extern int trace_set_buffer_size(const unsigned long int bufsize);
    /* extern int trace_set_processor_number(const int pnum); */
    extern int trace_fini(void);


    /* ----------- Trace Group Functions --------------- */

    /**
     * @brief Register a new top-level group for tracing.
     *
     * @param name The name of the group.
     * @param gid  The groupid assigned by the trace library.
     *
     * @return 0 on success, -1 otherwise.
     */
    extern int trace_register_group(const char *name, int *gid);

    /**
     * @brief Register a trace group as a child of an existing
     * group.
     *
     * @param name The name of the new group.
     * @param parent_id The group ID of the parent group.
     * @param gid The value of this group.
     */
    extern int trace_register_subgroup(const char *name, int parent_id, int *gid);

    /**
     * @brief Return the group id for a particular named group.
     *
     * @param name  The name of the group to find.
     *
     * @returns The group ID the group, if found; -1 otherwise.
     */
    extern int trace_get_gid(const char *name);

    /**
     * @brief Return the name of an existing group.
     *
     * @param gid  The group ID of the group to find.
     *
     * @returns The name of the group, if found; NULL otherwise.
     */
    extern const char *trace_get_gname(const int gid);


    /**
     * @brief Test to see if tracing is enabled.
     *
     * @param gid The groupID to check
     */
    extern int tracing_enabled(const int gid);

    /**
     * @brief Enable all groups listed in a comma-separated
     * list of group names.
     *
     * @param csl  The comma-separated list of groups to enable
     */
    extern int trace_enable(const char *csl);

    /**
     * @brief Enable the specified trace group and all its children.
     *
     * @param gid The group ID.
     *
     * @return 0 on success, -1 if the gid is invalid.
     */
    extern int trace_enable_gid(const int gid);

    /**
     * @brief Disable a specified trace group and all its children.
     *
     * @param gid The group ID of the trace group to disable.
     *
     * @return 0 on success, -1 if the gid is invalid.
     */
    extern int trace_disable(const int gid);

    /**
     * @brief Enable every registered trace group.
     *
     * This function enables all trace groups, thus turning on
     * all tracing.  It also sets the default enabled flag to
     * true so that any group registered after this call is
     * also enabled.
     */
    extern int trace_enable_all();

    /**
     * @brief Disable every registered trace group.
     *
     * This function disables all trace groups, thus disabling
     * all tracing.  It also sets the default enabled flag to
     * false so that any group registered after this call is
     * also disabled.
     */
    extern int trace_disable_all();

    /* ----------------- Trace Event Functions ------------ */


    /**
     * @brief Generate a generic trace event.
     *
     * @param gid       The trace group ID.
     * @param eventID   The ID of the trace.
     * @param pid       Process ID.
     * @param data      User-defined data passed in a character string.
     *
     * @return non-zero if successful
     * @return 0 if failure
     */
    extern int trace_event(
            const int gid,
            const int eventID,
            const int pid,
            const char *data);

    /* interval functions */



    /**
     * @brief Generate an interval start event.
     *
     * The \ref trace_start_interval "trace_start_interval()" function generates an
     * interval record that has information about the start, end, and duration of a
     * particular code fragment.
     *
     * We use the generic pablo trace and encode the additional information
     * we want in the data field.  The new data field will be, "interval:$name:duration".
     *
     * Pablo has its own interval records, but they are inadequate because it is
     * difficult to measure concurrent intervals (e.g., in threads).
     * A better (and more efficient) way to do this would be to create our own
     * Pablo record type, but this is a quick "hack" to address our needs.
     *
     * @param gid    The trace group ID.
     * @param pid    The ID of process/thread.
     *
     * @return non-zero if successful
     * @return 0 if failure
     */
    extern int trace_start_interval(
            const int gid,
            const int pid);

    /**
     * @brief Generate an end-interval event.
     *
     * The \ref trace_end_interval "trace_end_interval()" function generates an
     * interval record that has information about the start, end, and duration of a
     * particular code fragment.
     *
     * @param gid       The trace group ID.
     * @param event_id   The ID of the trace (could be the thread ID).
     * @param pid        Process ID.
     * @param data       User-defined data passed in a character
     *                             string.
     *
     * @return non-zero if successful
     * @return 0 if failure
     */
    extern int trace_end_interval(
            const int gid,
            const int event_id,
            const int pid,
            const char *data);


    /**
     * @brief Generate a throughput interval start event.
     *
     * The \ref trace_start_tput_interval "trace_start_tput_interval()" function generates an
     * throughput interval record that has information about the start, end, duration, and
     * number of elements processed during the interval. It is useful for measuring
     * bandwidth, ops/sec, ...
     *
     * We use the generic pablo trace and encode the additional information
     * we want in the data field.  The new data field will be, "interval:tput:$name:duration".
     *
     * Pablo has its own interval records, but they are inadequate because it is
     * difficult to measure concurrent intervals (e.g., in threads).
     * A better (and more efficient) way to do this would be to create our own
     * Pablo record type, but this is a quick "hack" to address our needs.
     *
     * @param interval_id   The ID of the interval (unique for each interval)
     * @param pid   The ID of the process/thread
     *
     * @return non-zero if successful
     * @return 0 if failure
     */
    extern int trace_start_tput_interval(const int gid, const int pid);

    /**
     * @brief Generate an end-interval event.
     *
     * The \ref trace_end_interval "trace_end_interval()" function generates an
     * interval record that has information about the start, end, and duration of a
     * particular code fragment.
     *
     * @param gid  The trace group ID.
     * @param event_id   The ID of the trace (could be the thread ID).
     * @param pid        Process ID.
     * @param data       User-defined data passed in a character
     *                             string.
     *
     * @return non-zero if successful.
     * @return 0 if failure
     */
    extern int trace_end_tput_interval(
            const int gid,
            const int event_id,
            const int pid,
            const long num_processed,
            const char *data);


    /* extern int trace_reset_interval(const int eventID); */

    /* count events */

    /**
     * @brief Get the current value of a trace counter.
     *
     * @param gid  The trace group ID.
     * @param cid   The ID of the counter.
     *
     * @return Value of the counter.
     */
    extern int trace_get_count(
            const int gid,
            const int cid);

    /**
     * @brief Generate an increment-count event.
     *
     * The trace_inc_count_event() function
     * increments a counter associated with the ID of counter and outputs
     * a count event to the tracefile.
     *
     * @param gid  The trace group ID.
     * @param cid  The ID of the counter.
     * @param pid  Process/thread ID.
     * @param data  User-defined data in a character string.
     *
     * @returns non-zero if successful.
     * @returns 0 if failure.
     */
    extern int trace_inc_count(
            const int gid,
            const int cid,
            const int pid,
            const char *data);


    extern int trace_dec_count(
            const int gid,
            const int cid,
            const int pid,
            const char *data);

    extern int trace_set_count(
            const int gid,
            const int cid,
            const int pid,
            const char *data,
            const int newval);

    extern int trace_reset_count(
            const int gid,
            const int cid,
            const int pid,
            const char *data);

    extern int trace_put_all_counts(
            const int gid,
            const int pid,
            const char *data);


#endif /* defined(__STDC__) || defined(__cplusplus) */


#ifdef __cplusplus
}
#endif

#endif /* _TRACE_H_ */
