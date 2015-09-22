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
/**  @file trace.c
 *
 *   @brief A simple library for generating timing traces.
 *
 *   @author Ron Oldfield (raoldfi@sandia.gov)
 *   @version $Revision: 406 $
 *   @date $Date: 2005-10-07 15:08:29 -0600 (Fri, 07 Oct 2005) $
 *
 */

#include "Trios_config.h"

#include "Trios_trace.h"


#include "TraceFile.h"

#ifndef TRIOS_ENABLE_TRACING
#include "trace_noop.c"
#else

#ifdef HAVE_TRIOS_PABLO
#include "SDDF.h"
#endif

#include <stack>
#include <map>

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include "Trios_trace.h"
#include "Trios_timer.h"
#include "Trios_logger.h"
#include "Trios_threads.h"

#include "TraceGroup.h"

#include <map>
#include <string>

/* The tracefile class */
static TraceFile *tracefile = 0;

/* hash tables for managing the trace groups */
static std::map<string, TraceGroup *> gname_map;
typedef std::map<string, TraceGroup *>::iterator gname_iterator_t;
typedef std::pair<string, TraceGroup *> gname_pair_t;
static nthread_lock_t gname_mutex;

static std::map<int, TraceGroup *> gid_map;
typedef std::map<int, TraceGroup *>::iterator gid_iterator_t;
typedef std::pair<int, TraceGroup *> gid_pair_t;
static nthread_lock_t gid_mutex;

/* --- Hashtable that stores the intervals --- */
static std::map<int, stack<double> > interval_map;
typedef std::map<int, stack<double> >::iterator interval_iterator_t;
typedef std::pair<int,stack<double> > interval_pair_t;
static nthread_lock_t interval_mutex;

/* --- Hashtable to store counter information --- */
static std::map<int, int> count_map;
typedef std::map<int, int>::iterator count_iterator_t;
typedef std::pair<int,int> count_pair_t;

static nthread_lock_t count_mutex;


static bool default_enable_flag = false;

#ifndef HAVE_TRIOS_PABLO
#warning Building trace library without Pablo
#endif


log_level trace_debug_level = LOG_UNDEFINED;

/* ------------------ PUBLIC API ----------------------- */


/*------------ Initialization and wrapup functions ------------*/



int trace_set_buffer_size(const unsigned long int bufsize)
{
    int rc = 1;
    if (tracefile) {
        rc = tracefile->set_buffer_size(bufsize);
    }
    return rc;
}


/**
 * @brief Register a new top-level group for tracing.
 */

int trace_register_group(const char *name, int *gid)
{
    static volatile int groupcount = 0;
    gname_iterator_t name_iter;

    /* Find out of the group already exists.  If so, return gid */
    nthread_lock(&gname_mutex);
    if ((name_iter = gname_map.find(name)) != gname_map.end()) {
        *gid = name_iter->second->gid;
        nthread_unlock(&gname_mutex);
        return 0;
    }

    /* increment the groupcount and assign the id of this new group */
    *gid = ++groupcount;

//    fprintf(logger_get_file(), "registering trace group %s (gid=%d)\n",
//            name, *gid);

    /* create the new group */
    TraceGroup *group = new TraceGroup(name, *gid, default_enable_flag);

    /* Add the group to the gname_map and the gid_map */
    gname_map.insert(pair<const char *, TraceGroup *>(group->name, group));
    nthread_unlock(&gname_mutex);


    nthread_lock(&gid_mutex);
    gid_map.insert(pair<int, TraceGroup *>(group->gid, group));
    nthread_unlock(&gid_mutex);

    return 0;
}

/**
 * @brief Register a trace group as a child of an existing
 * group.
 *
 * @param name The name of the new group.
 * @param parent_id The group ID of the parent group.
 * @param gid The value of this group.
 */
int trace_register_subgroup(
        const char *name,
        int parent_id,
        int *gid)
{
    int rc;
    gid_iterator_t parent_iter, child_iter;
    TraceGroup *child = NULL;

    rc = trace_register_group(name, gid);

    nthread_lock(&gid_mutex);
    child_iter = gid_map.find(*gid);
    child = child_iter->second;

    /* if the partent group exists, add this node as a child */
    parent_iter = gid_map.find(parent_id);
    if (parent_iter != gid_map.end()) {
        parent_iter->second->add_child(child);
    }
    nthread_unlock(&gid_mutex);

    return rc;
}


/**
 * @brief Return the group id for a particular named group.
 *
 * @param name  The name of the group to find.
 *
 * @returns The group ID the group, if found; -1 otherwise.
 */

int trace_get_gid(const char *name)
{
    int gid;
    gname_iterator_t name_iter;

    /* lookup the group by calling trace_register. If the
     * group is already registered it just returns the gid
     */
    trace_register_group(name, &gid);

    return gid;
}

/**
 * @brief Return true if the group is enabled.
 *
 * @param gid The groupID to check
 */

int tracing_enabled(const int gid)
{
    int result = false;

    /* first check that they initialized the trace environment */
    if (tracefile) {

        /* find the gid */
        nthread_lock(&gid_mutex);
        gid_iterator_t gid_iter = gid_map.find(gid);

        if (gid_iter == gid_map.end()) {
            result = false;
        }
        else {
            result = gid_iter->second->enabled;
        }
        nthread_unlock(&gid_mutex);
    }

    return result;
}

/**
 * @brief Return the name of an existing group.
 *
 * @param gid  The group ID of the group to find.
 *
 * @returns The name of the group, if found; NULL otherwise.
 */

const char *trace_get_gname(const int gid)
{
    gid_iterator_t gid_iter;

    /* lookup the group by gid */
    gid_iter = gid_map.find(gid);

    if (gid_iter == gid_map.end()) {
        return 0;
    }
    else {
        return gid_iter->second->name;
    }
}



/**
 * Enable a single TraceGroup that is part of a <int,TraceGroup *> pair.
 */
static void enable_group(pair< int, TraceGroup *> node)
{
    node.second->enabled = true;
}


/**
 * @brief Enable all groups listed in a comma-separated
 * list of group names.
 *
 * @param csl  The comma-separated list of groups to enable
 *
 * @return 0 if successfull, -1 on error.
 */
int trace_enable(const char *csl)
{
    int rc = 0;

    if (csl) {

        const char delim[] = ",";
        char *token;
        char *saveptr;

        token = strtok_r((char *)csl, delim, &saveptr);
        while (token != NULL) {
            int gid;

            if (strcasecmp(token, "all") == 0) {
                fprintf(logger_get_file(),
                        "enabling all registered trace groups\n" );
                trace_enable_all();
            }
            else {
                gid = trace_get_gid(token);
                fprintf(logger_get_file(), "enabling trace "
                        "group %s (gid=%d)\n",
                        token, gid);
                trace_enable_gid(gid);
            }

            token = strtok_r(NULL, delim, &saveptr);
        }
    }
    else {
        fprintf(logger_get_file(),"tracing disabled\n");
    }

    return rc;
}


/**
 * @brief Set the enable flag for a group and all its children.
 *
 * @param gid The group ID.
 *
 * @return 0 on success, -1 if the groupID is invalid.
 */
int trace_enable_gid(const int gid)
{
    int result = 0;
    gid_iterator_t gid_iter;

    /* lookup the group by its gid */
    nthread_lock(&gid_mutex);
    gid_iter = gid_map.find(gid);

    /* case for bad gid */
    if (gid_iter == gid_map.end()) {
        result = -1;
    }
    else {
        enable_group(*gid_iter);
    }
    nthread_unlock(&gid_mutex);


    return result;
}


/**
 * @brief Enable every registered trace group.
 */
int trace_enable_all()
{
    gid_iterator_t iter;
    /* Disable each trace group in the gid map */
    nthread_lock(&gid_mutex);

    for (iter=gid_map.begin(); iter!=gid_map.end(); iter++) {
        iter->second->enabled = true;
    }

    //for_each(gid_map.begin(), gid_map.end(), enable_group);

    /* set the default enable flag to false (for future groups) */
    default_enable_flag = true;
    nthread_unlock(&gid_mutex);

    return 0;
}



/**
 * @brief Disable a trace group.
 *
 * @param gid The group ID of the trace group to disable.
 *
 * @return 0 on success, -1 if the gid is invalid.
 */
int trace_disable(const int gid)
{
    gid_iterator_t gid_iter;
    int result = 0;

    /* lookup the group by its gid */
    nthread_lock(&gid_mutex);
    gid_iter = gid_map.find(gid);

    /* case for bad gid */
    if (gid_iter == gid_map.end()) {
        result = -1;
    }
    else {
        /* Call the enable() function on the TraceGroup */
        gid_iter->second->disable();
    }
    nthread_unlock(&gid_mutex);

    return result;
}

/**
 * Disable a single TraceGroup that is part of a <int,TraceGroup *> pair.
 */
static void disable_group(pair<int, TraceGroup *> node)
{
    node.second->enabled = false;
}


/**
 * @brief Disable all trace groups.
 */
int trace_disable_all()
{
    /* Disable each trace group in the gid map */
    nthread_lock(&gid_mutex);
    for_each(gid_map.begin(), gid_map.end(), disable_group);

    /* set the default enable flag to false (for future groups) */
    default_enable_flag = false;
    nthread_unlock(&gid_mutex);

    return 0;
}



/**
 * @brief Reset the tracing API.
 *
 * This means we close any existing files, reset all counters,
 * and open a new trace file.
 */
int trace_reset(const char *fname, const int ftype, const char *csl)
{
    int rc = 0;
    count_iterator_t iter;
    log_level debug_level = LOG_ALL;

    log_debug(debug_level, "trace_reset start: fname=%s, ftype=%d, enable=%s",
            fname, ftype, csl);

    /* Close the trace file */
    if (tracefile) {
        fprintf(logger_get_file(),"closing tracefile\n");
        delete tracefile;
        tracefile = 0;
    }

    /* Init the trace file */
    log_debug(trace_debug_level, "init tracing");
    rc = trace_init(fname, ftype);
    if (rc != 0) {
        return rc;
    }

    /* Reset all counters */
    log_debug(trace_debug_level, "reset counters");
    nthread_lock(&count_mutex);
    for (iter = count_map.begin(); iter != count_map.end(); iter++) {
        iter->second = 0;
    }
    nthread_unlock(&count_mutex);

    /* Disable all trace groups */
    log_debug(trace_debug_level, "disable trace groups");
    rc = trace_disable_all();
    if (rc != 0)
        return rc;

    /* Enable the groups in the csl */
    log_debug(trace_debug_level, "enable %s",csl);
    rc = trace_enable(csl);
    if (rc != 0) {
        return rc;
    }

    log_debug(trace_debug_level, "trace_reset exit");
    return rc;
}


/**
 * @brief Initialize the tracing environment.
 *
 * @param fname  The name of the trace file.
 * @param ftype  The type of tracing (currently unused).
 */
int trace_init(const char *fname, const int ftype)
{
    int rc = 0;
    log_level debug_level = trace_debug_level;

    log_debug(debug_level, "trace_init(fnam=%s)",fname);
    log_debug(debug_level, "trace_init(ftype=%d)",ftype);

    log_debug(trace_debug_level, "trace_init(%s,%d)",fname,ftype);

    nthread_lock_init(&gname_mutex);
    nthread_lock_init(&gid_mutex);
    nthread_lock_init(&interval_mutex);
    nthread_lock_init(&count_mutex);


    if (!tracefile && fname) {
        int fd;

        /* make sure we can even create the tracefile */
        fd = open(fname, O_CREAT | O_RDWR);
        if (fd < 0) {
            log_warn(trace_debug_level, "Could not open trace file %s... tracing disabled",
                    fname);
            /* do nothing (disable tracing) if we cannot create the file */
            return rc;
        }
        else {
            close (fd);
            remove(fname);
        }

        switch (ftype) {

#ifdef HAVE_TRIOS_PABLO

        case TRACE_SDDF_BINARY:
            log_debug(trace_debug_level, "SDDF output to \"%s\"",
                    fname);
            tracefile = new SDDF(fname, 0);
            break;

        case TRACE_SDDF_ASCII:
            tracefile = new SDDF(fname, 1);
            break;
#endif

        default:
            fprintf(stderr, "Trace file type %d not recognized\n",
                    ftype);
            rc = -1;
        }
    }

    log_debug(trace_debug_level, "trace_init exit");

    return rc;
}

static void delete_group(gname_pair_t entry)
{
    TraceGroup *grp = entry.second;
    delete grp;
}

/**
 * @brief Finalize the tracing environment.
 */
int trace_fini()
{
    int rc = 0;

    if (tracefile) {
        delete tracefile;
        tracefile = 0;
    }

    /* delete trace groups */
    nthread_lock(&gname_mutex);
    for_each(gname_map.begin(), gname_map.end(), delete_group);


    gname_map.clear();
    nthread_unlock(&gname_mutex);

    nthread_lock(&gid_mutex);
    gid_map.clear();
    nthread_unlock(&gid_mutex);

    nthread_lock(&interval_mutex);
    interval_map.clear();
    nthread_unlock(&interval_mutex);

    nthread_lock_fini(&gname_mutex);
    nthread_lock_fini(&gid_mutex);
    nthread_lock_fini(&interval_mutex);
    nthread_lock_fini(&count_mutex);

    return rc;
}



/**
 * @brief Generate a generic trace event.
 *
 * @param eventID @input_type  The ID of the trace.
 * @param pid     @input_type  Process ID.
 * @param data    @input_type  User-defined data passed in a character string.
 *
 * @return non-zero if successfull
 * @return 0 if failure
 */
int trace_event(
        const int gid,
        const int eventID,
        const int pid,
        const char *data)
{
    int rc = 0;

    if (tracefile && tracing_enabled(gid)) {
        rc = tracefile->output_generic_event(eventID,pid,data);
    }

    return rc;
}




/**
 * @brief Generate an interval start event.
 *
 * The \ref trace_start_interval "trace_start_interval()" function generates an
 * interval record that has information about the start, end, and duration of a
 * particular code fragment.
 *
 * Intervals are stored internally as a stack for each process ID.  When the
 * end_interval event arrives, we pop the starttime off the stack and use that
 * to calculate the duration.  Then we write an interval event to the
 * underlying trace framework.
 *
 * Pablo has its own interval records, but they are inadequate because it is
 * difficult to measure concurrent intervals (e.g., in threads).
 *
 * @param gid The trace group ID.
 * @param interval_id @output_type  The ID of the interval (unique for each interval)
 *
 * @return non-zero if successfull
 * @return 0 if failure
 */
int trace_start_interval(
        const int gid,
        const int pid)
{
    int rc = 0;

    log_debug(trace_debug_level, "trace_start_interval(gid=%d) start",
                    pid, gid);

    if (tracefile && tracing_enabled(gid)) {

        double time = trios_get_time();
        interval_iterator_t iter;

        /* see if an interval entry for this pid already exists */
        nthread_lock(&interval_mutex);
        iter = interval_map.find(pid);

        /* if the interval entry already exists, we push the starttime on the stack */
        if (iter != interval_map.end()) {
            iter->second.push(time);
        }

        /* register the start time of the interval */
        else {
            stack<double> s;
            s.push(time);
            interval_map.insert(pair< int, stack<double> >(pid, s));
        }
        nthread_unlock(&interval_mutex);
    }

    log_debug(trace_debug_level, "trace_start_interval(gid=%d) end", gid);

    return rc;
}


/**
 * @brief Generate an end-interval event.
 *
 * The \ref trace_interval_event "trace_interval_event()" function generates an
 * interval record that has information about the start, end, and duration of a
 * particular code fragment.
 *
 * @param intervalID @input_type The interval ID (unique for each interval).
 * @param eventID @input_type  The ID of the trace.
 * @param pid     @input_type  Process ID.
 * @param level   @input_type  The depth of the interval (0=root, 1=sub-interval, 2=sub-sub,...).
 * @param data    @input_type  User-defined data passed in a character
 *                             string.
 *
 * @return non-zero if successfull
 * @return 0 if failure
 */
int trace_end_interval(
        const int gid,
        const int eventID,
        const int pid,
        const char *data)
{
    int level = pid;
    int rc;
    log_level debug_level = trace_debug_level;

    log_debug(debug_level, "%d: trace_end_interval(gid=%d) start",
            pid, gid);

    if (tracefile && tracing_enabled(gid)) {
        double endtime = trios_get_time();
        double duration=0;
        double starttime=0;

        log_debug(debug_level, "%d: get interval lock", pid);

        /* get the interval event from the hashtable */
        nthread_lock(&interval_mutex);
        interval_iterator_t iter = interval_map.find(pid);

        /* case when interval is not found */
        if (iter == interval_map.end()) {
            log_info(trace_debug_level, "invalid gid %d\n", gid);
            nthread_unlock(&interval_mutex);
            return -1;
        }

        /* make sure the stack has at least one element */
        if (iter->second.size() == 0) {
            log_error(trace_debug_level, "Error empty interval stack for pid %d\n", pid);
            nthread_unlock(&interval_mutex);
            return -1;
        }

        /* The level is just the stack size */
        level = iter->second.size() - 1;

        starttime = iter->second.top();
        iter->second.pop();
        duration = endtime - starttime;

        log_debug(debug_level, "%d: release interval mutex", pid);
        nthread_unlock(&interval_mutex);

        /* erase the entry if the interval size is zero */
        //if (iter->second.size() == 0) {
        //    interval_map.erase(iter);
        //}


        /* output the interval event */
        rc = tracefile->output_interval_event(eventID, pid,
                level, data, duration);
    }

    log_debug(debug_level, "%d: trace_end_interval(gid=%d) end",
                pid,gid);

    return rc;
}

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
 * @param interval_id @input_type  The ID of the interval (unique for each interval)
 *
 * @return non-zero if successfull
 * @return 0 if failure
 */
int trace_start_tput_interval(
        const int gid,
        const int pid)
{
    return trace_start_interval(gid, pid);
}


/**
 * @brief Generate an end-interval event.
 *
 * The \ref trace_interval_event "trace_interval_event()" function generates an
 * interval record that has information about the start, end, and duration of a
 * particular code fragment.
 *
 * @param intervalID @input_type The interval ID (unique for each interval).
 * @param eventID @input_type  The ID of the trace.
 * @param pid     @input_type  Process ID.
 * @param level   @input_type  The depth of the interval (0=root, 1=sub-interval, 2=sub-sub,...).
 * @param data    @input_type  User-defined data passed in a character
 *                             string.
 *
 * @return non-zero if successfull
 * @return 0 if failure
 */
int trace_end_tput_interval(
        const int gid,
        const int eventID,
        const int pid,
        const long num_processed,
        const char *data)
{
    int level;
    int rc = 0;

    if (tracefile && tracing_enabled(gid)) {

        double endtime = trios_get_time();
        double duration=0;
        double starttime=0;

        /* get the interval event from the hashtable */
        nthread_lock(&interval_mutex);
        interval_iterator_t iter = interval_map.find(pid);

        /* case when interval is not found */
        if (iter == interval_map.end()) {
            fprintf(stderr, "Error: invalid gid\n");
            return -1;
        }

        /* make sure the stack has at least one element */
        if (iter->second.size() == 0) {
            fprintf(stderr, "Error empty interval stack for pid %d\n", pid);
            return -1;
        }

        starttime = iter->second.top();
        iter->second.pop();
        duration = endtime - starttime;
        nthread_unlock(&interval_mutex);

        /* output the interval event */
        rc = tracefile->output_tput_event(eventID, pid, level,
                data, duration, num_processed);
    }

    return rc;
}

int trace_get_count(
        const int gid,
        const int cid)
{
    int result = 0;

    if (tracing_enabled(gid)) {

        nthread_lock(&count_mutex);
        count_iterator_t iter = count_map.find(cid);

        if (iter == count_map.end()) {
            result = -1;
        }

        else {
            result = iter->second;
        }
        nthread_unlock(&count_mutex);
    }

    return result;
}

/**
 * @brief Generate an increment-count event.
 *
 * The \ref trace_inc_count_event "trace_inc_count_event()" function
 * increments a counter associated with the ID of counter and outputs
 * a count event to the tracefile.
 *
 * @param id @input_type The ID of the counter.
 * @param pid @input_type Process ID.
 * @param data @input_type User-defined data in a character string.
 *
 * @returns 0 on sucess, -1 if the cid is invalid.
 */
int trace_inc_count(
        const int gid,
        const int cid,
        const int pid,
        const char *data)
{
    int rc = 0;
    int val=0;


    if ((tracefile != NULL) && tracing_enabled(gid)) {

        log_debug(trace_debug_level, "%d: before count_mutex lock",nthread_self());
        nthread_lock(&count_mutex);
        log_debug(trace_debug_level, "%d: after count_mutex lock",nthread_self());

        count_iterator_t iter = count_map.find(cid);

        if (iter == count_map.end()) {
            trace_reset_count(gid, cid, pid, data);
            iter = count_map.find(cid);
        }

        val = ++(iter->second);
        nthread_unlock(&count_mutex);

        tracefile->output_count_event(cid, pid, data, val);
    }


    return rc;
}


/**
 * @brief Generate an decrement-count event.
 *
 * The \ref trace_inc_count_event "trace_inc_count_event()" function
 * increments a counter associated with the ID of counter and outputs
 * a count event to the tracefile.
 *
 * @param id @input_type The ID of the counter.
 * @param pid @input_type Process ID.
 * @param data @input_type User-defined data in a character string.
 *
 * @returns non-zero if successful.
 * @returns 0 if failure.
 */
int trace_dec_count(
        const int gid,
        const int cid,
        const int pid,
        const char *data)
{
    int rc = 0;
    int val;

    if (tracefile && tracing_enabled(gid)) {

        nthread_lock(&count_mutex);
        count_iterator_t iter = count_map.find(cid);

        /* return error if a counter does not exist */
        if (iter == count_map.end()) {
            log_warn(trace_debug_level, "%d: can't find cid=%d", cid);
            rc = -1;
        }
        else {
            val = --(iter->second);
            rc = tracefile->output_count_event(cid, pid, data, val);
        }
        nthread_unlock(&count_mutex);
    }

    return rc;
}

/**
 * @brief Generate a reset-count event.
 *
 * The \ref trace_reset_count "trace_reset_count()" function
 * sets the value of a counter associated with the ID to 0
 * and outputs a count event to the tracefile.
 *
 * @param id @input_type The ID of the counter.
 * @param pid @input_type Process ID.
 * @param data @input_type User-defined data in a character string.
 *
 * @returns non-zero if successful.
 * @returns 0 if failure.
 */
int trace_reset_count(
        const int gid,
        const int cid,
        const int pid,
        const char *data)
{
    return trace_set_count(gid, cid, pid, data, 0);
}

/**
 * @brief Generate an set-count event.
 *
 * The \ref trace_reset_count "trace_reset_count()" function
 * sets the value of a counter associated with the ID
 * and outputs a count event to the tracefile.
 *
 * @param id @input_type The ID of the counter.
 * @param pid @input_type Process ID.
 * @param data @input_type User-defined data in a character string.
 * @param new_count @input_type New value for the counter.
 *
 * @returns non-zero if successful.
 * @returns 0 if failure.
 */
int trace_set_count(
        const int gid,
        const int cid,
        const int pid,
        const char *data,
        const int new_count)
{
    int rc = 0;

    if (tracefile && tracing_enabled(gid)) {

        nthread_lock(&count_mutex);
        count_iterator_t iter = count_map.find(cid);

        /* if counter does not exist, create one */
        if (iter == count_map.end()) {
            count_map.insert(pair< int, int>(cid, new_count));
        }

        else {
            iter->second = new_count;
        }
        nthread_unlock(&count_mutex);

        rc = tracefile->output_count_event(cid, pid, data, new_count);
    }

    return rc;
}

struct CountWriter {

    const int pid;
    const char *data;

    CountWriter(const int p, const char *d): pid(p), data(d)
        { }

    void operator()(pair<int, int> p) {
        int val = p.second;
        int cid = p.first;
        if (tracefile) {
            tracefile->output_count_event(cid, pid, data, val);
        }
    }
};

int trace_put_all_counts(
        const int gid,
        const int pid,
        const char *data)
{
    int rc = 0;

    if (tracefile && tracing_enabled(gid)) {
        CountWriter writer(pid, data);

        /* we need to be able to iterate through all entries in the hashtable */
        nthread_lock(&count_mutex);
        for_each(count_map.begin(), count_map.end(), writer);
        nthread_unlock(&count_mutex);
    }

    return rc;
}

#endif  /* Trios_TRACIING_ENABLED */
