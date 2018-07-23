/*
 * Copyright(C) 2009-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h>

#if defined(__LIBCATAMOUNT__)
#include <sys/resource.h>
#include <sys/time.h>
#else
#include <sys/times.h>
#endif

#ifndef __USE_XOPEN
#define __USE_XOPEN
#endif
#include <stdio.h>

void add_to_log(const char *my_name, double elapsed)
{
#define LEN 512
  char time_string[LEN];
  char log_string[LEN];

  double         u_time, s_time;
  struct utsname sys_info;

  int    minutes;
  double seconds;

  char *access_dir = NULL;

  /* Don't log information if this environment variable is set */
  if (getenv("SEACAS_NO_LOGGING") != NULL) {
    fprintf(stderr, "SEACAS Audit logging disabled via SEACAS_NO_LOGGING setting.\n");
    return;
  }

  /* Now try to find the $ACCESS/etc/audit.log file */
  /* Don't need to try too hard since information is not critical; just useful */
  access_dir = getenv("ACCESS");
  if (access_dir != NULL) {
    char filename[LEN];
    snprintf(filename, LEN, "%s/etc/audit.log", access_dir);
    if (0 == access(filename, W_OK)) {
      FILE *audit = fopen(filename, "a");
      if (audit != NULL) {
        const char *codename = strrchr(my_name, '/');

        const char *username = getlogin();
        if (username == NULL) {
          username = getenv("LOGNAME");
        }
        if (username == NULL) {
          username = "UNKNOWN";
        }

        if (codename == NULL) {
          codename = my_name;
        }
        else {
          codename++;
        }

        {
          time_t     calendar_time = time(NULL);
          struct tm *local_time    = localtime(&calendar_time);
          strftime(time_string, LEN, "%a %b %d %H:%M:%S %Z %Y", local_time);
        }

        {
#if defined(__LIBCATAMOUNT__)
          struct rusage rusage;

          getrusage(RUSAGE_SELF, &rusage);
          /*pp
           * NOTE: Catamount seems to return the same values for user and system.
           *       To avoid double-counting cpu time, I only use the user time.
           *       and set the system time to 0.
           */
          u_time = rusage.ru_utime.tv_sec + rusage.ru_utime.tv_usec / 1.e6;
          s_time = 0.0;
#else
          int        ticks_per_second;
          struct tms time_buf;
          times(&time_buf);
          ticks_per_second = sysconf(_SC_CLK_TCK);
          u_time           = (double)(time_buf.tms_utime + time_buf.tms_cutime) / ticks_per_second;
          s_time           = (double)(time_buf.tms_stime + time_buf.tms_cstime) / ticks_per_second;
#endif
        }

        uname(&sys_info);

        minutes = (int)(elapsed / 60.0);
        seconds = elapsed - minutes * 60.0;

        snprintf(log_string, LEN, "%s %s %s %.3fu %.3fs %d:%5.2f 0.0%% 0+0k 0+0io 0pf+0w %s\n",
                 codename, username, time_string, u_time, s_time, minutes, seconds,
                 sys_info.nodename);

        fprintf(audit, "%s", log_string);
        fclose(audit);
      }
    }
  }
}
