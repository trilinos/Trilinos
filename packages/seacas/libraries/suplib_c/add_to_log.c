/*
 * Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#if !defined(WIN32) && !defined(__WIN32__) && !defined(_WIN32) && !defined(_MSC_VER) &&            \
    !defined(__MINGW32__) && !defined(_WIN64) && !defined(__MINGW64__)
/* Currently we just disable this functionality for windows-based systems... */
#include <stdlib.h>
#include <string.h>

#include <sys/times.h>
#include <sys/utsname.h>

#include <time.h>
#include <unistd.h>

#ifndef __USE_XOPEN
#define __USE_XOPEN
#endif
#include <stdio.h>
#endif

void add_to_log(const char *my_name, double elapsed)
{
#if !defined(WIN32) && !defined(__WIN32__) && !defined(_WIN32) && !defined(_MSC_VER) &&            \
    !defined(__MINGW32__) && !defined(_WIN64) && !defined(__MINGW64__)
#define LEN 512
  /* Don't log information if this environment variable is set */
  if (getenv("SEACAS_NO_LOGGING") != NULL) {
    fprintf(stderr, "SEACAS Audit logging disabled via SEACAS_NO_LOGGING setting.\n");
    return;
  }

  /* Now try to find the $ACCESS/etc/audit.log file */
  /* Don't need to try too hard since information is not critical; just useful */
  char *access_dir = getenv("ACCESS");
  if (access_dir != NULL) {
    char filename[LEN];
    snprintf(filename, LEN, "%s/etc/audit.log", access_dir);
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

      char       time_string[LEN];
      time_t     calendar_time = time(NULL);
      struct tm *local_time    = localtime(&calendar_time);
      strftime(time_string, LEN, "%a %b %d %H:%M:%S %Z %Y", local_time);

      int        ticks_per_second;
      struct tms time_buf;
      times(&time_buf);
      ticks_per_second = sysconf(_SC_CLK_TCK);
      double u_time    = (double)(time_buf.tms_utime + time_buf.tms_cutime) / ticks_per_second;
      double s_time    = (double)(time_buf.tms_stime + time_buf.tms_cstime) / ticks_per_second;

      struct utsname sys_info;
      uname(&sys_info);

      int    minutes = (int)(elapsed / 60.0);
      double seconds = elapsed - minutes * 60.0;

      char log_string[LEN];
      snprintf(log_string, LEN, "%s %s %s %.3fu %.3fs %d:%5.2f 0.0%% 0+0k 0+0io 0pf+0w %s\n",
               codename, username, time_string, u_time, s_time, minutes, seconds,
               sys_info.nodename) < 0
          ? abort()
          : (void)0;

      fprintf(audit, "%s", log_string);
      fclose(audit);
    }
  }
#else
  /* Try to avoid unused variable warning/error on windows-based system */
  (void)my_name;
  (void)elapsed;
#endif
}
