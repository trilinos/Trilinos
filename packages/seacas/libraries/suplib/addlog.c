/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <stdlib.h>
#include <string.h>
#if !defined(__CYGWIN__)
#include <sys/resource.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h>
#endif

#ifndef __USE_XOPEN
#define __USE_XOPEN
#endif
#include <stdio.h>

#if defined(ADDC_)
void addlog_(char *name, int len)
#else
void addlog(char *name, int len)
#endif
{
#if !defined(__CYGWIN__)
#define LEN 512
  char time_string[LEN];
  char log_string[LEN];
  char codename[LEN];

  double         u_time, s_time;
  struct utsname sys_info;
  const char *   username = NULL;

  /* Don't log information if this environment variable is set */
  if (getenv("SEACAS_NO_LOGGING") != NULL) {
    fprintf(stderr, "SEACAS Audit logging disabled via SEACAS_NO_LOGGING setting.\n");
    return;
  }

  username = getlogin();
  if (username == NULL) {
    username = getenv("LOGNAME");
  }
  if (username == NULL) {
    username = "UNKNOWN";
  }

  {
    int i;
    for (i = 0; i < len; i++)
      codename[i] = name[i];
    codename[len] = '\0';
  }

  {
    time_t     calendar_time = time(NULL);
    struct tm *local_time    = localtime(&calendar_time);
    strftime(time_string, LEN, "%a %b %d %H:%M:%S %Z %Y", local_time);
  }

  {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    u_time = rusage.ru_utime.tv_sec + rusage.ru_utime.tv_usec / 1.0e6;
    s_time = rusage.ru_stime.tv_sec + rusage.ru_stime.tv_usec / 1.0e6;
  }

  uname(&sys_info);

  snprintf(log_string, LEN, "%s %s %s %.3fu %.3fs 0:00.00 0.0%% 0+0k 0+0io 0pf+0w %s\n", codename,
           username, time_string, u_time, s_time, sys_info.nodename) < 0
      ? abort()
      : (void)0;

  /* Now try to find the $ACCESS/etc/audit.log file */
  /* Don't need to try too hard since information is not critical; just useful */
  {
    char *access_dir = getenv("ACCESS");
    if (access_dir != NULL) {
      char filename[LEN];
      snprintf(filename, LEN, "%s/etc/audit.log", access_dir);
      if (0 == access(filename, W_OK)) {
        FILE *audit = fopen(filename, "a");
        if (audit != NULL) {
          fprintf(audit, "%s", log_string);
          fclose(audit);
        }
      }
    }
  }
#endif
}
