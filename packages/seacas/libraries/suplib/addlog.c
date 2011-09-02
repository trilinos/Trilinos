/*
 * Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
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
 *     * Neither the name of Sandia Corporation nor the names of its
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
#include <string.h>
#include <stdlib.h>
#if !defined(__CYGWIN__)
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/resource.h>
#include <time.h>
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
#define LEN 256
  char time_string[LEN];
  char log_string[LEN];
  char codename[LEN];

  double u_time, s_time;
  struct utsname sys_info;
  char *username = NULL;

  /* Don't log information if this environment variable is set */
  if (getenv("SEACAS_NO_LOGGING") != NULL) {
    fprintf(stderr, "SEACAS Audit logging disabled via SEACAS_NO_LOGGING setting.\n");
    return;
  }
  
  username = getlogin();
  if (username == NULL) {
    username = getenv("LOGNAME");
  }
  
  {
    int i;
    for (i = 0; i < len; i++)
      codename[i] = name[i];
    codename[len] = '\0';
  }
  
  {
    time_t calendar_time = time(NULL);
    struct tm *local_time = localtime(&calendar_time);
    strftime(time_string, LEN, "%a %b %d %H:%M:%S %Z %Y", local_time);
  }

  {
    struct rusage rusage;
    getrusage(RUSAGE_SELF,&rusage);
    u_time = rusage.ru_utime.tv_sec + rusage.ru_utime.tv_usec / 1.0e6;
    s_time = rusage.ru_stime.tv_sec + rusage.ru_stime.tv_usec / 1.0e6;
  }
  
  uname(&sys_info);

  sprintf(log_string, "%s %s %s %.3fu %.3fs 0:00.00 0.0%% 0+0k 0+0io 0pf+0w %s\n",
	  codename, username, time_string, u_time, s_time, sys_info.nodename);

  /* Now try to find the $ACCESS/etc/audit.log file */
  /* Don't need to try too hard since information is not critical; just useful */
  {
    char *access_dir = getenv("ACCESS");
    if (access_dir != NULL) {
      char filename[LEN];
      sprintf(filename, "%s/etc/audit.log", access_dir);
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
