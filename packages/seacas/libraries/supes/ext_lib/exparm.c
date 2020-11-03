/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 */

/*
 *     DESCRIPTION:
 *     This routine defines various operating environment parameters. A
 *     character ID is supplied for both the processor hardware type and
 *     the operating system. The job processing mode, batch or
 *     interactive, is identified; for this purpose an interactive job
 *     is defined as one where the standard input device is attended by
 *     the user who can respond to unforeseen events. The number of
 *     character storage units and the number of numeric storage units
 *     in the smallest block of storage which contains an integral
 *     number of each are defined here. This routine further defines
 *     whether the record length of a direct access unformatted file is
 *     counted in character or numeric storage units.
 *
 *     This skeleton version returns blank strings for hardware and
 *     software IDs, zero to indicate batch mode, and unity for all
 *     other values.
 *
 *     FORMAL PARAMETERS:
 *     HARD      CHARACTER       System Hardware ID
 *     SOFT      CHARACTER       System Software ID
 *     MODE      INTEGER         Job Mode ( 0=batch , 1=interactive )
 *     KCSU      INTEGER         Number of Character Storage Units
 *     KNSU      INTEGER         Number of Numeric Storage Units
 *     IDAU      INTEGER         Unformatted Direct Access Units:
 *                                  0 = Character Storage Units
 *                                  1 = Numeric Storage Units
 *
 ************************************************************************
 *
 */

#include "fortranc.h"
#include <string.h>
#if defined(__NO_CYGWIN_OPTION__)
#include <windows.h>
#else
#include <sys/utsname.h>
#include <unistd.h> /* isatty  */
#endif
#include <stdio.h> /* sprintf */

#ifdef _MSC_VER
#include <io.h>
#include <sys/ioctl.h>
#define isatty _isatty
#endif

static char *copy_string(char *dest, char const *source, long int elements)
{
  char *d;
  for (d = dest; d + 1 < dest + elements && *source; d++, source++) {
    *d = *source;
  }
  *d = '\0';
  return d;
}

#define MAXCHAR 80
#define WORDLEN 8 /* Note that we *FORCE* the Fortran string */
                  /* length be 8 plus 1 for trailing null for the strings hard and soft. */
#if defined(ADDC_)
void exparm_(char *hard, char *soft, FTNINT *mode, FTNINT *kcsu, FTNINT *knsu, FTNINT *idau,
             FTNINT hlen, FTNINT slen)
#else
void exparm(char *hard, char *soft, FTNINT *mode, FTNINT *kcsu, FTNINT *knsu, FTNINT *idau,
            FTNINT hlen, FTNINT slen)
#endif
{
/********************************************************************/
#if defined(sgi)
  struct utsname SysInfo;
  char           hardname[MAXCHAR];
  char           softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL);
  *knsu = 1;

  uname(&SysInfo);

  sprintf(hardname, "SGI-%.4s", SysInfo.machine);
  sprintf(softname, "%.8s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* Silicon Graphics */
/********************************************************************/
#if defined(aix)

  struct utsname SysInfo;
  char           hardname[MAXCHAR];
  char           softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* IBM has 32 bit words */
  *knsu = 1;

  uname(&SysInfo);

  sprintf(hardname, "IBM %.4s", SysInfo.machine);
  sprintf(softname, "%03s %.4s", SysInfo.sysname, SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* IBM */
/********************************************************************/
#if defined(hpux)

  struct utsname SysInfo;
  char           hardname[MAXCHAR];
  char           softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* HP has 32 bit words */
  *knsu = 1;

  uname(&SysInfo);

  sprintf(hardname, "HP  %.4s", SysInfo.machine);
  sprintf(softname, "HP-UX%.3s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* HPUX */
/********************************************************************/
#if defined(__osf__)

  struct utsname SysInfo;
  char           hardname[MAXCHAR];
  char           softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* Chars/float */
  *knsu = 1;               /* Ditto */

  uname(&SysInfo);

  sprintf(hardname, "%.8s", SysInfo.machine);
  sprintf(softname, "%.4s%.4s", SysInfo.sysname, SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif
/********************************************************************/
#if defined(sun)
#if defined(SYSV) || defined(SVR4)

  struct utsname SysInfo;
  char           hardname[MAXCHAR];
  char           softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* See above */
  *knsu = 1;               /* Ditto */

  uname(&SysInfo);

  sprintf(hardname, "%.8s", SysInfo.machine);
  sprintf(softname, "SunOS%.3s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#else

  char  softname[MAXCHAR];
  char *darg;
  FILE *pipe;

  *idau = 0;
  *kcsu = 4; /* See above */
  *knsu = 1; /* Ditto */

#if defined(sparc)
  copy_string(hard, "Sun4", strlen("Sun4"));
#else  /* Then assume it's a SUN 3 */
  copy_string(hard, "Sun3", strlen("Sun3"));
#endif /* sparc */

  if ((darg = (char *)fgets(
           softname, MAXCHAR,
           pipe = popen("/usr/ucb/strings /vmunix|/usr/bin/fgrep Release", "r"))) == (char *)NULL) {
    perror("exparm: bad read from pipe");
    exit(1);
  }

  fclose(pipe);

  /*  pclose(pipe); */
  /* The previous system call to pclose() */
  /* is crapping out under SunView. */
  /* I think that this is due to the */
  /* event notification scheme that exists */
  /* under that window environment. */
  /* This occurs due to the wait() */
  /* function call in pclose(). */

  sscanf(softname, "%*s%*s%6s", soft + 2);

  strncat(hard, "    ", strlen("    ")); /* Another case of hardwiring the length of hard. */
  soft[0] = 'O';
  soft[1] = 'S';

#endif /* SYSV || SVR4 (Solaris 2.X) */
#endif

  /********************************************************************/

#if defined(__NO_CYGWIN_OPTION__)
  SYSTEM_INFO   SysInfo;
  OSVERSIONINFO OSInfo;
  char          hardname[MAXCHAR];
  char          softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* See above */
  *knsu = 1;               /* Ditto */

  GetSystemInfo(&SysInfo);

  sprintf(hardname, "%-8d", SysInfo.dwProcessorType);

  OSInfo.dwOSVersionInfoSize = sizeof(OSInfo);
  if (GetVersionEx(&OSInfo) > 0) {
    switch (OSInfo.dwPlatformId) {
    case 1: sprintf(softname, "Win98%d.%d", OSInfo.dwMajorVersion, OSInfo.dwMinorVersion); break;
    case 2: sprintf(softname, "WinNT%d.%d", OSInfo.dwMajorVersion, OSInfo.dwMinorVersion); break;
    default: sprintf(softname, "WinOS%d.%d", OSInfo.dwMajorVersion, OSInfo.dwMinorVersion); break;
    }
  }
  else {
    sprintf(softname, "Unknown OS");
  }
  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

/* cygwin native */
/********************************************************************/
#elif defined(__CYGWIN__)
  struct utsname SysInfo;
  char           hardname[MAXCHAR];
  char           softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* See above */
  *knsu = 1;               /* Ditto */

  uname(&SysInfo);

  sprintf(hardname, "%.8s", SysInfo.machine);
  sprintf(softname, "CW%.6s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif
/********************************************************************/
#if defined(__APPLE__)
  struct utsname SysInfo;
  char           hardname[MAXCHAR];
  char           softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* Darwin has 32 bit words */
  *knsu = 1;

  uname(&SysInfo);

  sprintf(hardname, "%.8s", SysInfo.machine);
  sprintf(softname, "OSX%.5s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* Darwin (Power Macintosh)*/

/********************************************************************/
#if defined(__linux__) || defined(interix)
  struct utsname SysInfo;
  char           hardname[MAXCHAR];
  char           softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* See above */
  *knsu = 1;               /* Ditto */

  uname(&SysInfo);

  sprintf(hardname, "%.8s", SysInfo.machine);
  sprintf(softname, "Lx%.6s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* Linux, Interix */
  /********************************************************************/

  if (isatty(0) != 0) /* Test stdin as to whether or not it's a terminal. */
    *mode = 1;        /* Indicates an interactive process. */
  else
    *mode = 0; /* Else running from a procedure. */
}
