/*
 * Copyright(C) 1999-2021, 2023 National Technology & Engineering Solutions
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
#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__) ||                             \
    defined(__NO_CYGWIN_OPTION__)
#define NOMINMAX
#include <io.h>
#include <windows.h>
#define isatty _isatty
#else
#include <sys/utsname.h>
#include <unistd.h>
#endif
#include <stdio.h>

static char *copy_string(char *dest, char const *source, long int elements)
{
  char *d;
  for (d = dest; d + 1 < dest + elements && *source; d++, source++) {
    *d = *source;
  }
  *d = '\0';
  return d;
}

#define SUPES_MAXCHAR 80
#define WORDLEN       8 /* Note that we *FORCE* the Fortran string */
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
  char           hardname[SUPES_MAXCHAR];
  char           softname[SUPES_MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL);
  *knsu = 1;

  uname(&SysInfo);

  snprintf(hardname, SUPES_MAXCHAR, "SGI-%.4s", SysInfo.machine);
  snprintf(softname, SUPES_MAXCHAR, "%.8s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* Silicon Graphics */
/********************************************************************/
#if defined(aix)

  struct utsname SysInfo;
  char           hardname[SUPES_MAXCHAR];
  char           softname[SUPES_MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* IBM has 32 bit words */
  *knsu = 1;

  uname(&SysInfo);

  snprintf(hardname, SUPES_MAXCHAR, "IBM %.4s", SysInfo.machine);
  snprintf(softname, SUPES_MAXCHAR, "%03s %.4s", SysInfo.sysname, SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* IBM */
/********************************************************************/
#if defined(hpux)

  struct utsname SysInfo;
  char           hardname[SUPES_MAXCHAR];
  char           softname[SUPES_MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* HP has 32 bit words */
  *knsu = 1;

  uname(&SysInfo);

  snprintf(hardname, SUPES_MAXCHAR, "HP  %.4s", SysInfo.machine);
  snprintf(softname, SUPES_MAXCHAR, "HP-UX%.3s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* HPUX */
/********************************************************************/
#if defined(__osf__)

  struct utsname SysInfo;
  char           hardname[SUPES_MAXCHAR];
  char           softname[SUPES_MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* Chars/float */
  *knsu = 1;               /* Ditto */

  uname(&SysInfo);

  snprintf(hardname, SUPES_MAXCHAR, "%.8s", SysInfo.machine);
  snprintf(softname, SUPES_MAXCHAR, "%.4s%.4s", SysInfo.sysname, SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif
  /********************************************************************/

#if defined(__NO_CYGWIN_OPTION__)
  SYSTEM_INFO   SysInfo;
  OSVERSIONINFO OSInfo;
  char          hardname[SUPES_MAXCHAR];
  char          softname[SUPES_MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* See above */
  *knsu = 1;               /* Ditto */

  GetSystemInfo(&SysInfo);

  snprintf(hardname, SUPES_MAXCHAR, "%-8d", SysInfo.dwProcessorType);

  OSInfo.dwOSVersionInfoSize = sizeof(OSInfo);
  if (GetVersionEx(&OSInfo) > 0) {
    switch (OSInfo.dwPlatformId) {
    case 1:
      snprintf(softname, SUPES_MAXCHAR, "Win98%d.%d", OSInfo.dwMajorVersion, OSInfo.dwMinorVersion);
      break;
    case 2:
      snprintf(softname, SUPES_MAXCHAR, "WinNT%d.%d", OSInfo.dwMajorVersion, OSInfo.dwMinorVersion);
      break;
    default:
      snprintf(softname, SUPES_MAXCHAR, "WinOS%d.%d", OSInfo.dwMajorVersion, OSInfo.dwMinorVersion);
      break;
    }
  }
  else {
    snprintf(softname, SUPES_MAXCHAR, "Unknown OS");
  }
  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

/* cygwin native */
/********************************************************************/
#elif defined(__CYGWIN__)
  struct utsname SysInfo;
  char           hardname[SUPES_MAXCHAR];
  char           softname[SUPES_MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* See above */
  *knsu = 1;               /* Ditto */

  uname(&SysInfo);

  snprintf(hardname, SUPES_MAXCHAR, "%.8s", SysInfo.machine);
  snprintf(softname, SUPES_MAXCHAR, "CW%.6s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif
/********************************************************************/
#if defined(__APPLE__)
  struct utsname SysInfo;
  char           hardname[SUPES_MAXCHAR];
  char           softname[SUPES_MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* Darwin has 32 bit words */
  *knsu = 1;

  uname(&SysInfo);

  snprintf(hardname, SUPES_MAXCHAR, "%.8s", SysInfo.machine);
  snprintf(softname, SUPES_MAXCHAR, "OSX%.5s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* Darwin (Power Macintosh)*/

/********************************************************************/
#if defined(__linux__) || defined(interix)
  struct utsname SysInfo;
  char           hardname[SUPES_MAXCHAR];
  char           softname[SUPES_MAXCHAR];

  *idau = 0;
  *kcsu = sizeof(FTNREAL); /* See above */
  *knsu = 1;               /* Ditto */

  uname(&SysInfo);

  snprintf(hardname, SUPES_MAXCHAR, "%.8s", SysInfo.machine);
  snprintf(softname, SUPES_MAXCHAR, "Lx%.6s", SysInfo.release);

  copy_string(hard, hardname, WORDLEN);
  copy_string(soft, softname, WORDLEN);

#endif /* Linux, Interix */
  /********************************************************************/

  if (isatty(0) != 0) { /* Test stdin as to whether or not it's a terminal. */
    *mode = 1;          /* Indicates an interactive process. */
  }
  else {
    *mode = 0; /* Else running from a procedure. */
  }
}
