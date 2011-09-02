/*
 * Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *           
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *                         
 * * Neither the name of Sandia Corporation nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
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
 * 
 */
/*
 * $Id: exparm.c,v 1.26 2008/03/14 13:22:37 gdsjaar Exp $
 */

/*
*     DESCRIPTION:
*     This routine defines various operating environment parameters. A
*     character ID is supplied for both the processor hardware type and
*     the operating system. The job processing mode, batch or
*     interactive, is identified; for this purpose an interactive job
*     is defined as one where the standard input device is attended by
*     the user who can respond to unforseen events. The number of
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
#include <unistd.h> /* isatty  */
#include <sys/utsname.h>
#endif
#include <stdio.h>  /* sprintf */

#define MAXCHAR 80
#define WORDLEN 8	/* Note that we *FORCE* the Fortran string */
			/* length be 8 for the strings hard and soft. */
#define HARD hardname
#define SOFT softname

#if defined(ADDC_)
void exparm_( char *hard, char *soft, FTNINT *mode,
	      FTNINT *kcsu, FTNINT *knsu, FTNINT *idau,
	      FTNINT hlen, FTNINT slen )
#else
void exparm( char *hard, char *soft, FTNINT *mode,
	     FTNINT *kcsu, FTNINT *knsu, FTNINT *idau,
	     FTNINT hlen, FTNINT slen )
#endif
{
/********************************************************************/
#if defined (sgi)
  struct utsname SysInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);
  *knsu = 1;

  uname( &SysInfo );

  sprintf( hardname, "SGI-%4s", SysInfo.machine );
  sprintf( softname, "%8s", SysInfo.release );

  strncpy( hard, HARD, WORDLEN );
  strncpy( soft, SOFT, WORDLEN );

#endif               /* Silicon Graphics */
/********************************************************************/
#if defined (aix)

  struct utsname SysInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);           /* IBM has 32 bit words */
  *knsu = 1;

  uname( &SysInfo );

  sprintf( hardname, "IBM %4s", SysInfo.machine );
  sprintf( softname, "%03s %4s", SysInfo.sysname, SysInfo.release );

  strncpy(  hard , HARD, WORDLEN );
  strncpy(  soft , SOFT, WORDLEN );

#endif			/* IBM */
/********************************************************************/
#if defined (hpux) 

  struct utsname SysInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);           /* HP has 32 bit words */
  *knsu = 1;

  uname( &SysInfo );

  sprintf( hardname, "HP  %4s", SysInfo.machine );
  sprintf( softname, "HP-UX%3s", SysInfo.release );

  strncpy(  hard , HARD, WORDLEN );
  strncpy(  soft , SOFT, WORDLEN );

#endif			/* HPUX */
/********************************************************************/
#if defined (paragon)
 
      struct utsname SysInfo;
      char hardname[MAXCHAR];
      char softname[MAXCHAR];
 
      *idau = 0;
      *kcsu = sizeof (FTNREAL);                /* 860 has 32 bit words */
      *knsu = 1;
 
      uname( &SysInfo );
 
      sprintf( hardname, "%8s", "i860 GP ");
      sprintf( softname, "OSF %2s.%1s", SysInfo.version, SysInfo.release );
 
      strncpy(  hard , HARD, WORDLEN );
      strncpy(  soft , SOFT, WORDLEN );
 
#endif                         
/********************************************************************/
#if defined (pumagon) || defined (p6)
 
      struct utsname SysInfo;
      char hardname[MAXCHAR];
      char softname[MAXCHAR];
 
      *idau = 0;
      *kcsu = sizeof (FTNREAL);                /* 860 has 32 bit words */
      *knsu = 1;
 
      sprintf( hardname, "%8s", "i860 GP ");
      sprintf( softname, "SUNMOS ");
 
      strncpy(  hard , HARD, WORDLEN );
      strncpy(  soft , SOFT, WORDLEN );
 
#endif                 
/********************************************************************/
#if defined (p6)
 
      struct utsname SysInfo;
      char hardname[MAXCHAR];
      char softname[MAXCHAR];
 
      *idau = 0;
      *kcsu = sizeof (FTNREAL);                /* P6 has 32 bit words */
      *knsu = 1;
 
      sprintf( hardname, "%8s", "P6      ");
      sprintf( softname, "Solari  ");
 
      strncpy(  hard , HARD, WORDLEN );
      strncpy(  soft , SOFT, WORDLEN );
#endif                       
/********************************************************************/
#if defined (cougar)
 
      struct utsname SysInfo;
      char hardname[MAXCHAR];
      char softname[MAXCHAR];
 
      *idau = 0;
      *kcsu = sizeof (FTNREAL);      /* p6 has 32 bit words */
      *knsu = 1;
 
      sprintf( hardname, "%8s", "PentPro ");
      sprintf( softname, "COUGAR  ");
 
      strncpy(  hard , HARD, WORDLEN );
      strncpy(  soft , SOFT, WORDLEN );
 
#endif           
/********************************************************************/
#if defined (__osf__)

  struct utsname SysInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);	/* Chars/float */
  *knsu = 1;			/* Ditto */

  uname( &SysInfo );

  sprintf( hardname, "%8s", SysInfo.machine );
  sprintf( softname, "%4s%4s", SysInfo.sysname, SysInfo.release );

  strncpy( hard, HARD, WORDLEN );
  strncpy( soft, SOFT, WORDLEN );

#endif
/********************************************************************/
#if defined (sun)
#if defined(SYSV) || defined(SVR4)

  struct utsname SysInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);	/* See above */
  *knsu = 1;			/* Ditto */

  uname( &SysInfo );

  sprintf( hardname, "%8s", SysInfo.machine );
  sprintf( softname, "SunOS%3s", SysInfo.release );

  strncpy(  hard , HARD, WORDLEN );
  strncpy(  soft , SOFT, WORDLEN );

#else

  char hardname[MAXCHAR];
  char softname[MAXCHAR];
  char *darg;
  FILE *pipe;

  *idau = 0;
  *kcsu = 4;			/* See above */
  *knsu = 1;			/* Ditto */

#if defined (sparc)
  strncpy( hard, "Sun4", strlen("Sun4") );
#else				/* Then assume it's a SUN 3 */
  strncpy( hard, "Sun3", strlen("Sun3") );
#endif				/* sparc */

  if ( ( darg = (char*)fgets( SOFT, MAXCHAR,
      pipe = popen("/usr/ucb/strings /vmunix|/usr/bin/fgrep Release", "r") ) ) 
      == (char *)NULL)
    {
      perror("exparm: bad read from pipe");
      exit(1);
    }

  fclose( pipe );

/*  pclose(pipe); */
				/* The previous system call to pclose() */
				/* is crapping out under SunView. */
				/* I think that this is due to the */
				/* event notification scheme that exists */
				/* under that window environment. */
				/* This occurs due to the wait() */
				/* function call in pclose(). */

  sscanf( SOFT, "%*s%*s%6s", soft+2 );

  strncat( hard, "    ", strlen("    ") );	/* Another case of hardwiring the length of hard. */
  soft[0]='O';
  soft[1]='S';

#endif                          /* SYSV || SVR4 (Solaris 2.X) */
#endif

/********************************************************************/

#if defined(__NO_CYGWIN_OPTION__)
  SYSTEM_INFO SysInfo;
  OSVERSIONINFO OSInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);	/* See above */
  *knsu = 1;			/* Ditto */

  GetSystemInfo(&SysInfo);

  sprintf( hardname, "%-8d", SysInfo.dwProcessorType );

  OSInfo.dwOSVersionInfoSize = sizeof(OSInfo);
  if (GetVersionEx(&OSInfo) > 0)
  {
    switch(OSInfo.dwPlatformId)
    {
      case 1: sprintf( softname, "Win98%d.%d",
                       OSInfo.dwMajorVersion, OSInfo.dwMinorVersion);
              break;
      case 2: sprintf( softname, "WinNT%d.%d",
                       OSInfo.dwMajorVersion, OSInfo.dwMinorVersion);
              break;
      default: sprintf( softname, "WinOS%d.%d",
                       OSInfo.dwMajorVersion, OSInfo.dwMinorVersion);
              break;
    }
  }
  else
  {
    sprintf( softname, "Unknown OS" );
  }
  strncpy(  hard , HARD, WORDLEN );
  strncpy(  soft , SOFT, WORDLEN );

  /* cygwin native */
/********************************************************************/
#elif defined(__CYGWIN__)
#undef linux
  struct utsname SysInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);	/* See above */
  *knsu = 1;			/* Ditto */

  uname( &SysInfo );

  sprintf( hardname, "%8s", SysInfo.machine );
  sprintf( softname, "CW%6s", SysInfo.release );

  strncpy(  hard , HARD, WORDLEN );
  strncpy(  soft , SOFT, WORDLEN );

#endif
/********************************************************************/
#if defined(__APPLE__)
#if defined linux
#undef linux
#endif
  struct utsname SysInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);           /* Darwin has 32 bit words */
  *knsu = 1;

  uname( &SysInfo );

  sprintf( hardname, "%8s", "PowerMac" );
  sprintf( softname, "Darw %3s", SysInfo.release );

  strncpy(  hard , HARD, WORDLEN );
  strncpy(  soft , SOFT, WORDLEN );

#endif/* Darwin (Power Macintosh)*/

/********************************************************************/
#if defined (linux) || defined (interix) 
  struct utsname SysInfo;
  char hardname[MAXCHAR];
  char softname[MAXCHAR];

  *idau = 0;
  *kcsu = sizeof (FTNREAL);	/* See above */
  *knsu = 1;			/* Ditto */

  uname( &SysInfo );

  sprintf( hardname, "%8s", SysInfo.machine );
  sprintf( softname, "Lx%6s", SysInfo.release );

  strncpy(  hard , HARD, WORDLEN );
  strncpy(  soft , SOFT, WORDLEN );

#endif                          /* Linux, Interix */
  /********************************************************************/

  if( isatty(0) != 0 )	/* Test stdin as to whether or not it's a terminal. */
    *mode = 1;			/* Indicates an interactive process. */
  else
    *mode = 0;			/* Else running from a procedure. */
}
