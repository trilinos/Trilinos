/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

/*
 *  Operating System
 */

#ifdef aix
#  include "os_aix.h"
#  define CONFIGURED
#endif

#ifdef bsd
#  include "os_bsd.h"
#  define CONFIGURED
#endif

#ifdef dec
#  include "os_dec.h"
#  define CONFIGURED
#endif

#ifdef sun
#  include "os_sun.h"
#  define CONFIGURED
#endif

#ifdef hpux
#  include "os_hpux.h"
#  define CONFIGURED
#endif

#ifdef sequent
#  include "os_dynix.h"
#  define CONFIGURED
#endif

#if defined (sgi) || defined(sgi10k)
# include "os_sgi.h"
# define CONFIGURED
#endif

#ifdef solaris
# include "os_solaris.h"
# define CONFIGURED
#endif

#ifdef ipsc
#  include "os_ipsc.h"
#  define CONFIGURED
#endif

#ifdef MSDOS
#  include "os_msdos.h"
#  define CONFIGURED
#endif

#ifdef NeXT
#  include "os_bsd.h"
#  define CONFIGURED
#endif

#ifdef THINK_C
#  include "os_mac.h"
#  define CONFIGURED
#endif

#ifdef alpha_osf
#  include "os_osf.h"
#  define CONFIGURED
#endif

#ifdef linux
#  include "os_linux.h"
#  define CONFIGURED
#endif

#ifdef freebsd
#  include "os_freebsd.h"
#  define CONFIGURED
#endif

#ifdef NMcplant
#  include "os_NMcplant.h"
#  define CONFIGURED
#endif

#ifdef mingw32
#  include "os_mingw32.h"
#  define CONFIGURED
#endif

#ifdef unix
#  include "os_unix.h"
#  define CONFIGURED
#endif

#ifdef ppc
#  include "os_linux.h"
#  define CONFIGURED
#endif

#ifndef CONFIGURED

error error error error

Operating system type unknown

error error error error

#endif
