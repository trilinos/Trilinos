#if 0 // Disabled!
// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_FloatingPointTrap.hpp"
#include "Teuchos_Assert.hpp"


//
// Implementation of floating point control
//


// ???


//
// We have floating point control!
//


static void ieee0(bool enableTrap);


void Teuchos::doFloatingPointTrap(bool enableTrap)
{
  ieee0(enableTrap);
}



//
// Stuff from uninit.c from f2c and them from Sacado!
//


//
// Up front stuff
//


#include <cstdio>
#include <cstring>

#define TYSHORT 2
#define TYLONG 3
#define TYREAL 4
#define TYDREAL 5
#define TYCOMPLEX 6
#define TYDCOMPLEX 7
#define TYINT1 11
#define TYQUAD 14

#ifndef Long
#define Long long
#endif // Long

#ifdef __mips
#define RNAN	0xffc00000
#define DNAN0	0xfff80000
#define DNAN1	0
#endif // __mips

#ifdef _PA_RISC1_1
#define RNAN	0xffc00000
#define DNAN0	0xfff80000
#define DNAN1	0
#endif // _PA_RISC1_1

#ifndef RNAN
#  define RNAN	0xff800001
#  ifdef IEEE_MC68k
#    define DNAN0	0xfff00000
#    define DNAN1	1
#  else
#    define DNAN0	1
#    define DNAN1	0xfff00000
#  endif
#endif /*RNAN*/


static unsigned Long rnan = RNAN, dnan0 = DNAN0, dnan1 = DNAN1;

double _0 = 0.;

#ifndef MSpc
#  ifdef MSDOS
#    define MSpc
#  else
#    ifdef _WIN32
#      define MSpc
#    endif
#  endif
#endif


//
// MSpc
//


#ifdef MSpc

#define IEEE0_done
#include "cfloat"
#include "csignal"

static void ieee0(bool enableTrap)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    enableTrap == false, std::logic_error,
    "Error, don't know how to turn off trap for MSpc!"
    );
#ifndef __alpha
  _control87(EM_DENORMAL | EM_UNDERFLOW | EM_INEXACT, MCW_EM);
#endif
  /* With MS VC++, compiling and linking with -Zi will permit */
  /* clicking to invoke the MS C++ debugger, which will show */
  /* the point of error -- provided SIGFPE is SIG_DFL. */
  signal(SIGFPE, SIG_DFL);
}

#endif // MSpc


//
// MIPS
//


#ifdef __mips	/* must link with -lfpe */

#define IEEE0_done
#include <cstdlib>
#include <cstdio>
#include "/usr/include/sigfpe.h"	/* full pathname for lcc -N */
#include "/usr/include/sys/fpu.h"

static void ieeeuserhand(unsigned std::exception[5], int val[2])
{
	fflush(stdout);
	fprintf(stderr,"ieee0() aborting because of ");
	if(std::exception[0]==_OVERFL) fprintf(stderr,"overflow\n");
	else if(std::exception[0]==_UNDERFL) fprintf(stderr,"underflow\n");
	else if(std::exception[0]==_DIVZERO) fprintf(stderr,"divide by 0\n");
	else if(std::exception[0]==_INVALID) fprintf(stderr,"invalid operation\n");
	else fprintf(stderr,"\tunknown reason\n");
	fflush(stderr);
	abort();
}

static void ieeeuserhand2(unsigned int **j)
{
	fprintf(stderr,"ieee0() aborting because of confusion\n");
	abort();
}

static void ieee0(bool enableTrap)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    enableTrap == false, std::logic_error,
    "Error, don't know how to turn off trap for MIPS!"
    );
	int i;
	for(i=1; i<=4; i++){
		sigfpe_[i].count = 1000;
		sigfpe_[i].trace = 1;
		sigfpe_[i].repls = _USER_DETERMINED;
		}
	sigfpe_[1].repls = _ZERO;	/* underflow */
	handle_sigfpes( _ON,
		_EN_UNDERFL|_EN_OVERFL|_EN_DIVZERO|_EN_INVALID,
		ieeeuserhand,_ABORT_ON_ERROR,ieeeuserhand2);
	}
#endif /* mips */


//
// Linux
//


#ifdef __linux__

#define IEEE0_done
#include "fpu_control.h"

#ifdef __alpha__
#  ifndef USE_setfpucw
#    define __setfpucw(x) __fpu_control = (x)
#  endif
#endif

#ifndef _FPU_SETCW
#  undef  Can_use__setfpucw
#  define Can_use__setfpucw
#endif

static void ieee0(bool enableTrap)
{

  TEUCHOS_TEST_FOR_EXCEPTION(
    enableTrap == false, std::logic_error,
    "Error, don't know how to turn off trap for LINUX!"
    );

#if (defined(__mc68000__) || defined(__mc68020__) || defined(mc68020) || defined (__mc68k__))

  /* Reported 20010705 by Alan Bain <alanb@chiark.greenend.org.uk> */
  /* Note that IEEE 754 IOP (illegal operation) */
  /* = Signaling NAN (SNAN) + operation error (OPERR). */

#ifdef Can_use__setfpucw /* Has __setfpucw gone missing from S.u.S.E. 6.3? */
	__setfpucw(_FPU_IEEE + _FPU_DOUBLE + _FPU_MASK_OPERR + _FPU_MASK_DZ + _FPU_MASK_SNAN + _FPU_MASK_OVFL);
#else
	__fpu_control = _FPU_IEEE + _FPU_DOUBLE + _FPU_MASK_OPERR + _FPU_MASK_DZ + _FPU_MASK_SNAN + _FPU_MASK_OVFL;
	_FPU_SETCW(__fpu_control);
#endif

#elif (defined(__powerpc__)||defined(_ARCH_PPC)||defined(_ARCH_PWR)) /* !__mc68k__ */
  /* Reported 20011109 by Alan Bain <alanb@chiark.greenend.org.uk> */

#ifdef Can_use__setfpucw

  /* The following is NOT a mistake -- the author of the fpu_control.h
     for the PPC has erroneously defined IEEE mode to turn on exceptions
     other than Inexact! Start from default then and turn on only the ones
     which we want*/

	__setfpucw(_FPU_DEFAULT +  _FPU_MASK_IM+_FPU_MASK_OM+_FPU_MASK_UM);

#else /* PPC && !Can_use__setfpucw */

	__fpu_control = _FPU_DEFAULT +_FPU_MASK_OM+_FPU_MASK_IM+_FPU_MASK_UM;
	_FPU_SETCW(__fpu_control);

#endif /*Can_use__setfpucw*/

#else /* !(mc68000||powerpc) */

#ifdef _FPU_IEEE
#  ifndef _FPU_EXTENDED /* e.g., ARM processor under Linux */
#  define _FPU_EXTENDED 0
#endif

#ifndef _FPU_DOUBLE
#  define _FPU_DOUBLE 0
#endif

#ifdef Can_use__setfpucw /* Has __setfpucw gone missing from S.u.S.E. 6.3? */
	__setfpucw(_FPU_IEEE - _FPU_EXTENDED + _FPU_DOUBLE - _FPU_MASK_IM - _FPU_MASK_ZM - _FPU_MASK_OM);
#else
	__fpu_control = _FPU_IEEE - _FPU_EXTENDED + _FPU_DOUBLE - _FPU_MASK_IM - _FPU_MASK_ZM - _FPU_MASK_OM;
	_FPU_SETCW(__fpu_control);
#endif

#else /* !_FPU_IEEE */

  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "Error, don't know how to trap floating-point errors on this Linux system!"
    );

#endif /* _FPU_IEEE */

#endif /* __mc68k__ */

} // ieee0()

#endif /* __linux__ */


//
// Alpha
//


#ifdef __alpha

#ifndef IEEE0_done

#define IEEE0_done
#include <machine/fpu.h>

static void ieee0(bool enableTrap)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    enableTrap == false, std::logic_error,
    "Error, don't know how to turn off trap for Alpha!"
    );
	ieee_set_fp_control(IEEE_TRAP_ENABLE_INV);
}

#endif /*IEEE0_done*/

#endif /*__alpha*/


//
// hpux
//



#ifdef __hpux

#define IEEE0_done
#define _INCLUDE_HPUX_SOURCE

#include <cmath>

#ifndef FP_X_INV
#  include <fenv.h>
#  define fpsetmask fesettrapenable
#  define FP_X_INV FE_INVALID
#endif

static void ieee0(bool enableTrap)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    enableTrap == false, std::logic_error,
    "Error, don't know how to turn off trap for HPUX!"
    );
	fpsetmask(FP_X_INV);
}

#endif /*__hpux*/


//
// AIX
//


#ifdef _AIX

#define IEEE0_done
#include <fptrap.h>

static void ieee0(bool enableTrap)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    enableTrap == false, std::logic_error,
    "Error, don't know how to turn off trap for AIX!"
    );
	fp_enable(TRP_INVALID);
	fp_trap(FP_TRAP_SYNC);
}

#endif /*_AIX*/


//
// SUN
//


#ifdef __sun

#define IEEE0_done
#include <ieeefp.h>

static void ieee0(bool enableTrap)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    enableTrap == false, std::logic_error,
    "Error, don't know how to turn off trap for SUN!"
    );
	fpsetmask(FP_X_INV);
}

#endif // __sun


//
// Default (none)
//


#ifndef IEEE0_done

static void ieee0(bool enableTrap)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "Error, Don't know how to implement floating-point traps on this platform!"
    );
}

#endif // IEEE0_done

#endif // 0 // Disabled!
