/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef MY_SIG_DH
#define MY_SIG_DH

#include "euclid_common.h"
#include <signal.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void sigRegister_dh();
extern void sigHandler_dh(int sig);

/* 
  list of signals the Euclid will handle
*/
#ifdef WIN32
int euclid_signals_len = 2;
int euclid_signals[] = { SIGSEGV, SIGFPE };
#else
int euclid_signals_len = 3;
int euclid_signals[] = { SIGSEGV, SIGFPE, SIGBUS };
#endif

/* 
   signal names and explanatory messages 
*/
static char *SIGNAME[] = {
    "Unknown signal",
    "HUP (Hangup detected on controlling terminal or death of controlling process)",
    "INT: Interrupt from keyboard",
    "QUIT: Quit from keyboard",
    "ILL: Illegal Instruction",
    "TRAP",
    "ABRT: Abort signal",
    "EMT",
    "FPE (Floating Point Exception)",
    "KILL: Kill signal",
    "BUS (Bus Error, possibly illegal memory access)",
    "SEGV (Segmentation Violation (memory access out of range?))",
    "SYS",
    "PIPE: Broken pipe: write to pipe with no readers",
    "ALRM: Timer signal",
    "TERM: Termination signal",
    "URG",
    "STOP",
    "TSTP",
    "CONT",
    "CHLD"
};
#ifdef __cplusplus
}
#endif
#endif

