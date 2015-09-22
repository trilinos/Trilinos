/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "sig_dh.h"
#include "Parser_dh.h"


#undef __FUNC__
#define __FUNC__ "sigHandler_dh"
void
sigHandler_dh (int sig)
{
  fprintf (stderr, "\n[%i] Euclid Signal Handler got: %s\n", myid_dh,
	   SIGNAME[sig]);
  fprintf (stderr,
	   "[%i] ========================================================\n",
	   myid_dh);
  fprintf (stderr,
	   "[%i] function calling sequence that led to the exception:\n",
	   myid_dh);
  fprintf (stderr,
	   "[%i] ========================================================\n",
	   myid_dh);
  printFunctionStack (stderr);
  fprintf (stderr, "\n\n");

  if (logFile != NULL)
    {
      fprintf (logFile, "\n[%i] Euclid Signal Handler got: %s\n", myid_dh,
	       SIGNAME[sig]);
      fprintf (logFile,
	       "[%i] ========================================================\n",
	       myid_dh);
      fprintf (logFile,
	       "[%i] function calling sequence that led to the exception:\n",
	       myid_dh);
      fprintf (logFile,
	       "[%i] ========================================================\n",
	       myid_dh);
      printFunctionStack (logFile);
      fprintf (logFile, "\n\n");
    }

  EUCLID_EXIT;
}

#undef __FUNC__
#define __FUNC__ "sigRegister_dh"
void
sigRegister_dh ()
{
  if (Parser_dhHasSwitch (parser_dh, "-sig_dh"))
    {
      int i;
      for (i = 0; i < euclid_signals_len; ++i)
	{
	  signal (euclid_signals[i], sigHandler_dh);
	}
    }
}
