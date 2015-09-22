/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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


#include <stdio.h>
#include <stdarg.h>

#include <az_aztec.h>

/** The following function pointers are initialized to 0 (NULL), and
   remain 0 unless a C++ user (such as the AztecOO class) calls the
   functions AZ_set_stream_out() and/or AZ_set_stream_err() which are
   in AZOO_printf.cpp.
*/
int (*azoo_printf_out)(const char* format, va_list aptr) = 0;
void (*azoo_flush_out)() = 0;
int (*azoo_printf_err)(const char* format, va_list aptr) = 0;




/** Most Aztec code (C code) calls the following function in place of
  printf. This function simply forwards to vprintf *UNLESS* the
  function AZOO_set_stream_out() has been used to set a C++ std::ostream
  as the target for Aztec's printf output. AZOO_set_stream_out()
  sets the function-pointer azoo_printf_out, which can then be
  called by this (AZ_printf_out) function.
*/
int AZ_printf_out(const char* format, ...)
{
  va_list aptr;
  va_start(aptr, format);

  if (azoo_printf_out != 0) {
    return( azoo_printf_out(format, aptr) );
  }

  return( vprintf(format, aptr) );
}

/** Most Aztec code (C code) calls the following function in place of
  fflush(stdout). This function simply forwards to fflush *UNLESS* the
  function AZOO_set_stream_out() has been used to set a C++ std::ostream
  as the target for Aztec's printf output. AZOO_set_stream_out()
  sets the function-pointer azoo_flush_out, which can then be
  called by this (AZ_flush_out) function.
*/
void AZ_flush_out()
{
  if (azoo_flush_out != 0) {
    azoo_flush_out();
    return;
  }

  fflush(stdout);
}

/** Most Aztec code (C code) calls the following function in place of
  fprintf(stderr,...). This function simply forwards to vfprintf *UNLESS* the
  function AZOO_set_stream_err() has been used to set a C++ std::ostream
  as the target for Aztec's stderr output. AZOO_set_stream_err()
  sets the function-pointer azoo_printf_err, which can then be
  called by this (AZ_printf_err) function.
*/
int AZ_printf_err(const char* format, ...)
{
  va_list aptr;
  va_start(aptr, format);

  if (azoo_printf_err != 0) {
    return( azoo_printf_err(format, aptr) );
  }

  return( vfprintf(stderr, format, aptr) );
}

