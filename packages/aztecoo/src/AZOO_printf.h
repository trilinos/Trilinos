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

#ifndef _AZOO_printf_h_
#define _AZOO_printf_h_

#if defined(AztecOO_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The AztecOO package is deprecated"
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*********************************************************************/
/* Declarations for functions that will be passed as function-pointers
   and called from C code in Aztec. These are not user-functions.
*/
extern int AZOO_printf_out(const char* format, ...);
extern void AZOO_flush_out();

extern int AZOO_printf_err(const char* format, ...);
/*********************************************************************/

#ifdef __cplusplus
}

#include <iostream>

/** C++ users can use this function to set a destination stream
   for the output that Aztec would normally send to stdout.
*/
void AZOO_set_stream_out(std::ostream& ostrm);

/** C++ users can use this function to set a destination stream
   for the output that Aztec would normally send to stderr.
*/
void AZOO_set_stream_err(std::ostream& ostrm);

/** C++ users can use this function to clear out and err streams
*/
void AZOO_clear_streams();
#endif

#endif
