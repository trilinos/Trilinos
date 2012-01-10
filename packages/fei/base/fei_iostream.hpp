/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_iostream_hpp_
#define _fei_iostream_hpp_

#include "fei_macros.hpp"

//
//The stuff in this file somewhat protects us from the fact that some
//platforms may not put stream-related stuff in the std namespace,
//even though most do.
//These days (2007) perhaps all platforms do put everything in std and
//perhaps no platforms still have iostream.h without having <iostream>, etc.
//But historically we've had to account for these possibilities and I see
//little to be gained from removing this flexibility at this point.
//
//The basic mechanism here is to use macros that are defined differently
//for certain situations. An alternative approach would be to import
//symbols into our namespace, but we choose not to do that since it is
//a well-known sin to perform namespace pollution from within a header.
//

#ifdef FEI_HAVE_IOSTREAM
#include <iostream>
#define FEI_OSTREAM std::ostream
#define FEI_ISTREAM std::istream
#define FEI_COUT std::cout
#define FEI_ENDL std::endl
#elif defined(FEI_HAVE_IOSTREAM_H)
#include <iostream.h>
#define FEI_OSTREAM ostream
#define FEI_ISTREAM istream
#define FEI_COUT cout
#define FEI_ENDL endl
#else
#error "must have <iostream> or <iostream.h>"
#endif
//endif for ifdef FEI_HAVE_IOSTREAM

#ifdef FEI_HAVE_IOMANIP
#include <iomanip>
#elif defined (FEI_HAVE_IOMANIP_H)
#include <iomanip.h>
#endif

#ifdef FEI_HAVE_STD_IOS_FMTFLAGS

#define IOS_FMTFLAGS std::ios_base::fmtflags
#define IOS_SCIENTIFIC std::ios_base::scientific
#define IOS_FLOATFIELD std::ios_base::floatfield
#define IOS_FIXED std::ios_base::fixed
#define IOS_APP std::ios_base::app
#define IOS_OUT std::ios_base::out

#else

#define IOS_FMTFLAGS long
#define IOS_SCIENTIFIC ios::scientific
#define IOS_FLOATFIELD ios::floatfield
#define IOS_FIXED ios::fixed
#define IOS_APP ios::app
#define IOS_OUT ios::out

#ifdef FEI_IOS_FMTFLAGS
#undef IOS_FMTFLAGS
#define IOS_FMTFLAGS ios::fmtflags
#endif

#endif
//endif for ifdef FEI_HAVE_STD_IOS_FMTFLAGS

#include <fei_console_ostream.hpp>
#endif
//endif for ifndef _fei_iostream_hpp_


