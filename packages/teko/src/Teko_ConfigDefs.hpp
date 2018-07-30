/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

#ifndef __Teko_ConfigDefs_hpp__
#define __Teko_ConfigDefs_hpp__

#include "Teko_Config.h"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map_decl.hpp"

namespace Teko {

typedef double ST;
typedef int LO;
#ifdef HAVE_Teko_USE_LONGLONG_GO
typedef long long GO;
#elif defined(HAVE_TPETRA_INST_INT_INT)
typedef int GO;
#elif defined(HAVE_TPETRA_INST_INT_LONG)
typedef long GO;
#elif defined(HAVE_TPETRA_INST_INT_LONG_LONG)
typedef long long GO;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED_LONG)
typedef unsigned long GO;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED)
typedef unsigned GO;
#else
#  error "Teko wants to use a GlobalOrdinal type which Tpetra does not enable.  None of the following types are enabled: int, long, long long, unsigned long, unsigned."
#endif
typedef Tpetra::Map<>::node_type NT;

}

#endif
