//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_OUTPUT_HPP
#define _FROSCH_OUTPUT_HPP

#include <ShyLU_DDFROSch_config.h>

// Width of indent before FROSch output
#ifndef FROSCH_OUTPUT_INDENT
    #define FROSCH_OUTPUT_INDENT 5
#endif

#ifndef FROSCH_ASSERT
    #define FROSCH_ASSERT(A,S) TEUCHOS_TEST_FOR_EXCEPTION(!(A),std::logic_error,S);
#endif

#ifndef FROSCH_WARNING
    #define FROSCH_WARNING(CLASS,VERBOSE,OUTPUT) if (VERBOSE) std::cerr << std::setw(FROSCH_OUTPUT_INDENT) << " " << CLASS << " : WARNING: " << OUTPUT << std::endl;
#endif

#ifndef FROSCH_NOTIFICATION
    #define FROSCH_NOTIFICATION(CLASS,VERBOSE,OUTPUT) if (VERBOSE) std::cout << std::setw(FROSCH_OUTPUT_INDENT) << " " << CLASS << " : NOTIFICATION: " << OUTPUT << std::endl;
#endif

#ifndef FROSCH_TEST_OUTPUT
    #define FROSCH_TEST_OUTPUT(COMM,VERBOSE,OUTPUT) COMM->barrier(); COMM->barrier(); COMM->barrier(); if (VERBOSE) std::cout << std::setw(FROSCH_OUTPUT_INDENT) << " " << OUTPUT << std::endl;
#endif

#endif
