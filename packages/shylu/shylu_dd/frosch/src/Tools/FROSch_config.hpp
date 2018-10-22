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

// ************************************************************************
//@HEADER

#ifndef FROSCH_CONFIG_HPP
#define FROSCH_CONFIG_HPP

#ifndef F77_BLAS_MANGLE
#cmakedefine F77_BLAS_MANGLE@F77_BLAS_MANGLE@
#endif

/* Options */

#cmakedefine HAVE_MPI

#cmakedefine HAVE_FROSCH_DEBUG

/* Optional Dependencies */

#cmakedefine HAVE_FROSCH_AMESOS

#cmakedefine HAVE_FROSCH_AMESOS2

#cmakedefine HAVE_FROSCH_AZTECOO

#cmakedefine HAVE_FROSCH_BELOS

#cmakedefine HAVE_FROSCH_EPETRA

#cmakedefine HAVE_FROSCH_EPETRAEXT

#cmakedefine HAVE_FROSCH_GALERI

#cmakedefine HAVE_FROSCH_TPETRA

/* Whether Tpetra is enabled with LocalOrdinal = int and GlobalOrdinal = int */

#cmakedefine HAVE_FROSCH_THYRA

#cmakedefine HAVE_FROSCH_STRATIMIKOS

#cmakedefine HAVE_MUELU_ZOLTAN

#cmakedefine HAVE_MUELU_ZOLTAN2





