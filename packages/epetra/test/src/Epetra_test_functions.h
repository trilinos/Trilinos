/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_TEST_FUNCTIONS_H
#define EPETRA_TEST_FUNCTIONS_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_ConfigDefs.h"

class Epetra_Comm;
class Epetra_CrsMatrix;

namespace epetra_test {

/** Check through a list of C-style string arguments searching for a specified flag
    on proc 0. Case-sensitive search.
   return true on all procs if flag occurs in strargs, false on all procs if not.
*/
EPETRA_LIB_DLL_EXPORT bool global_check_for_flag_on_proc_0(const char* flag,
                                     int numargs, char** strargs,
                                     const Epetra_Comm& comm);

/** If macro EPETRA_MPI is defined, call MPI_Init and then return new Epetra_MpiComm.
    Otherwise, return new Epetra_SerialComm.
*/
EPETRA_LIB_DLL_EXPORT Epetra_Comm* create_comm(int argc, char** argv);

/** Check whether the two CrsMatrix arguments have the same size, structure and coefs.
*/
EPETRA_LIB_DLL_EXPORT bool compare_matrices(const Epetra_CrsMatrix& A, const Epetra_CrsMatrix& B);

}//namespace epetra_test

#endif
