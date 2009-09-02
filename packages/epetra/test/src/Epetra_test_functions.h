/*
************************************************************************

              Epetra: Linear Algebra Services Package
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov)

************************************************************************
*/

#ifndef EPETRA_TEST_FUNCTIONS_H
#define EPETRA_TEST_FUNCTIONS_H

class Epetra_Comm;
class Epetra_CrsMatrix;

namespace epetra_test {

/** Check through a list of C-style string arguments searching for a specified flag
    on proc 0. Case-sensitive search.
   return true on all procs if flag occurs in strargs, false on all procs if not.
*/
bool global_check_for_flag_on_proc_0(const char* flag,
                                     int numargs, char** strargs,
                                     const Epetra_Comm& comm);

/** If macro EPETRA_MPI is defined, call MPI_Init and then return new Epetra_MpiComm.
    Otherwise, return new Epetra_SerialComm.
*/
Epetra_Comm* create_comm(int argc, char** argv);

/** Check whether the two CrsMatrix arguments have the same size, structure and coefs.
*/
bool compare_matrices(const Epetra_CrsMatrix& A, const Epetra_CrsMatrix& B);

}//namespace epetra_test

#endif

