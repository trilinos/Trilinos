/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "Epetra_Comm.h"
enum SparseSolverType { DSCPACK, UMFPACK, 
			SCALAPACK, MUMPS, TAUCS, PARDISO, CSS, PARAKLETE, KLU,
                        LAPACK, SUPERLUDIST, SUPERLU } ; 
enum AMESOS_MatrixType { AMESOS_Serial, AMESOS_Distributed } ; 

int Amesos_TestSolver( Epetra_Comm &Comm, char *matrix_file, 
		 SparseSolverType SparseSolver,
		 bool tranpose, int special, AMESOS_MatrixType MatrixType );

int Amesos_TestMultiSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special, AMESOS_MatrixType MatrixType );

int Amesos_TestMrhsSolver( Epetra_Comm &Comm, char *matrix_file, int numsolves,
		      SparseSolverType SparseSolver,
		      bool tranpose, int special, AMESOS_MatrixType MatrixType );


#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif
