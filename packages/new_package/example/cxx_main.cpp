//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

//
// hello_test
//
//  usage: 
//     hello_test
//
//  output:  
//     stdout
//
//  exits with 0 if test completed (does not imply that the test passed)
//  exits with -1 if command line options or file permissions are wrong 
//

//  #include "New_Package_config.h" - I suspect that Epatra_config.h has everything we need 
//  and hence obviates the need for New_Package_config.h
//  Or maybe not - we might still need HAVE_NEWP_SWAHILI 
//  
#include "Newp_Hello.h"
#ifdef HAVE_NEWP_SWAHILI
#include "Newp_Jambo.h"
#endif

main(int argc, char **argv)
{
  int exit_value = 0 ; 

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Newp_Hello Hello( Comm ) ; 
  Hello.Print( cout );
#ifdef HAVE_NEWP_SWAHILI
  Newp_Jambo Jambo( Comm ) ; 
  Jambo.Print( cout );
#endif

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  exit( exit_value ) ; 
}

  

