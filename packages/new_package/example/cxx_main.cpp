//@HEADER
// ***********************************************************************
//
//                     New_Package Example Package
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
//

// hello_test
//
//  usage: 
//     hello_test
//
//  output:  
//     prints a summary line and one line "Hello" for each process to standard out
//     If --enable-newp_swahili is set on the configure line:
//        prints a summary line and one line "Jambo" for each process to standard out
//
#include "Newp_Hello.h"
#ifdef HAVE_NEWP_SWAHILI
#include "Newp_Jambo.h"
#endif
#include "New_Package_Version.h"

int main(int argc, char **argv)
{
  //
  //  If --enable-mpi, an MPI communicator is used, otherwise a serial
  //  stub communicator is used.  
  //
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  //
  //  Print out a summary line followed by a "Hello" line from each process
  //

  if (Comm.MyPID()==0)
    cout << New_Package_Version() << endl << endl;

  Newp_Hello Hello( Comm ) ; 
  Hello.Print( cout );


  //
  //  If --enable-newp_swahili is set, HAVE_NEWP_SWAHILI is set in 
  //    New_Package_config.h which is included by Newp_Hello.h and hence:
  //      Print out a summary line followed by a "Jambo" line from each process
  //
#ifdef HAVE_NEWP_SWAHILI
  Newp_Jambo Jambo( Comm ) ; 
  Jambo.Print( cout );
#endif

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  return 0;
}

  

