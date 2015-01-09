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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "Epetra_SerialComm.h"
#include "Newp_Jambo.h"
#include <sstream>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[]){

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  using namespace std;

  // Get an Epetra_Comm
  Epetra_SerialComm epetra_serial_comm;
  Epetra_Comm * epetra_comm;
  epetra_comm = dynamic_cast<Epetra_Comm*>(&epetra_serial_comm);

  //Create an ostream that Newp_Hello can write to and that we can read from
  stringbuf string_buf;
  streambuf * stream_buf;
  stream_buf = dynamic_cast<streambuf*>(&string_buf);
  iostream io_stream(stream_buf);
  ostream * o_stream;
  o_stream = dynamic_cast<ostream*>(&io_stream);
  
  //Create a Newp_Jambo to test
  Newp_Jambo new_package_jambo(*epetra_comm);
  new_package_jambo.Print(*o_stream);

  //Read from the io_stream
  char temp[83];
  io_stream.getline(temp, 83, 0);
    
  char * expected = "This will print out one line for each of the 1 processes \n\nJambo.  I am process 0\n";

  if(strcmp(temp, expected) != 0){
    cout << "Test Failed!" << endl << "     Got::" << strlen(temp) << "::" << temp << "::" << endl << "Expected::" << strlen(expected) << "::" << expected << "::" << endl;
    return 1;
  }
  cout << "Test passed!" << endl;
  #ifdef EPETRA_MPI
  MPI_Finalize();
  #endif
  return 0;
}
