
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

#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Newp_Jambo.h"
#include <iostream>
#include <sstream>
#include <istream>
#include <ostream>
#include <string>

int main(){
  using namespace std;
  Epetra_SerialComm esc;
  Epetra_Comm * ec;
  ec = dynamic_cast<Epetra_Comm*>(&esc);
  stringbuf stringb;
  streambuf * sb;
  sb = dynamic_cast<streambuf*>(&stringb);
  iostream io(sb);
  ostream * os;
  os = dynamic_cast<ostream*>(&io);
  Newp_Jambo npj(*ec);
  npj.Print(*os);
  char results[100];
  io.readsome(results, 99);
  char * expected1 = "This will print out one line for each of the 1 processes";
  char * expected2 = "Jambo.  I am process 0";
  if(strstr(results, expected1) == NULL || strstr(results, expected2) == NULL){
    cout << "Test failed!" << endl;
    cout << "Expected:" << endl << expected1 << endl << expected2 << endl;
    cout << "Got:" << endl << results << endl;
    return 1;
  }
  cout << "Test passed!" << endl;
  return 0;
}
