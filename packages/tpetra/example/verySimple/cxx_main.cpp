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

#include "Tpetra_Object.hpp"
#include "Tpetra_SerialComm.hpp"
#include "Tpetra_SerialPlatform.hpp"

int main(int argc, char *argv[])
{
	cout << "*** Starting verySimple example..." << endl;

	// test Object
	cout << "*** Testing Object..." << endl;
  Tpetra::Object obj1;
	Tpetra::Object obj2("obj2");

	int temp1 = obj1.getTracebackMode();
	cout << obj1 << endl;

	// test SerialComm
	cout << "*** Testing SerialComm..." << endl;
	Tpetra::SerialComm<double, int> comm1;
	int temp2 = comm1.getNumImages();
	cout << comm1 << endl;

	// test Platform
	cout << "*** Testing SerialPlatform..." << endl;
	Tpetra::SerialPlatform<int, int> platform1;
	cout << platform1 << endl;

	cout << "*** Finished." << endl;

  return(0);
}

