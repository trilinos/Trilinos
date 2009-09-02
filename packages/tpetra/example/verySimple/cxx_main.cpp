//@HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

#include "Tpetra_OutputObject.hpp"
#include "Tpetra_SerialComm.hpp"
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_ElementSpace.hpp"

int main(int argc, char *argv[])
{
  cout << "*** Starting verySimple example..." << endl;

  cout << Tpetra::version() << endl;

  // test OutputObject
  cout << "*** Creating OutputObject..." << endl;
	Tpetra::OutputObject obj1;

	obj1.newPrint(cout); cout << endl;

	Teuchos::RCP<Tpetra::OutputManager> om = Teuchos::rcp (new Tpetra::OutputManager());
	om->setVerbosity(Tpetra::Signature+Tpetra::Summary);
	obj1.setOutputManager(om);
	obj1.newPrint(cout); cout << endl;

	// test Vector object
	int length = 10;
        int indexBase = 0;
	cout << "*** Creating Vector of length " << length << endl;
	const Tpetra::SerialPlatform <int, int> platformE;
	const Tpetra::SerialPlatform <int, double> platformV;
	Tpetra::ElementSpace<int> elementspace(length, indexBase, platformE);
	Tpetra::VectorSpace<int, double> vectorspace(elementspace, platformV);
	Tpetra::Vector<int, double> v1(vectorspace);
	cout << "Created v1, default constructor" << endl;

	cout << v1 << endl;

	v1.setOutputManager(om);
	cout << v1 << endl;

	cout << "*** Finished." << endl;

  return(0);
}

