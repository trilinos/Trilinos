
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

#include "EpetraExt_MatlabEngine.h"
#include "Epetra_Comm.h"

using namespace EpetraExt;
namespace EpetraExt {

//=============================================================================
MatlabEngine::MatlabEngine (const Epetra_Comm& Comm):Comm_(Comm) {

    // MATLAB engOpen, to construct the MATLAB engine

    MyPID_ = Comm_.MyPID();
    if (MyPID_ == 0) {
	cout << "Hello.  I am matlab " << endl ; 
	// Gentlemen, start your engines ...
	Engine_ = engOpen ((char *) 0) ;
    }

} 

//=============================================================================
MatlabEngine::~MatlabEngine (void) {

    // MATLAB engClose, to destruct the MATLAB engine

    if (MyPID_ == 0) {
	cout << "Goodbye.  I was matlab " << endl ; 

	int result = engClose (Engine_) ;

	if (result != 0)
	{
	    cout << "That was bad.  engClose failed." << endl ; 
	}

    }

}

//=======================================================================
void MatlabEngine::EvalString (char* command) const {

    // send a string command to the MATLAB engine
    if (MyPID_ == 0) {

	cout << "Sending command to matlab:" << command << endl ;

	int result = engEvalString (Engine_, command) ;

	if (result != 0)
	{
	    cout << "That was bad.  engEvalString failed." << endl ; 
	}

    }

}
} // namespace EpetraExt
