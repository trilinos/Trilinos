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

// Author: Ian Karlin ikarlin@sandia.gov 05-22-2008

#ifndef EPETRA_OSKIUTILS_H
#define EPETRA_OSKIUTILS_H

extern "C" {
#include "oski/oski.h"
}

//! Epetra_OskiUtils:  The Epetra OSKI Class to handle all operations that do not involve the use of a matrix, vector, error or permutation object.
/*! The Epetra_OskiUtils class is a helper class used to call OSKI functions that do not use matrix, vector, error or permutation objects.
	  It provides an interface to access the initialization and finalize routines of OSKI.

		All functions are public to allow access to methods needed by programs using OSKI.
		There are no data members of the class as all data is kept in the matrix, vector, 
		multi-vector, error and permutation classes.
*/

class Epetra_OskiUtils {
  public:
	//! @name Constructors/Destructor
        //@{        
	//! Default Constructor
	Epetra_OskiUtils();

	//! Destructor
	virtual ~Epetra_OskiUtils();
 	//@}

	//! @name Start/End
	//@{
        //! Initializes OSKI
	/*! Calls the OSKI routine to initialize the use of OSKI.  This routine is required before
	    OSKI can be used.
	*/	
  	void Init();

	//! Finalizes the use of OSKI
	/*! When done using OSKI this routine performs cleanup operations.  While not strictly required
	    it is highly recommended to be called when OSKI is no longer being used.
	*/
	void Close();
	//@}
};

#endif /* EPETRA_OSKIUTILS_H */
