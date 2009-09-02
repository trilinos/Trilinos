//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
/* The purpose of this file is to allow for 
informative error reporting that will not abort a test program after a single
error the way that an assert command does. */

/* The macro takes two arguments.  The first is the error code to be examined.  The second is an int that can be viewed either as an error count, or as a bool that simply indicates if any errors have occurred as of the time that the macro was run (zero -> no prior errors, non-zero -> prior errors). */

/* This is a rip off of the EPETRA_TEST_ERR function written by Mike Heroux for Epetra regression testing */

#ifndef _NOX_TEST_ERR_H_
#define _NOX_TEST_ERR_H_
#include "NOX_Common.H"
using namespace std;
// This function is to be used when first identifying an error.
#define NOX_TEST_ERR(a,b) { { int nox_testing_err = a; \
  if (nox_testing_err != 0) {\
    cerr << "Non zero error code " << nox_testing_err << \
       ", file: " << __FILE__ << ", line: " << __LINE__ << endl;\
 b+=1;\
  }\
  }\
}

#endif /*_NOX_TEST_ERR_H_ */
