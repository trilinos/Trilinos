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
                                                                                
#include "NOX.H"  // NOX headers
#include "NOX_LAPACK.H" // NOX LAPACK Interface headers
#include "NOX_TestError.H" // common file for testing 

#ifdef HAVE_MPI
#include <mpi.h>
#else 
#endif

int main(int argc, char *argv[]) {

  // Set up the printing utilities
  Teuchos::ParameterList noxParams;
  Teuchos::ParameterList& printParams = noxParams.sublist("Printing");
  printParams.set("Output Precision", 5);
  if (argc > 1) { 
    if (argv[1][0]=='-' && argv[1][1]=='v')
       printParams.set("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning +
			NOX::Utils::TestDetails);
    else
       printParams.set("Output Information", NOX::Utils::Error);
  }
  NOX::Utils printing(printParams);

  if (printing.isPrintType(NOX::Utils::TestDetails)) {
    cout << "Starting lapack/NOX_Group/NOX_Group.exe" << endl;
  }

  int status = 0;

  // Begin real testing here!
  if (status == 0) 
    cout << "Test passed!" << endl;
  else 
    cout << "Test failed!" << endl;

  // 0 is success
  return status;
}

/*
  end of file main.cc
*/
