// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX.H" // NOX library containing Parameter List class
#include "NOX_Petsc_Options.H" // class definition

using namespace NOX::Petsc;

Options::Options()
{
}

Options::Options(NOX::Parameter::List& params)
{
  cout << "Inside Options......will now fill params !!" << endl << endl;
  setOptions(params);
}

Options::~Options()
{
}


bool Options::setOptions(NOX::Parameter::List& nlParams)
{
  cout << "Inside Options::setOptions......will now fill params !!" << endl << endl;
 
  // First allow solution-type to be specified
  ierr = PetscOptionsHasName(PETSC_NULL,"-nox_trustregion_based",&flg);CHKERRQ(ierr);
  if(flg)
    nlParams.setParameter("Nonlinear Solver", "Trust Region Based");
  else // default
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

  // Now allow linesearch type to be specified
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  ierr = PetscOptionsGetString(PETSC_NULL,"-nox_linesearch_type",
               optionString, maxStringLength, &flg);CHKERRQ(ierr);

  if(flg)
  {
    cout << "linesearch optionString --> " << optionString << endl;
    if( !strcmp(optionString, "full_step") )
      searchParams.setParameter("Method", "Full Step");
    if( !strcmp(optionString, "polynomial") )
      searchParams.setParameter("Method", "Polynomial");
    if( !strcmp(optionString, "backtrack") )
      searchParams.setParameter("Method", "Backtrack");
    if( !strcmp(optionString, "more_thuente") )
      searchParams.setParameter("Method", "More'-Thuente");
    if( !strcmp(optionString, "more_thuente2") )
      searchParams.setParameter("Method", "More'-Thuente2");
#ifdef WITH_PRERELEASE
    if( !strcmp(optionString, "nonlinearcg") )
      searchParams.setParameter("Method", "NonlinearCG");
#endif
  }
  else // default
    searchParams.setParameter("Method", "Full Step");

  nlParams.print(cout);
  
  return true;
}
