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

Options::Options(NOX::Parameter::List& params, int rank_) :
  rank(rank_)
{
  setOptions(params);
}

Options::~Options()
{
}


bool Options::setOptions(NOX::Parameter::List& nlParams)
{

  // First allow solution-type to be specified
  ierr = PetscOptionsHasName(PETSC_NULL,"-nox_trustregion_based",&flg);
         CHKERRQ(ierr);
  if(flg)
    nlParams.setParameter("Nonlinear Solver", "Trust Region Based");
  else // default
    // This is done to allow PetscOptions to register that this option was used
    ierr = PetscOptionsHasName(PETSC_NULL,"-nox_linesearch_based",&flg);
           CHKERRQ(ierr);
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
#ifdef WITH_PRERELEASE
    if( !strcmp(optionString, "nonlinearcg") )
      searchParams.setParameter("Method", "NonlinearCG");
#endif
  }
  else // default
    searchParams.setParameter("Method", "Full Step");

  // Now allow direction type to be specified
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  ierr = PetscOptionsGetString(PETSC_NULL,"-nox_direction_type",
               optionString, maxStringLength, &flg);CHKERRQ(ierr);
  if(flg)
  {
    if( !strcmp(optionString, "newton") )
      dirParams.setParameter("Method", "Newton");
    if( !strcmp(optionString, "steepest_descent") )
    {
      dirParams.setParameter("Method", "Steepest Descent");

      // Check to see if any steepest_descent options are set
      PetscTruth lflg;
      ierr = PetscOptionsGetString(PETSC_NULL,"-nox_sd_scaling_type",
                   optionString, maxStringLength, &lflg);CHKERRQ(ierr);
      if(lflg)
      {
        NOX::Parameter::List& sdParams = dirParams.sublist("Steepest Descent");
        if( !strcmp(optionString, "none") )
          sdParams.setParameter("Scaling Type", "None");
        else if( !strcmp(optionString, "2norm") )
          sdParams.setParameter("Scaling Type", "2-Norm");
        else if( !strcmp(optionString, "quadratic_model_min") )
          sdParams.setParameter("Scaling Type", "Quadratic Model Min");
        else 
        {
          if(rank == 0) cout << "WARNING: Unsupported Steepest Descent "
                             << "Scaling Type --> " << optionString << endl;
          sdParams.setParameter("Scaling Type", "None"); // default
        }
      }
    } 
#ifdef WITH_PRERELEASE
    if( !strcmp(optionString, "nonlinearcg") )
      dirParams.setParameter("Method", "Nonlinear CG");
    // Need to make provision for the following
      //NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
        //nlcgParams.setParameter("Restart Frequency", 2000);
        //nlcgParams.setParameter("Precondition", "On");
        //nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
        //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");
#endif
  }
  else // default
    dirParams.setParameter("Method", "Newton");

  // Now set output parameters via the "Printing" sublist 
  // These are hard-coded for now
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", rank);
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  printParams.setParameter("Output Information",
                        NOX::Utils::OuterIteration +
                        NOX::Utils::OuterIterationStatusTest +
                        NOX::Utils::InnerIteration +
                        NOX::Utils::Parameters +
                        NOX::Utils::Details +
                        NOX::Utils::Warning);

  return true;
}
