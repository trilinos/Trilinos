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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX.H" // NOX library containing Parameter List class
#include "NOX_Petsc_Options.H" // class definition

using namespace NOX::Petsc;

Options::Options()
{
}

Options::Options(Teuchos::ParameterList& params, int rank_) :
  rank(rank_),
  testMaxIters(0),
  testNormF(0),
  testNormUpdate(0),
  testCombo(0)
{
  setOptions(params);
}

Options::~Options()
{
}


bool Options::setOptions(Teuchos::ParameterList& nlParams)
{

  // Set status tests if not already set
  if(!testCombo)
  {
    // Check for MaxIters option
    int maxIters;
    PetscTruth lflg;  // Needed to permit two ways of specification
    ierr = PetscOptionsGetInt(PETSC_NULL,"-snes_max_it",
                 &maxIters, &flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-nox_conv_maxiters",
                 &maxIters, &lflg);CHKERRQ(ierr);
    if(flg || lflg)
    {
      testMaxIters = new NOX::StatusTest::MaxIters(maxIters);
      if(!testCombo)
        testCombo = new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                               *testMaxIters);
      else
        testCombo->addStatusTest(*testMaxIters);
    }
   
    // Check for (absolute) residual norm (L2-norm) tolerance
    double absResNorm;
    PetscReal petscVal;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-snes_atol",
                 &petscVal, &flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-nox_conv_abs_res",
                 &petscVal, &lflg);CHKERRQ(ierr);
    if(flg || lflg)
    {
      absResNorm = (double) petscVal;
      testNormF = new NOX::StatusTest::NormF(absResNorm);
      if(!testCombo)
        testCombo = new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                               *testNormF);
      else
        testCombo->addStatusTest(*testNormF);
    }

    // Check for update norm (L2-norm) tolerance
    double absUpdateNorm;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-snes_stol",
                 &petscVal, &flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-nox_conv_update",
                 &petscVal, &lflg);CHKERRQ(ierr);
    if(flg || lflg)
    {
      absUpdateNorm = (double) petscVal;
      testNormUpdate = new NOX::StatusTest::NormUpdate(absUpdateNorm);
      if(!testCombo)
        testCombo = new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                               *testNormUpdate);
      else
        testCombo->addStatusTest(*testNormUpdate);
    }

    // Finally, provide a default test if none specified
    if(!testCombo) // No tests specified by the uesr
    {
      assert( testMaxIters == 0);
      testMaxIters = new NOX::StatusTest::MaxIters(20);
      assert( testNormF == 0);
      testNormF = new NOX::StatusTest::NormF(1.e-12);
      testCombo = new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
                                             *testMaxIters, *testNormF);
    }
    

  } // End of StatusTest construction


  // Allow solution-type to be specified
  ierr = PetscOptionsHasName(PETSC_NULL,"-nox_trustregion_based",&flg);
         CHKERRQ(ierr);
  if(flg)
    nlParams.set("Nonlinear Solver", "Trust Region Based");
  else // default
    // This is done to allow PetscOptions to register that this option was used
    ierr = PetscOptionsHasName(PETSC_NULL,"-nox_linesearch_based",&flg);
           CHKERRQ(ierr);
    nlParams.set("Nonlinear Solver", "Line Search Based");

  // Now allow linesearch type to be specified
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  ierr = PetscOptionsGetString(PETSC_NULL,"-nox_linesearch_type",
               optionString, maxStringLength, &flg);CHKERRQ(ierr);
  if(flg)
  {
    if( !strcmp(optionString, "full_step") )
      searchParams.set("Method", "Full Step");
    if( !strcmp(optionString, "polynomial") )
      searchParams.set("Method", "Polynomial");
    if( !strcmp(optionString, "backtrack") )
      searchParams.set("Method", "Backtrack");
    if( !strcmp(optionString, "more_thuente") )
      searchParams.set("Method", "More'-Thuente");
#ifdef WITH_PRERELEASE
    if( !strcmp(optionString, "nonlinearcg") )
      searchParams.set("Method", "NonlinearCG");
#endif
  }
  else // default
    searchParams.set("Method", "Full Step");

  // Now allow direction type to be specified
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  ierr = PetscOptionsGetString(PETSC_NULL,"-nox_direction_type",
               optionString, maxStringLength, &flg);CHKERRQ(ierr);
  if(flg)
  {
    if( !strcmp(optionString, "newton") )
      dirParams.set("Method", "Newton");
    if( !strcmp(optionString, "steepest_descent") )
    {
      dirParams.set("Method", "Steepest Descent");

      // Check to see if any steepest_descent options are set
      PetscTruth lflg;
      ierr = PetscOptionsGetString(PETSC_NULL,"-nox_sd_scaling_type",
                   optionString, maxStringLength, &lflg);CHKERRQ(ierr);
      if(lflg)
      {
        Teuchos::ParameterList& sdParams = dirParams.sublist("Steepest Descent");
        if( !strcmp(optionString, "none") )
          sdParams.set("Scaling Type", "None");
        else if( !strcmp(optionString, "2norm") )
          sdParams.set("Scaling Type", "2-Norm");
        else if( !strcmp(optionString, "quadratic_model_min") )
          sdParams.set("Scaling Type", "Quadratic Model Min");
        else 
        {
          if(rank == 0) cout << "WARNING: Unsupported Steepest Descent "
                             << "Scaling Type --> " << optionString << endl;
          sdParams.set("Scaling Type", "None"); // default
        }
      }
    } 
#ifdef WITH_PRERELEASE
    if( !strcmp(optionString, "nonlinearcg") )
      dirParams.set("Method", "Nonlinear CG");
    // Need to make provision for the following
      //Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
        //nlcgParams.set("Restart Frequency", 2000);
        //nlcgParams.set("Precondition", "On");
        //nlcgParams.set("Orthogonalize", "Polak-Ribiere");
        //nlcgParams.set("Orthogonalize", "Fletcher-Reeves");
#endif
  }
  else // default
    dirParams.set("Method", "Newton");

  // Now set output parameters via the "Printing" sublist 
  // These are hard-coded for now
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", rank);
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information",
                        NOX::Utils::OuterIteration +
                        NOX::Utils::OuterIterationStatusTest +
                        NOX::Utils::InnerIteration +
                        NOX::Utils::Parameters +
                        NOX::Utils::Details +
                        NOX::Utils::Warning);

  return true;
}

NOX::StatusTest::Combo& Options::getStatusTest()
{
  return *testCombo;
}
