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

#include "NOX.H" // NOX library containing Parameter List class
#include "NOX_Petsc_Options.H" // class definition

using namespace NOX::Petsc;

Options::Options()
{
}

Options::Options(Teuchos::ParameterList& params, int rank_) :
  rank(rank_)
{
  setOptions(params);
}

Options::~Options()
{
}


bool Options::setOptions(Teuchos::ParameterList& nlParams)
{

  // Set status tests if not already set
  if( Teuchos::is_null(testCombo) )
  {
    // Check for MaxIters option
    int maxIters;
    PetscTruth lflg;  // Needed to permit two ways of specification
    ierr = PetscOptionsGetInt(PETSC_NULL,"-snes_max_it", &maxIters, &flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-nox_conv_maxiters", &maxIters, &lflg);CHKERRQ(ierr);
    if(flg || lflg)
    {
      testMaxIters = Teuchos::rcp( new NOX::StatusTest::MaxIters(maxIters) );
      if( Teuchos::is_null(testCombo) )
        testCombo = Teuchos::rcp( new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, testMaxIters) );
      else
        testCombo->addStatusTest(testMaxIters);
    }
   
    // Check for (absolute) residual norm (L2-norm) tolerance
    double absResNorm;
    PetscReal petscVal;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-snes_atol", &petscVal, &flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-nox_conv_abs_res", &petscVal, &lflg);CHKERRQ(ierr);
    if(flg || lflg)
    {
      absResNorm = (double) petscVal;
      testNormF = Teuchos::rcp( new NOX::StatusTest::NormF(absResNorm) );
      if( Teuchos::is_null(testCombo) )
        testCombo = Teuchos::rcp( new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, testNormF) );
      else
        testCombo->addStatusTest(testNormF);
    }

    // Check for update norm (L2-norm) tolerance
    double absUpdateNorm;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-snes_stol", &petscVal, &flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL,"-nox_conv_update", &petscVal, &lflg);CHKERRQ(ierr);
    if(flg || lflg)
    {
      absUpdateNorm = (double) petscVal;
      testNormUpdate = Teuchos::rcp( new NOX::StatusTest::NormUpdate(absUpdateNorm) );
      if( Teuchos::is_null(testCombo) )
        testCombo = Teuchos::rcp( new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, testNormUpdate) );
      else
        testCombo->addStatusTest(testNormUpdate);
    }

    // Finally, provide a default test if none specified
    if( Teuchos::is_null(testCombo) ) // No tests specified by the uesr
    {
      assert( Teuchos::is_null(testMaxIters) );
      testMaxIters = Teuchos::rcp( new NOX::StatusTest::MaxIters(20) );
      assert( Teuchos::is_null(testNormF) );
      testNormF = Teuchos::rcp( new NOX::StatusTest::NormF(1.e-12) );
      testCombo = Teuchos::rcp( new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, testMaxIters, testNormF) );
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
          if(rank == 0) std::cout << "WARNING: Unsupported Steepest Descent "
                             << "Scaling Type --> " << optionString << std::endl;
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

Teuchos::RCP<NOX::StatusTest::Combo> & 
Options::getStatusTest()
{
  return testCombo;
}
