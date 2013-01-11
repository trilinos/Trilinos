/*
//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
*/

#ifndef CONTINUATION_MANAGER_H
#define CONTINUATION_MANAGER_H

// Trilinos headers
#include "Epetra_Time.h"
#include "LOCA_Stepper.H"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"
#include "LOCA_Epetra_Interface_xyzt.H"
#include "EpetraExt_MultiComm.h"

// DefCont headers
#include "ProblemLOCAPrototype.H"

class ContinuationManager {

  public:

    // Print or not a solution file
    enum SolutionFileAttribute { NoPrint, Print };

    // Constructor
    ContinuationManager( 
	//const Teuchos::RCP< Epetra_MpiComm > & aComm ,
	const Teuchos::RCP< Epetra_Comm > & aComm ,
	const std::string & taskFileName);

    // Destructor
    ~ContinuationManager();
      
    // Instantiante the LOCA stepper
    bool BuildLOCAStepper();

    // Instantiante the LOCA stepper for periodic (xyzt) problems
    bool BuildLOCAPeriodicStepper(const Teuchos::RCP<EpetraExt::MultiComm> globalComm);

    // Perform some checks on the LOCA stepper
    bool ValidateLOCAStepper();

    // Start the LOCA stepper
    bool RunLOCAStepper();

    // Print Statistics
    bool PrintLOCAStepperStatistics(); 

    // Returns the name of the actual solution file
    std::string GetSolutionFileName() const;

    // Returns the name of the continuation file
    std::string GetContinuationFileName() const;

    // Returns the stepper id 
    int GetStepID() const;

    // Getting the Task list it will 
    // be shared with a ProblemLOCAPrototype problem
    Teuchos::RCP <Teuchos::ParameterList> GetTaskList() const;

    // Determine whether a solution file is to be printed
    SolutionFileAttribute GetSolutionFileAttribute() const;

    // Set the ProblemLOCAPrototype problem 
    bool SetLOCAProblem( 
	const Teuchos::RCP <ProblemLOCAPrototype> & aProblem );

    // Set the ProblemLOCAPrototype problem 
    bool SetLOCAConstraint( 
	const Teuchos::RCP <LOCA::MultiContinuation::ConstraintInterface> &
	anInterfaceConstraint);

    // Setting the task list
    bool SetTaskList( Teuchos::RCP <Teuchos::ParameterList> aTaskList );

    // Makes a std::string from a given integer
    std::string StringifyInt( const int & intNumber , 
	const int & digits) const;

  private:

    // The communicator
    //Teuchos::RCP <Epetra_MpiComm> comm;
    Teuchos::RCP <Epetra_Comm> comm;

    // The continuation file name
    std::string continuationFileName;

    // The initial step label
    int iniStepLabel;

    // Intefrace to the eventual constraint
    Teuchos::RCP <LOCA::MultiContinuation::ConstraintInterface> interfaceConstraint;

//    // Is it a constrained problem?
//    bool isConstrainedProblem;

    // The LOCA global data
    Teuchos::RCP<LOCA::GlobalData> locaGlobalData;

    // The LOCA stepper
    Teuchos::RCP <LOCA::Stepper> locaStepper;

    // The LOCA stepper status
    LOCA::Abstract::Iterator::IteratorStatus locaStepperStatus;

    // Maximum number of steps allowed (iniStepLabel included)
    int maxAllowedSteps;

    // The output Directory
    std::string outputDir;

    // The ProblemLOCAPrototype problem
    Teuchos::RCP <ProblemLOCAPrototype> problem;

    // The extension for each solution (i.e. ".vtk")
    std::string solutionFilesExtension;

    // The flag to be assigned to each solution (i.e. "step_")
    std::string solutionFilesPrefix;

    // The task list
    Teuchos::RCP <Teuchos::ParameterList> taskList;

    // Time counter
    Epetra_Time timeCounter;

};

#endif
