// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
