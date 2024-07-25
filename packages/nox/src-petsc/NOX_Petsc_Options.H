// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_PETSC_OPTIONS_H
#define NOX_PETSC_OPTIONS_H

#include "petscversion.h"
#include "petscsys.h"
#include "NOX.H" // Needed for NOX::StatusTest classes

namespace NOX {
namespace Petsc {

  /*! \brief Class to convert solver options from command line (or command
   *  input file) into a form usable by NOX. */
  /*!  Details go here .....
   *   */
class Options {

public:

  //! Constructor
  Options();

  //! Constructor
  Options(Teuchos::ParameterList& params, int rank = 0);

  //! Destructor
  ~Options();

  //! Set parameter list based on options on command line or from file
  bool setOptions(Teuchos::ParameterList& params);

  //! Get reference to the final (combo) status test
  Teuchos::RCP<NOX::StatusTest::Combo> & getStatusTest();

private:

  enum {maxStringLength = 256};
  char optionString[maxStringLength];
#if  (PETSC_VERSION_MAJOR >= 3) || (PETSC_VERSION_MINOR >= 5)
  PetscBool flg;
#else
  PetscTruth flg;  // Boolean flag
#endif
  int ierr;   // Used in Petsc CHKERRQ macro
  int rank;   // proc ID used with output utils

  //! Supported convergence criteria
  Teuchos::RCP<NOX::StatusTest::MaxIters>    testMaxIters;
  Teuchos::RCP<NOX::StatusTest::NormF>       testNormF;
  Teuchos::RCP<NOX::StatusTest::NormUpdate>  testNormUpdate;
  Teuchos::RCP<NOX::StatusTest::Combo>       testCombo;

};
} // namespace Petsc
} // namespace NOX

#endif
