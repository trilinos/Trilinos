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

#ifndef NOX_PETSC_OPTIONS_H
#define NOX_PETSC_OPTIONS_H
  
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
  PetscTruth flg;  // Boolean flag
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
