// $Id$ 
// $Source$ 

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
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> & getStatusTest();

private:

  enum {maxStringLength = 256};
  char optionString[maxStringLength];
  PetscTruth flg;  // Boolean flag
  int ierr;   // Used in Petsc CHKERRQ macro
  int rank;   // proc ID used with output utils
 
  //! Supported convergence criteria
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters>    testMaxIters;
  Teuchos::RefCountPtr<NOX::StatusTest::NormF>       testNormF;
  Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate>  testNormUpdate;
  Teuchos::RefCountPtr<NOX::StatusTest::Combo>       testCombo;
  
};
} // namespace Petsc
} // namespace NOX

#endif
