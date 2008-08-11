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

#ifndef TEST_1DFEM_PREPOSTOPERATOR_H
#define TEST_1DFEM_PREPOSTOPERATOR_H

#include "NOX_Common.H"
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Utils.H"

class UserPrePostOperator : public NOX::Abstract::PrePostOperator {

public:

  //! Ctor.
  UserPrePostOperator(const NOX::Utils& u);

  //! Destructor.
  ~UserPrePostOperator();

  void runPreIterate(const NOX::Solver::Generic& solver);

  void runPostIterate(const NOX::Solver::Generic& solver);

  void runPreSolve(const NOX::Solver::Generic& solver);

  void runPostSolve(const NOX::Solver::Generic& solver);

  int getNumRunPreIterate() const { return numRunPreIterate; };

  int getNumRunPostIterate() const { return numRunPostIterate; };

  int getNumRunPreSolve() const { return numRunPreSolve; };

  int getNumRunPostSolve() const { return numRunPostSolve; };

protected:

  NOX::Utils utils;

  int numRunPreIterate;
  int numRunPostIterate;
  int numRunPreSolve;
  int numRunPostSolve;

};
#endif
