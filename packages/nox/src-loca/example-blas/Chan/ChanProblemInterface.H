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

#ifndef CHANPROBLEMINTERFACE_H
#define CHANPROBLEMINTERFACE_H

#include "LOCA_BLAS_Group.H"

class ChanProblemInterface : public LOCA::BLAS::Interface {

public:
 
  ChanProblemInterface(int N, double a, double b);

  ~ChanProblemInterface() {};

  const NOX::BLAS::Vector& getInitialGuess();

  bool computeF(NOX::BLAS::Vector& f, const NOX::BLAS::Vector &x);
  
  bool computeJacobian(NOX::BLAS::Matrix& J, const NOX::BLAS::Vector & x);

  bool setParams(const LOCA::ParameterVector& p);

private:

  double source_term(double x);

  double source_deriv(double x);

private:

  NOX::BLAS::Vector initialGuess;
  double alpha;
  double beta;
  int n;

};

#endif
