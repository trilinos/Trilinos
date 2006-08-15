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

#include "Anasazi_LOCA_Sort.H"
#include "LOCA_EigenvalueSort_Strategies.H"

Anasazi::LOCASort::LOCASort(
 const Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy>& strategy_)
  : strategy(strategy_)
{
}

Anasazi::LOCASort::~LOCASort()
{
}

Anasazi::ReturnType
Anasazi::LOCASort::sort(Anasazi::Eigensolver<double,
			                     Anasazi::LOCASort::MV,
			                     Anasazi::LOCASort::OP>* solver, 
			int n, double* evals, std::vector<int>* perm) const
{
  NOX::Abstract::Group::ReturnType res = strategy->sort(n, evals, perm);
  if (res == NOX::Abstract::Group::Ok)
    return Anasazi::Ok;
  else if (res == NOX::Abstract::Group::NotDefined)
    return Anasazi::Undefined;
  else
    return Anasazi::Failed;
}

Anasazi::ReturnType
Anasazi::LOCASort::sort(Anasazi::Eigensolver<double,
			                     Anasazi::LOCASort::MV,
			                     Anasazi::LOCASort::OP>* solver, 
			int n, double* r_evals, double* i_evals, 
			std::vector<int>* perm) const
{
  NOX::Abstract::Group::ReturnType res = 
    strategy->sort(n, r_evals, i_evals, perm);
  if (res == NOX::Abstract::Group::Ok)
    return Anasazi::Ok;
  else if (res == NOX::Abstract::Group::NotDefined)
    return Anasazi::Undefined;
  else
    return Anasazi::Failed;
}
