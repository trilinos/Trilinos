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
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

Anasazi::LOCASort::LOCASort(
 const Teuchos::RCP<LOCA::GlobalData>& global_data,
 const Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>& strategy_)
  : globalData(global_data),
    strategy(strategy_)
{
}

Anasazi::LOCASort::~LOCASort()
{
}

void
Anasazi::LOCASort::sort(std::vector<double>& evals, 
			                  Teuchos::RCP<std::vector<int> > perm, 
                        int n) const
{
  if (n == -1) {
    n = evals.size();
  }
  NOX::Abstract::Group::ReturnType res = strategy->sort(n, &evals[0], perm.get());
  globalData->locaErrorCheck->checkReturnType(res, "Anasazi::LOCASort::sort()");
}

void
Anasazi::LOCASort::sort(std::vector<double>& r_evals, 
			                  std::vector<double>& i_evals, 
			                  Teuchos::RCP<std::vector<int> > perm, 
                        int n) const
{
  if (n == -1) {
    n = r_evals.size();
  }
  NOX::Abstract::Group::ReturnType res = 
    strategy->sort(n, &r_evals[0], &i_evals[0], perm.get());
  globalData->locaErrorCheck->checkReturnType(res,
					      "Anasazi::LOCASort::sort()");
}
