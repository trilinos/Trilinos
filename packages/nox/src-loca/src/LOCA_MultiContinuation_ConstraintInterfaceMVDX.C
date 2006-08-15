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

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H"

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterfaceMVDX::multiplyDX(
		       double alpha, 
		       const NOX::Abstract::MultiVector& input_x,
		       NOX::Abstract::MultiVector::DenseMatrix& result_p) const
{
  

  if (!isDXZero()) {
    const NOX::Abstract::MultiVector* dgdx = getDX();
    input_x.multiply(alpha, *dgdx, result_p);
  }
  else
    result_p.putScalar(0.0);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::ConstraintInterfaceMVDX::addDX(
		              Teuchos::ETransp transb,
			      double alpha, 
		              const NOX::Abstract::MultiVector::DenseMatrix& b,
			      double beta,
			      NOX::Abstract::MultiVector& result_x) const
{
  if (!isDXZero()) {
    const NOX::Abstract::MultiVector* dgdx = getDX();
    result_x.update(transb, alpha, *dgdx, b, beta);
  }
  else
    result_x.scale(beta);

  return NOX::Abstract::Group::Ok;
}
