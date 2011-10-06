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

#include "NOX_Common.H"
#include "NOX_StatusTest_RelativeNormF.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

NOX::StatusTest::RelativeNormF::
RelativeNormF(double in_tolerance, bool in_scale_norms_by_vector_length, 
	      const NOX::Utils* u) :
  tolerance(in_tolerance),
  normF_0(0.0),
  normF(0.0),
  scale_norms_by_vector_length(in_scale_norms_by_vector_length)
{
  if (u != NULL)
    utils = *u;
}

NOX::StatusTest::StatusType NOX::StatusTest::RelativeNormF::
checkStatus(const NOX::Solver::Generic& problem,
	    NOX::StatusTest::CheckType checkType)
{ 
  // NOTE: This algorithm assumes a 2-norm!

  // On initial iteration, compute initial norm F
  if (problem.getNumIterations() == 0) {
    normF_0 = problem.getSolutionGroup().getF().norm(NOX::Abstract::Vector::TwoNorm);

    if (scale_norms_by_vector_length)
      normF_0 /= std::sqrt(Teuchos::as<double>(problem.getSolutionGroup().getF().length()));
  }

  if (checkType == NOX::StatusTest::None)
  {
    normF = -1.0;
    status = Unevaluated;
  }
  else
    {
      normF = problem.getSolutionGroup().getF().norm(NOX::Abstract::Vector::TwoNorm);

      if (scale_norms_by_vector_length)
	normF /= std::sqrt(Teuchos::as<double>(problem.getSolutionGroup().getF().length()));

      status = (normF < tolerance * normF_0) ? Converged : Unconverged;
  }

  return status;
}

NOX::StatusTest::StatusType NOX::StatusTest::RelativeNormF::getStatus() const
{
  return status;
}

ostream& NOX::StatusTest::RelativeNormF::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "F-Norm = " << Utils::sciformat(normF,3);
  stream << " < " << Utils::sciformat(tolerance * normF_0, 3);
  stream << " (" << Utils::sciformat(tolerance, 3);
  stream << " * " << Utils::sciformat(normF_0, 3);
  stream << ")" << endl;

  return stream;
}

