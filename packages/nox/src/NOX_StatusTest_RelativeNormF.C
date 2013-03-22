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

std::ostream& NOX::StatusTest::RelativeNormF::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "F-Norm = " << Utils::sciformat(normF,3);
  stream << " < " << Utils::sciformat(tolerance * normF_0, 3);
  stream << " (" << Utils::sciformat(tolerance, 3);
  stream << " * " << Utils::sciformat(normF_0, 3);
  stream << ")" << std::endl;

  return stream;
}

