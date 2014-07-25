//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

// ----------   Includes   ----------
#include <iostream>
#include "Problem_Interface.H"

// ----------   User Defined Includes   ----------
#include "Brusselator.H"
#include "Epetra_RowMatrix.h"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(Brusselator& Problem) :
  problem(Problem)
{ }

Problem_Interface::~Problem_Interface()
{ }

bool Problem_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& FVec,
                       NOX::Epetra::Interface::Required::FillType fillType)
{
  return problem.evaluate(F_ONLY, &x, &FVec, NULL);
}

bool Problem_Interface::computeJacobian(const Epetra_Vector& x,
                    Epetra_Operator& Jac)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,
              dynamic_cast<Epetra_RowMatrix*>(&Jac),
              1.0, 0.0);
}

void Problem_Interface::setParameters(const LOCA::ParameterVector& params)
{
  problem.setParameters(params.getValue("alpha"),
            params.getValue("beta"),
            params.getValue("D1"),
            params.getValue("D2"));
}

bool Problem_Interface::computeShiftedMatrix(double alpha, double beta,
                         const Epetra_Vector& x,
                         Epetra_Operator& A)
{
  return problem.evaluate(MATRIX_ONLY, &x, NULL,
              dynamic_cast<Epetra_RowMatrix*>(&A), alpha, beta);
}

//-----------------------------------------------------------------------------

