// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "PyTrilinos_LinearProblem.hpp"
#include "Epetra_MultiVector.h"

namespace PyTrilinos
{

//////////////////////////////////////////////////////////////////////

LinearProblem::LinearProblem()
{
}

//////////////////////////////////////////////////////////////////////

LinearProblem::LinearProblem(const Teuchos::RCP< Epetra_RowMatrix > matrix,
                             const Teuchos::RCP< Epetra_MultiVector > x,
                             const Teuchos::RCP< Epetra_MultiVector > b) :
  Epetra_LinearProblem(matrix.get(), x.get(), b.get()),
  _matrix(matrix),
  _x(x),
  _b(b)
{
  _operator = Teuchos::rcp_dynamic_cast< Epetra_Operator >(_matrix);
}

//////////////////////////////////////////////////////////////////////

LinearProblem::LinearProblem(const Teuchos::RCP< Epetra_Operator > op,
                             const Teuchos::RCP< Epetra_MultiVector > x,
                             const Teuchos::RCP< Epetra_MultiVector > b) :
  Epetra_LinearProblem(op.get(), x.get(), b.get()),
  _operator(op),
  _x(x),
  _b(b)
{
  _matrix = Teuchos::rcp_dynamic_cast< Epetra_RowMatrix >(op, false);
}

//////////////////////////////////////////////////////////////////////

LinearProblem::LinearProblem(const LinearProblem & source) :
  Epetra_LinearProblem(source),
  _matrix(  source._matrix  ),
  _operator(source._operator),
  _x(       source._x       ),
  _b(       source._b       )
{
}

//////////////////////////////////////////////////////////////////////

LinearProblem::LinearProblem(const Epetra_LinearProblem & source) :
  Epetra_LinearProblem()
{
  // Get the linear problem pointers
  Epetra_RowMatrix   * matrix = source.GetMatrix();
  Epetra_Operator    * op     = source.GetOperator();
  Epetra_MultiVector * lhs    = source.GetLHS();
  Epetra_MultiVector * rhs    = source.GetRHS();

  // Copy pointers to the base class
  Epetra_LinearProblem::SetOperator(matrix);
  Epetra_LinearProblem::SetOperator(op    );
  Epetra_LinearProblem::SetLHS(     lhs   );
  Epetra_LinearProblem::SetRHS(     rhs   );

  // Convert the pointers to RCPs
  if (matrix) _matrix   = Teuchos::rcp(matrix, false);
  if (op    ) _operator = Teuchos::rcp(op    , false);
  if (lhs   ) _x        = Teuchos::rcp(lhs   , false);
  if (rhs   ) _b        = Teuchos::rcp(rhs   , false);
}

//////////////////////////////////////////////////////////////////////

LinearProblem::~LinearProblem()
{
}

//////////////////////////////////////////////////////////////////////

void LinearProblem::SetOperator(Teuchos::RCP< Epetra_RowMatrix > & matrix)
{
  _matrix   = matrix;
  _operator = Teuchos::rcp_dynamic_cast< Epetra_Operator >(matrix);
  Epetra_LinearProblem::SetOperator(matrix.get());
}

//////////////////////////////////////////////////////////////////////

void LinearProblem::SetOperator(Teuchos::RCP< Epetra_Operator > & op)
{
  _matrix   = Teuchos::rcp_dynamic_cast< Epetra_RowMatrix >(op, false);
  _operator = op;
  Epetra_LinearProblem::SetOperator(op.get());
}

//////////////////////////////////////////////////////////////////////

void LinearProblem::SetLHS(Teuchos::RCP< Epetra_MultiVector > & x)
{
  _x = x;
  Epetra_LinearProblem::SetLHS(x.get());
}

//////////////////////////////////////////////////////////////////////

void LinearProblem::SetRHS(Teuchos::RCP< Epetra_MultiVector > & b)
{
  _b = b;
  Epetra_LinearProblem::SetRHS(b.get());
}

//////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_RowMatrix > LinearProblem::GetMatrix() const
{
  return _matrix;
}

//////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_Operator > LinearProblem::GetOperator() const
{
  return _operator;
}

//////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_MultiVector > LinearProblem::GetLHS() const
{
  return _x;
}

//////////////////////////////////////////////////////////////////////

Teuchos::RCP< Epetra_MultiVector > LinearProblem::GetRHS() const
{
  return _b;
}

}
