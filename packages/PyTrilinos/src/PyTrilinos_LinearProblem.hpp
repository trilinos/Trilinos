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
#ifndef PYTRILINOS_LINEARPROBLEM_HPP
#define PYTRILINOS_LINEARPROBLEM_HPP

#include "Epetra_LinearProblem.h"
#include "Teuchos_RCP.hpp"

namespace PyTrilinos
{
class LinearProblem : public Epetra_LinearProblem
{
private:
  Teuchos::RCP< Epetra_RowMatrix   > _matrix;
  Teuchos::RCP< Epetra_Operator    > _operator;
  Teuchos::RCP< Epetra_MultiVector > _x;
  Teuchos::RCP< Epetra_MultiVector > _b;

public:
  LinearProblem();

  LinearProblem(const Teuchos::RCP< Epetra_RowMatrix > matrix,
                const Teuchos::RCP< Epetra_MultiVector > x,
                const Teuchos::RCP< Epetra_MultiVector > b);

  LinearProblem(const Teuchos::RCP< Epetra_Operator > op,
                const Teuchos::RCP< Epetra_MultiVector > x,
                const Teuchos::RCP< Epetra_MultiVector > b);

  LinearProblem(const LinearProblem & source);

  LinearProblem(const Epetra_LinearProblem & source);

  virtual ~LinearProblem();

  using Epetra_LinearProblem::CheckInput;

  using Epetra_LinearProblem::AssertSymmetric;

  using Epetra_LinearProblem::SetPDL;

  void SetOperator(Teuchos::RCP< Epetra_RowMatrix > & matrix);

  void SetOperator(Teuchos::RCP< Epetra_Operator > & op);

  void SetLHS(Teuchos::RCP< Epetra_MultiVector > & x);

  void SetRHS(Teuchos::RCP< Epetra_MultiVector > & b);

  using Epetra_LinearProblem::LeftScale;

  using Epetra_LinearProblem::RightScale;

  Teuchos::RCP< Epetra_RowMatrix > GetMatrix() const;

  Teuchos::RCP< Epetra_Operator > GetOperator() const;

  Teuchos::RCP< Epetra_MultiVector > GetLHS() const;

  Teuchos::RCP< Epetra_MultiVector > GetRHS() const;

  using Epetra_LinearProblem::GetPDL;

  using Epetra_LinearProblem::IsOperatorSymmetric;

};
}

#endif

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

