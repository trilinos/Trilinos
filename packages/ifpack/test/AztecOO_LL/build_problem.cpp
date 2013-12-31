/*
//@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/


#include <string>
#include <sstream>
#include <iostream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

#include "ParameterHelper.hpp"
#include "read_matrix.hpp"

Teuchos::RCP<Epetra_LinearProblem> build_problem_mm(Teuchos::ParameterList& test_params, Epetra_CrsMatrix* A, Epetra_MultiVector* b)
{
  const Epetra_Map& rowmap = A->RowMap();

  Epetra_MultiVector* x = new Epetra_MultiVector(rowmap, 1);
  if (b == NULL) {
std::cout << "creating b = A*random" << std::endl;
    b = new Epetra_MultiVector(rowmap, 1);
    x->Random();

    A->Apply(*x, *b);
  }
  x->PutScalar(0);

  Teuchos::RCP<Epetra_LinearProblem> problem = Teuchos::rcp(new Epetra_LinearProblem(A,x,b));

  return problem;
}

Teuchos::RCP< Epetra_LinearProblem >
build_problem(Teuchos::ParameterList& test_params,
              const Epetra_Comm& comm)
{
  Teuchos::Time timer("build_problem");
  timer.start();

  Epetra_CrsMatrix* A;
  Epetra_Vector* b = NULL;

  std::string mm_file("not specified");
  std::string rhs_mm_file("not specified");
  helper::GetParameter(test_params, "mm_file", mm_file);
  helper::GetParameter(test_params, "rhs_mm_file", rhs_mm_file);
  std::string hb_file("not specified");
  helper::GetParameter(test_params, "hb_file", hb_file);

  if (mm_file != "not specified") {
    if (comm.MyPID() == 0) {
      std::cout << "Matrix-Market file: " << mm_file << std::endl;
    }
    A = read_matrix_mm(mm_file, comm);
    if (rhs_mm_file != "not specified") {
      if (comm.MyPID() == 0) {
        std::cout << "Matrix-Market file: " << rhs_mm_file << std::endl;
      }
      b = read_vector_mm(rhs_mm_file, comm);
    }
  }
  else if (hb_file != "not specified") {
    read_matrix_hb(hb_file, comm, A, b);
  }
  else {
    throw std::runtime_error("No matrix file specified.");
  }

  Teuchos::RCP<Epetra_LinearProblem> problem = build_problem_mm(test_params, A, b);
  timer.stop();
  if (comm.MyPID() == 0) {
    std::cout << "proc 0 time to read matrix & create problem: " << timer.totalElapsedTime()
      << std::endl;
  }

  return problem;
}

