/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef _build_solver_hpp_
#define _build_solver_hpp_

#include "Teuchos_RefCountPtr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"

template<class Scalar,class MV, class OP>
Teuchos::RCP<Belos::SolverManager<Scalar,MV,OP> >
build_solver(Teuchos::ParameterList& test_params,
             Teuchos::RCP<Belos::LinearProblem<Scalar,MV,OP> > problem)
{
  Teuchos::ParameterList bparams;
  if (test_params.isSublist("Belos")) {
    bparams = test_params.sublist("Belos");
  }
  Teuchos::RCP<Teuchos::ParameterList> rcpparams = Teuchos::rcp(&bparams,false);

  std::string solver_type("not specified");

  Ifpack2::getParameter(test_params, "solver_type", solver_type);

  if (solver_type == "not specified") {
    throw std::runtime_error("Error in build_solver: solver_type not specified.");
  }

  Teuchos::RCP<Belos::SolverManager<Scalar,MV,OP> > solver;
  Belos::SolverFactory<Scalar, MV, OP> factory;
  try {
    solver = factory.create (solver_type, rcpparams);
  } catch (std::exception& e) {
    std::cout << "*** FAILED: Belos::SolverFactory::create threw an exception: "
              << e.what () << std::endl;
    return Teuchos::null;
  }

  solver->setProblem( problem );

  return solver;
}
#endif

