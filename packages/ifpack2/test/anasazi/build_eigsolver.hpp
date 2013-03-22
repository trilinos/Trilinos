/*
// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

#ifndef _build_eigsolver_hpp_
#define _build_eigsolver_hpp_

#include <Teuchos_RefCountPtr.hpp>
#include <AnasaziEigenproblem.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>

template<class Scalar,class MV, class OP>
Teuchos::RCP<Anasazi::SolverManager<Scalar,MV,OP> >
build_eigsolver(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
             Teuchos::ParameterList& test_params,
             Teuchos::RCP<Anasazi::Eigenproblem<Scalar,MV,OP> > problem)
{
  typedef Anasazi::Eigenproblem<Scalar,MV,OP> AEigProb;
  Teuchos::RCP<Anasazi::SolverManager<Scalar,MV,OP> > solver;

  Teuchos::ParameterList aparams;
  if (test_params.isSublist("Anasazi")) {
    aparams = test_params.sublist("Anasazi");
  }

  std::string solver_type("not specified");
  Ifpack2::getParameter(test_params, "eigen_solver_type", solver_type);
  if (solver_type == "BlockKrylovSchur") {
    // if (comm->getRank() == 0) std::cout << aparams << std::endl;
    solver = Teuchos::rcp(new Anasazi::BlockKrylovSchurSolMgr<Scalar,MV,OP>(problem,aparams));
  }
  else if (solver_type == "not specified") {
    throw std::runtime_error("Error in build_eigsolver: solver_type not specified.");
  }
  else {
    std::ostringstream os;
    os << "Error in build_eigsolver: solver_type ("<<solver_type<<") not recognized.";
    os << "\nIfpack2's test-driver recognizes these solvers: PseudoBlockCG, PesudoBlockGmres, BlockGmres, TFQMR.";
    std::string str = os.str();
    throw std::runtime_error(str);
  }
  return solver;
}
#endif

