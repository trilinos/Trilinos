
// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  solver_def.hpp
    \brief Definition file for linear solvers for PDE-OPT.
*/

#ifndef ROL_PDEOPT_SOLVER_DEF_H
#define ROL_PDEOPT_SOLVER_DEF_H

#include "solver.hpp"

// Amesos2
#include "Amesos2.hpp"

// Ifpack2
#include "Ifpack2_Factory.hpp"

// MueLu
#include "MueLu.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"


template<class Real>
void Solver<Real>::setA(ROL::Ptr<Tpetra::CrsMatrix<> > &A) {
    if (useDirectSolver_) { // using Amesos2 direct solver
      if (firstSolve_) { // construct solver object
        try {
          solver_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >(directSolver_, A);
        }
        catch (std::invalid_argument e) {
          std::cout << e.what() << std::endl;
        }
        solver_->symbolicFactorization();
        firstSolve_ = false;
      }
      solver_->numericFactorization();
    }
    else { // construct MueLu or Ifpack2 preconditioner and Belos solver
      // Create transpose.
      Tpetra::RowMatrixTransposer<> transposer(A);
      A_trans_ = transposer.createTranspose();
      if (preconditioner_ == "Ifpack2") {
        Teuchos::ParameterList & parlistIfpack2 = parlist_.sublist("Ifpack2");
        // Create preconditioners.
        ifpack2Preconditioner_trans_ = Ifpack2::Factory::create<Tpetra::RowMatrix<Real,LO,GO,NO> > ("SCHWARZ", A_trans_);
        ifpack2Preconditioner_ = Ifpack2::Factory::create<Tpetra::RowMatrix<Real,LO,GO,NO> > ("SCHWARZ", A);
        ifpack2Preconditioner_trans_->setParameters(parlistIfpack2);
        ifpack2Preconditioner_->setParameters(parlistIfpack2);
        ifpack2Preconditioner_trans_->initialize();
        ifpack2Preconditioner_->initialize();
        ifpack2Preconditioner_trans_->compute();
        ifpack2Preconditioner_->compute();
      }
      else if (preconditioner_ == "MueLu") {
        Teuchos::ParameterList & parlistMuelu = parlist_.sublist("MueLu");
        // Create preconditioners.
        mueLuPreconditioner_trans_ = MueLu::CreateTpetraPreconditioner<Real,LO,GO,NO>(ROL::Ptr<OP>(A_trans_), parlistMuelu);
        mueLuPreconditioner_ = MueLu::CreateTpetraPreconditioner<Real,LO,GO,NO>(ROL::Ptr<OP>(A), parlistMuelu);
      }

      // Create Belos solver object and linear problem.
      if (firstSolve_) {
        Teuchos::RCP<Teuchos::ParameterList> parlistBelos = ROL::makePtrFromRef(parlist_.sublist("Belos"));
        // Transpose solver.
        problemBelos_trans_ = ROL::makePtr<Belos::LinearProblem<Real,MV,OP>>();
        problemBelos_trans_->setOperator(A_trans_);
        solverBelos_trans_ = ROL::makePtr<Belos::BlockGmresSolMgr<Real,MV,OP>>(problemBelos_trans_, parlistBelos);
        // Forward solver.
        problemBelos_ = ROL::makePtr<Belos::LinearProblem<Real,MV,OP>>();
        problemBelos_->setOperator(A);
        solverBelos_ = ROL::makePtr<Belos::BlockGmresSolMgr<Real,MV,OP>>(problemBelos_, parlistBelos);
        firstSolve_ = false;
      }
      else {
        problemBelos_trans_->setOperator(A_trans_);
      }
    }
}


template<class Real>
void Solver<Real>::solve(const ROL::Ptr<Tpetra::MultiVector<> > &x,
                         const ROL::Ptr<const Tpetra::MultiVector<> > &b,
                         const bool transpose) {

    if (useDirectSolver_) { // using Amesos2 direct solver
      //parlistAmesos2_->set("Transpose", transpose);
      parlistAmesos2_->sublist(directSolver_).set("Transpose", transpose);
      solver_->setParameters(parlistAmesos2_);
      solver_->setX(x);
      solver_->setB(b);
      solver_->solve();
    }
    else { // use Belos with MueLu or Ifpack2
      if ( transpose ) {
        problemBelos_trans_->setLHS(x);
        problemBelos_trans_->setRHS(b);
        if (preconditioner_ == "Ifpack2") {
          problemBelos_trans_->setRightPrec(ifpack2Preconditioner_trans_);
        }
        else if (preconditioner_ == "MueLu") {
          problemBelos_trans_->setRightPrec(mueLuPreconditioner_trans_);
        }
        problemBelos_trans_->setProblem();
        solverBelos_trans_->solve();
        }
      else {
        problemBelos_->setLHS(x);
        problemBelos_->setRHS(b);
        if (preconditioner_ == "Ifpack2") {
          problemBelos_->setRightPrec(ifpack2Preconditioner_);
        }
        else if (preconditioner_ == "MueLu") {
          problemBelos_->setRightPrec(mueLuPreconditioner_);
        }
        problemBelos_->setProblem();
        solverBelos_->solve();
      }
    }
  }

#endif
