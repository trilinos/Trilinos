
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

/*! \file  solver.hpp
    \brief Linear solvers for PDE-OPT.
*/

#ifndef ROL_PDEOPT_SOLVER_H
#define ROL_PDEOPT_SOLVER_H

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Amesos2
#include "Amesos2.hpp"

// MueLu
#include "MueLu.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosTpetraAdapter.hpp"

template<class Real>
class Solver {

  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NO;
  typedef Tpetra::MultiVector<Real,LO,GO,NO> MV;
  typedef Tpetra::Operator<Real,LO,GO,NO> OP;

private:

  // Linear solvers and preconditioners for Jacobian and adjoint Jacobian
  Teuchos::RCP<Amesos2::Solver< Tpetra::CrsMatrix<>, Tpetra::MultiVector<> > > solver_;
  Teuchos::RCP<MueLu::TpetraOperator<Real,LO,GO,NO> > mueLuPreconditioner_;
  Teuchos::RCP<MueLu::TpetraOperator<Real,LO,GO,NO> > mueLuPreconditioner_trans_;
  Teuchos::RCP<Belos::BlockGmresSolMgr<Real,MV,OP> > solverBelos_;
  Teuchos::RCP<Belos::BlockGmresSolMgr<Real,MV,OP> > solverBelos_trans_;
  Teuchos::RCP<Belos::LinearProblem<Real,MV,OP> > problemBelos_;
  Teuchos::RCP<Belos::LinearProblem<Real,MV,OP> > problemBelos_trans_;

  // Linear solver options.
  bool useDirectSolver_;
  std::string solverType_;

  // Matrix transpose.
  Teuchos::RCP<Tpetra::CrsMatrix<> > A_trans_;

  // Parameter list.
  Teuchos::ParameterList parlist_;
  Teuchos::RCP<Teuchos::ParameterList> parlistAmesos2_;

  // Construct solvers on first solve.
  bool firstSolve_;

public:

  virtual ~Solver() {}

  Solver(Teuchos::ParameterList & parlist) : parlist_(parlist), firstSolve_(true) {
    useDirectSolver_ = parlist.get("Use Direct Solver", true);
    solverType_ = parlist.sublist("Direct").get("Solver Type", "KLU2");
    parlistAmesos2_ = Teuchos::rcp(new Teuchos::ParameterList("Amesos2"));
  }

  void setA(Teuchos::RCP<Tpetra::CrsMatrix<> > &A) {
    if (useDirectSolver_) { // using Amesos2 direct solver
      if (firstSolve_) { // construct solver object
        try {
          solver_ = Amesos2::create< Tpetra::CrsMatrix<>,Tpetra::MultiVector<> >(solverType_, A);
        }
        catch (std::invalid_argument e) {
          std::cout << e.what() << std::endl;
        }
        solver_->symbolicFactorization();
        firstSolve_ = false;
      }
      solver_->numericFactorization();
    } // useDirectSolver_
    else { // construct MueLu preconditioner and Belos solver
      Teuchos::ParameterList & parlistMuelu = parlist_.sublist("MueLu");
      // Create transpose.
      Tpetra::RowMatrixTransposer<> transposer(A);
      A_trans_ = transposer.createTranspose();
      // Create preconditioners.
      mueLuPreconditioner_trans_ = MueLu::CreateTpetraPreconditioner<Real,LO,GO,NO>(Teuchos::RCP<OP>(A_trans_), parlistMuelu);
      mueLuPreconditioner_ = MueLu::CreateTpetraPreconditioner<Real,LO,GO,NO>(Teuchos::RCP<OP>(A), parlistMuelu);
      // Create Belos solver object and linear problem.
      if (firstSolve_) {
        Teuchos::RCP<Teuchos::ParameterList> parlistBelos = Teuchos::rcpFromRef(parlist_.sublist("Belos"));
        // Transpose solver.
        problemBelos_trans_ = Teuchos::rcp(new Belos::LinearProblem<Real,MV,OP>());
        problemBelos_trans_->setOperator(A_trans_);
        solverBelos_trans_ = Teuchos::rcp(new Belos::BlockGmresSolMgr<Real,MV,OP>(problemBelos_trans_, parlistBelos));
        // Forward solver.
        problemBelos_ = Teuchos::rcp(new Belos::LinearProblem<Real,MV,OP>());
        problemBelos_->setOperator(A);
        solverBelos_ = Teuchos::rcp(new Belos::BlockGmresSolMgr<Real,MV,OP>(problemBelos_, parlistBelos));
        firstSolve_ = false;
      }
    }
  }

  void solve(const Teuchos::RCP<Tpetra::MultiVector<> > &x,
             const Teuchos::RCP<const Tpetra::MultiVector<> > &b,
             const bool transpose = false) {
    if (useDirectSolver_) { // using Amesos2 direct solver
      //parlistAmesos2_->set("Transpose", transpose);
      parlistAmesos2_->sublist(solverType_).set("Transpose", transpose);
      solver_->setParameters(parlistAmesos2_);
      solver_->setX(x);
      solver_->setB(b);
      solver_->solve();
    }  // useDirectSolver_
    else { // use Belos
      if ( transpose ) {
        problemBelos_trans_->setLHS(x);
        problemBelos_trans_->setRHS(b);
        problemBelos_trans_->setRightPrec(mueLuPreconditioner_trans_);
        problemBelos_trans_->setProblem();
        solverBelos_trans_->solve();
      }
      else {
        problemBelos_->setLHS(x);
        problemBelos_->setRHS(b);
        problemBelos_->setRightPrec(mueLuPreconditioner_);
        problemBelos_->setProblem();
        solverBelos_->solve();
      }
    }
  }

};

#endif
