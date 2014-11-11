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

#ifndef ROL_NONLINEARCG_H
#define ROL_NONLINEARCG_H

/** \class ROL::NonlinearCG
    \brief Implementats nonlinear conjugate gradient methods.
*/

#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
struct NonlinearCGState {
  std::vector<Teuchos::RCP<Vector<Real> > > grad;   // Gradient Storage
  std::vector<Teuchos::RCP<Vector<Real> > > pstep;  // Step Storage
  int iter;                                         // Nonlinear-CG Iteration Counter
  int restart;                                      // Reinitialize every 'restart' iterations
  ENonlinearCG nlcg_type;                           // Nonlinear-CG Type
};

template<class Real>
class NonlinearCG {
private:

  Teuchos::RCP<NonlinearCGState<Real> > state_; // Nonlinear-CG State

  Teuchos::RCP<Vector<Real> > y_;
  Teuchos::RCP<Vector<Real> > yd_;

public:

  virtual ~NonlinearCG() {}

  // Constructor
  NonlinearCG(ENonlinearCG type, int restart = 100) {
    state_ = Teuchos::rcp( new NonlinearCGState<Real> ); 
    state_->iter = 0;
    state_->grad.resize(1);
    state_->pstep.resize(1);
    TEUCHOS_TEST_FOR_EXCEPTION(!(isValidNonlinearCG(type)),
                          std::invalid_argument,
                          ">>> ERROR (ROL_NonlinearCG.hpp): Invalid nonlinear CG type in constructor!");
    state_->nlcg_type = type;
    TEUCHOS_TEST_FOR_EXCEPTION((restart < 1),
                          std::invalid_argument,
                          ">>> ERROR (ROL_NonlinearCG.hpp): Non-positive restart integer in constructor!");
    state_->restart = restart;
  }

  Teuchos::RCP<NonlinearCGState<Real> >& get_state() { return this->state_; }

  // Run one step of nonlinear CG.
  virtual void run( Vector<Real> &s , const Vector<Real> &g, const Vector<Real> &x, Objective<Real> &obj ) {
    // Initialize vector storage
    if ( state_->iter == 0 ) {
      if ( state_->nlcg_type != NONLINEARCG_FLETCHER_REEVES && 
           state_->nlcg_type != NONLINEARCG_FLETCHER_CONJDESC ) {
        y_ = g.clone();
      }
      if ( state_->nlcg_type == NONLINEARCG_HAGAR_ZHANG ) {
        yd_ = g.clone();
      }
    }

    s.set(g.dual());

    if ((state_->iter % state_->restart) != 0) {
      Real beta = 0.0;
      switch(state_->nlcg_type) {

        case NONLINEARCG_HESTENES_STIEFEL: {
          y_->set(g);
          y_->axpy(-1.0, *(state_->grad[0]));
          beta =  - g.dot(*y_) / (state_->pstep[0]->dot(y_->dual()));
          beta = std::max(beta, 0.0);
          break;
          }

        case NONLINEARCG_FLETCHER_REEVES: {
          beta = g.dot(g) / (state_->grad[0])->dot(*(state_->grad[0]));
          break;
          }

        case NONLINEARCG_DANIEL: {
          Real htol = 0.0;
          obj.hessVec( *y_, *(state_->pstep[0]), x, htol );
          beta = - g.dot(*y_) / (state_->pstep[0])->dot(y_->dual());
          beta = std::max(beta, 0.0);
          break;
          }

        case NONLINEARCG_POLAK_RIBIERE: {
          y_->set(g);
          y_->axpy(-1.0, *(state_->grad[0]));
          beta = g.dot(*y_) / (state_->grad[0])->dot(*(state_->grad[0]));
          beta = std::max(beta, 0.0);
          break;
          }

        case NONLINEARCG_FLETCHER_CONJDESC: {
          beta =  g.dot(g) / (state_->pstep[0])->dot((state_->grad[0])->dual());
          break;
          }

        case NONLINEARCG_LIU_STOREY: {
          y_->set(g);
          y_->axpy(-1.0, *(state_->grad[0]));
          beta =  g.dot(*y_) / (state_->pstep[0])->dot((state_->grad[0])->dual());
          //beta = std::max(beta, 0.0); // Is this needed?  May need research.
          break;
          }

        case NONLINEARCG_DAI_YUAN: {
          y_->set(g);
          y_->axpy(-1.0, *(state_->grad[0]));
          beta =  - g.dot(g) / (state_->pstep[0])->dot(y_->dual());
          break;
          }

        case NONLINEARCG_HAGAR_ZHANG: {
          Real eta_0 = 1e-2; 
          y_->set(g);
          y_->axpy(-1.0, *(state_->grad[0]));
          yd_->set(*y_);
          Real mult = 2.0 * ( y_->dot(*y_) / (state_->pstep[0])->dot(y_->dual()) );
          yd_->axpy(-mult, (state_->pstep[0])->dual());
          beta = - yd_->dot(g) / (state_->pstep[0])->dot(y_->dual());
          Real eta = -1.0 / ((state_->pstep[0])->norm()*std::min(eta_0,(state_->grad[0])->norm()));
          beta = std::max(beta, eta);
          break;
          }

        default:
          TEUCHOS_TEST_FOR_EXCEPTION(!(isValidNonlinearCG(state_->nlcg_type)),
                          std::invalid_argument,
                          ">>> ERROR (ROL_NonlinearCG.hpp): Invalid nonlinear CG type in the 'run' method!");  
      }

      s.axpy(beta, *(state_->pstep[0]));
    }

    // Update storage.
    if (state_->iter == 0) {
      (state_->grad[0]) = g.clone();
      (state_->pstep[0]) = s.clone();
    }
    (state_->grad[0])->set(g);
    (state_->pstep[0])->set(s);
    state_->iter++;
  }

};

}

#endif
