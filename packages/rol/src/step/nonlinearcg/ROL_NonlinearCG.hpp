// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NONLINEARCG_H
#define ROL_NONLINEARCG_H

/** \class ROL::NonlinearCG
    \brief Implements nonlinear conjugate gradient methods.
    \detail Nonlinear CG methods have the update formulas
            \f[ x_{k+1} = x_k + \alpha_k d_k \f]
            \f[ d_{k+1} = -g_{k+1} + \beta_k d_k \f]
            where \f$\alpha_k\f$ is a step length determined by a linesearch method
            and \f$\beta_k\f$ is the parameter which distinguishes various CG methods.

            The standard notation is that \f$ y_k = g_{k+1}-g_k \$f, 
            \f$ s_k = \alpha_k d_k = x_{k+1}-x_k \f$

    Method                     | \f$\beta_k\f$                                  
    ---------------------------|------------------------------------------------------
    Hestenes-Stiefel           | \f$ \frac{g_{k+1}^\top y_k}{d_k^\top y_k } \f$
    Fletcher-Reeves            | \f$ \frac{\|g_{k+1}\|^2}{\|g_k\|^2} \f$
    Daniel (uses Hessian)      | \f$ \frac{g_{k+1}^\top \nabla^2 f(x_k) d_k}{d_k^\top \nabla^2 f(x_k) d_k} \f$
    Polak-Ribiere              | \f$ \frac{g_{k+1}^\top y_k}{\|g_k\|^2} \f$
    Fletcher Conjugate Descent | \f$ -\frac{\|g_{k+1}\|^2}{d_k^\top g_k} \f$
    Liu-Storey                 | \f$ -\frac{g_k^\top y_{k-1} }{d_{k-1}^\top g_{k-1} \f$
    Dai-Yuan                   | \f$ \frac{\|g_{k+1}\|^2}{d_k^\top y_k} \f$
    Hager-Zhang                | \f$ \frac{g_{k+1}^\top y_k}{d_k^\top y_k} - 2 \frac{\|y_k\|^2}{d_k^\top y_k} \frac{g_{k+1}^\top d_k}{d_k^\top y_k} \f$
    Oren-Luenberger            | \f$ \frac{g_{k+1}^\top y_k}{d_k^\top y_k} - \frac{\|y_k\|^2}{d_k^\top y_k} \frac{g_{k+1}^\top d_k}{d_k^\top y_k} \f$ 
            
*/

#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
struct NonlinearCGState {
  std::vector<ROL::Ptr<Vector<Real> > > grad;   // Gradient Storage
  std::vector<ROL::Ptr<Vector<Real> > > pstep;  // Step Storage
  int iter;                                         // Nonlinear-CG Iteration Counter
  int restart;                                      // Reinitialize every 'restart' iterations
  ENonlinearCG nlcg_type;                           // Nonlinear-CG Type
};

template<class Real>
class NonlinearCG {
private:

  ROL::Ptr<NonlinearCGState<Real> > state_; // Nonlinear-CG State

  ROL::Ptr<Vector<Real> > y_;
  ROL::Ptr<Vector<Real> > yd_;

public:

  virtual ~NonlinearCG() {}

  // Constructor
  NonlinearCG(ENonlinearCG type, int restart = 100) {
    state_ = ROL::makePtr<NonlinearCGState<Real>>(); 
    state_->iter = 0;
    state_->grad.resize(1);
    state_->pstep.resize(1);
    ROL_TEST_FOR_EXCEPTION(!(isValidNonlinearCG(type)),
                          std::invalid_argument,
                          ">>> ERROR (ROL_NonlinearCG.hpp): Invalid nonlinear CG type in constructor!");
    state_->nlcg_type = type;
    ROL_TEST_FOR_EXCEPTION((restart < 1),
                          std::invalid_argument,
                          ">>> ERROR (ROL_NonlinearCG.hpp): Non-positive restart integer in constructor!");
    state_->restart = restart;
  }

  ROL::Ptr<NonlinearCGState<Real> >& get_state() { return this->state_; }

  // Run one step of nonlinear CG.
  virtual void run( Vector<Real> &s , const Vector<Real> &g, const Vector<Real> &x, Objective<Real> &obj ) {
    Real one(1);
    // Initialize vector storage
    if ( state_->iter == 0 ) {
      if ( state_->nlcg_type != NONLINEARCG_FLETCHER_REEVES && 
           state_->nlcg_type != NONLINEARCG_FLETCHER_CONJDESC ) {
        y_ = g.clone();
      }
      if ( state_->nlcg_type == NONLINEARCG_HAGER_ZHANG ||
           state_->nlcg_type == NONLINEARCG_OREN_LUENBERGER ) {
        yd_ = g.clone();
      }
    }

    s.set(g.dual());

    if ((state_->iter % state_->restart) != 0) {
      Real beta(0), zero(0);
      switch(state_->nlcg_type) {

        case NONLINEARCG_HESTENES_STIEFEL: {
          y_->set(g);
          y_->axpy(-one, *(state_->grad[0]));
          beta =  - g.dot(*y_) / (state_->pstep[0]->dot(y_->dual()));
          beta = std::max(beta, zero);
          break;
          }

        case NONLINEARCG_FLETCHER_REEVES: {
          beta = g.dot(g) / (state_->grad[0])->dot(*(state_->grad[0]));
          break;
          }

        case NONLINEARCG_DANIEL: {
          Real htol(0);
          obj.hessVec( *y_, *(state_->pstep[0]), x, htol );
          beta = - g.dot(*y_) / (state_->pstep[0])->dot(y_->dual());
          beta = std::max(beta, zero);
          break;
          }

        case NONLINEARCG_POLAK_RIBIERE: {
          y_->set(g);
          y_->axpy(-one, *(state_->grad[0]));
          beta = g.dot(*y_) / (state_->grad[0])->dot(*(state_->grad[0]));
          beta = std::max(beta, zero);
          break;
          }

        case NONLINEARCG_FLETCHER_CONJDESC: {
          beta =  g.dot(g) / (state_->pstep[0])->dot((state_->grad[0])->dual());
          break;
          }

        case NONLINEARCG_LIU_STOREY: {
          y_->set(g);
          y_->axpy(-one, *(state_->grad[0]));
          beta =  g.dot(*y_) / (state_->pstep[0])->dot((state_->grad[0])->dual());
          //beta = std::max(beta, 0.0); // Is this needed?  May need research.
          break;
          }

        case NONLINEARCG_DAI_YUAN: {
          y_->set(g);
          y_->axpy(-one, *(state_->grad[0]));
          beta =  - g.dot(g) / (state_->pstep[0])->dot(y_->dual());
          break;
          }

        case NONLINEARCG_HAGER_ZHANG: {
          Real eta_0(1e-2), two(2); 
          y_->set(g);
          y_->axpy(-one, *(state_->grad[0]));
          yd_->set(*y_);
          Real mult = two * ( y_->dot(*y_) / (state_->pstep[0])->dot(y_->dual()) );
          yd_->axpy(-mult, (state_->pstep[0])->dual());
          beta = - yd_->dot(g) / (state_->pstep[0])->dot(y_->dual());
          Real eta = -one / ((state_->pstep[0])->norm()*std::min(eta_0,(state_->grad[0])->norm()));
          beta = std::max(beta, eta);
          break;
          }

        case NONLINEARCG_OREN_LUENBERGER: {
          Real eta_0(1e-2);
          y_->set(g);
          y_->axpy(-one, *(state_->grad[0]));
          yd_->set(*y_);
          Real mult = ( y_->dot(*y_) / (state_->pstep[0])->dot(y_->dual()) );
          yd_->axpy(-mult, (state_->pstep[0])->dual());
          beta = - yd_->dot(g) / (state_->pstep[0])->dot(y_->dual());
          Real eta = -one / ((state_->pstep[0])->norm()*std::min(eta_0,(state_->grad[0])->norm()));
          beta = std::max(beta, eta);
          break;
          }

        default:
          ROL_TEST_FOR_EXCEPTION(!(isValidNonlinearCG(state_->nlcg_type)),
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
