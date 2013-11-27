//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:  Denis Ridzal (dridzal@sandia.gov)
//                        Drew Kouri (dpkouri@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_NONLINEARCG_H
#define ROL_NONLINEARCG_H

/** \class ROL::NonlinearCG
    \brief Provides implementation of nonlinear conjugate gradient methods.
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

    s.set(g);

    if ((state_->iter % state_->restart) != 0) {
      Real beta = 0.0;
      switch(state_->nlcg_type) {

        case NONLINEARCG_HESTENES_STIEFEL: {
          Teuchos::RCP<Vector<Real> > y = g.clone();
          y->set(g);
          y->axpy(-1.0, *(state_->grad[0]));
          beta =  - g.dot(*y) / y->dot(*(state_->pstep[0]));
          beta = std::max(beta, 0.0);
          break;
          }

        case NONLINEARCG_FLETCHER_REEVES: {
          beta = g.dot(g) / (state_->grad[0])->dot(*(state_->grad[0]));
          break;
          }

        case NONLINEARCG_DANIEL: {
          Teuchos::RCP<Vector<Real> > y = s.clone();
          Real htol = 0.0;
          obj.hessVec( *y, *(state_->pstep[0]), x, htol );
          beta = - g.dot(*y) / y->dot(*(state_->pstep[0]));
          beta = std::max(beta, 0.0);
          break;
          }

        case NONLINEARCG_POLAK_RIBIERE: {
          Teuchos::RCP<Vector<Real> > y = g.clone();
          y->set(g);
          y->axpy(-1.0, *(state_->grad[0]));
          beta = g.dot(*y) / (state_->grad[0])->dot(*(state_->grad[0]));
          beta = std::max(beta, 0.0);
          break;
          }

        case NONLINEARCG_FLETCHER_CONJDESC: {
          beta =  g.dot(g) / (state_->grad[0])->dot(*(state_->pstep[0]));
          break;
          }

        case NONLINEARCG_LIU_STOREY: {
          Teuchos::RCP<Vector<Real> > y = g.clone();
          y->set(g);
          y->axpy(-1.0, *(state_->grad[0]));
          beta =  g.dot(*y) / (state_->grad[0])->dot(*(state_->pstep[0]));
          //beta = std::max(beta, 0.0); // Is this needed?  May need research.
          break;
          }

        case NONLINEARCG_DAI_YUAN: {
          Teuchos::RCP<Vector<Real> > y = g.clone();
          y->set(g);
          y->axpy(-1.0, *(state_->grad[0]));
          beta =  - g.dot(g) / y->dot(*(state_->pstep[0]));
          break;
          }

        case NONLINEARCG_HAGAR_ZHANG: {
          Real eta_0 = 1e-2; 
          Teuchos::RCP<Vector<Real> > y = g.clone();
          Teuchos::RCP<Vector<Real> > yd = g.clone();
          y->set(g);
          y->axpy(-1.0, *(state_->grad[0]));
          yd->set(*y);
          Real mult = 2.0 * ( y->dot(*y) / y->dot(*(state_->pstep[0])) );
          yd->axpy(-mult, *(state_->pstep[0]));
          beta = - yd->dot(g) / y->dot(*(state_->pstep[0]));
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
