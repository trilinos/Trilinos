// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for helper functions in ROL.
        \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_HELPERFUNCTIONS_HPP
#define ROL_HELPERFUNCTIONS_HPP

#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Secant.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

namespace ROL {

  template<class Real>
  Teuchos::SerialDenseMatrix<int, Real> computeDenseHessian(Objective<Real> &obj, const Vector<Real> &x) {

    Real tol = std::sqrt(ROL_EPSILON<Real>());

    int dim = x.dimension();
    Teuchos::SerialDenseMatrix<int, Real> H(dim, dim);

    ROL::Ptr<Vector<Real> > e = x.clone();
    ROL::Ptr<Vector<Real> > h = x.dual().clone();

    for (int i=0; i<dim; i++) {
      e = x.basis(i);
      obj.hessVec(*h, *e, x, tol);
      for (int j=0; j<dim; j++) {
        e = x.basis(j);
        //H(j,i) = e->dot(h->dual());
        H(j,i) = e->apply(*h);
      }
    }

    return H;

  }


  template<class Real>
  Teuchos::SerialDenseMatrix<int, Real> computeScaledDenseHessian(Objective<Real> &obj, const Vector<Real> &x) {

    Real tol = std::sqrt(ROL_EPSILON<Real>());

    int dim = x.dimension();
    Teuchos::SerialDenseMatrix<int, Real> H(dim, dim);

    ROL::Ptr<Vector<Real> > ei = x.clone();
    ROL::Ptr<Vector<Real> > ej = x.dual().clone();
    ROL::Ptr<Vector<Real> > h  = x.dual().clone();

    for (int i=0; i<dim; i++) {
      ei = ei->basis(i);
      obj.hessVec(*h, *ei, x, tol);
      for (int j=0; j<dim; j++) {
        ej = ej->basis(j);
        H(j,i) = ej->dot(*h);
      }
    }

    return H;

  }


  template<class Real>
  Teuchos::SerialDenseMatrix<int, Real> computeDotMatrix(const Vector<Real> &x) {

    int dim = x.dimension();
    Teuchos::SerialDenseMatrix<int, Real> M(dim, dim);

    ROL::Ptr<Vector<Real> > ei = x.clone();
    ROL::Ptr<Vector<Real> > ej = x.clone();

    for (int i=0; i<dim; i++) {
      ei = x.basis(i);
      for (int j=0; j<dim; j++) {
        ej = x.basis(j);
        M(j,i) = ej->dot(*ei);
      }
    }

    return M;

  }

  template<class Real>
  std::vector<std::vector<Real> > computeEigenvalues(const Teuchos::SerialDenseMatrix<int, Real> & mat) {

    Teuchos::LAPACK<int, Real> lapack;

    Teuchos::SerialDenseMatrix<int, Real> mymat(Teuchos::Copy, mat);

    char jobvl = 'N';
    char jobvr = 'N';

    int n = mat.numRows();

    std::vector<Real> real(n, 0);
    std::vector<Real> imag(n, 0);
    std::vector<std::vector<Real> > eigenvals;

    Real* vl = 0;
    Real* vr = 0;

    int ldvl = 1;
    int ldvr = 1;

    int lwork = 4*n;

    std::vector<Real> work(lwork, 0);

    int info = 0;

    lapack.GEEV(jobvl, jobvr, n, &mymat(0,0), n, &real[0], &imag[0], vl, ldvl, vr, ldvr, &work[0], lwork, &info);

    eigenvals.push_back(real);
    eigenvals.push_back(imag);

    return eigenvals;

  }


  template<class Real>
  std::vector<std::vector<Real> > computeGenEigenvalues(const Teuchos::SerialDenseMatrix<int, Real> & A,
                                                        const Teuchos::SerialDenseMatrix<int, Real> & B) {

    Teuchos::LAPACK<int, Real> lapack;

    Teuchos::SerialDenseMatrix<int, Real> myA(Teuchos::Copy, A);
    Teuchos::SerialDenseMatrix<int, Real> myB(Teuchos::Copy, B);

    char jobvl = 'N';
    char jobvr = 'N';

    int n = A.numRows();

    std::vector<Real> real(n, 0);
    std::vector<Real> imag(n, 0);
    std::vector<Real> beta(n, 0);
    std::vector<std::vector<Real> > eigenvals;

    Real* vl = 0;
    Real* vr = 0;

    int ldvl = 1;
    int ldvr = 1;

    int lwork = 10*n;

    std::vector<Real> work(lwork, 0);

    int info = 0;

    lapack.GGEV(jobvl, jobvr, n, &myA(0,0), n, &myB(0,0), n, &real[0], &imag[0], &beta[0], 
                vl, ldvl, vr, ldvr, &work[0], lwork, &info);

    for (int i=0; i<n; i++) {
      real[i] /= beta[i];
      imag[i] /= beta[i];
    }

    eigenvals.push_back(real);
    eigenvals.push_back(imag);

    return eigenvals;

  }


  template<class Real>
  Teuchos::SerialDenseMatrix<int, Real> computeInverse(const Teuchos::SerialDenseMatrix<int, Real> & mat) {

    Teuchos::LAPACK<int, Real> lapack;

    Teuchos::SerialDenseMatrix<int, Real> mymat(Teuchos::Copy, mat);

    int n = mat.numRows();

    std::vector<int> ipiv(n, 0);

    int lwork = 5*n;

    std::vector<Real> work(lwork, 0);

    int info = 0;

    lapack.GETRF(n, n, &mymat(0,0), n, &ipiv[0], &info);
    lapack.GETRI(n, &mymat(0,0), n, &ipiv[0], &work[0], lwork, &info);

    return mymat;

  }



  template<class Real>
  class ProjectedObjective : public Objective<Real> {
  private:
    ROL::Ptr<Objective<Real> >       obj_;
    ROL::Ptr<BoundConstraint<Real> > con_;
    ROL::Ptr<Secant<Real> >          secant_;

    ROL::Ptr<ROL::Vector<Real> > primalV_;
    ROL::Ptr<ROL::Vector<Real> > dualV_;
    bool isInitialized_;

    bool useSecantPrecond_;
    bool useSecantHessVec_;
    Real eps_;

  public:
    ProjectedObjective( Objective<Real> &obj, BoundConstraint<Real> &con,
                        ROL::Ptr<Secant<Real> > &secant,
                        bool useSecantPrecond = false,
                        bool useSecantHessVec = false,
                        Real eps = 0.0 )
      : isInitialized_(false), useSecantPrecond_(useSecantPrecond),
        useSecantHessVec_(useSecantHessVec), eps_(eps) {
      obj_    = ROL::makePtrFromRef(obj);
      con_    = ROL::makePtrFromRef(con);
      secant_ = secant;
    }

    void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
      obj_->update(x,flag,iter);
      con_->update(x,flag,iter);
    }

    Real value( const Vector<Real> &x, Real &tol ) {
      return obj_->value(x,tol);
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      obj_->gradient(g,x,tol);
    }

    Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
      return obj_->dirDeriv(x,d,tol);
    }

    void hessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      if ( useSecantHessVec_ ) {
        secant_->applyB( Hv, v );
      }
      else {
        obj_->hessVec( Hv, v, x, tol );
      }
    }

    void invHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      if ( useSecantHessVec_ ) {
        secant_->applyH(Hv,v);
      }
      else {
        obj_->invHessVec(Hv,v,x,tol);
      }
    }

    void precond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      if ( useSecantPrecond_ ) {
        secant_->applyH( Mv, v );
      }
      else {
        obj_->precond( Mv, v, x, tol );
      }
    }

    /** \brief Apply the reduced Hessian to a vector, v.
               The reduced Hessian first removes elements of v
               corresponding to the feasible indices from
               the point p in the direction -d.
                   Hv   the Hessian times a vector
                   v    input vector
                   p    starting point for tangent cone
                   d    negative of search direction
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &p,
                         const Vector<Real> &d, const Vector<Real> &x, Real &tol ) {
      if ( con_->isActivated() ) {
        if (!isInitialized_) {
          primalV_ = x.clone();
          dualV_   = x.dual().clone();
          isInitialized_ = true;
        }
        // Set vnew to v
        primalV_->set(v);
        // Remove elements of vnew corresponding to binding set
        con_->pruneActive(*primalV_,d,p,eps_);
        // Apply full Hessian to reduced vector
        hessVec(Hv,*primalV_,x,tol);
        // Remove elements of Hv corresponding to binding set
        con_->pruneActive(Hv,d,p,eps_);
        // Set vnew to v
        primalV_->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        con_->pruneInactive(*primalV_,d,p,eps_);
        dualV_->set(primalV_->dual());
        con_->pruneInactive(*dualV_,d,p,eps_);
        // Fill complement of binding set elements in Hp with v
        Hv.plus(*dualV_);
      }
      else {
        hessVec(Hv,v,x,tol);
      }
    }
 
    /** \brief Apply the reduced Hessian to a vector, v.
               The reduced Hessian first removes elements of v
               corresponding to the feasible indices from
               the point p.
                   Hv   the Hessian times a vector
                   v    input vector
                   p    starting point for tangent cone
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &p,
                         const Vector<Real> &x, Real &tol ) {
      if ( con_->isActivated() ) {
        if (!isInitialized_) {
          primalV_ = x.clone();
          dualV_   = x.dual().clone();
          isInitialized_ = true;
        }
        // Set vnew to v
        primalV_->set(v);
        // Remove elements of vnew corresponding to binding set
        con_->pruneActive(*primalV_,p,eps_);
        // Apply full Hessian to reduced vector
        hessVec(Hv,*primalV_,x,tol);
        // Remove elements of Hv corresponding to binding set
        con_->pruneActive(Hv,p,eps_);
        // Set vnew to v
        primalV_->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        con_->pruneInactive(*primalV_,p,eps_);
        dualV_->set(primalV_->dual());
        con_->pruneInactive(*dualV_,p,eps_);
        // Fill complement of binding set elements in Hp with v
        Hv.plus(*dualV_);
      }
      else {
        hessVec(Hv,v,x,tol);
      }
    }

    /** \brief Apply the reduced inverse Hessian to a vector, v.  
               The reduced inverse Hessian first removes elements 
               of v corresponding to the feasible indices from 
               the point p in the direction -d.
                   Hv   the inverse Hessian times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   d    negative of search direction
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedInvHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &p, 
                            const Vector<Real> &d, const Vector<Real> &x, Real &tol ) {
      if ( con_->isActivated() ) {
        if (!isInitialized_) {
          primalV_ = x.clone();
          dualV_   = x.dual().clone();
          isInitialized_ = true;
        }
        // Set vnew to v
        dualV_->set(v);
        // Remove elements of vnew corresponding to binding set
        con_->pruneActive(*dualV_,d,p,eps_);
        // Apply full Hessian to reduced vector
        invHessVec(Hv,*dualV_,x,tol);
        // Remove elements of Hv corresponding to binding set
        con_->pruneActive(Hv,d,p,eps_);
        // Set vnew to v
        dualV_->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        con_->pruneInactive(*dualV_,d,p,eps_);
        primalV_->set(dualV_->dual());
        con_->pruneInactive(*primalV_,d,p,eps_);
        // Fill complement of binding set elements in Hv with v
        Hv.plus(*primalV_);
      }
      else {
        invHessVec(Hv,v,x,tol);
      }
    }

    /** \brief Apply the reduced inverse Hessian to a vector, v.  
               The reduced inverse Hessian first removes elements 
               of v corresponding to the feasible indices from 
               the point p.
                   Hv   the inverse Hessian times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedInvHessVec( Vector<Real> &Hv, const Vector<Real> &v, const Vector<Real> &p, 
                            const Vector<Real> &x, Real &tol ) {
      if ( con_->isActivated() ) {
        if (!isInitialized_) {
          primalV_ = x.clone();
          dualV_   = x.dual().clone();
          isInitialized_ = true;
        }
        // Set vnew to v
        dualV_->set(v);
        // Remove elements of vnew corresponding to binding set
        con_->pruneActive(*dualV_,p,eps_);
        // Apply full Hessian to reduced vector
        invHessVec(Hv,*dualV_,x,tol);
        // Remove elements of Hv corresponding to binding set
        con_->pruneActive(Hv,p,eps_);
        // Set vnew to v
        dualV_->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        con_->pruneInactive(*dualV_,p,eps_);
        primalV_->set(dualV_->dual());
        con_->pruneInactive(*primalV_,p,eps_);
        // Fill complement of binding set elements in Hv with v
        Hv.plus(*primalV_);
      }
      else {
        invHessVec(Hv,v,x,tol);
      }
    }

    /** \brief Apply the reduced preconditioner to a vector, v.  
               The reduced preconditioner first removes elements 
               of v corresponding to the feasible indices from 
               the point p in the direction -d.
                   Hv   the preconditioner times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   d    negative of search direction
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedPrecond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &p, 
                         const Vector<Real> &d, const Vector<Real> &x, Real &tol ) {
      if ( con_->isActivated() ) {
        if (!isInitialized_) {
          primalV_ = x.clone();
          dualV_   = x.dual().clone();
          isInitialized_ = true;
        }
        // Set vnew to v
        dualV_->set(v);
        // Remove elements of vnew corresponding to binding set
        con_->pruneActive(*dualV_,d,p,eps_);
        // Apply full Hessian to reduced vector
        precond(Mv,*dualV_,x,tol);
        // Remove elements of Mv corresponding to binding set
        con_->pruneActive(Mv,d,p,eps_);
        // Set vnew to v
        dualV_->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        con_->pruneInactive(*dualV_,d,p,eps_);
        primalV_->set(dualV_->dual());
        con_->pruneInactive(*primalV_,d,p,eps_);
        // Fill complement of binding set elements in Mv with v
        Mv.plus(*primalV_);
      }
      else {
        precond(Mv,v,x,tol);
      }
    }

    /** \brief Apply the reduced preconditioner to a vector, v.  
               The reduced preconditioner first removes elements 
               of v corresponding to the feasible indices from 
               the point p.
                   Hv   the preconditioner times a vector
                   v    input vector 
                   p    starting point for tangent cone
                   x    current iteration vector
                   tol  objective function tolerance
    */
    void reducedPrecond( Vector<Real> &Mv, const Vector<Real> &v, const Vector<Real> &p, 
                         const Vector<Real> &x, Real &tol ) {
      if ( con_->isActivated() ) {
        if (!isInitialized_) {
          primalV_ = x.clone();
          dualV_   = x.dual().clone();
          isInitialized_ = true;
        }
        // Set vnew to v
        dualV_->set(v);
        // Remove elements of vnew corresponding to binding set
        con_->pruneActive(*dualV_,p,eps_);
        // Apply full Hessian to reduced vector
        precond(Mv,*dualV_,x,tol);
        // Remove elements of Mv corresponding to binding set
        con_->pruneActive(Mv,p,eps_);
        // Set vnew to v
        dualV_->set(v);
        // Remove Elements of vnew corresponding to complement of binding set
        con_->pruneInactive(*dualV_,p,eps_);
        primalV_->set(dualV_->dual());
        con_->pruneInactive(*primalV_,p,eps_);
        // Fill complement of binding set elements in Mv with v
        Mv.plus(*primalV_);
      }
      else {
        precond(Mv,v,x,tol);
      }
    }

    void project( Vector<Real> &x ) {
      con_->project(x);
    } 

    void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x ) {
      con_->pruneActive(v,g,x,eps_);
    }

    void pruneActive( Vector<Real> &v, const Vector<Real> &x ) {
      con_->pruneActive(v,x,eps_);
    }

    void pruneInactive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x ) {
      con_->pruneInactive(v,g,x,eps_);
    }

    void pruneInactive( Vector<Real> &v, const Vector<Real> &x ) {
      con_->pruneInactive(v,x,eps_);
    }

    bool isFeasible( const Vector<Real> &v ) {
      return con_->isFeasible(v);
    }

    bool isConActivated(void) {
      return con_->isActivated();
    }

    void computeProjectedStep( Vector<Real> &v, const Vector<Real> &x ) {
      con_->computeProjectedStep(v,x);
    } 
  }; 

} // namespace ROL

#endif 
