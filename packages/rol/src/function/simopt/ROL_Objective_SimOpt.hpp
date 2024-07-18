// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OBJECTIVE_SIMOPT_H
#define ROL_OBJECTIVE_SIMOPT_H

#include "ROL_Objective.hpp"
#include "ROL_Vector_SimOpt.hpp"

/** @ingroup func_group
    \class ROL::Objective_SimOpt
    \brief Provides the interface to evaluate simulation-based objective functions.
*/


namespace ROL {

template <class Real>
class Objective_SimOpt : public Objective<Real> {
public:

  /** \brief Update objective function.  
                u is an iterate, 
                z is an iterate, 
                flag = true if the iterate has changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    this->update(*(xs.get_1()),*(xs.get_2()),flag,iter);
  }

  virtual void update( const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter = -1 ) {}

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    this->update(*(xs.get_1()),*(xs.get_2()),type,iter);
  }


  /** \brief Compute value.
  */
  virtual Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) = 0;

  Real value( const Vector<Real> &x, Real &tol ) {
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    return this->value(*(xs.get_1()),*(xs.get_2()),tol);
  }


  /** \brief Compute gradient with respect to first component.
  */
  virtual void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    Real ftol  = std::sqrt(ROL_EPSILON<Real>());
    Real h     = 0.0;
    this->update(u,z,UpdateType::Temp);
    Real v     = this->value(u,z,ftol);
    Real deriv = 0.0;
    ROL::Ptr<Vector<Real> > unew = u.clone();
    g.zero();
    for (int i = 0; i < g.dimension(); i++) {
      h = u.dot(*u.basis(i))*tol;
      unew->set(u);
      unew->axpy(h,*(u.basis(i)));
      this->update(*unew,z,UpdateType::Temp);
      deriv = (this->value(*unew,z,ftol) - v)/h;
      g.axpy(deriv,*(g.basis(i)));
    }
    this->update(u,z,UpdateType::Temp);
  }
  /** \brief Compute gradient with respect to second component.
  */
  virtual void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    Real ftol  = std::sqrt(ROL_EPSILON<Real>());
    Real h     = 0.0;
    this->update(u,z,UpdateType::Temp);
    Real v     = this->value(u,z,ftol);
    Real deriv = 0.0;
    ROL::Ptr<Vector<Real> > znew = z.clone();
    g.zero();
    for (int i = 0; i < g.dimension(); i++) {
      h = z.dot(*z.basis(i))*tol;
      znew->set(z);
      znew->axpy(h,*(z.basis(i)));
      this->update(u,*znew,UpdateType::Temp);
      deriv = (this->value(u,*znew,ftol) - v)/h;
      g.axpy(deriv,*(g.basis(i)));
    }
    this->update(u,z,UpdateType::Temp);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    ROL::Vector_SimOpt<Real> &gs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(g));
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    ROL::Ptr<Vector<Real> > g1 = gs.get_1()->clone();
    ROL::Ptr<Vector<Real> > g2 = gs.get_2()->clone();
    this->gradient_1(*g1,*(xs.get_1()),*(xs.get_2()),tol);
    this->gradient_2(*g2,*(xs.get_1()),*(xs.get_2()),tol);
    gs.set_1(*g1);
    gs.set_2(*g2);
  }


  /** \brief Apply Hessian approximation to vector.
  */
  virtual void hessVec_11( Vector<Real> &hv, const Vector<Real> &v, 
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate gradient of first component at (u+hv,z)
    ROL::Ptr<Vector<Real> > unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    this->update(*unew,z,UpdateType::Temp);
    hv.zero();
    this->gradient_1(hv,*unew,z,gtol);
    // Evaluate gradient of first component at (u,z)
    ROL::Ptr<Vector<Real> > g = hv.clone();
    this->update(u,z,UpdateType::Temp);
    this->gradient_1(*g,u,z,gtol);
    // Compute Newton quotient
    hv.axpy(-1.0,*g);
    hv.scale(1.0/h);
  }

  virtual void hessVec_12( Vector<Real> &hv, const Vector<Real> &v, 
                           const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate gradient of first component at (u,z+hv)
    ROL::Ptr<Vector<Real> > znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    this->update(u,*znew,UpdateType::Temp);
    hv.zero();
    this->gradient_1(hv,u,*znew,gtol);
    // Evaluate gradient of first component at (u,z)
    ROL::Ptr<Vector<Real> > g = hv.clone();
    this->update(u,z,UpdateType::Temp);
    this->gradient_1(*g,u,z,gtol);
    // Compute Newton quotient
    hv.axpy(-1.0,*g);
    hv.scale(1.0/h);
  }

  virtual void hessVec_21( Vector<Real> &hv, const Vector<Real> &v, 
                           const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate gradient of first component at (u+hv,z)
    ROL::Ptr<Vector<Real> > unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    this->update(*unew,z,UpdateType::Temp);
    hv.zero();
    this->gradient_2(hv,*unew,z,gtol);
    // Evaluate gradient of first component at (u,z)
    ROL::Ptr<Vector<Real> > g = hv.clone();
    this->update(u,z,UpdateType::Temp);
    this->gradient_2(*g,u,z,gtol);
    // Compute Newton quotient
    hv.axpy(-1.0,*g);
    hv.scale(1.0/h);
  }

  virtual void hessVec_22( Vector<Real> &hv, const Vector<Real> &v, 
                     const Vector<Real> &u,  const Vector<Real> &z, Real &tol ) {
    Real gtol = std::sqrt(ROL_EPSILON<Real>());
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON<Real>())) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate gradient of first component at (u,z+hv)
    ROL::Ptr<Vector<Real> > znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    this->update(u,*znew,UpdateType::Temp);
    hv.zero();
    this->gradient_2(hv,u,*znew,gtol);
    // Evaluate gradient of first component at (u,z)
    ROL::Ptr<Vector<Real> > g = hv.clone();
    this->update(u,z,UpdateType::Temp);
    this->gradient_2(*g,u,z,gtol);
    // Compute Newton quotient
    hv.axpy(-1.0,*g);
    hv.scale(1.0/h);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    ROL::Vector_SimOpt<Real> &hvs = dynamic_cast<ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<ROL::Vector<Real>&>(hv));
    const ROL::Vector_SimOpt<Real> &vs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(v));
    const ROL::Vector_SimOpt<Real> &xs = dynamic_cast<const ROL::Vector_SimOpt<Real>&>(
      dynamic_cast<const ROL::Vector<Real>&>(x));
    ROL::Ptr<Vector<Real> > h11 = (hvs.get_1())->clone();
    this->hessVec_11(*h11,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    ROL::Ptr<Vector<Real> > h12 = (hvs.get_1())->clone();
    this->hessVec_12(*h12,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    ROL::Ptr<Vector<Real> > h21 = (hvs.get_2())->clone();
    this->hessVec_21(*h21,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    ROL::Ptr<Vector<Real> > h22 = (hvs.get_2())->clone();
    this->hessVec_22(*h22,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    h11->plus(*h12);
    hvs.set_1(*h11);
    h22->plus(*h21);
    hvs.set_2(*h22);
  }

  std::vector<std::vector<Real> > checkGradient_1( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &d,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                   const int order = 1 ) {
    return checkGradient_1(u, z, u.dual(), d, printToStream, outStream, numSteps, order);
  }

  std::vector<std::vector<Real> > checkGradient_1( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &g,
                                                   const Vector<Real> &d,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int numSteps,
                                                   const int order ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
  
    return checkGradient_1(u,z,g,d,steps,printToStream,outStream,order);
  } // checkGradient_1

  std::vector<std::vector<Real> > checkGradient_1( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &g,
                                                   const Vector<Real> &d,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int order ) {
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument, 
                                "Error: finite difference order must be 1,2,3, or 4" );
  
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real> > gCheck(numSteps, tmp);
  
    // Save the format state of the original outStream.
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Evaluate objective value at x.
    this->update(u,z,UpdateType::Temp);
    Real val = this->value(u,z,tol);
  
    // Compute gradient at x.
    ROL::Ptr<Vector<Real> > gtmp = g.clone();
    this->gradient_1(*gtmp, u, z, tol);
    //Real dtg = d.dot(gtmp->dual());
    Real dtg = d.apply(*gtmp);
  
    // Temporary vectors.
    ROL::Ptr<Vector<Real> > unew = u.clone();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      unew->set(u);
  
      // Compute gradient, finite-difference gradient, and absolute error.
      gCheck[i][0] = eta;
      gCheck[i][1] = dtg;
  
      gCheck[i][2] = weights[order-1][0] * val;
  
      for(int j=0; j<order; ++j) {
        // Evaluate at x <- x+eta*c_i*d.
        unew->axpy(eta*shifts[order-1][j], d);
  
        // Only evaluate at shifts where the weight is nonzero  
        if( weights[order-1][j+1] != 0 ) {
          this->update(*unew,z,UpdateType::Temp);
          gCheck[i][2] += weights[order-1][j+1] * this->value(*unew,z,tol);
        }
      }
      gCheck[i][2] /= eta;
  
      gCheck[i][3] = std::abs(gCheck[i][2] - gCheck[i][1]);
  
      if (printToStream) {
        if (i==0) {
          outStream << std::right
                    << std::setw(20) << "Step size"
                    << std::setw(20) << "grad'*dir"
                    << std::setw(20) << "FD approx"
                    << std::setw(20) << "abs error"
                    << "\n"
                    << std::setw(20) << "---------"
                    << std::setw(20) << "---------"
                    << std::setw(20) << "---------"
                    << std::setw(20) << "---------"
                    << "\n";
        }
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << gCheck[i][0]
                  << std::setw(20) << gCheck[i][1]
                  << std::setw(20) << gCheck[i][2]
                  << std::setw(20) << gCheck[i][3]
                  << "\n";
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return gCheck;
  } // checkGradient_1


  std::vector<std::vector<Real> > checkGradient_2( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &d,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                   const int order = 1 ) {
    return checkGradient_2(u, z, z.dual(), d, printToStream, outStream, numSteps, order);
  }

  std::vector<std::vector<Real> > checkGradient_2( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &g,
                                                   const Vector<Real> &d,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int numSteps,
                                                   const int order ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
  
    return checkGradient_2(u,z,g,d,steps,printToStream,outStream,order);
  } // checkGradient_2

  std::vector<std::vector<Real> > checkGradient_2( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &g,
                                                   const Vector<Real> &d,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int order ) {
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument, 
                                "Error: finite difference order must be 1,2,3, or 4" );
  
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real> > gCheck(numSteps, tmp);
  
    // Save the format state of the original outStream.
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Evaluate objective value at x.
    this->update(u,z,UpdateType::Temp);
    Real val = this->value(u,z,tol);
  
    // Compute gradient at x.
    ROL::Ptr<Vector<Real> > gtmp = g.clone();
    this->gradient_2(*gtmp, u, z, tol);
    //Real dtg = d.dot(gtmp->dual());
    Real dtg = d.apply(*gtmp);
  
    // Temporary vectors.
    ROL::Ptr<Vector<Real> > znew = z.clone();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      znew->set(z);
  
      // Compute gradient, finite-difference gradient, and absolute error.
      gCheck[i][0] = eta;
      gCheck[i][1] = dtg;
  
      gCheck[i][2] = weights[order-1][0] * val;
  
      for(int j=0; j<order; ++j) {
        // Evaluate at x <- x+eta*c_i*d.
        znew->axpy(eta*shifts[order-1][j], d);
  
        // Only evaluate at shifts where the weight is nonzero  
        if( weights[order-1][j+1] != 0 ) {
          this->update(u,*znew,UpdateType::Temp);
          gCheck[i][2] += weights[order-1][j+1] * this->value(u,*znew,tol);
        }
      }
      gCheck[i][2] /= eta;
  
      gCheck[i][3] = std::abs(gCheck[i][2] - gCheck[i][1]);
  
      if (printToStream) {
        if (i==0) {
          outStream << std::right
                    << std::setw(20) << "Step size"
                    << std::setw(20) << "grad'*dir"
                    << std::setw(20) << "FD approx"
                    << std::setw(20) << "abs error"
                    << "\n"
                    << std::setw(20) << "---------"
                    << std::setw(20) << "---------"
                    << std::setw(20) << "---------"
                    << std::setw(20) << "---------"
                    << "\n";
        }
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << gCheck[i][0]
                  << std::setw(20) << gCheck[i][1]
                  << std::setw(20) << gCheck[i][2]
                  << std::setw(20) << gCheck[i][3]
                  << "\n";
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return gCheck;
  } // checkGradient_2


  std::vector<std::vector<Real> > checkHessVec_11( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &v,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                   const int order = 1 ) {

    return checkHessVec_11(u, z, u.dual(), v, printToStream, outStream, numSteps, order);

  }

  std::vector<std::vector<Real> > checkHessVec_11( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &v,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int order = 1 ) {

    return checkHessVec_11(u, z, u.dual(), v, steps, printToStream, outStream, order);
  }


  std::vector<std::vector<Real> > checkHessVec_11( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &hv,
                                                   const Vector<Real> &v,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int numSteps,
                                                   const int order ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
  
    return checkHessVec_11(u,z,hv,v,steps,printToStream,outStream,order);
  } // checkHessVec_11


  std::vector<std::vector<Real> > checkHessVec_11( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &hv,
                                                   const Vector<Real> &v,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int order ) {
  
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument,
                                "Error: finite difference order must be 1,2,3, or 4" );
  
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real> > hvCheck(numSteps, tmp);
  
    // Save the format state of the original outStream.
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Compute gradient at x.
    ROL::Ptr<Vector<Real> > g = hv.clone();
    this->update(u,z,UpdateType::Temp);
    this->gradient_1(*g, u, z, tol);
  
    // Compute (Hessian at x) times (vector v).
    ROL::Ptr<Vector<Real> > Hv = hv.clone();
    this->hessVec_11(*Hv, v, u, z, tol);
    Real normHv = Hv->norm();
  
    // Temporary vectors.
    ROL::Ptr<Vector<Real> > gdif = hv.clone();
    ROL::Ptr<Vector<Real> > gnew = hv.clone();
    ROL::Ptr<Vector<Real> > unew = u.clone();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      // Evaluate objective value at x+eta*d.
      unew->set(u);
  
      gdif->set(*g);
      gdif->scale(weights[order-1][0]);
  
      for(int j=0; j<order; ++j) {
  
          // Evaluate at x <- x+eta*c_i*d.
          unew->axpy(eta*shifts[order-1][j], v);
  
          // Only evaluate at shifts where the weight is nonzero  
          if( weights[order-1][j+1] != 0 ) {
              this->update(*unew,z,UpdateType::Temp);
              this->gradient_1(*gnew, *unew, z, tol);
              gdif->axpy(weights[order-1][j+1],*gnew);
          }
  
      }
  
      gdif->scale(1.0/eta);
  
      // Compute norms of hessvec, finite-difference hessvec, and error.
      hvCheck[i][0] = eta;
      hvCheck[i][1] = normHv;
      hvCheck[i][2] = gdif->norm();
      gdif->axpy(-1.0, *Hv);
      hvCheck[i][3] = gdif->norm();
  
      if (printToStream) {
        if (i==0) {
        outStream << std::right
                  << std::setw(20) << "Step size"
                  << std::setw(20) << "norm(Hess*vec)"
                  << std::setw(20) << "norm(FD approx)"
                  << std::setw(20) << "norm(abs error)"
                  << "\n"
                  << std::setw(20) << "---------"
                  << std::setw(20) << "--------------"
                  << std::setw(20) << "---------------"
                  << std::setw(20) << "---------------"
                  << "\n";
        }
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << hvCheck[i][0]
                  << std::setw(20) << hvCheck[i][1]
                  << std::setw(20) << hvCheck[i][2]
                  << std::setw(20) << hvCheck[i][3]
                  << "\n";
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return hvCheck;
  } // checkHessVec_11


  std::vector<std::vector<Real> > checkHessVec_12( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &v,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                   const int order = 1 ) {
    return checkHessVec_12(u, z, u.dual(), v, printToStream, outStream, numSteps, order);
  }

  std::vector<std::vector<Real> > checkHessVec_12( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &v,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int order = 1 ) {
    return checkHessVec_12(u, z, u.dual(), v, steps, printToStream, outStream, order);
  }


  std::vector<std::vector<Real> > checkHessVec_12( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &hv,
                                                   const Vector<Real> &v,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int numSteps,
                                                   const int order ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
  
    return checkHessVec_12(u,z,hv,v,steps,printToStream,outStream,order);
  } // checkHessVec_12


  std::vector<std::vector<Real> > checkHessVec_12( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &hv,
                                                   const Vector<Real> &v,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int order ) {
  
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument,
                                "Error: finite difference order must be 1,2,3, or 4" );
  
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real> > hvCheck(numSteps, tmp);
  
    // Save the format state of the original outStream.
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Compute gradient at x.
    ROL::Ptr<Vector<Real> > g = hv.clone();
    this->update(u,z,UpdateType::Temp);
    this->gradient_1(*g, u, z, tol);
  
    // Compute (Hessian at x) times (vector v).
    ROL::Ptr<Vector<Real> > Hv = hv.clone();
    this->hessVec_12(*Hv, v, u, z, tol);
    Real normHv = Hv->norm();
  
    // Temporary vectors.
    ROL::Ptr<Vector<Real> > gdif = hv.clone();
    ROL::Ptr<Vector<Real> > gnew = hv.clone();
    ROL::Ptr<Vector<Real> > znew = z.clone();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      // Evaluate objective value at x+eta*d.
      znew->set(z);
  
      gdif->set(*g);
      gdif->scale(weights[order-1][0]);
  
      for(int j=0; j<order; ++j) {
  
          // Evaluate at x <- x+eta*c_i*d.
          znew->axpy(eta*shifts[order-1][j], v);
  
          // Only evaluate at shifts where the weight is nonzero  
          if( weights[order-1][j+1] != 0 ) {
              this->update(u,*znew,UpdateType::Temp);
              this->gradient_1(*gnew, u, *znew, tol);
              gdif->axpy(weights[order-1][j+1],*gnew);
          }
  
      }
  
      gdif->scale(1.0/eta);
  
      // Compute norms of hessvec, finite-difference hessvec, and error.
      hvCheck[i][0] = eta;
      hvCheck[i][1] = normHv;
      hvCheck[i][2] = gdif->norm();
      gdif->axpy(-1.0, *Hv);
      hvCheck[i][3] = gdif->norm();
  
      if (printToStream) {
        if (i==0) {
        outStream << std::right
                  << std::setw(20) << "Step size"
                  << std::setw(20) << "norm(Hess*vec)"
                  << std::setw(20) << "norm(FD approx)"
                  << std::setw(20) << "norm(abs error)"
                  << "\n"
                  << std::setw(20) << "---------"
                  << std::setw(20) << "--------------"
                  << std::setw(20) << "---------------"
                  << std::setw(20) << "---------------"
                  << "\n";
        }
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << hvCheck[i][0]
                  << std::setw(20) << hvCheck[i][1]
                  << std::setw(20) << hvCheck[i][2]
                  << std::setw(20) << hvCheck[i][3]
                  << "\n";
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return hvCheck;
  } // checkHessVec_12


  std::vector<std::vector<Real> > checkHessVec_21( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &v,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                   const int order = 1 ) {

    return checkHessVec_21(u, z, z.dual(), v, printToStream, outStream, numSteps, order);

  }

  std::vector<std::vector<Real> > checkHessVec_21( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &v,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int order = 1 ) {

    return checkHessVec_21(u, z, z.dual(), v, steps, printToStream, outStream, order);
  }


  std::vector<std::vector<Real> > checkHessVec_21( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &hv,
                                                   const Vector<Real> &v,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int numSteps,
                                                   const int order ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
  
    return checkHessVec_21(u,z,hv,v,steps,printToStream,outStream,order);
  } // checkHessVec_21


  std::vector<std::vector<Real> > checkHessVec_21( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &hv,
                                                   const Vector<Real> &v,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int order ) {
  
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument,
                                "Error: finite difference order must be 1,2,3, or 4" );
  
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real> > hvCheck(numSteps, tmp);
  
    // Save the format state of the original outStream.
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Compute gradient at x.
    ROL::Ptr<Vector<Real> > g = hv.clone();
    this->update(u,z,UpdateType::Temp);
    this->gradient_2(*g, u, z, tol);
  
    // Compute (Hessian at x) times (vector v).
    ROL::Ptr<Vector<Real> > Hv = hv.clone();
    this->hessVec_21(*Hv, v, u, z, tol);
    Real normHv = Hv->norm();
  
    // Temporary vectors.
    ROL::Ptr<Vector<Real> > gdif = hv.clone();
    ROL::Ptr<Vector<Real> > gnew = hv.clone();
    ROL::Ptr<Vector<Real> > unew = u.clone();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      // Evaluate objective value at x+eta*d.
      unew->set(u);
  
      gdif->set(*g);
      gdif->scale(weights[order-1][0]);
  
      for(int j=0; j<order; ++j) {
  
          // Evaluate at x <- x+eta*c_i*d.
          unew->axpy(eta*shifts[order-1][j], v);
  
          // Only evaluate at shifts where the weight is nonzero  
          if( weights[order-1][j+1] != 0 ) {
              this->update(*unew,z,UpdateType::Temp);
              this->gradient_2(*gnew, *unew, z, tol);
              gdif->axpy(weights[order-1][j+1],*gnew);
          }
  
      }
  
      gdif->scale(1.0/eta);
  
      // Compute norms of hessvec, finite-difference hessvec, and error.
      hvCheck[i][0] = eta;
      hvCheck[i][1] = normHv;
      hvCheck[i][2] = gdif->norm();
      gdif->axpy(-1.0, *Hv);
      hvCheck[i][3] = gdif->norm();
  
      if (printToStream) {
        if (i==0) {
        outStream << std::right
                  << std::setw(20) << "Step size"
                  << std::setw(20) << "norm(Hess*vec)"
                  << std::setw(20) << "norm(FD approx)"
                  << std::setw(20) << "norm(abs error)"
                  << "\n"
                  << std::setw(20) << "---------"
                  << std::setw(20) << "--------------"
                  << std::setw(20) << "---------------"
                  << std::setw(20) << "---------------"
                  << "\n";
        }
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << hvCheck[i][0]
                  << std::setw(20) << hvCheck[i][1]
                  << std::setw(20) << hvCheck[i][2]
                  << std::setw(20) << hvCheck[i][3]
                  << "\n";
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return hvCheck;
  } // checkHessVec_21


  std::vector<std::vector<Real> > checkHessVec_22( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &v,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                   const int order = 1 ) {

    return checkHessVec_22(u, z, z.dual(), v, printToStream, outStream, numSteps, order);

  }

  std::vector<std::vector<Real> > checkHessVec_22( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &v,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream = true,
                                                   std::ostream & outStream = std::cout,
                                                   const int order = 1 ) {

    return checkHessVec_22(u, z, z.dual(), v, steps, printToStream, outStream, order);
  }


  std::vector<std::vector<Real> > checkHessVec_22( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &hv,
                                                   const Vector<Real> &v,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int numSteps,
                                                   const int order ) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
  
    return checkHessVec_22(u,z,hv,v,steps,printToStream,outStream,order);
  } // checkHessVec_22


  std::vector<std::vector<Real> > checkHessVec_22( const Vector<Real> &u,
                                                   const Vector<Real> &z,
                                                   const Vector<Real> &hv,
                                                   const Vector<Real> &v,
                                                   const std::vector<Real> &steps,
                                                   const bool printToStream,
                                                   std::ostream & outStream,
                                                   const int order ) {
  
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument,
                                "Error: finite difference order must be 1,2,3, or 4" );
  
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
  
  
    Real tol = std::sqrt(ROL_EPSILON<Real>());
  
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real> > hvCheck(numSteps, tmp);
  
    // Save the format state of the original outStream.
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
  
    // Compute gradient at x.
    ROL::Ptr<Vector<Real> > g = hv.clone();
    this->update(u,z,UpdateType::Temp);
    this->gradient_2(*g, u, z, tol);
  
    // Compute (Hessian at x) times (vector v).
    ROL::Ptr<Vector<Real> > Hv = hv.clone();
    this->hessVec_22(*Hv, v, u, z, tol);
    Real normHv = Hv->norm();
  
    // Temporary vectors.
    ROL::Ptr<Vector<Real> > gdif = hv.clone();
    ROL::Ptr<Vector<Real> > gnew = hv.clone();
    ROL::Ptr<Vector<Real> > znew = z.clone();
  
    for (int i=0; i<numSteps; i++) {
  
      Real eta = steps[i];
  
      // Evaluate objective value at x+eta*d.
      znew->set(z);
  
      gdif->set(*g);
      gdif->scale(weights[order-1][0]);
  
      for(int j=0; j<order; ++j) {
  
          // Evaluate at x <- x+eta*c_i*d.
          znew->axpy(eta*shifts[order-1][j], v);
  
          // Only evaluate at shifts where the weight is nonzero  
          if( weights[order-1][j+1] != 0 ) {
              this->update(u,*znew,UpdateType::Temp);
              this->gradient_2(*gnew, u, *znew, tol);
              gdif->axpy(weights[order-1][j+1],*gnew);
          }
  
      }
  
      gdif->scale(1.0/eta);
  
      // Compute norms of hessvec, finite-difference hessvec, and error.
      hvCheck[i][0] = eta;
      hvCheck[i][1] = normHv;
      hvCheck[i][2] = gdif->norm();
      gdif->axpy(-1.0, *Hv);
      hvCheck[i][3] = gdif->norm();
  
      if (printToStream) {
        if (i==0) {
        outStream << std::right
                  << std::setw(20) << "Step size"
                  << std::setw(20) << "norm(Hess*vec)"
                  << std::setw(20) << "norm(FD approx)"
                  << std::setw(20) << "norm(abs error)"
                  << "\n"
                  << std::setw(20) << "---------"
                  << std::setw(20) << "--------------"
                  << std::setw(20) << "---------------"
                  << std::setw(20) << "---------------"
                  << "\n";
        }
        outStream << std::scientific << std::setprecision(11) << std::right
                  << std::setw(20) << hvCheck[i][0]
                  << std::setw(20) << hvCheck[i][1]
                  << std::setw(20) << hvCheck[i][2]
                  << std::setw(20) << hvCheck[i][3]
                  << "\n";
      }
  
    }
  
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
  
    return hvCheck;
  } // checkHessVec_22

}; // class Objective_SimOpt

} // namespace ROL

#endif
