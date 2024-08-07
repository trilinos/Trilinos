// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for the equality constrained NLP
            from Nocedal/Wright, 2nd edition, page 574, example 18.2;
            note the typo in reversing the initial guess and the solution.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_SIMPLEEQCONSTRAINED_HPP
#define ROL_SIMPLEEQCONSTRAINED_HPP

#include "ROL_TestProblem.hpp"
#include "ROL_StdVector.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

namespace ROL {
namespace ZOO {

  /** \brief Objective function:
             f(x) = exp(x1*x2*x3*x4*x5) + 0.5*(x1^3+x2^3+1)^2
   */
  template< class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real> >
  class Objective_SimpleEqConstrained : public Objective<Real> {

  typedef std::vector<Real> vector;
  typedef Vector<Real>      V;

  typedef typename vector::size_type uint;
   

  private:

    template<class VectorType>
    ROL::Ptr<const vector> getVector( const V& x ) {
      
      return dynamic_cast<const VectorType&>(x).getVector();
    }

    template<class VectorType>
    ROL::Ptr<vector> getVector( V& x ) {
      
      return dynamic_cast<VectorType&>(x).getVector();
    }

  public:
    Objective_SimpleEqConstrained() {}

    Real value( const Vector<Real> &x, Real &tol ) {
 
     
     ROL::Ptr<const vector> xp = getVector<XPrim>(x); 

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, objective value): "
                                                                   "Primal vector x must be of length 5.");

      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];
      Real x3 = (*xp)[2];
      Real x4 = (*xp)[3];
      Real x5 = (*xp)[4];

      Real arg = x1*x2*x3*x4*x5;
      Real val = -0.5*pow(pow(x1,3)+pow(x2,3)+1.0,2);
      if (std::abs(arg) < 1e5) val += exp(x1*x2*x3*x4*x5);
      else if (arg > 1e5) val += 1e10; 
      //Real val = exp(x1*x2*x3*x4*x5) - 0.5 * pow( (pow(x1,3)+pow(x2,3)+1.0), 2);

      return val;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<vector> gp = getVector<XDual>(g); 

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, objective gradient): "
                                                                   " Primal vector x must be of length 5."); 

      n = gp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, objective gradient): "
                                                                   "Gradient vector g must be of length 5."); 

      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];
      Real x3 = (*xp)[2];
      Real x4 = (*xp)[3];
      Real x5 = (*xp)[4];

      Real expxi = exp(x1*x2*x3*x4*x5);

      (*gp)[0] = x2*x3*x4*x5 * expxi - 3*pow(x1,2) * (pow(x1,3) + pow(x2,3) + 1);
      (*gp)[1] = x1*x3*x4*x5 * expxi - 3*pow(x2,2) * (pow(x1,3) + pow(x2,3) + 1);
      (*gp)[2] = x1*x2*x4*x5 * expxi;
      (*gp)[3] = x1*x2*x3*x5 * expxi;
      (*gp)[4] = x1*x2*x3*x4 * expxi;
    }

    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<const vector> vp = getVector<XPrim>(v);
      ROL::Ptr<vector> hvp = getVector<XDual>(hv);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, objective hessVec): "
                                                                   "Primal vector x must be of length 5."); 

      n = vp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, objective hessVec): "
                                                                   "Input vector v must be of length 5."); 

      n = hvp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, objective hessVec): "
                                                                   "Output vector hv must be of length 5."); 

      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];
      Real x3 = (*xp)[2];
      Real x4 = (*xp)[3];
      Real x5 = (*xp)[4];

      Real v1 = (*vp)[0];
      Real v2 = (*vp)[1];
      Real v3 = (*vp)[2];
      Real v4 = (*vp)[3];
      Real v5 = (*vp)[4];

      Real expxi = exp(x1*x2*x3*x4*x5);

      (*hvp)[0] = ( pow(x2,2)*pow(x3,2)*pow(x4,2)*pow(x5,2)*expxi-9.0*pow(x1,4)-6.0*(pow(x1,3)+pow(x2,3)+1.0)*x1 ) * v1  +
                  ( x3*x4*x5*expxi+x2*pow(x3,2)*pow(x4,2)*pow(x5,2)*x1*expxi-9.0*pow(x2,2)*pow(x1,2) ) * v2 +
                  ( x2*x4*x5*expxi+pow(x2,2)*x3*pow(x4,2)*pow(x5,2)*x1*expxi ) * v3 +
                  ( x2*x3*x5*expxi+pow(x2,2)*pow(x3,2)*x4*pow(x5,2)*x1*expxi ) * v4 +
                  ( x2*x3*x4*expxi+pow(x2,2)*pow(x3,2)*pow(x4,2)*x5*x1*expxi ) * v5;

      (*hvp)[1] = ( x3*x4*x5*expxi+x2*pow(x3,2)*pow(x4,2)*pow(x5,2)*x1*expxi-9.0*pow(x2,2)*pow(x1,2) ) * v1  +
                  ( pow(x1,2)*pow(x3,2)*pow(x4,2)*pow(x5,2)*expxi-9.0*pow(x2,4)-6.0*(pow(x1,3)+pow(x2,3)+1.0)*x2 ) * v2  +
                  ( x1*x4*x5*expxi+pow(x1,2)*x3*pow(x4,2)*pow(x5,2)*x2*expxi ) * v3  +
                  ( x1*x3*x5*expxi+pow(x1,2)*pow(x3,2)*x4*pow(x5,2)*x2*expxi ) * v4  +
                  ( x1*x3*x4*expxi+pow(x1,2)*pow(x3,2)*pow(x4,2)*x5*x2*expxi ) * v5;

      (*hvp)[2] = ( x2*x4*x5*expxi+pow(x2,2)*x3*pow(x4,2)*pow(x5,2)*x1*expxi ) * v1  +
                  ( x1*x4*x5*expxi+pow(x1,2)*x3*pow(x4,2)*pow(x5,2)*x2*expxi ) * v2  +
                  ( pow(x1,2)*pow(x2,2)*pow(x4,2)*pow(x5,2)*expxi ) * v3  +
                  ( x1*x2*x5*expxi+pow(x1,2)*pow(x2,2)*x4*pow(x5,2)*x3*expxi ) * v4  +
                  ( x1*x2*x4*expxi+pow(x1,2)*pow(x2,2)*pow(x4,2)*x5*x3*expxi ) * v5;

      (*hvp)[3] = ( x2*x3*x5*expxi+pow(x2,2)*pow(x3,2)*x4*pow(x5,2)*x1*expxi ) * v1  +
                  ( x1*x3*x5*expxi+pow(x1,2)*pow(x3,2)*x4*pow(x5,2)*x2*expxi ) * v2  +
                  ( x1*x2*x5*expxi+pow(x1,2)*pow(x2,2)*x4*pow(x5,2)*x3*expxi ) * v3  +
                  ( pow(x1,2)*pow(x2,2)*pow(x3,2)*pow(x5,2)*expxi ) * v4  +
                  ( x1*x2*x3*expxi+pow(x1,2)*pow(x2,2)*pow(x3,2)*x5*x4*expxi ) * v5;

      (*hvp)[4] = ( x2*x3*x4*expxi+pow(x2,2)*pow(x3,2)*pow(x4,2)*x5*x1*expxi ) * v1  +
                  ( x1*x3*x4*expxi+pow(x1,2)*pow(x3,2)*pow(x4,2)*x5*x2*expxi ) * v2  +
                  ( x1*x2*x4*expxi+pow(x1,2)*pow(x2,2)*pow(x4,2)*x5*x3*expxi ) * v3  +
                  ( x1*x2*x3*expxi+pow(x1,2)*pow(x2,2)*pow(x3,2)*x5*x4*expxi ) * v4  +
                  ( pow(x1,2)*pow(x2,2)*pow(x3,2)*pow(x4,2)*expxi ) * v5;
    }

  };


  /** \brief Equality constraints c_i(x) = 0, where:
             c1(x) = x1^2+x2^2+x3^2+x4^2+x5^2 - 10
             c2(x) = x2*x3-5*x4*x5
             c3(x) = x1^3 + x2^3 + 1
   */
  template<class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real>, class CPrim=StdVector<Real>, class CDual=StdVector<Real> >
  class EqualityConstraint_SimpleEqConstrained : public Constraint<Real> {

    typedef std::vector<Real> vector;
    typedef Vector<Real>      V;

    typedef typename vector::size_type uint;

  private:
    template<class VectorType>
    ROL::Ptr<const vector> getVector( const V& x ) {
      
      return dynamic_cast<const VectorType&>(x).getVector();
    }

    template<class VectorType> 
    ROL::Ptr<vector> getVector( V& x ) {
      
      return dynamic_cast<VectorType&>(x).getVector(); 
    }

  public:
    EqualityConstraint_SimpleEqConstrained() {}

    void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<vector> cp = getVector<CPrim>(c);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint value): "
                                                                   "Primal vector x must be of length 5.");

      uint m = cp->size();
      ROL_TEST_FOR_EXCEPTION( (m != 3), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint value): "
                                                                   "Constraint vector c must be of length 3.");

      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];
      Real x3 = (*xp)[2];
      Real x4 = (*xp)[3];
      Real x5 = (*xp)[4];

      (*cp)[0] = x1*x1+x2*x2+x3*x3+x4*x4+x5*x5 - 10.0;
      (*cp)[1] = x2*x3 - 5.0*x4*x5;
      (*cp)[2] = x1*x1*x1 + x2*x2*x2 + 1.0;
    }
  
    void applyJacobian( Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<const vector> vp = getVector<XPrim>(v);
      ROL::Ptr<vector> jvp = getVector<CPrim>(jv);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyJacobian): "
                                                                   "Primal vector x must be of length 5.");

      uint d = vp->size();
      ROL_TEST_FOR_EXCEPTION( (d != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyJacobian): "
                                                                   "Input vector v must be of length 5.");
      d = jvp->size();
      ROL_TEST_FOR_EXCEPTION( (d != 3), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyJacobian): "
                                                                   "Output vector jv must be of length 3.");
      
      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];
      Real x3 = (*xp)[2];
      Real x4 = (*xp)[3];
      Real x5 = (*xp)[4];

      Real v1 = (*vp)[0];
      Real v2 = (*vp)[1];
      Real v3 = (*vp)[2];
      Real v4 = (*vp)[3];
      Real v5 = (*vp)[4];

      (*jvp)[0] = 2.0*(x1*v1+x2*v2+x3*v3+x4*v4+x5*v5);
      (*jvp)[1] = x3*v2+x2*v3-5.0*x5*v4-5.0*x4*v5;
      (*jvp)[2] = 3.0*x1*x1*v1+3.0*x2*x2*v2;

    } //applyJacobian

    void applyAdjointJacobian( Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<const vector> vp = getVector<CDual>(v);
      ROL::Ptr<vector> ajvp = getVector<XDual>(ajv);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyAdjointJacobian): "
                                                                   "Primal vector x must be of length 5.");

      uint d = vp->size();
      ROL_TEST_FOR_EXCEPTION( (d != 3), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyAdjointJacobian): "
                                                                   "Input vector v must be of length 3.");

      d = ajvp->size();
      ROL_TEST_FOR_EXCEPTION( (d != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyAdjointJacobian): "
                                                                   "Output vector ajv must be of length 5.");
      
      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];
      Real x3 = (*xp)[2];
      Real x4 = (*xp)[3];
      Real x5 = (*xp)[4];

      Real v1 = (*vp)[0];
      Real v2 = (*vp)[1];
      Real v3 = (*vp)[2];

      (*ajvp)[0] = 2.0*x1*v1+3.0*x1*x1*v3;
      (*ajvp)[1] = 2.0*x2*v1+x3*v2+3.0*x2*x2*v3;
      (*ajvp)[2] = 2.0*x3*v1+x2*v2;
      (*ajvp)[3] = 2.0*x4*v1-5.0*x5*v2;
      (*ajvp)[4] = 2.0*x5*v1-5.0*x4*v2;

    } //applyAdjointJacobian

    void applyAdjointHessian( Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      
      ROL::Ptr<const vector> xp = getVector<XPrim>(x);
      ROL::Ptr<const vector> up = getVector<CDual>(u);
      ROL::Ptr<const vector> vp = getVector<XPrim>(v);
      ROL::Ptr<vector> ahuvp = getVector<XDual>(ahuv);

      uint n = xp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyAdjointHessian): "
                                                                   "Primal vector x must be of length 5.");

      n = vp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyAdjointHessian): "
                                                                   "Direction vector v must be of length 5.");

      n = ahuvp->size();
      ROL_TEST_FOR_EXCEPTION( (n != 5), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyAdjointHessian): "
                                                                   "Output vector ahuv must be of length 5.");
      uint d = up->size();
      ROL_TEST_FOR_EXCEPTION( (d != 3), std::invalid_argument, ">>> ERROR (ROL_SimpleEqConstrained, constraint applyAdjointHessian): "
                                                                   "Dual constraint vector u must be of length 3.");
      
      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];

      Real v1 = (*vp)[0];
      Real v2 = (*vp)[1];
      Real v3 = (*vp)[2];
      Real v4 = (*vp)[3];
      Real v5 = (*vp)[4];

      Real u1 = (*up)[0];
      Real u2 = (*up)[1];
      Real u3 = (*up)[2];

      (*ahuvp)[0] = 2.0*u1*v1 +             6.0*u3*x1*v1;
      (*ahuvp)[1] = 2.0*u1*v2 +     u2*v3 + 6.0*u3*x2*v2;
      (*ahuvp)[2] = 2.0*u1*v3 +     u2*v2;
      (*ahuvp)[3] = 2.0*u1*v4 - 5.0*u2*v5;
      (*ahuvp)[4] = 2.0*u1*v5 - 5.0*u2*v4;

    } //applyAdjointHessian

    /*std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2, const Vector<Real> &b1, const Vector<Real> &b2, const Vector<Real> &x, Real &tol) {
      ROL::Ptr<std::vector<Real> > v1p =
        ROL::constPtrCast<std::vector<Real> >((dynamic_cast<XPrim&>(v1)).getVector());    
      ROL::Ptr<std::vector<Real> > v2p =
        ROL::constPtrCast<std::vector<Real> >((dynamic_cast<CDual&>(v2)).getVector());
      ROL::Ptr<const std::vector<Real> > b1p =
        (dynamic_cast<XDual>(const_cast<Vector<Real> &&>(b1))).getVector();
      ROL::Ptr<const std::vector<Real> > b2p =
        (dynamic_cast<CPrim>(const_cast<Vector<Real> &&>(b2))).getVector();
      ROL::Ptr<const std::vector<Real> > xp =
        (dynamic_cast<XPrim>(const_cast<Vector<Real> &&>(x))).getVector();

      Real x1 = (*xp)[0];
      Real x2 = (*xp)[1];
      Real x3 = (*xp)[2];
      Real x4 = (*xp)[3];
      Real x5 = (*xp)[4];

      int i = 0;

      // Assemble augmented system.
      Teuchos::SerialDenseMatrix<int, Real> augmat(8,8);
      for (i=0; i<5; i++) {
        augmat(i,i) = 1.0;
      }
      augmat(5,0) = 2.0*x1;    augmat(5,1) = 2.0*x2;    augmat(5,2) = 2.0*x3; augmat(5,3) = 2.0*x4;  augmat(5,4) = 2.0*x5;
                               augmat(6,1) = x3;        augmat(6,2) = x2;     augmat(6,3) = -5.0*x5; augmat(6,4) = -5.0*x4;
      augmat(7,0) = 3.0*x1*x1; augmat(7,1) = 3.0*x2*x2;
      augmat(0,5) = 2.0*x1;    augmat(1,5) = 2.0*x2;    augmat(2,5) = 2.0*x3; augmat(3,5) = 2.0*x4;  augmat(4,5) = 2.0*x5;
                               augmat(1,6) = x3;        augmat(2,6) = x2;     augmat(3,6) = -5.0*x5; augmat(4,6) = -5.0*x4;
      augmat(0,7) = 3.0*x1*x1; augmat(1,7) = 3.0*x2*x2;
      Teuchos::SerialDenseVector<int, Real> lhs(8);
      Teuchos::SerialDenseVector<int, Real> rhs(8);
      for (i=0; i<5; i++) {
        rhs(i) = (*b1p)[i];
      }
      for (i=5; i<8; i++) {
        rhs(i) = (*b2p)[i-5];
      }

      // Solve augmented system.
      Teuchos::SerialDenseSolver<int, Real> augsolver;
      augsolver.setMatrix(&augmat, false);
      augsolver.setVectors(&lhs, false), Teuchos::&rhs, false;
      augsolver.solve();

      // Retrieve solution.
      for (i=0; i<5; i++) {
        (*v1p)[i] = lhs(i);
      }
      for (i=0; i<3; i++) {
        (*v2p)[i] = lhs(i+5);
      }

      return std::vector<Real>(0);
        
    }*/ //solveAugmentedSystem

  };


  template<class Real, class XPrim=StdVector<Real>, class XDual=StdVector<Real>, class CPrim=StdVector<Real>, class CDual=StdVector<Real> >
  class getSimpleEqConstrained : public TestProblem<Real> {
    typedef std::vector<Real> vector;
    typedef typename vector::size_type uint;
  public:
    getSimpleEqConstrained(void) {}

    Ptr<Objective<Real>> getObjective(void) const {
      // Instantiate objective function.
      return ROL::makePtr<Objective_SimpleEqConstrained<Real,XPrim,XDual>>();
    }

    Ptr<Vector<Real>> getInitialGuess(void) const { 
      uint n = 5;
      // Get initial guess.
      Ptr<vector> x0p = makePtr<vector>(n,0);
      (*x0p)[0] = -1.8;
      (*x0p)[1] = 1.7;
      (*x0p)[2] = 1.9;
      (*x0p)[3] = -0.8;
      (*x0p)[4] = -0.8;
      return makePtr<XPrim>(x0p);
    }

    Ptr<Vector<Real>> getSolution(const int i = 0) const { 
      uint n = 5;
      // Get solution.
      Ptr<vector> solp = makePtr<vector>(n,0);
      (*solp)[0] = -1.717143570394391e+00;
      (*solp)[1] =  1.595709690183565e+00;
      (*solp)[2] =  1.827245752927178e+00;
      (*solp)[3] = -7.636430781841294e-01;
      (*solp)[4] = -7.636430781841294e-01;
      return makePtr<XPrim>(solp);
    }

    Ptr<Constraint<Real>> getEqualityConstraint(void) const {
      // Instantiate constraints.
      return ROL::makePtr<EqualityConstraint_SimpleEqConstrained<Real,XPrim,XDual,CPrim,CDual>>();
    }

    Ptr<Vector<Real>> getEqualityMultiplier(void) const {
      Ptr<vector> lp = makePtr<vector>(3,0);
      return makePtr<CDual>(lp);
    }
  };

} // End ZOO Namespace
} // End ROL Namespace

#endif
