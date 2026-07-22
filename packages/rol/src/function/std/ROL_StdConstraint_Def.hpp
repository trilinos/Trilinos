// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STDEQUALITY_CONSTRAINT_DEF_H
#define ROL_STDEQUALITY_CONSTRAINT_DEF_H

namespace ROL {

template<typename Real>
void StdConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  update(*(xs.getVector()),flag,iter);
}

template<typename Real>
void StdConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  update(*(xs.getVector()),type,iter);
}

template<typename Real>
void StdConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
  StdVector<Real> cs = dynamic_cast<StdVector<Real>&>(c);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  value(*(cs.getVector()),*(xs.getVector()),tol);
}

template<typename Real>
void StdConstraint<Real>::applyJacobian(Vector<Real> &jv,
                                  const Vector<Real> &v, 
                                  const Vector<Real> &x, Real &tol) {
  StdVector<Real> jvs = dynamic_cast<StdVector<Real>&>(jv);
  const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  try {
    applyJacobian(*(jvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);      
  } 
  catch (std::exception &e ){
    Constraint<Real>::applyJacobian(jv,v,x,tol);
  }
}

template<typename Real>
void StdConstraint<Real>::applyJacobian( std::vector<Real> &jv,
                                   const std::vector<Real> &v, 
                                   const std::vector<Real> &x, Real &tol ) {
  ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
    ">>> ERROR (ROL::StdConstraint): applyJacobian not implemented!");
}

template<typename Real>
void StdConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv,
                                         const Vector<Real> &v,
                                         const Vector<Real> &x, Real &tol) {
  StdVector<Real> ajvs = dynamic_cast<StdVector<Real>&>(ajv);
  const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  try {
     applyAdjointJacobian(*(ajvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);      
  } 
  catch (std::exception &e ){
    Constraint<Real>::applyAdjointJacobian(ajv,v,x,tol);
  }
}

template<typename Real>
void StdConstraint<Real>::applyAdjointJacobian( std::vector<Real> &ajv,
                                          const std::vector<Real> &v, 
                                          const std::vector<Real> &x, Real &tol ) {
  ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
    ">>> ERROR (ROL::StdConstraint): applyAdjointJacobian not implemented!");
}

template<typename Real>
void StdConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv,
                                        const Vector<Real> &u,
                                        const Vector<Real> &v,
                                        const Vector<Real> &x, Real &tol) {
  StdVector<Real> ahuvs = dynamic_cast<StdVector<Real>&>(ahuv);
  const StdVector<Real> us = dynamic_cast<const StdVector<Real>&>(u);
  const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  try {
    applyAdjointHessian( *(ahuvs.getVector()), *(us.getVector()), *(vs.getVector()), 
                         *(xs.getVector()), tol );
  }
  catch (std::exception &e) {
    Constraint<Real>::applyAdjointHessian(ahuv,u,v,x,tol);
  }   
}

template<typename Real>
void StdConstraint<Real>::applyAdjointHessian( std::vector<Real> &ahuv,
                                         const std::vector<Real> &u,
                                         const std::vector<Real> &v,
                                         const std::vector<Real> &x, Real &tol ) {
  ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
    ">>> ERROR (ROL::StdConstraint) : applyAdjointHessian not implemented!");
}

template<typename Real>
std::vector<Real> StdConstraint<Real>::solveAugmentedSystem(Vector<Real> &v1,
                                                            Vector<Real> &v2,
                                                      const Vector<Real> &b1,
                                                      const Vector<Real> &b2,
                                                      const Vector<Real> &x, Real &tol) {
  StdVector<Real> v1s = dynamic_cast<StdVector<Real>&>(v1);
  StdVector<Real> v2s = dynamic_cast<StdVector<Real>&>(v2);
  const StdVector<Real> b1s = dynamic_cast<const StdVector<Real>&>(b1);
  const StdVector<Real> b2s = dynamic_cast<const StdVector<Real>&>(b2);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  try {
    return solveAugmentedSystem( *(v1s.getVector()), *(v2s.getVector()), *(b1s.getVector()),
                          *(b2s.getVector()), *(xs.getVector()), tol );
  }
  catch (std::exception &e) {
    return Constraint<Real>::solveAugmentedSystem(v1,v2,b1,b2,x,tol); 
  }
}

template<typename Real>
std::vector<Real> StdConstraint<Real>::solveAugmentedSystem( std::vector<Real> &v1,
                                                             std::vector<Real> &v2,
                                                       const std::vector<Real> &b1,
                                                       const std::vector<Real> &b2,
                                                       const std::vector<Real> &x, Real tol ) {
  ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
    ">>> ERROR (ROL::StdConstraint) : solveAugmentedSystem not implemented!");
  return std::vector<Real>();
}

template<typename Real>
void StdConstraint<Real>::applyPreconditioner(Vector<Real> &pv,
                                        const Vector<Real> &v,
                                        const Vector<Real> &x,
                                        const Vector<Real> &g, Real &tol) {
  StdVector<Real> pvs = dynamic_cast<StdVector<Real>&>(pv);
  const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
  const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
  const StdVector<Real> gs = dynamic_cast<const StdVector<Real>&>(g);
  try {
    applyPreconditioner( *(pvs.getVector()), *(vs.getVector()), *(xs.getVector()),
                         *(gs.getVector()), tol );
  }
  catch (std::exception &e) {
    Constraint<Real>::applyPreconditioner(pv,v,x,g,tol);
  }
}

template<typename Real>
void StdConstraint<Real>::applyPreconditioner( std::vector<Real> &pv,
                                         const std::vector<Real> &v,
                                         const std::vector<Real> &x,
                                         const std::vector<Real> &g, Real &tol ) {
  ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument, 
    ">>> ERROR (ROL::StdConstraint) : applyPreconditioner not implemented!");
}

} // namespace ROL

#endif
