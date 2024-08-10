// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_ALGORITHM_H
#define ROL_TYPEB_ALGORITHM_H

#include "ROL_CombinedStatusTest.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_PolyhedralProjectionFactory.hpp"
#include "ROL_Problem.hpp"

/** \class ROL::TypeB::Algorithm
    \brief Provides an interface to run bound constrained optimization algorithms.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
struct AlgorithmState : public ROL::AlgorithmState<Real> {
  Real searchSize;
  Ptr<Vector<Real>> stepVec;
  Ptr<Vector<Real>> gradientVec;
  int nproj;

  AlgorithmState()
    : searchSize(1),
      stepVec(nullPtr),
      gradientVec(nullPtr),
      nproj(0) {}

  void reset() {
    ROL::AlgorithmState<Real>::reset();
    searchSize = static_cast<Real>(1);
    if (stepVec != nullPtr) stepVec->zero();
    if (gradientVec != nullPtr) gradientVec->zero();
    nproj = 0;
  }
};

template<typename Real>
class Algorithm {
protected:
  const Ptr<CombinedStatusTest<Real>> status_;
  const Ptr<AlgorithmState<Real>>     state_;
  Ptr<PolyhedralProjection<Real>>     proj_;

  void initialize(const Vector<Real> &x, const Vector<Real> &g); 

  Real optimalityCriterion(const Vector<Real> &x, const Vector<Real> &g, Vector<Real> &primal,
                           std::ostream &outStream = std::cout) const;

public:

  virtual ~Algorithm() {}

  /** \brief Constructor, given a step and a status test.
  */
  Algorithm();

  void setStatusTest(const Ptr<StatusTest<Real>> &status,
                     const bool combineStatus = false);

  /** \brief Run algorithm on bound constrained problems (Type-B).
             This is the primary Type-B interface.
  */
  virtual void run( Problem<Real> &problem,
                    std::ostream  &outStream = std::cout );

  /** \brief Run algorithm on bound constrained problems (Type-B).
             This is the primary Type-B interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on bound constrained problems (Type-B).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    std::ostream          &outStream = std::cout) = 0;

  /** \brief Run algorithm on bound constrained problems with explicit
             linear constraints (Type-B).
	     This is the primary Type-B with explicit linear
             constraints interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on bound constrained problems with explicit
             linear constraints (Type-B).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    const Vector<Real>    &linear_eres,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on bound constrained problems with explicit
             linear constraints (Type-B).
	     This is the primary Type-B with explicit linear
             constraints interface.
  */
  virtual void run( Vector<Real>          &x,
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &bnd,
                                        Constraint<Real>      &linear_icon,
                                        Vector<Real>          &linear_imul,
                                        BoundConstraint<Real> &linear_ibnd,
                                        std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on bound constrained problems with explicit
             linear constraints (Type-B).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual void run( Vector<Real>          &x,
                                        const Vector<Real>    &g,
                                        Objective<Real>       &obj,
                                        BoundConstraint<Real> &bnd,
                                        Constraint<Real>      &linear_icon,
                                        Vector<Real>          &linear_imul,
                                        BoundConstraint<Real> &linear_ibnd,
                                        const Vector<Real>    &linear_ires,
                                        std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on bound constrained problems with explicit
             linear constraints (Type-B).
	     This is the primary Type-B with explicit linear
             constraints interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    Constraint<Real>      &linear_icon,
                    Vector<Real>          &linear_imul,
                    BoundConstraint<Real> &linear_ibnd,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on bound constrained problems with explicit
             linear constraints (Type-B).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    const Vector<Real>    &linear_eres,
                    Constraint<Real>      &linear_icon,
                    Vector<Real>          &linear_imul,
                    BoundConstraint<Real> &linear_ibnd,
                    const Vector<Real>    &linear_ires,
                    std::ostream          &outStream = std::cout );

  /** \brief Print iterate header.
  */
  virtual void writeHeader( std::ostream& os ) const;

  /** \brief Print step name.
  */
  virtual void writeName( std::ostream& os ) const;

  /** \brief Print iterate status.
  */
  virtual void writeOutput( std::ostream& os, const bool write_header = false ) const;

  virtual void writeExitStatus( std::ostream& os ) const;

  //Ptr<const AlgorithmState<Real>>& getState() const;
  Ptr<const AlgorithmState<Real>> getState() const;

  void reset();

}; // class ROL::Type::Algorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_Algorithm_Def.hpp"

#endif
