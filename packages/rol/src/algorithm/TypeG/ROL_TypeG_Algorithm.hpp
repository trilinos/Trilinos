// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEG_ALGORITHM_H
#define ROL_TYPEG_ALGORITHM_H

#include "ROL_CombinedStatusTest.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_PolyhedralProjectionFactory.hpp"
#include "ROL_Problem.hpp"

/** \class ROL::TypeG::Algorithm
    \brief Provides an interface to run general constrained optimization algorithms.
*/

namespace ROL {
namespace TypeG {

template<typename Real>
struct AlgorithmState : public ROL::AlgorithmState<Real> {
  Real searchSize;
  Ptr<Vector<Real>> stepVec;
  Ptr<Vector<Real>> gradientVec;
  Ptr<Vector<Real>> constraintVec;

  AlgorithmState()
    : searchSize(1),
      stepVec(nullPtr),
      gradientVec(nullPtr),
      constraintVec(nullPtr) {}

  void reset() {
    ROL::AlgorithmState<Real>::reset();
    searchSize = static_cast<Real>(1);
    if (stepVec != nullPtr) {
      stepVec->zero();
    }
    if (gradientVec != nullPtr) {
      gradientVec->zero();
    }
    if (constraintVec != nullPtr) {
      constraintVec->zero();
    }
  }
};

template<typename Real>
class Algorithm {
protected:
  const Ptr<CombinedStatusTest<Real>> status_;
  const Ptr<AlgorithmState<Real>>     state_;
  Ptr<PolyhedralProjection<Real>>     proj_;

  void initialize( const Vector<Real> &x,
                   const Vector<Real> &g,
                   const Vector<Real> &mul,
                   const Vector<Real> &c);

public:

  virtual ~Algorithm() {}

  /** \brief Constructor, given a step and a status test.
  */
  Algorithm();

  void setStatusTest( const Ptr<StatusTest<Real>> &status,
                      bool combineStatus = false);

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Problem<Real> &problem,
                    std::ostream &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    std::ostream          &outStream = std::cout );


  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    const Vector<Real>    &eres,
                    std::ostream          &outStream = std::cout ) = 0;

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    const Vector<Real>    &ires,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    const Vector<Real>    &ires,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    const Vector<Real>    &eres,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    const Vector<Real>    &ires,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    const Vector<Real>    &eres,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    const Vector<Real>    &ires,
                    std::ostream          &outStream = std::cout );




  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    std::ostream          &outStream = std::cout );




  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    const Vector<Real>    &eres,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    const Vector<Real>    &linear_eres,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    const Vector<Real>    &ires,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    const Vector<Real>    &linear_eres,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    const Vector<Real>    &ires,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    const Vector<Real>    &linear_eres,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    const Vector<Real>    &eres,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    const Vector<Real>    &ires,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    const Vector<Real>    &linear_eres,
                    std::ostream          &outStream = std::cout );

  /** \brief Run algorithm on general constrained problems (Type-G).
             This is the primary Type-G interface.
  */
  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    BoundConstraint<Real> &bnd,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    const Vector<Real>    &eres,
                    Constraint<Real>      &icon,
                    Vector<Real>          &imul,
                    BoundConstraint<Real> &ibnd,
                    const Vector<Real>    &ires,
                    Constraint<Real>      &linear_econ,
                    Vector<Real>          &linear_emul,
                    const Vector<Real>    &linear_eres,
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

}; // class ROL::Algorithm

} // namespace TypeG
} // namespace ROL

#include "ROL_TypeG_Algorithm_Def.hpp"

#endif
