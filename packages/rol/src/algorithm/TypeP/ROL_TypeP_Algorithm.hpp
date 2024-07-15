// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEP_ALGORITHM_H
#define ROL_TYPEP_ALGORITHM_H

#include "ROL_CombinedStatusTest.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Problem.hpp"

/** \class ROL::TypeP::Algorithm
    \brief Provides an interface to run optimization algorithms to minimize composite optimization problems f+phi.
*/

namespace ROL {
namespace TypeP { 

template<typename Real>
struct AlgorithmState : public ROL::AlgorithmState<Real> {
  Real searchSize, svalue, nvalue;
  Ptr<Vector<Real>> stepVec, gradientVec;
  int nprox, nsval, nnval;

  AlgorithmState()
    : searchSize(1),
      svalue(ROL_INF<Real>()),
      nvalue(ROL_INF<Real>()),
      stepVec(nullPtr),
      gradientVec(nullPtr),
      nprox(0),
      nsval(0),
      nnval(0) {}

  void reset() {
    ROL::AlgorithmState<Real>::reset();
    searchSize = static_cast<Real>(1);
    svalue     = ROL_INF<Real>(); 
    nvalue     = ROL_INF<Real>(); 
    if (stepVec != nullPtr)
      stepVec->zero();
    if (gradientVec != nullPtr)
      gradientVec->zero();
    nprox = 0;
    nsval = 0; 
    nnval = 0; 
  }
};

template<typename Real>
class Algorithm {
protected:
  const Ptr<CombinedStatusTest<Real>> status_;
  const Ptr<AlgorithmState<Real>>     state_;

  void initialize(const Vector<Real> &x, const Vector<Real> &g); 
  void pgstep(Vector<Real>       &pgiter, // pgiter = Prox(x - t dg)
              Vector<Real>       &pgstep, // pgstep = pgiter - x
              Objective<Real>    &nobj,   // nobj   = nonsmooth objective
              const Vector<Real> &x,      // x      = current iterate
              const Vector<Real> &dg,     // dg     = dual of current gradient
              Real                t,      // t      = prox parameter
              Real               &tol) const;

public:

  virtual ~Algorithm() {}

  /** \brief Constructor, given a step and a status test.
  */
  Algorithm();

  void setStatusTest(const Ptr<StatusTest<Real>> &status,
                     bool combineStatus = false);

  /** \brief Run algorithm on unconstrained problems (Type-U).
             This is the primary Type-U interface.
  */
  virtual void run( Problem<Real> &problem,
                    std::ostream  &outStream = std::cout );

  /** \brief Run algorithm on unconstrained problems (Type-U).
             This is the primary Type-U interface.
  */
  virtual void run( Vector<Real>    &x,
                    Objective<Real> &sobj,
                    Objective<Real> &nobj,
                    std::ostream    &outStream = std::cout );

   /** \brief Run algorithm on unconstrained problems (Type-U).
             This general interface supports the use of dual optimization vector spaces,
             where the user does not define the dual() method.
  */
  virtual void run( Vector<Real>       &x,
                    const Vector<Real> &g, 
                    Objective<Real>    &sobj,
                    Objective<Real>    &nobj,
                    std::ostream       &outStream = std::cout) = 0;

  /** \brief Print iterate header.
  */
  virtual void writeHeader( std::ostream& os ) const;

  /** \brief Print step name.
  */
  virtual void writeName( std::ostream& os ) const;

  /** \brief Print iterate status.
  */
  virtual void writeOutput( std::ostream& os, bool write_header = false ) const;

  virtual void writeExitStatus( std::ostream& os ) const;

  //Ptr<const AlgorithmState<Real>>& getState() const;
  Ptr<const AlgorithmState<Real>> getState() const;

  void reset();

}; // class ROL::TypeP::Algorithm
} // namespace TypeP
} // namespace ROL

#include "ROL_TypeP_Algorithm_Def.hpp"

#endif
