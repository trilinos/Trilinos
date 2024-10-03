// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_QUASINEWTONALGORITHM_HPP
#define ROL_TYPEB_QUASINEWTONALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"
#include "ROL_SecantFactory.hpp"

/** \class ROL::TypeB::QuasiNewtonAlgorithm
    \brief Provides an interface to run the projected secant algorithm.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class QuasiNewtonAlgorithm : public TypeB::Algorithm<Real> {
private:
  Ptr<Secant<Real>> secant_; ///< Secant object (used for quasi-Newton)
  ESecant esec_;             ///< Secant type
  std::string secantName_;   ///< Secant name

  int maxit_;         ///< Maximum number of line search steps (default: 20)
  Real rhodec_;       ///< Backtracking rate (default: 0.5)
  Real c1_;           ///< Sufficient Decrease Parameter (default: 1e-4)
  Real sigma1_;       ///< Lower safeguard for quadratic line search (default: 0.1)
  Real sigma2_;       ///< Upper safeguard for quadratic line search (default: 0.9)
  Real sp_tol1_;
  Real sp_tol2_;
  Real sp_tol_min_;
  std::string algoName_;

  ParameterList list_;

  bool hasLEC_;
  int ls_nfval_, spgIter_;
  int verbosity_;
  bool writeHeader_;

  using TypeB::Algorithm<Real>::status_;
  using TypeB::Algorithm<Real>::state_;
  using TypeB::Algorithm<Real>::proj_;

  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &obj,
                  BoundConstraint<Real> &bnd,
                  std::ostream &outStream = std::cout); 

public:

  QuasiNewtonAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using TypeB::Algorithm<Real>::run;

  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream &os, const bool write_header = false ) const override;

}; // class ROL::TypeB::QuasiNewtonAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_QuasiNewtonAlgorithm_Def.hpp"

#endif
