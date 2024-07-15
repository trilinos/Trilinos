// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEP_QUASINEWTONALGORITHM_HPP
#define ROL_TYPEP_QUASINEWTONALGORITHM_HPP

#include "ROL_TypeP_Algorithm.hpp"
#include "ROL_SecantFactory.hpp"

/** \class ROL::TypeP::QuasiNewtonAlgorithm
    \brief Provides an interface to run the projected secant algorithm.
*/

namespace ROL {
namespace TypeP {

template<typename Real>
class QuasiNewtonAlgorithm : public TypeP::Algorithm<Real> {
private:
  Ptr<Secant<Real>> secant_; ///< Secant object (used for quasi-Newton)
  ESecant esec_;             ///< Secant type
  std::string secantName_;   ///< Secant name

  int t0_;
  bool initProx_;

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

  using TypeP::Algorithm<Real>::pgstep;
  using TypeP::Algorithm<Real>::status_;
  using TypeP::Algorithm<Real>::state_;

  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &sobj,
                  Objective<Real>       &nobj,
                  Vector<Real>          &dg,
                  std::ostream &outStream = std::cout); 

public:

  QuasiNewtonAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using TypeP::Algorithm<Real>::run;

  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &sobj,
            Objective<Real>       &nobj,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream &os, bool write_header = false ) const override;

}; // class ROL::TypeP::QuasiNewtonAlgorithm

} // namespace TypeP
} // namespace ROL

#include "ROL_TypeP_QuasiNewtonAlgorithm_Def.hpp"

#endif
