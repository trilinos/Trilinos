// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_QUASINEWTON_U_H
#define ROL_QUASINEWTON_U_H

#include "ROL_DescentDirection_U.hpp"

#include "ROL_Types.hpp"
#include "ROL_SecantFactory.hpp"

/** @ingroup step_group
    \class ROL::QuasiNewton_U
    \brief Provides the interface to compute optimization steps
           with a secant method.
*/

namespace ROL {

template<typename Real>
class QuasiNewton_U : public DescentDirection_U<Real> {
private:

  Ptr<Secant<Real>> secant_; ///< Secant object (used for quasi-Newton)
  ESecant esec_;             ///< Secant type
  std::string secantName_;   ///< Secant name

public:

  /** \brief Constructor.

      Constructor to build a QuasiNewton object with a user-defined 
      secant object.  Algorithmic specifications are passed in through 
      a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
      @param[in]     secant     is a user-defined secant object
  */
  QuasiNewton_U( ParameterList &parlist,
           const Ptr<Secant<Real>> &secant = nullPtr)
    : secant_(secant), esec_(SECANT_USERDEFINED) {
    // Initialize secant object
    if ( secant == nullPtr ) {
      secantName_ = parlist.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS");
      esec_ = StringToESecant(secantName_);
      secant_ = SecantFactory<Real>(parlist);
    }
    else {
      secantName_ = parlist.sublist("General").sublist("Secant").get("User Defined Secant Name",
                                                                     "Unspecified User Defined Secant Method"); 
    }
  }

  void compute( Vector<Real> &s, Real &snorm, Real &sdotg, int &iter, int &flag,
          const Vector<Real> &x, const Vector<Real> &g, Objective<Real> &obj) override {
    secant_->applyH(s,g);
    //sdotg = -s.dot(g.dual());
    sdotg = -s.apply(g);
    if (sdotg >= static_cast<Real>(0)) {
      s.set(g.dual());
      //sdotg = -s.dot(g.dual());
      sdotg = -s.apply(g);
    }
    s.scale(static_cast<Real>(-1));
    snorm = s.norm();
    iter  = 0;
    flag  = 0;
  }

  void update(const Vector<Real> &x, const Vector<Real> &s,
              const Vector<Real> &gold, const Vector<Real> &gnew,
              const Real snorm, const int iter) override {
    // Update Secant Information
    secant_->updateStorage(x,gnew,gold,s,snorm,iter+1);
  }

  std::string printName(void) const override {
    std::stringstream name;
    name << "Quasi-Newton Method with " << secantName_;
    return name.str();
  }
}; // class ROL::QuasiNewton_U

} // namespace ROL

#endif
