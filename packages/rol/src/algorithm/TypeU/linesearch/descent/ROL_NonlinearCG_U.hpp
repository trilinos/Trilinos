// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NONLINEARCG_U_H
#define ROL_NONLINEARCG_U_H

#include "ROL_DescentDirection_U.hpp"
#include "ROL_Types.hpp"
#include "ROL_NonlinearCG.hpp"

/** @ingroup step_group
    \class ROL::NonlinearCG_U
    \brief Provides the interface to compute optimization steps
           with nonlinear CG.
*/

#include "ROL_DescentDirection_U.hpp"

namespace ROL {

template<typename Real>
class NonlinearCG_U : public DescentDirection_U<Real> {
private:

  Ptr<NonlinearCG<Real>> nlcg_; ///< NonlinearCG object (used for quasi-Newton)
  ENonlinearCG enlcg_;
  std::string ncgName_;

public:

  /** \brief Constructor.

      Constructor to build a NonlinearCG object with a user-defined 
      nonlinear CG object.  Algorithmic specifications are passed in through 
      a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
      @param[in]     nlcg       is a user-defined NonlinearCG object
  */
  NonlinearCG_U( ParameterList &parlist,
                 const Ptr<NonlinearCG<Real>> &nlcg = nullPtr)
    : nlcg_(nlcg), enlcg_(NONLINEARCG_USERDEFINED) {
    // Initialize secant object
    ParameterList& Llist = parlist.sublist("Step").sublist("Line Search");
    if ( nlcg == nullPtr ) {
      ncgName_ = Llist.sublist("Descent Method").get("Nonlinear CG Type","Oren-Luenberger");
      enlcg_
        = StringToENonlinearCG(ncgName_);
      nlcg_ = ROL::makePtr<NonlinearCG<Real>>(enlcg_);
    }
    else {
      ncgName_ = Llist.sublist("Descent Method").get("User Defined Nonlinear CG Name",
                                                     "Unspecified User Define Nonlinear CG Method");
    }
  }

  void compute( Vector<Real> &s, Real &snorm, Real &sdotg, int &iter, int &flag,
          const Vector<Real> &x, const Vector<Real> &g, Objective<Real> &obj) override {
    nlcg_->run(s,g,x,obj);
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

  std::string printName(void) const override {
    std::stringstream name;
    name << ncgName_ << " Nonlinear CG";
    return name.str();
  }
}; // class ROL::NonlinearCG_U

} // namespace ROL

#endif
