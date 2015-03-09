// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions of test objective functions.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_TESTOBJECTIVES_HPP
#define ROL_TESTOBJECTIVES_HPP

#include "ROL_Rosenbrock.hpp"
#include "ROL_FreudensteinRoth.hpp"
#include "ROL_Beale.hpp"
#include "ROL_Powell.hpp"
#include "ROL_SumOfSquares.hpp"
#include "ROL_LeastSquares.hpp"
#include "ROL_PoissonControl.hpp"
#include "ROL_PoissonInversion.hpp"
#include "ROL_Zakharov.hpp"

#include "ROL_HS1.hpp"
#include "ROL_HS2.hpp"
#include "ROL_HS3.hpp"
#include "ROL_HS4.hpp"
#include "ROL_HS5.hpp"
#include "ROL_HS25.hpp"
#include "ROL_HS38.hpp"
#include "ROL_HS45.hpp"
#include "ROL_BVP.hpp"

#include "ROL_Types.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL {

  template<class Real>
  void getTestObjectives( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x, 
                          const ETestObjectives test ) {
    switch (test) {
      case TESTOBJECTIVES_ROSENBROCK:          ZOO::getRosenbrock<Real,StdVector<Real>,StdVector<Real> > (obj,x0,x);       break;
      case TESTOBJECTIVES_FREUDENSTEINANDROTH: ZOO::getFreudensteinRoth(obj,x0,x); break;
      case TESTOBJECTIVES_BEALE:               ZOO::getBeale(obj,x0,x);            break;
      case TESTOBJECTIVES_POWELL:              ZOO::getPowell(obj,x0,x);           break;
      case TESTOBJECTIVES_SUMOFSQUARES:        ZOO::getSumOfSquares(obj,x0,x);     break;
      case TESTOBJECTIVES_LEASTSQUARES:        ZOO::getLeastSquares(obj,x0,x);     break; 
      case TESTOBJECTIVES_POISSONCONTROL:      ZOO::getPoissonControl(obj,x0,x);   break;
      case TESTOBJECTIVES_POISSONINVERSION:    ZOO::getPoissonInversion(obj,x0,x); break;
      case TESTOBJECTIVES_ZAKHAROV:            ZOO::getZakharov(obj,x0,x);         break; 
      case TESTOBJECTIVES_LAST:                break;
    }
  }

  template<class Real>
  void getTestObjectives( Teuchos::RCP<Objective<Real> > &obj, Teuchos::RCP<BoundConstraint<Real> > &con, 
                          Vector<Real> &x0, Vector<Real> &x, 
                          const ETestOptProblem test ) {
    switch (test) {
      case TESTOPTPROBLEM_HS1:  ZOO::getHS1(obj,con,x0,x);  break;
      case TESTOPTPROBLEM_HS2:  ZOO::getHS2(obj,con,x0,x);  break;
      case TESTOPTPROBLEM_HS3:  ZOO::getHS3(obj,con,x0,x);  break;
      case TESTOPTPROBLEM_HS4:  ZOO::getHS4(obj,con,x0,x);  break;
      case TESTOPTPROBLEM_HS5:  ZOO::getHS5(obj,con,x0,x);  break;
      case TESTOPTPROBLEM_HS25: ZOO::getHS25(obj,con,x0,x); break;
      case TESTOPTPROBLEM_HS38: ZOO::getHS38(obj,con,x0,x); break;
      case TESTOPTPROBLEM_HS45: ZOO::getHS45(obj,con,x0,x); break;
      case TESTOPTPROBLEM_BVP:  ZOO::getBVP(obj,con,x0,x);  break;
      case TESTOPTPROBLEM_LAST: break;
    }
  }
} // namespace ROL

#endif
