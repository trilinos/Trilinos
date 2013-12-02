// @HEADER
// ************************************************************************
//
// Questions? Contact Denis Ridzal (dridzal@sandia.gov), or
//                    Drew Kouri   (dpkouri@sandia.gov)
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

#include "ROL_Types.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {

  template<class Real>
  void getTestObjectives( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x, 
                          const ETestObjectives test ) {
    switch (test) {
      case TESTOBJECTIVES_ROSENBROCK:          getRosenbrock(obj,x0,x);       break;
      case TESTOBJECTIVES_FREUDENSTEINANDROTH: getFreudensteinRoth(obj,x0,x); break;
      case TESTOBJECTIVES_BEALE:               getBeale(obj,x0,x);            break;
      case TESTOBJECTIVES_POWELL:              getPowell(obj,x0,x);           break;
      case TESTOBJECTIVES_SUMOFSQUARES:        getSumOfSquares(obj,x0,x);     break;
      case TESTOBJECTIVES_LEASTSQUARES:        getLeastSquares(obj,x0,x);     break; 
      case TESTOBJECTIVES_POISSONCONTROL:      getPoissonControl(obj,x0,x);   break;
      case TESTOBJECTIVES_LAST:                break;
    }
  }

} // namespace ROL

#endif
