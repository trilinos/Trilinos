// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "ROL_HS9.hpp"
#include "ROL_HS14.hpp"
#include "ROL_HS21.hpp"
#include "ROL_HS24.hpp"
#include "ROL_HS25.hpp"
#include "ROL_HS28.hpp"
#include "ROL_HS29.hpp"
#include "ROL_HS32.hpp"
#include "ROL_HS38.hpp"
#include "ROL_HS39.hpp"
#include "ROL_HS41.hpp"
#include "ROL_HS42.hpp"
#include "ROL_HS45.hpp"
#include "ROL_HS48.hpp"
#include "ROL_HS49.hpp"
#include "ROL_HS50.hpp"
#include "ROL_HS51.hpp"
#include "ROL_HS52.hpp"
#include "ROL_HS53.hpp"
#include "ROL_HS55.hpp"
#include "ROL_HS63.hpp"
#include "ROL_BVP.hpp"
#include "ROL_ParaboloidCircle.hpp"
#include "ROL_SimpleEqConstrained.hpp"
#include "ROL_CantileverBeam.hpp"
#include "ROL_Cubic.hpp"
#include "ROL_Quartic.hpp"
#include "ROL_CylinderHead.hpp"
#include "ROL_Cantilever.hpp"
#include "ROL_Minimax1.hpp"
#include "ROL_Minimax2.hpp"
#include "ROL_Minimax3.hpp"


namespace ROL {

  /** \enum   ROL::ETestOptProblem
      \brief  Enumeration of test optimization problems.

      \arg    ROSENBROCK                describe
      \arg    FREUDENSTEINANDROTH       describe
      \arg    BEALE                     describe
      \arg    POWELL                    describe
      \arg    SUMOFSQUARES              describe
      \arg    LEASTSQUARES              describe
      \arg    POISSONCONTROL            describe
      \arg    POISSONINVERSION          describe
      \arg    ZAKHAROV                  describe
      \arg    HS1                       describe
      \arg    HS2                       describe
      \arg    HS3                       describe
      \arg    HS4                       describe
      \arg    HS5                       describe
      \arg    HS9                       describe
      \arg    HS14                      describe
      \arg    HS21                      describe
      \arg    HS24                      describe
      \arg    HS25                      describe
      \arg    HS28                      describe
      \arg    HS29                      describe
      \arg    HS32                      describe
      \arg    HS38                      describe
      \arg    HS39                      describe
      \arg    HS41                      describe
      \arg    HS42                      describe
      \arg    HS45                      describe
      \arg    HS48                      describe
      \arg    HS49                      describe
      \arg    HS50                      describe
      \arg    HS51                      describe
      \arg    HS52                      describe
      \arg    HS53                      describe
      \arg    HS55                      describe
      \arg    HS63                      describe
      \arg    BVP                       describe
      \arg    PARABOLOIDCIRCLE          describe
      \arg    SIMPLEEQCONSTRAINED       describe
      \arg    CANTILEVERBEAM            describe
      \arg    CUBIC                     describe
      \arg    QUARTIC                   describe
      \arg    CYLINDERHEAD              describe
      \arg    CANTILEVER                describe
      \arg    MINIMAX1                  describe
      \arg    MINIMAX2                  describe
      \arg    MINIMAX3                  describe
   */
  enum ETestOptProblem {
    TESTOPTPROBLEM_ROSENBROCK = 0,
    TESTOPTPROBLEM_FREUDENSTEINANDROTH,
    TESTOPTPROBLEM_BEALE,
    TESTOPTPROBLEM_POWELL,
    TESTOPTPROBLEM_SUMOFSQUARES,
    TESTOPTPROBLEM_LEASTSQUARES,
    TESTOPTPROBLEM_POISSONCONTROL,
    TESTOPTPROBLEM_POISSONINVERSION,
    TESTOPTPROBLEM_ZAKHAROV,
    TESTOPTPROBLEM_HS1,
    TESTOPTPROBLEM_HS2,
    TESTOPTPROBLEM_HS3,
    TESTOPTPROBLEM_HS4,
    TESTOPTPROBLEM_HS5,
    TESTOPTPROBLEM_HS9,
    TESTOPTPROBLEM_HS14,
    TESTOPTPROBLEM_HS21,
    TESTOPTPROBLEM_HS24,
    TESTOPTPROBLEM_HS25,
    TESTOPTPROBLEM_HS28,
    TESTOPTPROBLEM_HS29,
    TESTOPTPROBLEM_HS32,
    TESTOPTPROBLEM_HS38,
    TESTOPTPROBLEM_HS39,
    TESTOPTPROBLEM_HS41,
    TESTOPTPROBLEM_HS42,
    TESTOPTPROBLEM_HS45,
    TESTOPTPROBLEM_HS48,
    TESTOPTPROBLEM_HS49,
    TESTOPTPROBLEM_HS50,
    TESTOPTPROBLEM_HS51,
    TESTOPTPROBLEM_HS52,
    TESTOPTPROBLEM_HS53,
    TESTOPTPROBLEM_HS55,
    TESTOPTPROBLEM_HS63,
    TESTOPTPROBLEM_BVP,
    TESTOPTPROBLEM_PARABOLOIDCIRCLE,
    TESTOPTPROBLEM_SIMPLEEQCONSTRAINED,
    TESTOPTPROBLEM_CANTILEVERBEAM,
    TESTOPTPROBLEM_CUBIC,
    TESTOPTPROBLEM_QUARTIC,
    TESTOPTPROBLEM_CYLINDERHEAD,
    TESTOPTPROBLEM_CANTILEVER,
    TESTOPTPROBLEM_MINIMAX1,
    TESTOPTPROBLEM_MINIMAX2,
    TESTOPTPROBLEM_MINIMAX3,
    TESTOPTPROBLEM_LAST
  };

  inline std::string ETestOptProblemToString(ETestOptProblem to) {
    std::string retString;
    switch(to) {
      case TESTOPTPROBLEM_ROSENBROCK:          retString = "Rosenbrock's Function";                  break;
      case TESTOPTPROBLEM_FREUDENSTEINANDROTH: retString = "Freudenstein and Roth's Function";       break;
      case TESTOPTPROBLEM_BEALE:               retString = "Beale's Function";                       break;
      case TESTOPTPROBLEM_POWELL:              retString = "Powell's Badly Scaled Function";         break;
      case TESTOPTPROBLEM_SUMOFSQUARES:        retString = "Sum of Squares Function";                break;
      case TESTOPTPROBLEM_LEASTSQUARES:        retString = "Least Squares Function";                 break;
      case TESTOPTPROBLEM_POISSONCONTROL:      retString = "Poisson Optimal Control";                break;
      case TESTOPTPROBLEM_POISSONINVERSION:    retString = "Poisson Inversion Problem";              break;
      case TESTOPTPROBLEM_ZAKHAROV:            retString = "Zakharov's Function";                    break;
      case TESTOPTPROBLEM_HS1:                 retString = "Hock and Schittkowski Test Problem #1";  break;
      case TESTOPTPROBLEM_HS2:                 retString = "Hock and Schittkowski Test Problem #2";  break;
      case TESTOPTPROBLEM_HS3:                 retString = "Hock and Schittkowski Test Problem #3";  break;
      case TESTOPTPROBLEM_HS4:                 retString = "Hock and Schittkowski Test Problem #4";  break;
      case TESTOPTPROBLEM_HS5:                 retString = "Hock and Schittkowski Test Problem #5";  break;
      case TESTOPTPROBLEM_HS9:                 retString = "Hock and Schittkowski Test Problem #9";  break;
      case TESTOPTPROBLEM_HS14:                retString = "Hock and Schittkowski Test Problem #14"; break;
      case TESTOPTPROBLEM_HS21:                retString = "Hock and Schittkowski Test Problem #21"; break;
      case TESTOPTPROBLEM_HS24:                retString = "Hock and Schittkowski Test Problem #24"; break;
      case TESTOPTPROBLEM_HS25:                retString = "Hock and Schittkowski Test Problem #25"; break;
      case TESTOPTPROBLEM_HS28:                retString = "Hock and Schittkowski Test Problem #28"; break;
      case TESTOPTPROBLEM_HS29:                retString = "Hock and Schittkowski Test Problem #29"; break;
      case TESTOPTPROBLEM_HS32:                retString = "Hock and Schittkowski Test Problem #32"; break;
      case TESTOPTPROBLEM_HS38:                retString = "Hock and Schittkowski Test Problem #38"; break;
      case TESTOPTPROBLEM_HS39:                retString = "Hock and Schittkowski Test Problem #39"; break;
      case TESTOPTPROBLEM_HS41:                retString = "Hock and Schittkowski Test Problem #41"; break;
      case TESTOPTPROBLEM_HS42:                retString = "Hock and Schittkowski Test Problem #42"; break;
      case TESTOPTPROBLEM_HS45:                retString = "Hock and Schittkowski Test Problem #45"; break;
      case TESTOPTPROBLEM_HS48:                retString = "Hock and Schittkowski Test Problem #48"; break;
      case TESTOPTPROBLEM_HS49:                retString = "Hock and Schittkowski Test Problem #49"; break;
      case TESTOPTPROBLEM_HS50:                retString = "Hock and Schittkowski Test Problem #50"; break;
      case TESTOPTPROBLEM_HS51:                retString = "Hock and Schittkowski Test Problem #51"; break;
      case TESTOPTPROBLEM_HS52:                retString = "Hock and Schittkowski Test Problem #52"; break;
      case TESTOPTPROBLEM_HS53:                retString = "Hock and Schittkowski Test Problem #53"; break;
      case TESTOPTPROBLEM_HS55:                retString = "Hock and Schittkowski Test Problem #55"; break;
      case TESTOPTPROBLEM_HS63:                retString = "Hock and Schittkowski Test Problem #63"; break;
      case TESTOPTPROBLEM_BVP:                 retString = "Boundary Value Problem";                 break;
      case TESTOPTPROBLEM_PARABOLOIDCIRCLE:    retString = "Paraboloid Circle";                      break;
      case TESTOPTPROBLEM_SIMPLEEQCONSTRAINED: retString = "Simple Equality Constrained";            break;
      case TESTOPTPROBLEM_CANTILEVERBEAM:      retString = "Cantilever Beam";                        break;
      case TESTOPTPROBLEM_CUBIC:               retString = "Cubic";                                  break;
      case TESTOPTPROBLEM_QUARTIC:             retString = "Quartic";                                break;
      case TESTOPTPROBLEM_CYLINDERHEAD:        retString = "Cylinder Head";                          break;
      case TESTOPTPROBLEM_CANTILEVER:          retString = "Cantilever";                             break;
      case TESTOPTPROBLEM_MINIMAX1:            retString = "Minimax #1";                             break;
      case TESTOPTPROBLEM_MINIMAX2:            retString = "Minimax #2";                             break;
      case TESTOPTPROBLEM_MINIMAX3:            retString = "Minimax #3";                             break;
      case TESTOPTPROBLEM_LAST:                retString = "Last Type (Dummy)";                      break;
      default:                                 retString = "INVALID ETestOptProblem";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a TestOptProblem enum.
    
      \param  to  [in]  - enum of the TestOptProblem
      \return 1 if the argument is a valid TestOptProblem; 0 otherwise.
    */
  inline int isValidTestOptProblem(ETestOptProblem to){
    return( (to == TESTOPTPROBLEM_ROSENBROCK)          ||
            (to == TESTOPTPROBLEM_FREUDENSTEINANDROTH) ||
            (to == TESTOPTPROBLEM_BEALE)               ||
            (to == TESTOPTPROBLEM_POWELL)              ||
            (to == TESTOPTPROBLEM_SUMOFSQUARES)        ||
            (to == TESTOPTPROBLEM_LEASTSQUARES)        ||
            (to == TESTOPTPROBLEM_POISSONCONTROL)      ||
            (to == TESTOPTPROBLEM_POISSONINVERSION)    ||
            (to == TESTOPTPROBLEM_ZAKHAROV)            ||
            (to == TESTOPTPROBLEM_HS1)                 ||
            (to == TESTOPTPROBLEM_HS2)                 ||
            (to == TESTOPTPROBLEM_HS3)                 ||
            (to == TESTOPTPROBLEM_HS4)                 ||
            (to == TESTOPTPROBLEM_HS5)                 ||
            (to == TESTOPTPROBLEM_HS9)                 ||
            (to == TESTOPTPROBLEM_HS14)                ||
            (to == TESTOPTPROBLEM_HS21)                ||
            (to == TESTOPTPROBLEM_HS24)                ||
            (to == TESTOPTPROBLEM_HS25)                ||
            (to == TESTOPTPROBLEM_HS28)                ||
            (to == TESTOPTPROBLEM_HS29)                ||
            (to == TESTOPTPROBLEM_HS32)                ||
            (to == TESTOPTPROBLEM_HS38)                ||
            (to == TESTOPTPROBLEM_HS39)                ||
	    (to == TESTOPTPROBLEM_HS41)                ||
	    (to == TESTOPTPROBLEM_HS42)                ||
            (to == TESTOPTPROBLEM_HS45)                ||
            (to == TESTOPTPROBLEM_HS48)                ||
            (to == TESTOPTPROBLEM_HS49)                ||
            (to == TESTOPTPROBLEM_HS50)                ||
            (to == TESTOPTPROBLEM_HS51)                ||
            (to == TESTOPTPROBLEM_HS52)                ||
            (to == TESTOPTPROBLEM_HS53)                ||
            (to == TESTOPTPROBLEM_HS55)                ||
            (to == TESTOPTPROBLEM_HS63)                ||
            (to == TESTOPTPROBLEM_BVP)                 ||
            (to == TESTOPTPROBLEM_PARABOLOIDCIRCLE)    ||
            (to == TESTOPTPROBLEM_SIMPLEEQCONSTRAINED) ||
            (to == TESTOPTPROBLEM_CANTILEVERBEAM)      ||
            (to == TESTOPTPROBLEM_CUBIC)               ||
            (to == TESTOPTPROBLEM_QUARTIC)             ||
            (to == TESTOPTPROBLEM_CYLINDERHEAD)        ||
            (to == TESTOPTPROBLEM_CANTILEVER)          ||
            (to == TESTOPTPROBLEM_MINIMAX1)            ||
            (to == TESTOPTPROBLEM_MINIMAX2)            ||
            (to == TESTOPTPROBLEM_MINIMAX3) );
  }

  inline ETestOptProblem & operator++(ETestOptProblem &type) {
    return type = static_cast<ETestOptProblem>(type+1);
  }

  inline ETestOptProblem operator++(ETestOptProblem &type, int) {
    ETestOptProblem oldval = type;
    ++type;
    return oldval;
  }

  inline ETestOptProblem & operator--(ETestOptProblem &type) {
    return type = static_cast<ETestOptProblem>(type-1);
  }

  inline ETestOptProblem operator--(ETestOptProblem &type, int) {
    ETestOptProblem oldval = type;
    --type;
    return oldval;
  }

  inline ETestOptProblem StringToETestOptProblem(std::string s) {
    s = removeStringFormat(s);
    for ( ETestOptProblem to = TESTOPTPROBLEM_ROSENBROCK; to < TESTOPTPROBLEM_LAST; to++ ) {
      if ( !s.compare(removeStringFormat(ETestOptProblemToString(to))) ) {
        return to;
      }
    }
    return TESTOPTPROBLEM_ROSENBROCK;
  }

  template<class Real>
  void GetTestProblem( Ptr<OptimizationProblem<Real> > &problem,
                       Ptr<Vector<Real> > &x0,
                       std::vector<Ptr<Vector<Real> > > &x, 
                       const ETestOptProblem test ) {
    Ptr<TestProblem<Real>> tp;
    switch (test) {
      case TESTOPTPROBLEM_ROSENBROCK:          tp = makePtr<ZOO::getRosenbrock<Real>>();          break;
      case TESTOPTPROBLEM_FREUDENSTEINANDROTH: tp = makePtr<ZOO::getFreudensteinRoth<Real>>();    break;
      case TESTOPTPROBLEM_BEALE:               tp = makePtr<ZOO::getBeale<Real>>();               break;
      case TESTOPTPROBLEM_POWELL:              tp = makePtr<ZOO::getPowell<Real>>();              break;
      case TESTOPTPROBLEM_SUMOFSQUARES:        tp = makePtr<ZOO::getSumOfSquares<Real>>();        break;
      case TESTOPTPROBLEM_LEASTSQUARES:        tp = makePtr<ZOO::getLeastSquares<Real>>();        break; 
      case TESTOPTPROBLEM_POISSONCONTROL:      tp = makePtr<ZOO::getPoissonControl<Real>>();      break;
      case TESTOPTPROBLEM_POISSONINVERSION:    tp = makePtr<ZOO::getPoissonInversion<Real>>();    break;
      case TESTOPTPROBLEM_ZAKHAROV:            tp = makePtr<ZOO::getZakharov<Real>>();            break; 
      case TESTOPTPROBLEM_HS1:                 tp = makePtr<ZOO::getHS1<Real>>();                 break;
      case TESTOPTPROBLEM_HS2:                 tp = makePtr<ZOO::getHS2<Real>>();                 break;
      case TESTOPTPROBLEM_HS3:                 tp = makePtr<ZOO::getHS3<Real>>();                 break;
      case TESTOPTPROBLEM_HS4:                 tp = makePtr<ZOO::getHS4<Real>>();                 break;
      case TESTOPTPROBLEM_HS5:                 tp = makePtr<ZOO::getHS5<Real>>();                 break;
      case TESTOPTPROBLEM_HS9:                 tp = makePtr<ZOO::getHS9<Real>>();                 break;
      case TESTOPTPROBLEM_HS14:                tp = makePtr<ZOO::getHS14<Real>>();                break;
      case TESTOPTPROBLEM_HS21:                tp = makePtr<ZOO::getHS21<Real>>();                break;
      case TESTOPTPROBLEM_HS24:                tp = makePtr<ZOO::getHS24<Real>>();                break;
      case TESTOPTPROBLEM_HS25:                tp = makePtr<ZOO::getHS25<Real>>();                break;
      case TESTOPTPROBLEM_HS28:                tp = makePtr<ZOO::getHS28<Real>>();                break;
      case TESTOPTPROBLEM_HS29:                tp = makePtr<ZOO::getHS29<Real>>();                break;
      case TESTOPTPROBLEM_HS32:                tp = makePtr<ZOO::getHS32<Real>>();                break;
      case TESTOPTPROBLEM_HS38:                tp = makePtr<ZOO::getHS38<Real>>();                break;
      case TESTOPTPROBLEM_HS39:                tp = makePtr<ZOO::getHS39<Real>>();                break;
      case TESTOPTPROBLEM_HS41:                tp = makePtr<ZOO::getHS41<Real>>();                break;
      case TESTOPTPROBLEM_HS42:                tp = makePtr<ZOO::getHS42<Real>>();                break;
      case TESTOPTPROBLEM_HS45:                tp = makePtr<ZOO::getHS45<Real>>();                break;
      case TESTOPTPROBLEM_HS48:                tp = makePtr<ZOO::getHS48<Real>>();                break;
      case TESTOPTPROBLEM_HS49:                tp = makePtr<ZOO::getHS49<Real>>();                break;
      case TESTOPTPROBLEM_HS50:                tp = makePtr<ZOO::getHS50<Real>>();                break;
      case TESTOPTPROBLEM_HS51:                tp = makePtr<ZOO::getHS51<Real>>();                break;
      case TESTOPTPROBLEM_HS52:                tp = makePtr<ZOO::getHS52<Real>>();                break;
      case TESTOPTPROBLEM_HS53:                tp = makePtr<ZOO::getHS53<Real>>();                break;
      case TESTOPTPROBLEM_HS55:                tp = makePtr<ZOO::getHS55<Real>>();                break;
      case TESTOPTPROBLEM_HS63:                tp = makePtr<ZOO::getHS63<Real>>();                break;
      case TESTOPTPROBLEM_BVP:                 tp = makePtr<ZOO::getBVP<Real>>();                 break;
      case TESTOPTPROBLEM_PARABOLOIDCIRCLE:    tp = makePtr<ZOO::getParaboloidCircle<Real>>();    break;
      case TESTOPTPROBLEM_SIMPLEEQCONSTRAINED: tp = makePtr<ZOO::getSimpleEqConstrained<Real>>(); break;
      case TESTOPTPROBLEM_CANTILEVERBEAM:      tp = makePtr<ZOO::getCantileverBeam<Real>>();      break;
      case TESTOPTPROBLEM_CUBIC:               tp = makePtr<ZOO::getCubic<Real>>(2);              break;
      case TESTOPTPROBLEM_QUARTIC:             tp = makePtr<ZOO::getQuartic<Real>>();             break;
      case TESTOPTPROBLEM_CYLINDERHEAD:        tp = makePtr<ZOO::getCylinderHead<Real>>();        break;
      case TESTOPTPROBLEM_CANTILEVER:          tp = makePtr<ZOO::getCantilever<Real>>();          break;
      case TESTOPTPROBLEM_MINIMAX1:            tp = makePtr<ZOO::getMinimax1<Real>>();            break;
      case TESTOPTPROBLEM_MINIMAX2:            tp = makePtr<ZOO::getMinimax2<Real>>();            break;
      case TESTOPTPROBLEM_MINIMAX3:            tp = makePtr<ZOO::getMinimax3<Real>>();            break;
      case TESTOPTPROBLEM_LAST:                                                                   break;
    }
    tp->get(problem,x0,x);
  }
} // namespace ROL

#endif
