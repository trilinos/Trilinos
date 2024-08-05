// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HS_PROBLEMFACTORY_HPP
#define HS_PROBLEMFACTORY_HPP

#include "HS_Problem_001.hpp"
#include "HS_Problem_002.hpp"
#include "HS_Problem_003.hpp"
#include "HS_Problem_004.hpp"
#include "HS_Problem_005.hpp"
#include "HS_Problem_006.hpp"
#include "HS_Problem_007.hpp"
#include "HS_Problem_008.hpp"
#include "HS_Problem_009.hpp"
#include "HS_Problem_010.hpp"
#include "HS_Problem_011.hpp"
#include "HS_Problem_012.hpp"
#include "HS_Problem_013.hpp"
#include "HS_Problem_014.hpp"
#include "HS_Problem_015.hpp"
#include "HS_Problem_016.hpp"
#include "HS_Problem_017.hpp"
#include "HS_Problem_018.hpp"
#include "HS_Problem_019.hpp"
#include "HS_Problem_020.hpp"
#include "HS_Problem_021.hpp"
#include "HS_Problem_022.hpp"
#include "HS_Problem_023.hpp"
#include "HS_Problem_024.hpp"
#include "HS_Problem_025.hpp"
#include "HS_Problem_026.hpp"
#include "HS_Problem_027.hpp"
#include "HS_Problem_028.hpp"
#include "HS_Problem_029.hpp"
#include "HS_Problem_030.hpp"
#include "HS_Problem_031.hpp"
#include "HS_Problem_032.hpp"
#include "HS_Problem_033.hpp"
#include "HS_Problem_034.hpp"
#include "HS_Problem_035.hpp"
#include "HS_Problem_036.hpp"
#include "HS_Problem_037.hpp"
#include "HS_Problem_038.hpp"
#include "HS_Problem_039.hpp"
#include "HS_Problem_040.hpp"
#include "HS_Problem_041.hpp"

namespace HS {
template<class Real> 
class ProblemFactory {
public:
  ROL::Ptr<ROL::NonlinearProgram<Real> > getProblem(int n) {
    ROL::Ptr<ROL::NonlinearProgram<Real> > np;
    switch(n) {
      case   1: np = ROL::makePtr<Problem_001<Real>>(); break;
      case   2: np = ROL::makePtr<Problem_002<Real>>(); break;
      case   3: np = ROL::makePtr<Problem_003<Real>>(); break;
      case   4: np = ROL::makePtr<Problem_004<Real>>(); break;
      case   5: np = ROL::makePtr<Problem_005<Real>>(); break;
      case   6: np = ROL::makePtr<Problem_006<Real>>(); break;
      case   7: np = ROL::makePtr<Problem_007<Real>>(); break;
      case   8: np = ROL::makePtr<Problem_008<Real>>(); break;
      case   9: np = ROL::makePtr<Problem_009<Real>>(); break;
      case  10: np = ROL::makePtr<Problem_010<Real>>(); break;
      case  11: np = ROL::makePtr<Problem_011<Real>>(); break;
      case  12: np = ROL::makePtr<Problem_012<Real>>(); break;
      case  13: np = ROL::makePtr<Problem_013<Real>>(); break;
      case  14: np = ROL::makePtr<Problem_014<Real>>(); break;
      case  15: np = ROL::makePtr<Problem_015<Real>>(); break;
      case  16: np = ROL::makePtr<Problem_016<Real>>(); break;
      case  17: np = ROL::makePtr<Problem_017<Real>>(); break;
      case  18: np = ROL::makePtr<Problem_018<Real>>(); break;
      case  19: np = ROL::makePtr<Problem_019<Real>>(); break;
      case  20: np = ROL::makePtr<Problem_020<Real>>(); break;
      case  21: np = ROL::makePtr<Problem_021<Real>>(); break;
      case  22: np = ROL::makePtr<Problem_022<Real>>(); break;
      case  23: np = ROL::makePtr<Problem_023<Real>>(); break;
      case  24: np = ROL::makePtr<Problem_024<Real>>(); break;
      case  25: np = ROL::makePtr<Problem_025<Real>>(); break;
      case  26: np = ROL::makePtr<Problem_026<Real>>(); break;
      case  27: np = ROL::makePtr<Problem_027<Real>>(); break;
      case  28: np = ROL::makePtr<Problem_028<Real>>(); break;
      case  29: np = ROL::makePtr<Problem_029<Real>>(); break;
      case  30: np = ROL::makePtr<Problem_030<Real>>(); break;
      case  31: np = ROL::makePtr<Problem_031<Real>>(); break;
      case  32: np = ROL::makePtr<Problem_032<Real>>(); break;
      case  33: np = ROL::makePtr<Problem_033<Real>>(); break;
      case  34: np = ROL::makePtr<Problem_034<Real>>(); break;
      case  35: np = ROL::makePtr<Problem_035<Real>>(); break;
      case  36: np = ROL::makePtr<Problem_036<Real>>(); break;
      case  37: np = ROL::makePtr<Problem_037<Real>>(); break;
      case  38: np = ROL::makePtr<Problem_038<Real>>(); break;
      case  39: np = ROL::makePtr<Problem_039<Real>>(); break;
      case  40: np = ROL::makePtr<Problem_040<Real>>(); break;
      case  41: np = ROL::makePtr<Problem_041<Real>>(); break;


      default:
        ROL_TEST_FOR_EXCEPTION(true,std::logic_error,"Unknown problem number.");
      break;
    }
    return np;
  }
};
}
 
#endif // HS_PROBLEMFACTORY_HPP

