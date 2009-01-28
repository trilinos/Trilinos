//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#include "Rythmos_RKButcherTableau.hpp"

namespace Rythmos {

  Teuchos::Array<std::string> getS_RKButcherTableauMethodNames()
  {
    Teuchos::Array<std::string> a;
    a.push_back(RKBT_ForwardEuler_name());
    a.push_back(RKBT_BackwardEuler_name());
    a.push_back(Explicit4Stage_name());
    a.push_back(Explicit3_8Rule_name());
    a.push_back(Explicit2Stage2ndOrderRunge_name());
    a.push_back(Explicit3Stage3rdOrderHeun_name());
    a.push_back(Explicit3Stage3rdOrder_name());
    a.push_back(Explicit4Stage3rdOrderRunge_name());
    a.push_back(Implicit1Stage2ndOrderGauss_name());
    a.push_back(Implicit2Stage4thOrderGauss_name());
    a.push_back(Implicit3Stage6thOrderGauss_name());
    a.push_back(Implicit1Stage1stOrderRadauA_name());
    a.push_back(Implicit2Stage3rdOrderRadauA_name());
    a.push_back(Implicit3Stage5thOrderRadauA_name());
    a.push_back(Implicit1Stage1stOrderRadauB_name());
    a.push_back(Implicit2Stage3rdOrderRadauB_name());
    a.push_back(Implicit3Stage5thOrderRadauB_name());
    a.push_back(Implicit2Stage2ndOrderLobattoA_name());
    a.push_back(Implicit3Stage4thOrderLobattoA_name());
    a.push_back(Implicit4Stage6thOrderLobattoA_name());
    a.push_back(Implicit2Stage2ndOrderLobattoB_name());
    a.push_back(Implicit3Stage4thOrderLobattoB_name());
    a.push_back(Implicit4Stage6thOrderLobattoB_name());
    a.push_back(Implicit2Stage2ndOrderLobattoC_name());
    a.push_back(Implicit3Stage4thOrderLobattoC_name());
    a.push_back(Implicit4Stage6thOrderLobattoC_name());
    a.push_back(Implicit2Stage4thOrderHammerHollingsworth_name());
    a.push_back(Implicit3Stage6thOrderKuntzmannButcher_name());
    //a.push_back(Implicit4Stage8thOrderKuntzmannButcher_name()); // tscoffe 11/10/08 not passing convergence testing yet
    a.push_back(DIRK2Stage3rdOrder_name());
    a.push_back(SDIRK2Stage3rdOrder_name());
    a.push_back(SDIRK5Stage5thOrder_name());
    a.push_back(SDIRK5Stage4thOrder_name());
    a.push_back(SDIRK3Stage4thOrder_name());
    return a;
  }
   
} // namespace Rythmos

