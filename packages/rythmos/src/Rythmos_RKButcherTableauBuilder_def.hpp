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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_RK_BUTCHER_TABLEAU_BUILDER_DEF_HPP
#define RYTHMOS_RK_BUTCHER_TABLEAU_BUILDER_DEF_HPP

#include "Rythmos_Types.hpp"

#include "Rythmos_RKButcherTableauBuilder_decl.hpp"
#include "Rythmos_RKButcherTableau.hpp"

namespace Rythmos {

// Nonmember constructor
template<class Scalar>
RCP<RKButcherTableauBuilder<Scalar> > rKButcherTableauBuilder()
{
  RCP<RKButcherTableauBuilder<Scalar> > rkbtfn = rcp(new RKButcherTableauBuilder<Scalar>() );
  return rkbtfn;
}
// Nonmember helper function
template<class Scalar>
RCP<RKButcherTableauBase<Scalar> > createRKBT(const std::string& rkbt_name)
{
  RCP<RKButcherTableauBuilder<Scalar> > rkbtfn = rKButcherTableauBuilder<Scalar>();
  RCP<RKButcherTableauBase<Scalar> > rkbt = rkbtfn->create(rkbt_name);
  return rkbt;
}

template<class Scalar>
RKButcherTableauBuilder<Scalar>::RKButcherTableauBuilder()
{
  this->initializeDefaults_();
}

template<class Scalar>
void RKButcherTableauBuilder<Scalar>::setRKButcherTableauFactory(
    const RCP<const Teuchos::AbstractFactory<RKButcherTableauBase<Scalar> > > &rkbtFactory,
    const std::string &rkbtFactoryName
    )
{
  builder_.setObjectFactory(rkbtFactory, rkbtFactoryName);
}

template<class Scalar>
void RKButcherTableauBuilder<Scalar>::setParameterList(
  RCP<Teuchos::ParameterList> const& paramList
  )
{
  builder_.setParameterList(paramList);
}

template<class Scalar>
RCP<Teuchos::ParameterList>
RKButcherTableauBuilder<Scalar>::getNonconstParameterList()
{
  return builder_.getNonconstParameterList();
}


template<class Scalar>
RCP<Teuchos::ParameterList>
RKButcherTableauBuilder<Scalar>::unsetParameterList()
{
  return builder_.unsetParameterList();
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
RKButcherTableauBuilder<Scalar>::getParameterList() const
{
  return builder_.getParameterList();
}


template<class Scalar>
RCP<const Teuchos::ParameterList>
RKButcherTableauBuilder<Scalar>::getValidParameters() const
{
  return builder_.getValidParameters();
}

template<class Scalar>
RCP<RKButcherTableauBase<Scalar> >
RKButcherTableauBuilder<Scalar>::create(
    const std::string &rkbt_name
    ) const
{
  return builder_.create(rkbt_name);
}

template<class Scalar>
void RKButcherTableauBuilder<Scalar>::initializeDefaults_()
{

  using Teuchos::abstractFactoryStd;

  builder_.setObjectName("Rythmos::RKButcherTableau");
  builder_.setObjectTypeName("Runge Kutta Butcher Tableau Type");

  //
  // RK Butcher Tableaus:
  //

  // Explicit
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          ForwardEuler_RKBT<Scalar> >(),
      RKBT_ForwardEuler_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Explicit2Stage2ndOrderRunge_RKBT<Scalar> >(),
      Explicit2Stage2ndOrderRunge_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          ExplicitTrapezoidal_RKBT<Scalar> >(),
      ExplicitTrapezoidal_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Explicit3Stage3rdOrder_RKBT<Scalar> >(),
      Explicit3Stage3rdOrder_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Explicit3Stage3rdOrderHeun_RKBT<Scalar> >(),
      Explicit3Stage3rdOrderHeun_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Explicit3Stage3rdOrderTVD_RKBT<Scalar> >(),
      Explicit3Stage3rdOrderTVD_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Explicit4Stage3rdOrderRunge_RKBT<Scalar> >(),
      Explicit4Stage3rdOrderRunge_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Explicit4Stage4thOrder_RKBT<Scalar> >(),
      Explicit4Stage_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Explicit3_8Rule_RKBT<Scalar> >(),
      Explicit3_8Rule_name());

  // Implicit
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          BackwardEuler_RKBT<Scalar> >(),
      RKBT_BackwardEuler_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          IRK1StageTheta_RKBT<Scalar> >(),
      IRK1StageTheta_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          IRK2StageTheta_RKBT<Scalar> >(),
      IRK2StageTheta_name());

  // SDIRK
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          SDIRK2Stage2ndOrder_RKBT<Scalar> >(),
      SDIRK2Stage2ndOrder_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          SDIRK2Stage3rdOrder_RKBT<Scalar> >(),
      SDIRK2Stage3rdOrder_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          SDIRK3Stage4thOrder_RKBT<Scalar> >(),
      SDIRK3Stage4thOrder_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          SDIRK5Stage4thOrder_RKBT<Scalar> >(),
      SDIRK5Stage4thOrder_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          SDIRK5Stage5thOrder_RKBT<Scalar> >(),
      SDIRK5Stage5thOrder_name());

  // DIRK
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          DIRK2Stage3rdOrder_RKBT<Scalar> >(),
      DIRK2Stage3rdOrder_name());

  // IRK
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit1Stage2ndOrderGauss_RKBT<Scalar> >(),
      Implicit1Stage2ndOrderGauss_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit2Stage4thOrderGauss_RKBT<Scalar> >(),
      Implicit2Stage4thOrderGauss_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit3Stage6thOrderGauss_RKBT<Scalar> >(),
      Implicit3Stage6thOrderGauss_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit2Stage4thOrderHammerHollingsworth_RKBT<Scalar> >(),
      Implicit2Stage4thOrderHammerHollingsworth_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit3Stage6thOrderKuntzmannButcher_RKBT<Scalar> >(),
      Implicit3Stage6thOrderKuntzmannButcher_name());

  //  This RKBT does not pass convergence testing, so we're disbaling it for now.
//  builder_.setObjectFactory(
//      abstractFactoryStd< RKButcherTableauBase<Scalar>, Implicit4Stage8thOrderKuntzmannButcher_RKBT<Scalar> >(),
//      Implicit4Stage8thOrderKuntzmannButcher_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit1Stage1stOrderRadauA_RKBT<Scalar> >(),
      Implicit1Stage1stOrderRadauA_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit2Stage3rdOrderRadauA_RKBT<Scalar> >(),
      Implicit2Stage3rdOrderRadauA_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit3Stage5thOrderRadauA_RKBT<Scalar> >(),
      Implicit3Stage5thOrderRadauA_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit1Stage1stOrderRadauB_RKBT<Scalar> >(),
      Implicit1Stage1stOrderRadauB_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit2Stage3rdOrderRadauB_RKBT<Scalar> >(),
      Implicit2Stage3rdOrderRadauB_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit3Stage5thOrderRadauB_RKBT<Scalar> >(),
      Implicit3Stage5thOrderRadauB_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit2Stage2ndOrderLobattoA_RKBT<Scalar> >(),
      Implicit2Stage2ndOrderLobattoA_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit3Stage4thOrderLobattoA_RKBT<Scalar> >(),
      Implicit3Stage4thOrderLobattoA_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit4Stage6thOrderLobattoA_RKBT<Scalar> >(),
      Implicit4Stage6thOrderLobattoA_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit2Stage2ndOrderLobattoB_RKBT<Scalar> >(),
      Implicit2Stage2ndOrderLobattoB_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit3Stage4thOrderLobattoB_RKBT<Scalar> >(),
      Implicit3Stage4thOrderLobattoB_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit4Stage6thOrderLobattoB_RKBT<Scalar> >(),
      Implicit4Stage6thOrderLobattoB_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit2Stage2ndOrderLobattoC_RKBT<Scalar> >(),
      Implicit2Stage2ndOrderLobattoC_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit3Stage4thOrderLobattoC_RKBT<Scalar> >(),
      Implicit3Stage4thOrderLobattoC_name());

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableauBase<Scalar>,
                          Implicit4Stage6thOrderLobattoC_RKBT<Scalar> >(),
      Implicit4Stage6thOrderLobattoC_name());

  builder_.setDefaultObject("None");

}

//
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_RK_BUTCHER_TABLEAU_BUILDER_INSTANT(SCALAR) \
  \
  template class RKButcherTableauBuilder< SCALAR >; \
  \
  template RCP<RKButcherTableauBuilder< SCALAR > > rKButcherTableauBuilder(); \
  \
  template RCP<RKButcherTableauBase< SCALAR > > createRKBT(const std::string& rkbt_name);


} // namespace Rythmos


#endif // RYTHMOS_RK_BUTCHER_TABLEAU_BUILDER_DEF_HPP
