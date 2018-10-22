// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_RKButcherTableauBuilder_impl_hpp
#define Tempus_RKButcherTableauBuilder_impl_hpp

#include "Tempus_RKButcherTableauBuilder_decl.hpp"
#include "Tempus_RKButcherTableau.hpp"

namespace Tempus {

// Nonmember constructor
template<class Scalar>
Teuchos::RCP<RKButcherTableauBuilder<Scalar> > rKButcherTableauBuilder()
{
  Teuchos::RCP<RKButcherTableauBuilder<Scalar> >
    rkbtfn = rcp(new RKButcherTableauBuilder<Scalar>() );
  return rkbtfn;
}
// Nonmember helper function
template<class Scalar>
Teuchos::RCP<RKButcherTableau<Scalar> > createRKBT(
  const std::string& rkbt_name, Teuchos::RCP<Teuchos::ParameterList> pl)
{
  Teuchos::RCP<RKButcherTableauBuilder<Scalar> >
    rkbtfn = rKButcherTableauBuilder<Scalar>();
  Teuchos::RCP<RKButcherTableau<Scalar> > rkbt = rkbtfn->create(rkbt_name);
  rkbt->setParameterList(pl);
  return rkbt;
}

template<class Scalar>
RKButcherTableauBuilder<Scalar>::RKButcherTableauBuilder()
{
  this->initializeDefaults_();
}

template<class Scalar>
void RKButcherTableauBuilder<Scalar>::setRKButcherTableauFactory(
    const Teuchos::RCP<const Teuchos::AbstractFactory<RKButcherTableau<Scalar> > >
      &rkbtFactory,
    const std::string &rkbtFactoryName
    )
{
  builder_.setObjectFactory(rkbtFactory, rkbtFactoryName);
}

template<class Scalar>
void RKButcherTableauBuilder<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  builder_.setParameterList(paramList);
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
RKButcherTableauBuilder<Scalar>::getNonconstParameterList()
{
  return builder_.getNonconstParameterList();
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
RKButcherTableauBuilder<Scalar>::unsetParameterList()
{
  return builder_.unsetParameterList();
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
RKButcherTableauBuilder<Scalar>::getParameterList() const
{
  return builder_.getParameterList();
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
RKButcherTableauBuilder<Scalar>::getValidParameters() const
{
  return builder_.getValidParameters();
}

template<class Scalar>
Teuchos::RCP<RKButcherTableau<Scalar> >
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

  builder_.setObjectName("Tempus::RKButcherTableau");
  builder_.setObjectTypeName("Runge Kutta Butcher Tableau Type");

  //
  // RK Butcher Tableaus:
  //

  // Explicit
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          GeneralExplicit_RKBT<Scalar> >(),
      "General ERK");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          ForwardEuler_RKBT<Scalar> >(),
      "RK Forward Euler");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Explicit2Stage2ndOrderRunge_RKBT<Scalar> >(),
      "RK Explicit 2 Stage 2nd order by Runge");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          ExplicitTrapezoidal_RKBT<Scalar> >(),
      "RK Explicit Trapezoidal");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Explicit3Stage3rdOrder_RKBT<Scalar> >(),
      "RK Explicit 3 Stage 3rd order");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Explicit3Stage3rdOrderHeun_RKBT<Scalar> >(),
      "RK Explicit 3 Stage 3rd order by Heun");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Explicit3Stage3rdOrderTVD_RKBT<Scalar> >(),
      "RK Explicit 3 Stage 3rd order TVD");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Explicit4Stage3rdOrderRunge_RKBT<Scalar> >(),
      "RK Explicit 4 Stage 3rd order by Runge");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Explicit5Stage3rdOrderKandG_RKBT<Scalar> >(),
      "RK Explicit 5 Stage 3rd order by Kinnmark and Gray");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Explicit4Stage4thOrder_RKBT<Scalar> >(),
      "RK Explicit 4 Stage");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Explicit3_8Rule_RKBT<Scalar> >(),
      "RK Explicit 3/8 Rule");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          ExplicitBogackiShampine32_RKBT<Scalar> >(),
      "Bogacki-Shampine 3(2) Pair");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          ExplicitMerson45_RKBT<Scalar> >(),
      "Merson 4(5) Pair");

  // Implicit
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          BackwardEuler_RKBT<Scalar> >(),
      "RK Backward Euler");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          IRK1StageTheta_RKBT<Scalar> >(),
      "IRK 1 Stage Theta Method");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          IRK1StageTheta_RKBT<Scalar> >(),
      "Implicit Midpoint");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          EDIRK2StageTheta_RKBT<Scalar> >(),
      "EDIRK 2 Stage Theta Method");

  // SDIRK
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          GeneralDIRK_RKBT<Scalar> >(),
      "General DIRK");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          SDIRK1Stage1stOrder_RKBT<Scalar> >(),
      "SDIRK 1 Stage 1st order");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          SDIRK2Stage2ndOrder_RKBT<Scalar> >(),
      "SDIRK 2 Stage 2nd order");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          SDIRK2Stage3rdOrder_RKBT<Scalar> >(),
      "SDIRK 2 Stage 3rd order");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          SDIRK3Stage4thOrder_RKBT<Scalar> >(),
      "SDIRK 3 Stage 4th order");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          SDIRK5Stage4thOrder_RKBT<Scalar> >(),
      "SDIRK 5 Stage 4th order");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          SDIRK5Stage5thOrder_RKBT<Scalar> >(),
      "SDIRK 5 Stage 5th order");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          SDIRK21_RKBT<Scalar> >(),
      "SDIRK 2(1) Pair");

  // EDIRK
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          EDIRK2Stage3rdOrder_RKBT<Scalar> >(),
      "EDIRK 2 Stage 3rd order");

  // IRK
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit1Stage2ndOrderGauss_RKBT<Scalar> >(),
      "RK Implicit 1 Stage 2nd order Gauss");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit2Stage4thOrderGauss_RKBT<Scalar> >(),
      "RK Implicit 2 Stage 4th order Gauss");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit3Stage6thOrderGauss_RKBT<Scalar> >(),
      "RK Implicit 3 Stage 6th order Gauss");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit2Stage4thOrderHammerHollingsworth_RKBT<Scalar> >(),
      "RK Implicit 2 Stage 4th Order Hammer & Hollingsworth");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit3Stage6thOrderKuntzmannButcher_RKBT<Scalar> >(),
      "RK Implicit 3 Stage 6th Order Kuntzmann & Butcher");

  //  This RKBT does not pass convergence testing, so we're disbaling it for now.
//  builder_.setObjectFactory(
//      abstractFactoryStd< RKButcherTableau<Scalar>, Implicit4Stage8thOrderKuntzmannButcher_RKBT<Scalar> >(),
//      "RK Implicit 4 Stage 8th Order Kuntzmann & Butcher");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit1Stage1stOrderRadauA_RKBT<Scalar> >(),
      "RK Implicit 1 Stage 1st order Radau left");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit2Stage3rdOrderRadauA_RKBT<Scalar> >(),
      "RK Implicit 2 Stage 3rd order Radau left");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit3Stage5thOrderRadauA_RKBT<Scalar> >(),
      "RK Implicit 3 Stage 5th order Radau left");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit1Stage1stOrderRadauB_RKBT<Scalar> >(),
      "RK Implicit 1 Stage 1st order Radau right");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit2Stage3rdOrderRadauB_RKBT<Scalar> >(),
      "RK Implicit 2 Stage 3rd order Radau right");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit3Stage5thOrderRadauB_RKBT<Scalar> >(),
      "RK Implicit 3 Stage 5th order Radau right");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit2Stage2ndOrderLobattoA_RKBT<Scalar> >(),
      "RK Implicit 2 Stage 2nd order Lobatto A");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit3Stage4thOrderLobattoA_RKBT<Scalar> >(),
      "RK Implicit 3 Stage 4th order Lobatto A");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit4Stage6thOrderLobattoA_RKBT<Scalar> >(),
      "RK Implicit 4 Stage 6th order Lobatto A");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit2Stage2ndOrderLobattoB_RKBT<Scalar> >(),
      "RK Implicit 2 Stage 2nd order Lobatto B");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit3Stage4thOrderLobattoB_RKBT<Scalar> >(),
      "RK Implicit 3 Stage 4th order Lobatto B");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit4Stage6thOrderLobattoB_RKBT<Scalar> >(),
      "RK Implicit 4 Stage 6th order Lobatto B");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit2Stage2ndOrderLobattoC_RKBT<Scalar> >(),
      "RK Implicit 2 Stage 2nd order Lobatto C");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit3Stage4thOrderLobattoC_RKBT<Scalar> >(),
      "RK Implicit 3 Stage 4th order Lobatto C");

  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          Implicit4Stage6thOrderLobattoC_RKBT<Scalar> >(),
      "RK Implicit 4 Stage 6th order Lobatto C");

  builder_.setDefaultObject("None");

}

} // namespace Tempus


#endif // Tempus_RKButcherTableauBuilder_impl_hpp
