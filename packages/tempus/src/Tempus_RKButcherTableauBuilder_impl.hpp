#ifndef Tempus_RKButcherTableauBuilder_impl_hpp
#define Tempus_RKButcherTableauBuilder_impl_hpp

#include "Tempus_RKButcherTableauBuilder_decl.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Tempus_String_Utilities.hpp"

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
Teuchos::RCP<RKButcherTableau<Scalar> > createRKBT(const std::string& rkbt_name)
{
  Teuchos::RCP<RKButcherTableauBuilder<Scalar> >
    rkbtfn = rKButcherTableauBuilder<Scalar>();
  Teuchos::RCP<RKButcherTableau<Scalar> > rkbt = rkbtfn->create(rkbt_name);
  return rkbt;
}

template<class Scalar>
Teuchos::RCP<RKButcherTableau<Scalar> > parseRKTableau(const Teuchos::ParameterList & pl)
{
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::string error_msg = "The format of the Butcher Tableau parameter list is\n"
                          "  <Parameter name=\"A\" type=\"string\" value=\"# # # ;\n"
                          "                                                # # # ;\n"
                          "                                                # # #\"/>\n"
                          "  <Parameter name=\"b\" type=\"string\" value=\"# # #\"/>\n"
                          "  <Parameter name=\"c\" type=\"string\" value=\"# # #\"/>\n\n"
                          "Note the number of stages is implicit in the number of entries.\n" 
                          "The number of stages must be consistent.\n";

  std::string embedded_error_msg = "The format of the Butcher Tableau parameter list is\n"
                          "  <Parameter name=\"A\" type=\"string\" value=\"# # # ;\n"
                          "                                                # # # ;\n"
                          "                                                # # #\"/>\n"
                          "  <Parameter name=\"b\" type=\"string\" value=\"# # #\"/>\n"
                          "  <Parameter name=\"bhat\" type=\"string\" value=\"# # #\"/>\n"
                          "  <Parameter name=\"c\" type=\"string\" value=\"# # #\"/>\n\n"
                          "Note the number of stages is implicit in the number of entries.\n" 
                          "The number of stages must be consistent.\n";

  // do some error checking
    bool has_A = pl.isType<std::string>("A");
    bool has_b = pl.isType<std::string>("b");
    bool has_c = pl.isType<std::string>("c");
    bool has_order = pl.isType<int>("order");

    TEUCHOS_TEST_FOR_EXCEPTION(!has_A,std::runtime_error,error_msg);
    TEUCHOS_TEST_FOR_EXCEPTION(!has_b,std::runtime_error,error_msg);
    TEUCHOS_TEST_FOR_EXCEPTION(!has_c,std::runtime_error,error_msg);

  std::size_t numStages = -1;
  int order =  has_order ? pl.get<int>("order") : 1;
  Teuchos::SerialDenseMatrix<int,double> A_tbl;
  Teuchos::SerialDenseVector<int,double> b_tbl;
  Teuchos::SerialDenseVector<int,double> bhat_tbl;
  Teuchos::SerialDenseVector<int,double> c_tbl;


  // read in the A matrix
  {
    std::vector<std::string> A_row_tokens; 
    Tempus::StringTokenizer(A_row_tokens,pl.get<std::string>("A"),";",true);

    // this is the only place where numStages is set
    numStages = A_row_tokens.size();

    // allocate the matrix
    A_tbl.shape(as<int>(numStages),as<int>(numStages));

    // fill the rows
    for(std::size_t r=0;r<numStages;r++) {
      // parse the row (tokenize on space)
      std::vector<std::string> tokens;
      Tempus::StringTokenizer(tokens,A_row_tokens[r]," ",true);

      std::vector<double> values;
      Tempus::TokensToDoubles(values,tokens);

      TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
                                 "Error parsing A matrix, wrong number of stages in row " << r << "\n" + error_msg);

      for(std::size_t c=0;c<numStages;c++)
        A_tbl(r,c) = values[c]; 
    }
  }

  // size b and c vectors
  b_tbl.size(as<int>(numStages));
  c_tbl.size(as<int>(numStages));


  // read in the b vector
  {
    std::vector<std::string> tokens; 
    Tempus::StringTokenizer(tokens,pl.get<std::string>("b")," ",true);
    std::vector<double> values;
    Tempus::TokensToDoubles(values,tokens);

    TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
                               "Error parsing b vector, wrong number of stages.\n" + error_msg);

    for(std::size_t i=0;i<numStages;i++)
      b_tbl(i) = values[i]; 
  }

  // read in the c vector
  {
    std::vector<std::string> tokens; 
    Tempus::StringTokenizer(tokens,pl.get<std::string>("c")," ",true);
    std::vector<double> values;
    Tempus::TokensToDoubles(values,tokens);

    TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
                               "Error parsing c vector, wrong number of stages.\n" + error_msg);

    for(std::size_t i=0;i<numStages;i++)
      c_tbl(i) = values[i]; 
  }


  RCP<RKButcherTableau<double> > tableau = rcp(new RKButcherTableau<double>);

  tableau->initialize(A_tbl,b_tbl,c_tbl,order, "Built from Tempus::parseRKTableau");

  return tableau;
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
                          IRK2StageTheta_RKBT<Scalar> >(),
      "IRK 2 Stage Theta Method");

  // SDIRK
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

  // DIRK
  builder_.setObjectFactory(
      abstractFactoryStd< RKButcherTableau<Scalar>,
                          DIRK2Stage3rdOrder_RKBT<Scalar> >(),
      "Diagonal IRK 2 Stage 3rd order");

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
