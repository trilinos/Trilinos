// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperFactory_hpp
#define Tempus_StepperFactory_hpp

#include "Teuchos_ParameterList.hpp"
#include "Tempus_StepperForwardEuler.hpp"
#include "Tempus_StepperBackwardEuler.hpp"
#include "Tempus_StepperExplicitRK.hpp"
#include "Tempus_StepperBDF2.hpp"
#include "Tempus_StepperNewmarkImplicitAForm.hpp"
#include "Tempus_StepperNewmarkImplicitDForm.hpp"
#include "Tempus_StepperNewmarkExplicitAForm.hpp"
#include "Tempus_StepperHHTAlpha.hpp"
#include "Tempus_StepperRKButcherTableau.hpp"
#include "Tempus_StepperIMEX_RK.hpp"
#include "Tempus_StepperIMEX_RK_Partition.hpp"
#include "Tempus_StepperLeapfrog.hpp"
#include "Tempus_StepperOperatorSplit.hpp"
#include "Tempus_StepperTrapezoidal.hpp"
#include "Tempus_StepperSubcycling.hpp"

#include "NOX_Thyra.H"

namespace Tempus {

/** \brief Stepper factory.
 *
 * <b>Adding Steppers</b>
 *    -#
 */
template<class Scalar>
class StepperFactory
{
public:

  /// Constructor
  StepperFactory(){}

  /// Destructor
  virtual ~StepperFactory() {}

  /// Create default stepper from stepper type (e.g., "Forward Euler").
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    std::string stepperType = "Forward Euler",
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
      model = Teuchos::null)
  {
    if (stepperType == "") stepperType = "Forward Euler";
    return this->createStepper(model, stepperType, Teuchos::null);
  }

  /// Create stepper from ParameterList with its details.
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
      model = Teuchos::null)
  {
    std::string stepperType = "Forward Euler";
    if (stepperPL != Teuchos::null)
      stepperType = stepperPL->get<std::string>("Stepper Type","Forward Euler");
    return this->createStepper(model, stepperType, stepperPL);
  }

  /// Create stepper from ParameterList with its details.
  Teuchos::RCP<Stepper<Scalar> > createMultiSteppers(
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models)
  {
    std::string stepperType = "Operator Split";
    if (stepperPL != Teuchos::null)
      stepperType = stepperPL->get<std::string>("Stepper Type");
    return this->createStepper(models, stepperType, stepperPL);
  }


  // ---------------------------------------------------------------------------

  /// Set Stepper member data from the ParameterList.
  void setStepperValues(
    Teuchos::RCP<Stepper<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    if (stepperPL != Teuchos::null) {
      auto stepperType =
        stepperPL->get<std::string>("Stepper Type", stepper->getStepperType());
      TEUCHOS_TEST_FOR_EXCEPTION(
        stepperType != stepper->getStepperType() ,std::runtime_error,
        "  ParameterList 'Stepper Type' (='" + stepperType +"')\n"
        "  does not match type for stepper Stepper (='"
        + stepper->getStepperType() + "').");
      stepper->setStepperType(stepperType);

      stepper->setUseFSAL(
        stepperPL->get<bool>("Use FSAL", stepper->getUseFSALDefault()));

      stepper->setICConsistency(
        stepperPL->get<std::string>("Initial Condition Consistency",
                                    stepper->getICConsistencyDefault()));

      stepper->setICConsistencyCheck(
        stepperPL->get<bool>("Initial Condition Consistency Check",
                             stepper->getICConsistencyCheckDefault()));
    }
  }


  /// Create a tableau from the ParameterList.
  Teuchos::RCP<RKButcherTableau<Scalar> > createTableau(
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp_const_cast;
    using Teuchos::ParameterList;

    TEUCHOS_TEST_FOR_EXCEPTION(stepperPL == Teuchos::null,std::runtime_error,
      "Error parsing general tableau.  ParameterList is null.\n");

    Teuchos::RCP<RKButcherTableau<Scalar> > tableau;

    RCP<ParameterList> tableauPL = sublist(stepperPL,"Tableau",true);
    std::size_t numStages = 0;
    int order = tableauPL->get<int>("order");
    Teuchos::SerialDenseMatrix<int,Scalar> A;
    Teuchos::SerialDenseVector<int,Scalar> b;
    Teuchos::SerialDenseVector<int,Scalar> c;
    Teuchos::SerialDenseVector<int,Scalar> bstar;

    // read in the A matrix
    {
      std::vector<std::string> A_row_tokens;
      Tempus::StringTokenizer(A_row_tokens, tableauPL->get<std::string>("A"),
                              ";",true);

      // this is the only place where numStages is set
      numStages = A_row_tokens.size();

      // allocate the matrix
      A.shape(as<int>(numStages),as<int>(numStages));

      // fill the rows
      for(std::size_t row=0;row<numStages;row++) {
        // parse the row (tokenize on space)
        std::vector<std::string> tokens;
        Tempus::StringTokenizer(tokens,A_row_tokens[row]," ",true);

        std::vector<double> values;
        Tempus::TokensToDoubles(values,tokens);

        TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
          "Error parsing A matrix, wrong number of stages in row "
          << row << ".\n");

        for(std::size_t col=0;col<numStages;col++)
          A(row,col) = values[col];
      }
    }

    // size b and c vectors
    b.size(as<int>(numStages));
    c.size(as<int>(numStages));

    // read in the b vector
    {
      std::vector<std::string> tokens;
      Tempus::StringTokenizer(tokens,tableauPL->get<std::string>("b")," ",true);
      std::vector<double> values;
      Tempus::TokensToDoubles(values,tokens);

      TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
        "Error parsing b vector, wrong number of stages.\n");

      for(std::size_t i=0;i<numStages;i++)
        b(i) = values[i];
    }

    // read in the c vector
    {
      std::vector<std::string> tokens;
      Tempus::StringTokenizer(tokens,tableauPL->get<std::string>("c")," ",true);
      std::vector<double> values;
      Tempus::TokensToDoubles(values,tokens);

      TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
        "Error parsing c vector, wrong number of stages.\n");

      for(std::size_t i=0;i<numStages;i++)
        c(i) = values[i];
    }

    if (tableauPL->isParameter("bstar") and
        tableauPL->get<std::string>("bstar") != "") {
      bstar.size(as<int>(numStages));
      // read in the bstar vector
      {
        std::vector<std::string> tokens;
        Tempus::StringTokenizer(
          tokens, tableauPL->get<std::string>("bstar"), " ", true);
        std::vector<double> values;
        Tempus::TokensToDoubles(values,tokens);

        TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
          "Error parsing bstar vector, wrong number of stages.\n"
          "      Number of RK stages    = " << numStages << "\n"
          "      Number of bstar values = " << values.size() << "\n");

        for(std::size_t i=0;i<numStages;i++)
          bstar(i) = values[i];
      }
      tableau = rcp(new RKButcherTableau<Scalar>(
        "From ParameterList",A,b,c,order,order,order,false,-4.0,bstar));
    } else {
      tableau = rcp(new RKButcherTableau<Scalar>(
        "From ParameterList",A,b,c,order,order,order,false,-4.0));
    }
    return tableau;
  }


  /// Set StepperExplicit member data from the ParameterList.
  void setStepperExplicitValues(
    Teuchos::RCP<StepperExplicit<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    if (stepperPL != Teuchos::null) {
      stepperPL->validateParametersAndSetDefaults(
                                              *stepper->getValidParameters());
      setStepperValues(stepper, stepperPL);
    }
  }

  /// Set StepperImplicit member data from the ParameterList.
  void setStepperImplicitValues(
    Teuchos::RCP<StepperImplicit<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto solver = rcp(new Thyra::NOXNonlinearSolver());
    auto solverPL = Tempus::defaultSolverParameters();
    auto subPL = sublist(solverPL, "NOX");
    solver->setParameterList(subPL);
    if (stepperPL != Teuchos::null) {
      // Can not validate because of optional Parameters, e.g., 'Solver Name'.
      //stepperPL->validateParametersAndSetDefaults(
      //                                      *stepper->getValidParameters());
      setStepperValues(stepper, stepperPL);
      stepper->setZeroInitialGuess(
        stepperPL->get<bool>("Zero Initial Guess", false));
    }
  }

  /// Set StepperRK member data from the ParameterList.
  void setStepperRKValues(
    Teuchos::RCP<StepperExplicitRK<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    if (stepperPL != Teuchos::null) {
      stepperPL->validateParametersAndSetDefaults(
                                              *stepper->getValidParameters());
      setStepperValues(stepper, stepperPL);
      stepper->setUseEmbedded(
        stepperPL->get<bool>("Use Embedded",stepper->getUseEmbeddedDefault()));
    }
  }

  /// Set solver from ParameterList.
  void setStepperSolverValues(
    Teuchos::RCP<StepperImplicit<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto solver = rcp(new Thyra::NOXNonlinearSolver());
    auto defaultSolverPL = Tempus::defaultSolverParameters();
    auto subPL = sublist(defaultSolverPL, "NOX");
    solver->setParameterList(subPL);
    if (stepperPL != Teuchos::null) {
      std::string solverName = stepperPL->get<std::string>("Solver Name");
      if ( stepperPL->isSublist(solverName) ) {
        auto solverPL = Teuchos::parameterList();
        solverPL = Teuchos::sublist(stepperPL, solverName);
        Teuchos::RCP<Teuchos::ParameterList> noxPL =
          Teuchos::sublist(solverPL,"NOX",true);
        solver->setParameterList(noxPL);
      }
    }
    stepper->setSolver(solver);
  }

  /// Set StepperDIRK member data from the ParameterList.
  void setStepperDIRKValues(
    Teuchos::RCP<StepperDIRK<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto solver = rcp(new Thyra::NOXNonlinearSolver());
    auto solverPL = Tempus::defaultSolverParameters();
    auto subPL = sublist(solverPL, "NOX");
    solver->setParameterList(subPL);
    if (stepperPL != Teuchos::null) {
      // Can not validate because of optional Parameters, e.g., 'Solver Name'.
      //stepperPL->validateParametersAndSetDefaults(
      //                                      *stepper->getValidParameters());
      setStepperValues(stepper, stepperPL);
      stepper->setUseEmbedded(
        stepperPL->get<bool>("Use Embedded",stepper->getUseEmbeddedDefault()));
      stepper->setZeroInitialGuess(
        stepperPL->get<bool>("Zero Initial Guess", false));
    }
  }

  void setTableausPartition(
    Teuchos::RCP<StepperIMEX_RK_Partition<Scalar> > stepper,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::string stepperType)
  {
    using Teuchos::RCP;
    if (stepperType == "") {
      if (stepperPL == Teuchos::null)
        stepperType = "Partitioned IMEX RK SSP2";
      else
        stepperType = stepperPL->get<std::string>("Stepper Type",
                                                  "Partitioned IMEX RK SSP2");
    }

    if (stepperType != "General Partitioned IMEX RK") {
      stepper->setTableaus(stepperType);
    } else {
      if (stepperPL != Teuchos::null) {
        Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau;
        Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau;
        if (stepperPL->isSublist("IMEX-RK Explicit Stepper")) {
          RCP<Teuchos::ParameterList> explicitPL = Teuchos::rcp(
              new Teuchos::ParameterList(
                  stepperPL->sublist("IMEX-RK Explicit Stepper")));
          auto stepperTemp = createStepper(explicitPL, Teuchos::null);
          auto stepperERK = Teuchos::rcp_dynamic_cast<StepperExplicitRK<Scalar> > (
                          stepperTemp,true);
          TEUCHOS_TEST_FOR_EXCEPTION(stepperERK == Teuchos::null, std::logic_error,
              "Error - The explicit component of a general partitioned IMEX RK stepper was not specified as an ExplicitRK stepper");
          explicitTableau = stepperERK->getTableau();
        }

        if (stepperPL->isSublist("IMEX-RK Implicit Stepper")) {
          RCP<Teuchos::ParameterList> implicitPL = Teuchos::rcp(
              new Teuchos::ParameterList(
                  stepperPL->sublist("IMEX-RK Implicit Stepper")));
          auto stepperTemp = createStepper(implicitPL, Teuchos::null);
          auto stepperDIRK = Teuchos::rcp_dynamic_cast<StepperDIRK<Scalar> > (
              stepperTemp,true);
          TEUCHOS_TEST_FOR_EXCEPTION(stepperDIRK == Teuchos::null, std::logic_error,
              "Error - The implicit component of a general partitioned IMEX RK stepper was not specified as an DIRK stepper");
          implicitTableau = stepperDIRK->getTableau();
        }

        TEUCHOS_TEST_FOR_EXCEPTION(
          !(explicitTableau!=Teuchos::null && implicitTableau!=Teuchos::null), std::logic_error,
          "Error - A parameter list was used to setup a general partitioned IMEX RK stepper, but did not "
          "specify both an explicit and an implicit tableau!\n");

        stepper->setTableaus(stepperType, explicitTableau, implicitTableau);

        stepper->setOrder(stepperPL->get<int>("overall order", 1));
      }
    }
  }


  void setTableaus(Teuchos::RCP<StepperIMEX_RK<Scalar> > stepper,
                   Teuchos::RCP<Teuchos::ParameterList> stepperPL,
                   std::string stepperType)
  {
    using Teuchos::RCP;
    if (stepperType == "") {
      if (stepperPL == Teuchos::null)
        stepperType = "IMEX RK SSP2";
      else
        stepperType = stepperPL->get<std::string>("Stepper Type",
                                                  "IMEX RK SSP2");
    }

    if (stepperType != "General IMEX RK") {
      stepper->setTableaus(stepperType);
    } else {
      if (stepperPL != Teuchos::null) {
        Teuchos::RCP<const RKButcherTableau<Scalar> > explicitTableau;
        Teuchos::RCP<const RKButcherTableau<Scalar> > implicitTableau;
        if (stepperPL->isSublist("IMEX-RK Explicit Stepper")) {
          RCP<Teuchos::ParameterList> explicitPL = Teuchos::rcp(
              new Teuchos::ParameterList(
                  stepperPL->sublist("IMEX-RK Explicit Stepper")));
          auto stepperTemp = createStepper(explicitPL, Teuchos::null);
          auto stepperERK = Teuchos::rcp_dynamic_cast<StepperExplicitRK<Scalar> > (
                          stepperTemp,true);
          TEUCHOS_TEST_FOR_EXCEPTION(stepperERK == Teuchos::null, std::logic_error,
              "Error - The explicit component of a general IMEX RK stepper was not specified as an ExplicitRK stepper");
          explicitTableau = stepperERK->getTableau();
        }

        if (stepperPL->isSublist("IMEX-RK Implicit Stepper")) {
          RCP<Teuchos::ParameterList> implicitPL = Teuchos::rcp(
              new Teuchos::ParameterList(
                  stepperPL->sublist("IMEX-RK Implicit Stepper")));
          auto stepperTemp = createStepper(implicitPL, Teuchos::null);
          auto stepperDIRK = Teuchos::rcp_dynamic_cast<StepperDIRK<Scalar> > (
              stepperTemp,true);
          TEUCHOS_TEST_FOR_EXCEPTION(stepperDIRK == Teuchos::null, std::logic_error,
              "Error - The implicit component of a general IMEX RK stepper was not specified as an DIRK stepper");
          implicitTableau = stepperDIRK->getTableau();
        }

        TEUCHOS_TEST_FOR_EXCEPTION(
          !(explicitTableau!=Teuchos::null && implicitTableau!=Teuchos::null), std::logic_error,
          "Error - A parameter list was used to setup a general IMEX RK stepper, but did not "
          "specify both an explicit and an implicit tableau!\n");

        stepper->setTableaus(stepperType, explicitTableau, implicitTableau);
        stepper->setOrder(stepperPL->get<int>("overall order", 1));
      }
    }
  }


  // ---------------------------------------------------------------------------
  // Create individual Steppers.

  Teuchos::RCP<StepperSubcycling<Scalar> >
  createStepperSubcycling(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSubcycling<Scalar>());

    TEUCHOS_TEST_FOR_EXCEPTION(stepperPL != Teuchos::null, std::logic_error,
      "Error - Construction of StepperSubcycling with a ParameterList\n"
      "is not implemented yet!\n");

    if (stepperPL != Teuchos::null) {
      setStepperValues(stepper, stepperPL);
      //setStepperSubcyclingValues(stepper, stepperPL);
    }
    //else {
    //  integrator->setTempusParameterList(Teuchos::null);
    //}

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperIMEX_RK_Partition<Scalar> >
  createStepperIMEX_RK_Partition(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperIMEX_RK_Partition<Scalar>());
    stepper->setStepperType(stepperType);
    setStepperImplicitValues(stepper, stepperPL);
    setTableausPartition(stepper, stepperPL, stepperType);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperIMEX_RK<Scalar> >
  createStepperIMEX_RK(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperIMEX_RK<Scalar>());
    stepper->setStepperType(stepperType);
    setStepperImplicitValues(stepper, stepperPL);
    setTableaus(stepper, stepperPL, stepperType);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperHHTAlpha<Scalar> >
  createStepperHHTAlpha(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperHHTAlpha<Scalar>());
    setStepperImplicitValues(stepper, stepperPL);
    if (stepperPL != Teuchos::null) {
      if (stepperPL->isSublist("HHT-Alpha Parameters")) {
        auto hhtalphaPL = stepperPL->sublist("HHT-Alpha Parameters", true);
        std::string schemeName =
          hhtalphaPL.get<std::string>("Scheme Name", "Newmark Beta Average Acceleration");
        stepper->setSchemeName(schemeName);
        if (schemeName == "Newmark Beta User Defined") {
          stepper->setBeta (hhtalphaPL.get<double>("Beta",  0.25));
          stepper->setGamma(hhtalphaPL.get<double>("Gamma", 0.5 ));
        }
        stepper->setAlphaF(hhtalphaPL.get<double>("Alpha_f",  0.0));
        stepper->setAlphaM(hhtalphaPL.get<double>("Alpha_m",  0.0));
      } else {
        stepper->setSchemeName("Newmark Beta Average Acceleration");
        stepper->setAlphaF(0.0);
        stepper->setAlphaM(0.0);
      }
    }

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperNewmarkImplicitDForm<Scalar> >
  createStepperNewmarkImplicitDForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperNewmarkImplicitDForm<Scalar>());
    setStepperImplicitValues(stepper, stepperPL);
    if (stepperPL != Teuchos::null) {
      if (stepperPL->isSublist("Newmark Parameters")) {
        auto newmarkPL = stepperPL->sublist("Newmark Parameters", true);
        std::string schemeName =
          newmarkPL.get<std::string>("Scheme Name", "Average Acceleration");
        stepper->setSchemeName(schemeName);
        if (schemeName == "User Defined") {
          stepper->setBeta (newmarkPL.get<double>("Beta",  0.25));
          stepper->setGamma(newmarkPL.get<double>("Gamma", 0.5 ));
        }
      } else {
        stepper->setSchemeName("Average Acceleration");
      }
    }

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperNewmarkImplicitAForm<Scalar> >
  createStepperNewmarkImplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperNewmarkImplicitAForm<Scalar>());
    setStepperImplicitValues(stepper, stepperPL);
    if (stepperPL != Teuchos::null) {
      if (stepperPL->isSublist("Newmark Parameters")) {
        auto newmarkPL = stepperPL->sublist("Newmark Parameters", true);
        std::string schemeName =
          newmarkPL.get<std::string>("Scheme Name", "Average Acceleration");
        stepper->setSchemeName(schemeName);
        if (schemeName == "User Defined") {
          stepper->setBeta (newmarkPL.get<double>("Beta",  0.25));
          stepper->setGamma(newmarkPL.get<double>("Gamma", 0.5 ));
        }
      } else {
        stepper->setSchemeName("Average Acceleration");
      }
    }

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperBDF2<Scalar> >
  createStepperBDF2(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperBDF2<Scalar>());
    setStepperImplicitValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      std::string startUpStepper = "DIRK 1 Stage Theta Method";
      if (stepperPL != Teuchos::null) startUpStepper =
        stepperPL->get<std::string>("Start Up Stepper Type", startUpStepper);
      stepper->setStartUpStepper(startUpStepper);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperTrapezoidal<Scalar> >
  createStepperTrapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperTrapezoidal<Scalar>());
    setStepperImplicitValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperBackwardEuler<Scalar> >
  createStepperBackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperBackwardEuler<Scalar>());
    setStepperImplicitValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);

      std::string predictor = "None";
      if (stepperPL != Teuchos::null) predictor =
        stepperPL->get<std::string>("Predictor Stepper Type", predictor);
      stepper->setPredictor(predictor);

      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperNewmarkExplicitAForm<Scalar> >
  createStepperNewmarkExplicitAForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperNewmarkExplicitAForm<Scalar>());
    setStepperExplicitValues(stepper, stepperPL);
    if (stepperPL != Teuchos::null) {
      Scalar gamma = stepperPL->sublist("Newmark Explicit Parameters")
                                   .template get<double>("Gamma", 0.5);
      stepper->setGamma(gamma);
    }

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperLeapfrog<Scalar> >
  createStepperLeapfrog(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperLeapfrog<Scalar>());
    setStepperExplicitValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperForwardEuler<Scalar> >
  createStepperForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperForwardEuler<Scalar>());
    setStepperExplicitValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_General<Scalar> >
  createStepperERK_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_General<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (stepperPL != Teuchos::null) {
      if (stepperPL->isParameter("Tableau")) {
        auto t = createTableau(stepperPL);
        stepper->setTableau( t->A(),t->b(),t->c(),
                             t->order(),t->orderMin(),t->orderMax(),
                             false,-4.0,
                             t->bstar() );
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      stepper->getTableau()->isImplicit() == true, std::logic_error,
      "Error - General ERK received an implicit Butcher Tableau!\n");

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_ForwardEuler<Scalar> >
  createStepperERK_ForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_ForwardEuler<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_4Stage4thOrder<Scalar> >
  createStepperERK_4Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_4Stage4thOrder<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_3_8Rule<Scalar> >
  createStepperERK_3_8Rule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_3_8Rule<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_4Stage3rdOrderRunge<Scalar> >
  createStepperERK_4Stage3rdOrderRunge(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_4Stage3rdOrderRunge<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_5Stage3rdOrderKandG<Scalar> >
  createStepperERK_5Stage3rdOrderKandG(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_5Stage3rdOrderKandG<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_3Stage3rdOrder<Scalar> >
  createStepperERK_3Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_3Stage3rdOrder<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_3Stage3rdOrderTVD<Scalar> >
  createStepperERK_3Stage3rdOrderTVD(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    const std::string stepperType)
  {
    auto stepper = Teuchos::rcp(new StepperERK_3Stage3rdOrderTVD<Scalar>());
    stepper->setStepperType(stepperType);
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_SSPERK54<Scalar> >
  createStepperERK_SSPERK54(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_SSPERK54<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_3Stage3rdOrderHeun<Scalar> >
  createStepperERK_3Stage3rdOrderHeun(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_3Stage3rdOrderHeun<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_Midpoint<Scalar> >
  createStepperERK_Midpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_Midpoint<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_Trapezoidal<Scalar> >
  createStepperERK_Trapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::string stepperType)
  {
    auto stepper = Teuchos::rcp(new StepperERK_Trapezoidal<Scalar>());
    stepper->setStepperType(stepperType);
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_Ralston<Scalar> >
  createStepperERK_Ralston(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL,
    std::string stepperType)
  {
    auto stepper = Teuchos::rcp(new StepperERK_Ralston<Scalar>());
    stepper->setStepperType(stepperType);
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_BogackiShampine32<Scalar> >
  createStepperERK_BogackiShampine32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_BogackiShampine32<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperERK_Merson45<Scalar> >
  createStepperERK_Merson45(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperERK_Merson45<Scalar>());
    setStepperRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      stepper->initialize();
    }

    return stepper;
  }


  Teuchos::RCP<StepperDIRK_General<Scalar> >
  createStepperDIRK_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperDIRK_General<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (stepperPL != Teuchos::null) {
      if (stepperPL->isParameter("Tableau")) {
        auto t = createTableau(stepperPL);
        stepper->setTableau( t->A(),t->b(),t->c(),
                             t->order(),t->orderMin(),t->orderMax(),
                             false, -4.0,
                             t->bstar() );
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      stepper->getTableau()->isDIRK() != true, std::logic_error,
      "Error - General DIRK did not receive a DIRK Butcher Tableau!\n");

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperDIRK_BackwardEuler<Scalar> >
  createStepperDIRK_BackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperDIRK_BackwardEuler<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_2Stage2ndOrder<Scalar> >
  createStepperSDIRK_2Stage2ndOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_2Stage2ndOrder<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);
    if (stepperPL != Teuchos::null)
      stepper->setGamma(stepperPL->get<double>("gamma", 0.2928932188134524));

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_3Stage2ndOrder<Scalar> >
  createStepperSDIRK_3Stage2ndOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_3Stage2ndOrder<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }


  Teuchos::RCP<StepperSDIRK_SSPDIRK22<Scalar> >
  createStepperSDIRK_SSPDIRK22(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_SSPDIRK22<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);
    //if (stepperPL != Teuchos::null)
      //stepper->setGamma(stepperPL->get<double>("gamma", 0.2928932188134524));

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_SSPDIRK32<Scalar> >
  createStepperSDIRK_SSPDIRK32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_SSPDIRK32<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);
    //if (stepperPL != Teuchos::null)
      //stepper->setGamma(stepperPL->get<double>("gamma", 0.2928932188134524));

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_SSPDIRK23<Scalar> >
  createStepperSDIRK_SSPDIRK23(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_SSPDIRK23<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);
    //if (stepperPL != Teuchos::null)
      //stepper->setGamma(stepperPL->get<double>("gamma", 0.2928932188134524));

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_SSPDIRK33<Scalar> >
  createStepperSDIRK_SSPDIRK33(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_SSPDIRK33<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);
    //if (stepperPL != Teuchos::null)
      //stepper->setGamma(stepperPL->get<double>("gamma", 0.2928932188134524));

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }


  Teuchos::RCP<StepperSDIRK_2Stage3rdOrder<Scalar> >
  createStepperSDIRK_2Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_2Stage3rdOrder<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);
    if (stepperPL != Teuchos::null) {
      stepper->setGammaType(
        stepperPL->get<std::string>("Gamma Type", "3rd Order A-stable"));
      stepper->setGamma(stepperPL->get<double>("gamma", 0.7886751345948128));
    }

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperEDIRK_2Stage3rdOrder<Scalar> >
  createStepperEDIRK_2Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperEDIRK_2Stage3rdOrder<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperDIRK_1StageTheta<Scalar> >
  createStepperDIRK_1StageTheta(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperDIRK_1StageTheta<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);
    if (stepperPL != Teuchos::null) {
      stepper->setTheta(stepperPL->get<double>("theta", 0.5));
    }

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperEDIRK_2StageTheta<Scalar> >
  createStepperEDIRK_2StageTheta(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperEDIRK_2StageTheta<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);
    if (stepperPL != Teuchos::null) {
      stepper->setTheta(stepperPL->get<double>("theta", 0.5));
    }

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperEDIRK_TrapezoidalRule<Scalar> >
  createStepperEDIRK_TrapezoidalRule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    // This stepper goes by various names (e.g., 'RK Trapezoidal Rule'
    // and 'RK Crank-Nicolson').  Make sure it is set to the default name.
    if (stepperPL != Teuchos::null)
      stepperPL->set<std::string>("Stepper Type", "RK Trapezoidal Rule");

    auto stepper = Teuchos::rcp(new StepperEDIRK_TrapezoidalRule<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_ImplicitMidpoint<Scalar> >
  createStepperSDIRK_ImplicitMidpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_ImplicitMidpoint<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperDIRK_1Stage1stOrderRadauIA<Scalar> >
  createStepperDIRK_1Stage1stOrderRadauIA(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperDIRK_1Stage1stOrderRadauIA<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperDIRK_2Stage2ndOrderLobattoIIIB<Scalar> >
  createStepperDIRK_2Stage2ndOrderLobattoIIIB(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperDIRK_2Stage2ndOrderLobattoIIIB<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_5Stage4thOrder<Scalar> >
  createStepperSDIRK_5Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_5Stage4thOrder<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_3Stage4thOrder<Scalar> >
  createStepperSDIRK_3Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_3Stage4thOrder<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_5Stage5thOrder<Scalar> >
  createStepperSDIRK_5Stage5thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_5Stage5thOrder<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  Teuchos::RCP<StepperSDIRK_21Pair<Scalar> >
  createStepperSDIRK_21Pair(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperSDIRK_21Pair<Scalar>());
    setStepperDIRKValues(stepper, stepperPL);

    if (model != Teuchos::null) {
      stepper->setModel(model);
      setStepperSolverValues(stepper, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

  // ---------------------------------------------------------------------------

private:

  /// Very simple factory method
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    using Teuchos::rcp;
    if (stepperType == "Forward Euler")
      return createStepperForwardEuler(model, stepperPL);
    else if (stepperType == "Backward Euler")
      return createStepperBackwardEuler(model, stepperPL);
    else if (stepperType == "Trapezoidal Method")
      return createStepperTrapezoidal(model, stepperPL);
    else if (stepperType == "BDF2")
      return createStepperBDF2(model, stepperPL);
    else if (stepperType == "Newmark Implicit a-Form")
      return createStepperNewmarkImplicitAForm(model, stepperPL);
    else if (stepperType == "Newmark Implicit d-Form")
      return createStepperNewmarkImplicitDForm(model, stepperPL);
    else if (stepperType == "Newmark Explicit a-Form")
      return createStepperNewmarkExplicitAForm(model, stepperPL);
    else if (stepperType == "HHT-Alpha")
      return createStepperHHTAlpha(model, stepperPL);
    else if (stepperType == "General ERK" )
      return createStepperERK_General(model, stepperPL);
    else if (stepperType == "RK Forward Euler" || stepperType == "RK1" )
      return createStepperERK_ForwardEuler(model, stepperPL);
    else if (stepperType == "RK Explicit 4 Stage" )
      return createStepperERK_4Stage4thOrder(model, stepperPL);
    else if (stepperType == "RK Explicit 3/8 Rule" )
      return createStepperERK_3_8Rule(model, stepperPL);
    else if (stepperType == "RK Explicit 4 Stage 3rd order by Runge" )
      return createStepperERK_4Stage3rdOrderRunge(model, stepperPL);
    else if (stepperType == "RK Explicit 5 Stage 3rd order by Kinnmark and Gray" )
      return createStepperERK_5Stage3rdOrderKandG(model, stepperPL);
    else if (stepperType == "RK Explicit 3 Stage 3rd order" )
      return createStepperERK_3Stage3rdOrder(model, stepperPL);
    else if (stepperType == "RK Explicit 3 Stage 3rd order TVD" ||
             stepperType == "SSPERK33" || stepperType == "SSPRK3" )
      return createStepperERK_3Stage3rdOrderTVD(model, stepperPL, stepperType);
    else if (stepperType == "RK Explicit 3 Stage 3rd order by Heun" )
      return createStepperERK_3Stage3rdOrderHeun(model, stepperPL);
    else if (stepperType == "RK Explicit Midpoint" )
      return createStepperERK_Midpoint(model, stepperPL);
    else if (stepperType == "RK Explicit Trapezoidal" ||
             stepperType == "Heuns Method" || stepperType == "SSPERK22" ||
             stepperType == "SSPRK2" )
      return createStepperERK_Trapezoidal(model, stepperPL, stepperType);
    else if (stepperType == "RK Explicit Ralston" || stepperType == "RK2" )
      return createStepperERK_Ralston(model, stepperPL, stepperType);
    else if (stepperType == "SSPERK54" )
      return createStepperERK_SSPERK54(model, stepperPL);
    else if (stepperType == "Bogacki-Shampine 3(2) Pair" )
      return createStepperERK_BogackiShampine32(model, stepperPL);
    else if (stepperType == "Merson 4(5) Pair" )
      return createStepperERK_Merson45(model, stepperPL);
    else if (stepperType == "General DIRK" )
      return createStepperDIRK_General(model, stepperPL);
    else if (stepperType == "RK Backward Euler" )
      return createStepperDIRK_BackwardEuler(model, stepperPL);
    else if (stepperType == "SDIRK 2 Stage 2nd order" )
      return createStepperSDIRK_2Stage2ndOrder(model, stepperPL);
    else if (stepperType == "SSPDIRK22" )
      return createStepperSDIRK_SSPDIRK22(model, stepperPL);
    else if (stepperType == "SDIRK 3 Stage 2nd order" )
      return createStepperSDIRK_3Stage2ndOrder(model, stepperPL);
    else if (stepperType == "SSPDIRK32" )
      return createStepperSDIRK_SSPDIRK32(model, stepperPL);
    else if (stepperType == "SSPDIRK23" )
      return createStepperSDIRK_SSPDIRK23(model, stepperPL);
    else if (stepperType == "SSPDIRK33" )
      return createStepperSDIRK_SSPDIRK33(model, stepperPL);
    else if (stepperType == "SDIRK 2 Stage 3rd order" )
      return createStepperSDIRK_2Stage3rdOrder(model, stepperPL);
    else if (stepperType == "EDIRK 2 Stage 3rd order" )
      return createStepperEDIRK_2Stage3rdOrder(model, stepperPL);
    else if (stepperType == "DIRK 1 Stage Theta Method" )
      return createStepperDIRK_1StageTheta(model, stepperPL);
    else if (stepperType == "EDIRK 2 Stage Theta Method" )
      return createStepperEDIRK_2StageTheta(model, stepperPL);
    else if (stepperType == "RK Trapezoidal Rule" ||
             stepperType == "RK Crank-Nicolson" )
      return createStepperEDIRK_TrapezoidalRule(model, stepperPL);
    else if (stepperType == "RK Implicit Midpoint" )
      return createStepperSDIRK_ImplicitMidpoint(model, stepperPL);
    else if (stepperType == "RK Implicit 1 Stage 1st order Radau IA" )
      return createStepperDIRK_1Stage1stOrderRadauIA(model, stepperPL);
    else if (stepperType == "RK Implicit 2 Stage 2nd order Lobatto IIIB" )
      return createStepperDIRK_2Stage2ndOrderLobattoIIIB(model, stepperPL);
    else if (stepperType == "SDIRK 5 Stage 4th order" )
      return createStepperSDIRK_5Stage4thOrder(model, stepperPL);
    else if (stepperType == "SDIRK 3 Stage 4th order" )
      return createStepperSDIRK_3Stage4thOrder(model, stepperPL);
    else if (stepperType == "SDIRK 5 Stage 5th order" )
      return createStepperSDIRK_5Stage5thOrder(model, stepperPL);
    else if ( stepperType == "SDIRK 2(1) Pair" )
      return createStepperSDIRK_21Pair(model, stepperPL);
    else if (
      stepperType == "IMEX RK 1st order" ||
      stepperType == "SSP1_111"          ||
      stepperType == "IMEX RK SSP2"      ||
      stepperType == "IMEX RK SSP3"      ||
      stepperType == "SSP3_332"          ||
      stepperType == "SSP2_222"          ||
      stepperType == "SSP2_222_L"        ||
      stepperType == "SSP2_222_A"        ||
      stepperType == "IMEX RK ARS 233"   ||
      stepperType == "ARS 233"           ||
      stepperType == "General IMEX RK" )
      return createStepperIMEX_RK(model, stepperType, stepperPL);
    else if (
      stepperType == "Partitioned IMEX RK 1st order" ||
      stepperType == "Partitioned IMEX RK SSP2"      ||
      stepperType == "Partitioned IMEX RK ARS 233"   ||
      stepperType == "General Partitioned IMEX RK" )
      return createStepperIMEX_RK_Partition(model, stepperType, stepperPL);
    else if (stepperType == "Leapfrog")
      return createStepperLeapfrog(model, stepperPL);
    else if (stepperType == "Subcycling")
      return createStepperSubcycling(model, stepperPL);
    else {
      Teuchos::RCP<Teuchos::FancyOStream> out =
        Teuchos::VerboseObjectBase::getDefaultOStream();
      Teuchos::OSTab ostab(out,1,"StepperFactory::createStepper");
      *out
      << "Unknown Stepper Type!  ('"+stepperType+"').\n"
      << "Here is a list of available Steppers.\n"
      << "  One-Step Methods:\n"
      << "    'Forward Euler'\n"
      << "    'Backward Euler'\n"
      << "    'Trapezoidal Method'\n"
      << "  Multi-Step Methods:\n"
      << "    'BDF2'\n"
      << "  Second-order PDE Methods:\n"
      << "    'Leapfrog'\n"
      << "    'Newmark Implicit a-Form'\n"
      << "    'Newmark Implicit d-Form'\n"
      << "    'Newmark Explicit a-Form'\n"
      << "    'HHT-Alpha'\n"
      << "  Explicit Runge-Kutta Methods:\n"
      << "    'RK Forward Euler (RK1)'\n"
      << "    'RK Explicit 4 Stage'\n"
      << "    'RK Explicit 3/8 Rule'\n"
      << "    'RK Explicit 4 Stage 3rd order by Runge'\n"
      << "    'RK Explicit 5 Stage 3rd order by Kinnmark and Gray'\n"
      << "    'RK Explicit 3 Stage 3rd order'\n"
      << "    'RK Explicit 3 Stage 3rd order TVD'\n"
      << "    'RK Explicit 3 Stage 3rd order by Heun'\n"
      << "    'RK Explicit Midpoint'\n"
      << "    'RK Explicit Trapezoidal' or 'Heuns Method'\n"
      << "    'Bogacki-Shampine 3(2) Pair'\n"
      << "    'SSPERK22 (SSPRK2)'\n"
      << "    'SSPERK33 (SSPRK3)'\n"
      << "    'SSPERK54'\n"
      << "    'General ERK'\n"
      << "  Implicit Runge-Kutta Methods:\n"
      << "    'RK Backward Euler'\n"
      << "    'DIRK 1 Stage Theta Method'\n"
      << "    'RK Implicit Midpoint'\n"
      << "    'SDIRK 1 Stage 1st order'\n"
      << "    'SDIRK 2 Stage 2nd order'\n"
      << "    'SDIRK 2 Stage 3rd order'\n"
      << "    'EDIRK 2 Stage 3rd order'\n"
      << "    'EDIRK 2 Stage Theta Method'\n"
      << "    'SDIRK 3 Stage 4th order'\n"
      << "    'SDIRK 5 Stage 4th order'\n"
      << "    'SDIRK 5 Stage 5th order'\n"
      << "    'SSPDIRK22'\n"
      << "    'SSPDIRK32'\n"
      << "    'SSPDIRK23'\n"
      << "    'SSPDIRK33'\n"
      << "    'SDIRK 2(1) Pair'\n"
      << "    'SDIRK 3 Stage 2nd order'\n"
      << "    'RK Trapezoidal Rule' or 'RK Crank-Nicolson'\n"
      << "    'General DIRK'\n"
      << "  Implicit-Explicit (IMEX) Methods:\n"
      << "    'IMEX RK 1st order'\n"
      << "    'SSP1_111'\n"
      << "    'IMEX RK SSP2 (SSP2_222)'\n"
      << "    'IMEX RK SSP3 (SSP3_332)'\n"
      << "    'IMEX RK ARS 233'\n"
      << "    'General IMEX RK'\n"
      << "    'Partitioned IMEX RK 1st order'\n"
      << "    'Partitioned IMEX RK SSP2'\n"
      << "    'Partitioned IMEX RK ARS 233'\n"
      << "    'General Partitioned IMEX RK'\n"
      << "  Steppers with subSteppers:\n"
      << "    'Operator Split'\n"
      << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Unknown 'Stepper Type' = " << stepperType);
    }
  }

// ==========================================================================

  void createSubSteppers(
    Teuchos::RCP<StepperOperatorSplit<Scalar> > stepper,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    if (stepperPL != Teuchos::null) {
      using Teuchos::RCP;
      using Teuchos::ParameterList;

      // Parse Stepper List String
      std::vector<std::string> stepperListStr;
      stepperListStr.clear();
      std::string str = stepperPL->get<std::string>("Stepper List");
      std::string delimiters(",");
      // Skip delimiters at the beginning
      std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
      // Find the first delimiter
      std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
      while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
        std::string token = str.substr(lastPos,pos-lastPos);
        // Strip single quotes
        std::string::size_type beg = token.find_first_of("'") + 1;
        std::string::size_type end = token.find_last_of ("'");
        stepperListStr.push_back(token.substr(beg,end-beg));

        lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters
        pos = str.find_first_of(delimiters, lastPos);     // Find next delimiter
      }

      TEUCHOS_TEST_FOR_EXCEPTION(stepperListStr.size() != appModels.size(),
        std::logic_error, "Error - Number of models and Steppers do not match!\n"
        << "  There are " << appModels.size() << " models.\n"
        << "  There are " << stepperListStr.size() << " steppers.\n"
        << "    " << str << "\n");

      typename
        std::vector<RCP<const Thyra::ModelEvaluator<Scalar> > >::iterator
          aMI = appModels.begin();
      typename std::vector<std::string>::iterator sLSI = stepperListStr.begin();

      for (; aMI<appModels.end() || sLSI<stepperListStr.end(); aMI++, sLSI++) {
        RCP<ParameterList> subStepperPL = Teuchos::sublist(stepperPL,*sLSI,true);
        bool useFSAL = subStepperPL->template get<bool>("Use FSAL",false);
        auto subStepper = createStepper(subStepperPL, *aMI);
        if (useFSAL) {
          Teuchos::RCP<Teuchos::FancyOStream> out =
            Teuchos::VerboseObjectBase::getDefaultOStream();
          Teuchos::OSTab ostab(out,1,"StepperFactory::createSubSteppers()");
          *out << "Warning -- subStepper = '"
               << subStepper->getStepperType() << "' has \n"
               << "  subStepper->getUseFSAL() = " << useFSAL << ".\n"
               << "  subSteppers usually can not use the FSAL priniciple with\n"
               << "  operator splitting.  Proceeding with it set to true.\n"
               << std::endl;
        }
        stepper->addStepper(subStepper, useFSAL);
      }
    }
  }

  Teuchos::RCP<StepperOperatorSplit<Scalar> >
  createStepperOperatorSplit(
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    auto stepper = Teuchos::rcp(new StepperOperatorSplit<Scalar>());

    if (stepperPL != Teuchos::null) {
      setStepperValues(stepper, stepperPL);
      stepper->setOrderMin(stepperPL->get<int>("Minimum Order", 1));
      stepper->setOrder   (stepperPL->get<int>("Order", 1));
      stepper->setOrderMax(stepperPL->get<int>("Maximum Order", 1));
    }

    if ( !(appModels.empty()) ) {
      createSubSteppers(stepper, appModels, stepperPL);
      stepper->initialize();
    }

    return stepper;
  }

// --------------------------------------------------------------------------
  Teuchos::RCP<Stepper<Scalar> > createStepper(
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models,
    std::string stepperType,
    Teuchos::RCP<Teuchos::ParameterList> stepperPL)
  {
    if (stepperType == "Operator Split")
      return createStepperOperatorSplit(models, stepperPL);
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Unknown 'Stepper Type' = " << stepperType);
    }
  }

};


} // namespace Tempus
#endif // Tempus_StepperFactory_hpp
