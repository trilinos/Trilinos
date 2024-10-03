//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperRKBase_hpp
#define Tempus_StepperRKBase_hpp

#include "Thyra_VectorBase.hpp"

#include "Tempus_config.hpp"
#include "Tempus_Stepper.hpp"
#include "Tempus_RKButcherTableau.hpp"
#include "Tempus_StepperRKAppAction.hpp"
#include "Tempus_StepperRKModifierDefault.hpp"
#include "Tempus_Stepper_ErrorNorm.hpp"

namespace Tempus {

/** \brief Base class for Runge-Kutta methods, ExplicitRK, DIRK and IMEX.
 *
 *  Only common RK methods should be implemented in StepperRKBase.  All
 *  other Stepper methods should be implemented through Stepper,
 *  StepperExplicit or StepperImplicit.
 */
template <class Scalar>
class StepperRKBase : virtual public Tempus::Stepper<Scalar> {
 public:
  virtual Teuchos::RCP<const RKButcherTableau<Scalar>> getTableau() const
  {
    return tableau_;
  }

  virtual Scalar getOrder() const { return getTableau()->order(); }
  virtual Scalar getOrderMin() const { return getTableau()->orderMin(); }
  virtual Scalar getOrderMax() const { return getTableau()->orderMax(); }
  virtual int getNumberOfStages() const { return getTableau()->numStages(); }

  virtual int getStageNumber() const { return stageNumber_; }
  virtual void setStageNumber(int s) { stageNumber_ = s; }

  virtual void setUseEmbedded(bool a)
  {
    useEmbedded_ = a;
    this->setEmbeddedMemory();
    this->isInitialized_ = false;
  }

  virtual bool getUseEmbedded() const { return useEmbedded_; }

  virtual void setErrorNorm(const Teuchos::RCP<Stepper_ErrorNorm<Scalar>>
                                &errCalculator = Teuchos::null)
  {
    if (errCalculator != Teuchos::null) {
      stepperErrorNormCalculator_ = errCalculator;
    }
    else {
      auto er                     = Teuchos::rcp(new Stepper_ErrorNorm<Scalar>());
      stepperErrorNormCalculator_ = er;
    }
  }

  virtual void setAppAction(Teuchos::RCP<StepperRKAppAction<Scalar>> appAction)
  {
    if (appAction == Teuchos::null) {
      // Create default appAction
      stepperRKAppAction_ =
          Teuchos::rcp(new StepperRKModifierDefault<Scalar>());
    }
    else {
      stepperRKAppAction_ = appAction;
    }
    this->isInitialized_ = false;
  }

  virtual Teuchos::RCP<StepperRKAppAction<Scalar>> getAppAction() const
  {
    return stepperRKAppAction_;
  }

  /// Set StepperRK member data from the ParameterList.
  virtual void setStepperRKValues(Teuchos::RCP<Teuchos::ParameterList> pl)
  {
    if (pl != Teuchos::null) {
      pl->validateParametersAndSetDefaults(*this->getValidParameters());
      this->setStepperValues(pl);
      if (pl->isParameter("Use Embedded"))
        this->setUseEmbedded(pl->get<bool>("Use Embedded"));
    }
  }

  virtual Teuchos::RCP<RKButcherTableau<Scalar>> createTableau(
      Teuchos::RCP<Teuchos::ParameterList> pl)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        pl == Teuchos::null, std::runtime_error,
        "Error parsing general tableau.  ParameterList is null.\n");

    Teuchos::RCP<RKButcherTableau<Scalar>> tableau;

    Teuchos::RCP<Teuchos::ParameterList> tableauPL =
        sublist(pl, "Tableau", true);
    std::size_t numStages = 0;
    int order             = tableauPL->get<int>("order");
    Teuchos::SerialDenseMatrix<int, Scalar> A;
    Teuchos::SerialDenseVector<int, Scalar> b;
    Teuchos::SerialDenseVector<int, Scalar> c;
    Teuchos::SerialDenseVector<int, Scalar> bstar;

    // read in the A matrix
    {
      std::vector<std::string> A_row_tokens;
      Tempus::StringTokenizer(A_row_tokens, tableauPL->get<std::string>("A"),
                              ";", true);

      // this is the only place where numStages is set
      numStages = A_row_tokens.size();

      // allocate the matrix
      A.shape(Teuchos::as<int>(numStages), Teuchos::as<int>(numStages));

      // fill the rows
      for (std::size_t row = 0; row < numStages; row++) {
        // parse the row (tokenize on space)
        std::vector<std::string> tokens;
        Tempus::StringTokenizer(tokens, A_row_tokens[row], " ", true);

        std::vector<double> values;
        Tempus::TokensToDoubles(values, tokens);

        TEUCHOS_TEST_FOR_EXCEPTION(
            values.size() != numStages, std::runtime_error,
            "Error parsing A matrix, wrong number of stages in row "
                << row << ".\n");

        for (std::size_t col = 0; col < numStages; col++)
          A(row, col) = values[col];
      }
    }

    // size b and c vectors
    b.size(Teuchos::as<int>(numStages));
    c.size(Teuchos::as<int>(numStages));

    // read in the b vector
    {
      std::vector<std::string> tokens;
      Tempus::StringTokenizer(tokens, tableauPL->get<std::string>("b"), " ",
                              true);
      std::vector<double> values;
      Tempus::TokensToDoubles(values, tokens);

      TEUCHOS_TEST_FOR_EXCEPTION(
          values.size() != numStages, std::runtime_error,
          "Error parsing b vector, wrong number of stages.\n");

      for (std::size_t i = 0; i < numStages; i++) b(i) = values[i];
    }

    // read in the c vector
    {
      std::vector<std::string> tokens;
      Tempus::StringTokenizer(tokens, tableauPL->get<std::string>("c"), " ",
                              true);
      std::vector<double> values;
      Tempus::TokensToDoubles(values, tokens);

      TEUCHOS_TEST_FOR_EXCEPTION(
          values.size() != numStages, std::runtime_error,
          "Error parsing c vector, wrong number of stages.\n");

      for (std::size_t i = 0; i < numStages; i++) c(i) = values[i];
    }

    if (tableauPL->isParameter("bstar") &&
        tableauPL->get<std::string>("bstar") != "") {
      bstar.size(Teuchos::as<int>(numStages));
      // read in the bstar vector
      {
        std::vector<std::string> tokens;
        Tempus::StringTokenizer(tokens, tableauPL->get<std::string>("bstar"),
                                " ", true);
        std::vector<double> values;
        Tempus::TokensToDoubles(values, tokens);

        TEUCHOS_TEST_FOR_EXCEPTION(
            values.size() != numStages, std::runtime_error,
            "Error parsing bstar vector, wrong number of stages.\n"
                << "      Number of RK stages    = " << numStages
                << "\n      Number of bstar values = "
                << values.size() << "\n");

        for (std::size_t i = 0; i < numStages; i++) bstar(i) = values[i];
      }
      tableau = rcp(new RKButcherTableau<Scalar>("From ParameterList", A, b, c,
                                                 order, order, order, bstar));
    }
    else {
      tableau = rcp(new RKButcherTableau<Scalar>("From ParameterList", A, b, c,
                                                 order, order, order));
    }
    return tableau;
  }

 protected:
  virtual void setEmbeddedMemory() {}

  Teuchos::RCP<RKButcherTableau<Scalar>> tableau_;

  // For Embedded RK
  bool useEmbedded_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> ee_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> abs_u0;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> abs_u;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> sc;

  Teuchos::RCP<Stepper_ErrorNorm<Scalar>> stepperErrorNormCalculator_;

  /// The current Runge-Kutta stage number, {0,...,s-1}.  -1 indicates outside
  /// stage loop.
  int stageNumber_;
  Teuchos::RCP<StepperRKAppAction<Scalar>> stepperRKAppAction_;
};

}  // namespace Tempus

#endif  // Tempus_StepperRKBase_hpp
