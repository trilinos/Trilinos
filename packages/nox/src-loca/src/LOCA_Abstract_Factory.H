// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_ABSTRACT_FACTORY_H
#define LOCA_ABSTRACT_FACTORY_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"

#include <vector>

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace MultiPredictor {
    class AbstractStrategy;
  }
  namespace MultiContinuation {
    class AbstractStrategy;
    class AbstractGroup;
  }
  namespace StepSize {
    class AbstractStrategy;
  }
  namespace BorderedSolver {
    class AbstractStrategy;
  }
  namespace Eigensolver {
    class AbstractStrategy;
  }
  namespace EigenvalueSort {
    class AbstractStrategy;
  }
  namespace SaveEigenData {
    class AbstractStrategy;
  }
  namespace AnasaziOperator {
    class AbstractStrategy;
  }
  namespace TurningPoint {
    namespace MooreSpence {
      class SolverStrategy;
    }
  }
  namespace Pitchfork {
    namespace MooreSpence {
      class SolverStrategy;
    }
  }
  namespace Hopf {
    namespace MooreSpence {
      class SolverStrategy;
    }
  }
}

namespace LOCA {

  namespace Abstract {

    //! %Abstract interface for providing a user-defined factory
    /*!
     * %LOCA::Abstract::Factory provides an abstract interface for providing
     * user-defined factories to the LOCA::Factory.  The LOCA::Factory provides
     * a mechanism for instantiating different strategies based on
     * parameter list choices.  This class allows additional strategies to
     * be instantiated by the factory without modifying the factory itself.
     * This is done by deriving a user-defined factory from this interface,
     * implementing any of the create methods for the user-defined strategies,
     * and passing an instance of the derived factory to the LOCA::Factory
     * object.  Any derived class must implement the
     * \c init() method to set the global data object which the
     * factory can then pass to any instantiated strategies.
     */
    class Factory {

    public:

      //! Constructor.
      Factory() {}

      //! Destructor
      virtual ~Factory() {}

      //! Initialize factory
      /*!
       * The LOCA::Factory will call this method to initialize the
       * user provided factory.  The user-provided factory should perform
       * any needed initialization here that cannot occur at construction.
       */
      virtual void
      init(const Teuchos::RCP<LOCA::GlobalData>& global_data) = 0;

      /*!
       * @name Strategy create methods
       */
      //@{

      //! Create predictor strategy
      virtual bool
      createPredictorStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& predictorParams,
       Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& strategy);

      //! Create continuation strategy
      virtual bool
      createContinuationStrategy(
      const std::string& strategyName,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& stepperParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs,
      Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>& strategy);

      //! Create bifurcation strategy
      virtual bool
      createBifurcationStrategy(
      const std::string& strategyName,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& bifurcationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& strategy);

      //! Create step size strategy
      virtual bool
      createStepSizeStrategy(
         const std::string& strategyName,
         const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams,
     Teuchos::RCP<LOCA::StepSize::AbstractStrategy>& strategy);

      //! Create bordered system solver strategy
      virtual bool
      createBorderedSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>& strategy);

      //! Create eigensolver strategy
      virtual bool
      createEigensolverStrategy(
         const std::string& strategyName,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
     Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>& strategy);

      //! Create eigenvalue sorting strategy
      virtual bool
      createEigenvalueSortStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
       Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>& strategy);

      //! Create strategy to save eigenvector/value data
      virtual bool
      createSaveEigenDataStrategy(
        const std::string& strategyName,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>& strategy);

      //! Create Anasazi operator
      virtual bool
      createAnasaziOperatorStrategy(
      const std::string& strategyName,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
      const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
      const Teuchos::RCP<NOX::Abstract::Group>& grp,
      Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>& strategy);

      //! Create Moore-Spence turning point solver strategy
      virtual bool
      createMooreSpenceTurningPointSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy>& strategy);

      //! Create Moore-Spence pitchfork solver strategy
      virtual bool
      createMooreSpencePitchforkSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::Pitchfork::MooreSpence::SolverStrategy>& strategy);

      //! Create Moore-Spence Hopf solver strategy
      virtual bool
      createMooreSpenceHopfSolverStrategy(
       const std::string& strategyName,
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy>& strategy);

      //@}

    }; // Class Factory

  } // Namespace Abstract

} // Namespace LOCA

#endif
