//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef RYTHMOS_CHARON_SOLVER_HPP
#define RYTHMOS_CHARON_SOLVER_HPP

// General Trilinos includes:
#include "Thyra_ModelEvaluatorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"

// Rythmos includes:
#include "Rythmos_StepperSupportTypes.hpp"

namespace RythmosCharon {
  class CharonIntegrationControlAndObserver;
}
namespace Thyra {
  template<class Scalar> class NonlinearSolverBase;
  class EpetraModelEvaluator;
}
namespace Rythmos {
  template<class Scalar> class SolverAcceptingStepperBase;
  template<class Scalar> class IntegratorBase;
}

class RythmosCharonSolver
  : virtual public Teuchos::VerboseObject<RythmosCharonSolver>,
    virtual public Teuchos::ParameterListAcceptor
{
public:
  static Teuchos::RCP<const Teuchos::ParameterList> getStaticValidParameters();
  static Teuchos::RCP<Thyra::NonlinearSolverBase<double> >
  buildTimeStepNonlinearSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &timeStepNonlinearSolverSublist,
    Teuchos::FancyOStream &out
    );
  static Teuchos::RCP<Rythmos::SolverAcceptingStepperBase<double> >
  buildRythmosTimeStepper(
    const Teuchos::RCP<Teuchos::ParameterList> &rythmosStepperSelectionPL,
    Teuchos::FancyOStream &out
    );
  static Teuchos::RCP<RythmosCharon::CharonIntegrationControlAndObserver>
  buildIntegrationControlAndObserverStrategy(
    const Teuchos::RCP<Teuchos::ParameterList> &charonStepControlAndObservationSettingsPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<double> > &epetraThyraModel,
    Teuchos::FancyOStream &out
    );
  Thyra::ModelEvaluatorBase::InArgs<double> getStateInitialCondition();
  static void updateCharonState(
    const Thyra::VectorBase<double> &x_dot,
    const Thyra::VectorBase<double> &x,
    const double currentTime,
    const double timeStep,
    const bool guaranteeUpToDateAuxiliaryData,
    Teuchos::FancyOStream &out
    );
  static void writeOutput(
    const double currentTime,
    const double timeStep,
    const int outIter,
    Teuchos::FancyOStream &out
    );
  RythmosCharonSolver();
  void setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const& paramList
    );
  Teuchos::RCP<Teuchos::ParameterList> getParameterList()
  { return getNonconstParameterList(); }
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  void setup();
  bool solve();
private:
  Teuchos::RCP<Teuchos::ParameterList> paramList_;
  Teuchos::RCP<Rythmos::SolverAcceptingStepperBase<double> > rythmosStepper_;
  Teuchos::RCP<RythmosCharon::CharonIntegrationControlAndObserver>
    ryhmosCharonIntegrationControlAndObserver_;
  Teuchos::RCP<Rythmos::IntegratorBase<double> > rythmosStateIntegrator_;
  Teuchos::RCP<Thyra::ModelEvaluator<double> > epetraThyraModel_;
};


#endif // RYTHMOS_CHARON_SOLVER_HPP

