// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_ADAPTIVESOLUTIONMANAGER
#define THYRA_ADAPTIVESOLUTIONMANAGER

#include "LOCA.H"
#include "LOCA_Thyra.H"
#include "LOCA_Thyra_GroupWrapper.H"
#include "Thyra_MultiVectorBase.hpp"


namespace Thyra {

class AdaptiveStateBase
{
public:

  AdaptiveStateBase(const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model);
  virtual ~AdaptiveStateBase () {}

  virtual void buildSolutionGroup() = 0;

  Teuchos::RCP< ::Thyra::ModelEvaluator<double> > getModel(){ return model_; }

protected:

  Teuchos::RCP< ::Thyra::ModelEvaluator<double> > model_;

};

class LOCAAdaptiveState : public AdaptiveStateBase
{
public:

  LOCAAdaptiveState(const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
                  const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy,
                  const Teuchos::RCP<LOCA::GlobalData>& global_data,
                  const Teuchos::RCP<LOCA::ParameterVector>& p,
                  int p_index);

  virtual ~LOCAAdaptiveState () {}

  //! Build the LOCA solution group
  void buildSolutionGroup();

  //! Accessor for the LOCA solution group
  Teuchos::RCP<LOCA::Thyra::Group>
       getSolutionGroup(){ return grp_; }


private:

  Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> saveDataStrategy_;
  Teuchos::RCP<LOCA::GlobalData> globalData_;
  Teuchos::RCP<LOCA::ParameterVector> paramVector_;
  int p_index_;

  //! The solution group
  Teuchos::RCP<LOCA::Thyra::GroupWrapper> grp_;


};

class TransAdaptiveState : public AdaptiveStateBase
{
public:

  TransAdaptiveState(const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model);

  virtual ~TransAdaptiveState () {}

  void buildSolutionGroup(){}

private:

};

class AdaptiveSolutionManager
{
public:

  AdaptiveSolutionManager ();

  virtual ~AdaptiveSolutionManager () {}

  void initialize(const Teuchos::RCP<Thyra::AdaptiveStateBase>& state){
        base = state;
  }

  Teuchos::RCP<Thyra::AdaptiveStateBase> getState(){ return base; }

  //! Remap "old" solution into new data structures
  virtual void projectCurrentSolution() = 0;

  //! Return current (adapted and remapped) solution from discretization
  virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > getCurrentSolution() = 0;

  //! Accessor functions

  virtual Teuchos::RCP<Teuchos::ParameterList> getAdaptParamsNonConst() { return adaptParams_; }

  virtual Teuchos::RCP<const Teuchos::ParameterList> getAdaptParams() const { return adaptParams_; }

  bool isAdaptive(){ return adaptiveMesh_; }

  //! Track the time state of the mesh
  void setTime(const double time){ time_ = time; }

  //! Track the time state of the mesh
  void setIteration(const int iter){ iter_ = iter; }

  //! Method called by Piro NOXSolver to determine if the mesh needs adapting
  // A return type of true means that the mesh should be adapted
  virtual bool queryAdaptationCriteria() = 0;

  //! Apply adaptation method to mesh and problem. Returns true if adaptation is performed successfully.
  virtual bool adaptProblem() = 0;

protected:

  //! The adaptation parameter list
  // These are set by the top-level adaptive solution manager
  Teuchos::RCP<Teuchos::ParameterList> adaptParams_;

  //! The parent of the solution parameter list
  // These are set by the top-level adaptive solution manager
  Teuchos::RCP<Teuchos::ParameterList> piroParams_;

  bool adaptiveMesh_;
  double time_;
  int iter_;

  Teuchos::RCP<Thyra::AdaptiveStateBase> base;

};

}

#endif //THYRA_ADAPTIVESOLUTIONMANAGER
