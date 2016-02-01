//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2012) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
// Glen Hansen (gahanse@sandia.gov), Sandia National Laboratories.
//
//@HEADER

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
