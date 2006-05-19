// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_STATE_ELIMINATION_MODEL_EVALUATOR_HPP
#define THYRA_DEFAULT_STATE_ELIMINATION_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {

/** \brief This class wraps any ModelEvaluator object along with a NonlinearSolverBase object
 * and eliminates the steady-state equations f(x,...)=0
 *
 * ToDo: Finish documentation!
 *
 */
template<class Scalar>
class DefaultStateEliminationModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultStateEliminationModelEvaluator();

  /** \brief . */
  DefaultStateEliminationModelEvaluator(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
    ,const Teuchos::RefCountPtr<NonlinearSolverBase<Scalar> >           &stateSolver
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
    ,const Teuchos::RefCountPtr<NonlinearSolverBase<Scalar> >           &stateSolver
    );

  /** \brief . */
  void uninitialize(
    Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 *thyraModel  = NULL
    ,Teuchos::RefCountPtr<NonlinearSolverBase<Scalar> >           *stateSolver = NULL
    );

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DfDp_op(int l) const;
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> > create_DgDx_op(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                thyraModel_;
  Teuchos::RefCountPtr<NonlinearSolverBase<Scalar> >           stateSolver_;

  Teuchos::RefCountPtr<DefaultNominalBoundsOverrideModelEvaluator<Scalar> >  wrappedThyraModel_;

  mutable Teuchos::RefCountPtr<VectorBase<Scalar> >            x_guess_solu_;
  
};

// /////////////////////////////////
// Implementations

// Constructors/initializers/accessors/utilities

template<class Scalar>
DefaultStateEliminationModelEvaluator<Scalar>::DefaultStateEliminationModelEvaluator()
{}

template<class Scalar>
DefaultStateEliminationModelEvaluator<Scalar>::DefaultStateEliminationModelEvaluator(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
  ,const Teuchos::RefCountPtr<NonlinearSolverBase<Scalar> >           &stateSolver
  )
{
  initialize(thyraModel,stateSolver);
}

template<class Scalar>
void DefaultStateEliminationModelEvaluator<Scalar>::initialize(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 &thyraModel
  ,const Teuchos::RefCountPtr<NonlinearSolverBase<Scalar> >           &stateSolver
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  TEST_FOR_EXCEPT(!stateSolver.get());
  stateSolver_ = stateSolver;
  const ModelEvaluatorBase::InArgs<Scalar>
    nominalValues = thyraModel->getNominalValues();
  if(nominalValues.get_x().get()) {
    x_guess_solu_ = nominalValues.get_x()->clone_v();
  }
  else {
    x_guess_solu_ = createMember(thyraModel->get_x_space());
    assign(&*x_guess_solu_,Scalar(0.0));
  }
  wrappedThyraModel_ = rcp(
    new DefaultNominalBoundsOverrideModelEvaluator<Scalar>(
      Teuchos::rcp_const_cast<ModelEvaluator<Scalar> >(thyraModel)
      ,Teuchos::null
      )
    );
  stateSolver_->setModel(wrappedThyraModel_);
}

template<class Scalar>
void DefaultStateEliminationModelEvaluator<Scalar>::uninitialize(
  Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                 *thyraModel
  ,Teuchos::RefCountPtr<NonlinearSolverBase<Scalar> >           *stateSolver
  )
{
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  if(stateSolver) *stateSolver = stateSolver_;
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
  stateSolver_ = Teuchos::null;
  wrappedThyraModel_ = Teuchos::null;
  x_guess_solu_ = Teuchos::null;
}

// Overridden from ModelEvaulator.

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultStateEliminationModelEvaluator<Scalar>::get_x_space() const
{
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >
DefaultStateEliminationModelEvaluator<Scalar>::get_f_space() const
{
  return Teuchos::null;
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultStateEliminationModelEvaluator<Scalar>::getNominalValues() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  MEB::InArgsSetup<Scalar> nominalValues(thyraModel->getNominalValues());
  nominalValues.setModelEvalDescription(this->description());
  nominalValues.setUnsupportsAndRelated(MEB::IN_ARG_x); // Wipe out x, x_dot ...
  return nominalValues;
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultStateEliminationModelEvaluator<Scalar>::getLowerBounds() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  MEB::InArgsSetup<Scalar> lowerBounds(thyraModel->getLowerBounds());
  lowerBounds.setModelEvalDescription(this->description());
  lowerBounds.setUnsupportsAndRelated(MEB::IN_ARG_x); // Wipe out x, x_dot ...
  return lowerBounds;
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultStateEliminationModelEvaluator<Scalar>::getUpperBounds() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  MEB::InArgsSetup<Scalar> upperBounds(thyraModel->getUpperBounds());
  upperBounds.setModelEvalDescription(this->description());
  upperBounds.setUnsupportsAndRelated(MEB::IN_ARG_x); // Wipe out x, x_dot ...
  return upperBounds;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
DefaultStateEliminationModelEvaluator<Scalar>::create_W() const
{
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultStateEliminationModelEvaluator<Scalar>::create_W_op() const
{
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultStateEliminationModelEvaluator<Scalar>::create_DfDp_op(int l) const
{
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultStateEliminationModelEvaluator<Scalar>::create_DgDx_op(int j) const
{
  return Teuchos::null;
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultStateEliminationModelEvaluator<Scalar>::createInArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::InArgs<Scalar> wrappedInArgs = thyraModel->createInArgs();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(wrappedInArgs.Np());
  inArgs.setSupports(wrappedInArgs);
  inArgs.setUnsupportsAndRelated(MEB::IN_ARG_x); // Wipe out x, x_dot ...
  return inArgs;
}

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultStateEliminationModelEvaluator<Scalar>::createOutArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  const int Np = wrappedOutArgs.Np(), Ng = wrappedOutArgs.Ng();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np,Ng);
  outArgs.setSupports(wrappedOutArgs);
  outArgs.setUnsupportsAndRelated(MEB::IN_ARG_x); // wipe out DgDx ...
  outArgs.setUnsupportsAndRelated(MEB::OUT_ARG_f); // wipe out f, W, DfDp ...
  return outArgs;
}

template<class Scalar>
void DefaultStateEliminationModelEvaluator<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar>     &inArgs
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &outArgs
  ) const
{
  typedef ModelEvaluatorBase MEB;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const Teuchos::RefCountPtr<Teuchos::FancyOStream> out       = this->getOStream();
  const Teuchos::EVerbosityLevel                    verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::DefaultStateEliminationModelEvaluator<Scalar>::evalModel(...) ...\n";

  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();

  const int Np = outArgs.Np(), Ng = outArgs.Ng();

  // Reset the nominal values
  MEB::InArgs<Scalar> wrappedNominalValues = thyraModel->getNominalValues();
  wrappedNominalValues.setArgs(inArgs,true);
  wrappedNominalValues.set_x(x_guess_solu_);
  
  typedef Teuchos::VerboseObjectTempState<ModelEvaluatorBase> VOTSME;
  //VOTSME thyraModel_outputTempState(rcp(&wrappedThyraModel,false),out,verbLevel);

  typedef Teuchos::VerboseObjectTempState<NonlinearSolverBase<Scalar> > VOTSNSB;
  VOTSNSB statSolver_outputTempState(
    stateSolver_,out
    ,static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW) ? Teuchos::VERB_LOW : Teuchos::VERB_NONE 
    );

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\ninArgs =\n" << Teuchos::describe(inArgs,verbLevel)
      << "\noutArgs on input =\n" << Teuchos::describe(outArgs,Teuchos::VERB_LOW);

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nSolving f(x,...) for x ...\n";

  wrappedThyraModel_->setNominalValues(
    rcp(new MEB::InArgs<Scalar>(wrappedNominalValues))
    );
  
  SolveStatus<Scalar> solveStatus = stateSolver_->solve(&*x_guess_solu_,NULL);

  if( solveStatus.solveStatus == SOLVE_STATUS_CONVERGED ) {
    
    if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
      *out << "\nComputing the output functions at the solved state solution ...\n";

    MEB::InArgs<Scalar>   wrappedInArgs  = thyraModel->createInArgs();
    MEB::OutArgs<Scalar>  wrappedOutArgs = thyraModel->createOutArgs();
    wrappedInArgs.setArgs(inArgs,true);
    wrappedInArgs.set_x(x_guess_solu_);
    wrappedOutArgs.setArgs(outArgs,true);
    
    for( int l = 0; l < Np; ++l ) {
      for( int j = 0; j < Ng; ++j ) {
        if(
          outArgs.supports(MEB::OUT_ARG_DgDp,j,l).none()==false
          && outArgs.get_DgDp(j,l).isEmpty()==false
          )
        {
          // Set DfDp(l) and DgDx(j) to be computed!
          //wrappedOutArgs.set_DfDp(l,...);
          //wrappedOutArgs.set_DgDx(j,...);
          TEST_FOR_EXCEPT(true);
        }
      }
    }
    
    thyraModel->evalModel(wrappedInArgs,wrappedOutArgs);

    //
    // Compute DgDp(j,l) using direct sensitivties
    //
    for( int l = 0; l < Np; ++l ) {
      if(
        wrappedOutArgs.supports(MEB::OUT_ARG_DfDp,l).none()==false
        && wrappedOutArgs.get_DfDp(l).isEmpty()==false
        )
      {
        //
        // Compute:  D(l) = -inv(DfDx)*DfDp(l)
        //
        TEST_FOR_EXCEPT(true);
        for( int j = 0; j < Ng; ++j ) {
          if(
            outArgs.supports(MEB::OUT_ARG_DgDp,j,l).none()==false
            && outArgs.get_DgDp(j,l).isEmpty()==false
            )
          {
            //
            // Compute:  DgDp(j,l) = DgDp(j,l) + DgDx(j)*D
            //
            TEST_FOR_EXCEPT(true);
          }
        }
      }
    }
    // ToDo: Add a mode to compute DgDp(l) using adjoint sensitivities?
    
  }
  else {
    
    if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
      *out << "\nFailed to converge, returning NaNs ...\n";
    outArgs.setFailed();
    
  }
  
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\noutArgs on output =\n" << Teuchos::describe(outArgs,verbLevel);

  totalTimer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nTotal evaluation time = "<<totalTimer.totalElapsedTime()<<" sec\n"
      << "\nLeaving Thyra::DefaultStateEliminationModelEvaluator<Scalar>::evalModel(...) ...\n";
  
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultStateEliminationModelEvaluator<Scalar>::description() const
{
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultStateEliminationModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << ",stateSolver=";
  if(stateSolver_.get())
    oss << "\'"<<stateSolver_->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}

} // namespace Thyra

#endif // THYRA_DEFAULT_STATE_ELIMINATION_MODEL_EVALUATOR_HPP
