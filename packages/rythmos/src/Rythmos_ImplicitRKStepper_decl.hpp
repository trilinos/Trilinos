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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_IMPLICIT_RK_STEPPER_DECL_H
#define Rythmos_IMPLICIT_RK_STEPPER_DECL_H

#include "Rythmos_Types.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_DataStore.hpp"
#include "Rythmos_SolverAcceptingStepperBase.hpp"
#include "Rythmos_RKButcherTableauAcceptingStepperBase.hpp"
#include "Rythmos_RKButcherTableauBase.hpp"

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_NonlinearSolverBase.hpp"

namespace Rythmos {


/** \brief . */
template<class Scalar>
class ImplicitRKStepper : 
  virtual public SolverAcceptingStepperBase<Scalar>,
  virtual public RKButcherTableauAcceptingStepperBase<Scalar>
{
public:
  
  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;
  
  /** \name Constructors, intializers, Misc. */
  //@{

  /** \brief . */
  ImplicitRKStepper();
  
  void set_W_factory( const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &irk_W_factory );

  /** \brief . */
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  
  /** \name Overridden from RKButcherTableauAcceptingStepperBase */
  //@{
  
  /** \brief . */
  void setRKButcherTableau( const RCP<const RKButcherTableauBase<Scalar> > &rkButcherTableau );

  /** \brief . */
  RCP<const RKButcherTableauBase<Scalar> > getRKButcherTableau() const;

  //@}

  /** \brief . */
  // This function is mostly for testing purposes to explicitely over-ride the
  // internal RKBT detection to allow testing of 1-stage RKBTs as both fully
  // implicit RK and as DIRK methods.
  void setDirk(bool isDirk);

  //@}
  
  /** \name Overridden from SolverAcceptingStepperBase */
  //@{

  /** \brief . */
  void setSolver(
    const RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    );

  /** \brief . */
  RCP<Thyra::NonlinearSolverBase<Scalar> >
  getNonconstSolver();

  /** \brief . */
  RCP<const Thyra::NonlinearSolverBase<Scalar> >
  getSolver() const;

  //@}

  /** \name Overridden from StepperBase */
  //@{

  /** \brief Returns true. */
  bool isImplicit() const;
 
  /** \brief Returns true. */
  bool supportsCloning() const;

  /** \brief . */
  RCP<StepperBase<Scalar> > cloneStepperAlgorithm() const;

  /** \brief . */
  void setModel(const RCP<const Thyra::ModelEvaluator<Scalar> >& model);

  /** \brief . */
  void setNonconstModel(const RCP<Thyra::ModelEvaluator<Scalar> >& model);
  
  /** \brief . */
  RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const;

  /** \brief . */
  RCP<Thyra::ModelEvaluator<Scalar> > getNonconstModel();

  /** \brief . */
  void setInitialCondition(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
    );

  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getInitialCondition() const;

  /** \brief . */
  Scalar takeStep(Scalar dt, StepSizeType flag);
  
  /** \brief . */
  const StepStatus<Scalar> getStepStatus() const;
  
  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> >
  get_x_space() const;

  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    );
  
  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;
  
  /** \brief . */
  void getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;
  
  /** \brief . */
  void getNodes(Array<Scalar>* time_vec) const;
  
  /** \brief . */
  void removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}
  
  /** \name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  
  /** \brief. */
  RCP<const ParameterList> getValidParameters() const;
 
  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{
  
  /** \brief . */
  void describe(
    FancyOStream  &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

private:

  // ///////////////////////
  // Private date members

  bool isInitialized_;
  RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  RCP<Thyra::NonlinearSolverBase<Scalar> > solver_;
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > irk_W_factory_;
  RCP<ParameterList> paramList_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
  RCP<Thyra::VectorBase<Scalar> > x_;
  RCP<Thyra::VectorBase<Scalar> > x_old_;
  RCP<Thyra::VectorBase<Scalar> > x_dot_;

  TimeRange<Scalar> timeRange_;

  RCP<Thyra::ModelEvaluator<Scalar> > irkModel_;
  RCP<const RKButcherTableauBase<Scalar> > irkButcherTableau_;

  bool isDirk_; // Used for Diagonal Implicit RK 

  int numSteps_;

  bool haveInitialCondition_;

  // Cache
  RCP<Thyra::ProductVectorBase<Scalar> > x_stage_bar_;

  // //////////////////////////
  // Private member functions

  void defaultInitializeAll_();
  void initialize_();

};


/** \brief Nonmember constructor.
 *
 * \relates ImplicitRKStepper
 */
template<class Scalar>
RCP<ImplicitRKStepper<Scalar> >
implicitRKStepper();

template<class Scalar>
RCP<ImplicitRKStepper<Scalar> >
implicitRKStepper(
  const RCP<const Thyra::ModelEvaluator<Scalar> >& model,
  const RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> >& irk_W_factory,
  const RCP<const RKButcherTableauBase<Scalar> >& irkbt
  );


} // namespace Rythmos

#endif //Rythmos_IMPLICIT_RK_STEPPER_DECL_H
