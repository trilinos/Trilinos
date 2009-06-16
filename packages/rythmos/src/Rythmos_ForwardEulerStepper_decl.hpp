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

#ifndef Rythmos_FORWARDEULER_STEPPER_DECL_H
#define Rythmos_FORWARDEULER_STEPPER_DECL_H

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_Types.hpp"
#include "Thyra_ModelEvaluator.hpp"

namespace Rythmos {

/** \brief . */
template<class Scalar>
class ForwardEulerStepper : virtual public StepperBase<Scalar>
{
  public:

    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief . */
    ForwardEulerStepper();

    /** \brief . */
    bool supportsCloning() const;

    /** \brief . */
    RCP<StepperBase<Scalar> > cloneStepperAlgorithm() const;

    /** \brief . */
    void setModel(const RCP<const Thyra::ModelEvaluator<Scalar> > &model);

    /** \brief . */
    RCP<const Thyra::ModelEvaluator<Scalar> >
    getModel() const;

    /** \brief . */
    void setInitialCondition(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
      );

    /** \brief . */
    RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
    
    /** \brief . */
    ~ForwardEulerStepper();

    /** \brief . */
    Scalar takeStep(Scalar dt, StepSizeType flag);

    /** \brief . */
    const StepStatus<Scalar> getStepStatus() const;

    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(
      Teuchos::FancyOStream            &out,
      const Teuchos::EVerbosityLevel   verbLevel
      ) const;
    
    /// Redefined from InterpolationBufferBase 
    /// Add points to buffer
    void addPoints(
      const Array<Scalar>& time_vec
      ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec
      ,const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
      );
    
    /// Get values from buffer
    void getPoints(
      const Array<Scalar>& time_vec
      ,Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec
      ,Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec
      ,Array<ScalarMag>* accuracy_vec
      ) const;

    /// Fill data in from another interpolation buffer
    void setRange(
      const TimeRange<Scalar>& range,
      const InterpolationBufferBase<Scalar> & IB
      );

    /** \brief . */
    TimeRange<Scalar> getTimeRange() const;

    /// Get interpolation nodes
    void getNodes(Array<Scalar>* time_vec) const;

    /// Remove interpolation nodes
    void removeNodes(Array<Scalar>& time_vec);

    /// Get order of interpolation
    int getOrder() const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(RCP<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    RCP<Teuchos::ParameterList> getNonconstParameterList();

    /** \brief . */
    RCP<Teuchos::ParameterList> unsetParameterList();

    /** \brief . */
    RCP<const Teuchos::ParameterList> getValidParameters() const;


  private:

    RCP<const Thyra::ModelEvaluator<Scalar> > model_;
    RCP<Thyra::VectorBase<Scalar> > solution_vector_;
    RCP<Thyra::VectorBase<Scalar> > residual_vector_;
    Scalar t_;
    Scalar dt_;
    Scalar t_old_;
    RCP<Thyra::VectorBase<Scalar> > solution_vector_old_;
    Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;
    int numSteps_;
    bool haveInitialCondition_;

    RCP<Teuchos::ParameterList> parameterList_;
    bool isInitialized_;

    // Private member functions:
    void defaultInitializAll_();
    void initialize_();

};

// Nonmember constructor
template<class Scalar>
RCP<ForwardEulerStepper<Scalar> > forwardEulerStepper();

// Nonmember constructor
template<class Scalar>
RCP<ForwardEulerStepper<Scalar> > forwardEulerStepper(const RCP<const Thyra::ModelEvaluator<Scalar> > &model);

} // namespace Rythmos

#endif //Rythmos_FORWARDEULER_STEPPER_DECL_H
