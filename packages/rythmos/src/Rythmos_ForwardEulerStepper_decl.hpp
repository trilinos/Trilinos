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
#include "Rythmos_MomentoBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Rythmos_StateSerializerStrategy.hpp"

namespace Rythmos {

/** \brief Concrete momento class for the ForwardEulerStepper.
 * 
 * Note:  The model is not contained in the momento and must be set on the stepper in addition to the momento.
 */
  /*
template<class Scalar>
  class ForwardEulerStepperMomento :
    virtual public MomentoBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
    ForwardEulerStepperMomento() {}
    virtual ~ForwardEulerStepperMomento() {}

    void serialize(
        const StateSerializerStrategy<Scalar>& stateSerializer,
        std::ostream& oStream
        ) const
    {
      using Teuchos::is_null;
      TEUCHOS_ASSERT( !is_null(model_) );
      if (is_null(solution_vector_)) {
        solution_vector_ = Thyra::createMember(model_->get_x_space());
      }
      if (is_null(residual_vector_)) {
        residual_vector_ = Thyra::createMember(model_->get_f_space());
      }
      if (is_null(solution_vector_old_)) {
        solution_vector_old_ = Thyra::createMember(model_->get_x_space());
      }
      stateSerializer.serializeVectorBase(*solution_vector_,oStream);
      stateSerializer.serializeVectorBase(*residual_vector_,oStream);
      stateSerializer.serializeVectorBase(*solution_vector_old_,oStream);
      stateSerializer.serializeScalar(t_,oStream);
      stateSerializer.serializeScalar(t_old_,oStream);
      stateSerializer.serializeScalar(dt_,oStream);
      stateSerializer.serializeInt(numSteps_,oStream);
      stateSerializer.serializeBool(isInitialized_,oStream);
      stateSerializer.serializeBool(haveInitialCondition_,oStream);
      RCP<ParameterList> pl = parameterList_;
      if (Teuchos::is_null(pl)) {
        pl = Teuchos::parameterList();
      }
      stateSerializer.serializeParameterList(*pl,oStream);
    }

    void deSerialize(
        const StateSerializerStrategy<Scalar>& stateSerializer,
        std::istream& iStream
        )
    {
      using Teuchos::outArg;
      using Teuchos::is_null;
      TEUCHOS_ASSERT( !is_null(model_) );
      if (is_null(solution_vector_)) {
        solution_vector_ = Thyra::createMember(*model_->get_x_space());
      }
      if (is_null(residual_vector_)) {
        residual_vector_ = Thyra::createMember(*model_->get_f_space());
      }
      if (is_null(solution_vector_old_)) {
        solution_vector_old_ = Thyra::createMember(*model_->get_x_space());
      }
      stateSerializer.deSerializeVectorBase(outArg(*solution_vector_),iStream);
      stateSerializer.deSerializeVectorBase(outArg(*residual_vector_),iStream);
      stateSerializer.deSerializeVectorBase(outArg(*solution_vector_old_),iStream);
      stateSerializer.deSerializeScalar(outArg(t_),iStream);
      stateSerializer.deSerializeScalar(outArg(t_old_),iStream);
      stateSerializer.deSerializeScalar(outArg(dt_),iStream);
      stateSerializer.deSerializeInt(outArg(numSteps_),iStream);
      stateSerializer.deSerializeBool(outArg(isInitialized_),iStream);
      stateSerializer.deSerializeBool(outArg(haveInitialCondition_),iStream);
      if (is_null(parameterList_)) {
        parameterList_ = Teuchos::parameterList();
      }
      stateSerializer.deSerializeParameterList(outArg(*parameterList_),iStream);
    }

    RCP<MomentoBase<Scalar> > clone() const
    {
      RCP<ForwardEulerStepperMomento<Scalar> > m = rcp(new ForwardEulerStepperMomento<Scalar>());
      m->set_solution_vector(solution_vector_);
      m->set_residual_vector(residual_vector_);
      m->set_solution_vector_old(solution_vector_old_);
      m->set_t(t_);
      m->set_t_old(t_old_);
      m->set_dt(dt_);
      m->set_numSteps(numSteps_);
      m->set_isInitialized(isInitialized_);
      m->set_haveInitialCondition(haveInitialCondition_);
      m->set_parameterList(parameterList_);
      if (!Teuchos::is_null(this->getMyParamList())) {
        m->setParameterList(Teuchos::parameterList(*(this->getMyParamList())));
      }
      m->setModel(model_);
      m->setBasePoint(basePoint_);
      // How do I copy the VerboseObject data?  
      // 07/10/09 tscoffe:  Its not set up in Teuchos to do this yet
      return m;
    }

    void set_solution_vector(const RCP<const VectorBase<Scalar> >& solution_vector )
    { 
      solution_vector_ = Teuchos::null;
      if (!Teuchos::is_null(solution_vector)) {
        solution_vector_ = solution_vector->clone_v(); 
      }
    }
    RCP<VectorBase<Scalar> > get_solution_vector() const
    { return solution_vector_; }

    void set_residual_vector(const RCP<const VectorBase<Scalar> >& residual_vector )
    { 
      residual_vector_ = Teuchos::null;
      if (!Teuchos::is_null(residual_vector)) {
        residual_vector_ = residual_vector->clone_v(); 
      }
    }
    RCP<VectorBase<Scalar> > get_residual_vector() const
    { return residual_vector_; }

    void set_solution_vector_old(const RCP<const VectorBase<Scalar> >& solution_vector_old )
    { 
      solution_vector_old_ = Teuchos::null;
      if (!Teuchos::is_null(solution_vector_old)) {
        solution_vector_old_ = solution_vector_old->clone_v(); 
      }
    }
    RCP<VectorBase<Scalar> > get_solution_vector_old() const
    { return solution_vector_old_; }

    void set_t(const Scalar & t)
    { t_ = t; }
    Scalar get_t() const
    { return t_; }

    void set_t_old(const Scalar & t_old)
    { t_old_ = t_old; }
    Scalar get_t_old() const
    { return t_old_; }

    void set_dt(const Scalar & dt)
    { dt_ = dt; }
    Scalar get_dt() const
    { return dt_; }

    void set_numSteps(const int & numSteps)
    { numSteps_ = numSteps; }
    int get_numSteps() const
    { return numSteps_; }

    void set_isInitialized(const bool & isInitialized)
    { isInitialized_ = isInitialized; }
    bool get_isInitialized() const
    { return isInitialized_; }

    void set_haveInitialCondition(const bool & haveInitialCondition)
    { haveInitialCondition_ = haveInitialCondition; }
    bool get_haveInitialCondition() const
    { return haveInitialCondition_; }

    void set_parameterList(const RCP<const ParameterList>& pl)
    { 
      parameterList_ = Teuchos::null;
      if (!Teuchos::is_null(pl)) {
        parameterList_ = Teuchos::parameterList(*pl); 
      }
    }
    RCP<ParameterList> get_parameterList() const
    { return parameterList_; }

    void setParameterList(const RCP<ParameterList>& paramList)
    { this->setMyParamList(paramList); }
    RCP<const ParameterList> getValidParameters() const
    { return Teuchos::null; }

    void set_model(const RCP<const Thyra::ModelEvaluator<Scalar> >& model)
    { model_ = model; }
    RCP<const Thyra::ModelEvaluator<Scalar> > get_model() const
    { return model_; }

    void set_basePoint(const RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> >& basePoint)
    { basePoint_ = basePoint; }
    RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> > get_basePoint() const
    { return basePoint_; }

  private:
    RCP<Thyra::VectorBase<Scalar> > solution_vector_;
    RCP<Thyra::VectorBase<Scalar> > residual_vector_;
    RCP<Thyra::VectorBase<Scalar> > solution_vector_old_;
    Scalar t_;
    Scalar t_old_;
    Scalar dt_;
    int numSteps_;
    bool isInitialized_;
    bool haveInitialCondition_;
    RCP<ParameterList> parameterList_;
     
    // Objects that must be set prior to serialization and deSerialization:
    RCP<const Thyra::ModelEvaluator<Scalar> > model_;
    // Objects that must be set prior to calling ForwardEulerStepper::setMomento: 
    RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> > basePoint_;
};
*/
template<class Scalar>
  class ForwardEulerStepperMomento :
    virtual public MomentoBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
    ForwardEulerStepperMomento();
    virtual ~ForwardEulerStepperMomento();

    void serialize(
        const StateSerializerStrategy<Scalar>& stateSerializer,
        std::ostream& oStream
        ) const;

    void deSerialize(
        const StateSerializerStrategy<Scalar>& stateSerializer,
        std::istream& iStream
        );

    RCP<MomentoBase<Scalar> > clone() const;

    void set_solution_vector(const RCP<const VectorBase<Scalar> >& solution_vector );
    RCP<VectorBase<Scalar> > get_solution_vector() const;

    void set_residual_vector(const RCP<const VectorBase<Scalar> >& residual_vector );
    RCP<VectorBase<Scalar> > get_residual_vector() const;

    void set_solution_vector_old(const RCP<const VectorBase<Scalar> >& solution_vector_old );
    RCP<VectorBase<Scalar> > get_solution_vector_old() const;

    void set_t(const Scalar & t);
    Scalar get_t() const;

    void set_t_old(const Scalar & t_old);
    Scalar get_t_old() const;

    void set_dt(const Scalar & dt);
    Scalar get_dt() const;

    void set_numSteps(const int & numSteps);
    int get_numSteps() const;

    void set_isInitialized(const bool & isInitialized);
    bool get_isInitialized() const;

    void set_haveInitialCondition(const bool & haveInitialCondition);
    bool get_haveInitialCondition() const;

    void set_parameterList(const RCP<const ParameterList>& pl);
    RCP<ParameterList> get_parameterList() const;

    void setParameterList(const RCP<ParameterList>& paramList);
    RCP<const ParameterList> getValidParameters() const;

    void set_model(const RCP<Thyra::ModelEvaluator<Scalar> >& model);
    RCP<Thyra::ModelEvaluator<Scalar> > get_model() const;

    void set_basePoint(const RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> >& basePoint);
    RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> > get_basePoint() const;

  private:
    RCP<Thyra::VectorBase<Scalar> > solution_vector_;
    RCP<Thyra::VectorBase<Scalar> > residual_vector_;
    RCP<Thyra::VectorBase<Scalar> > solution_vector_old_;
    Scalar t_;
    Scalar t_old_;
    Scalar dt_;
    int numSteps_;
    bool isInitialized_;
    bool haveInitialCondition_;
    RCP<ParameterList> parameterList_;
     
    // Objects that must be set prior to serialization and deSerialization:
    RCP<Thyra::ModelEvaluator<Scalar> > model_;
    // Objects that must be set prior to calling ForwardEulerStepper::setMomento: 
    RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar> > basePoint_;
};


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
    void setModel(const RCP<Thyra::ModelEvaluator<Scalar> >& model);

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

    /** \brief Get momento object for use in restarts
    *
    */
    RCP<const MomentoBase<Scalar> > getMomento() const;

    /** \brief Set momento object for use in restarts
    *
    */
    void setMomento( const Ptr<const MomentoBase<Scalar> >& momentoPtr );

  private:

    RCP<Thyra::ModelEvaluator<Scalar> > model_;
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
    void checkConsistentState_();

};

// Nonmember constructor
template<class Scalar>
RCP<ForwardEulerStepper<Scalar> > forwardEulerStepper();

// Nonmember constructor
template<class Scalar>
RCP<ForwardEulerStepper<Scalar> > forwardEulerStepper(const RCP<Thyra::ModelEvaluator<Scalar> >& model);

} // namespace Rythmos

#endif //Rythmos_FORWARDEULER_STEPPER_DECL_H
