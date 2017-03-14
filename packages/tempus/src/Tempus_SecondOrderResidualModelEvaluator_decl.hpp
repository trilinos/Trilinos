// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_SecondOrderResidualModelEvaluator_decl_hpp
#define Tempus_SecondOrderResidualModelEvaluator_decl_hpp

#include <functional>
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

/** \brief A ModelEvaluator for residual evaluations given a state.
 *  This ModelEvaluator takes a state, x, and determines its residual,
 *  \f$ g(x) \f$, which is suitable for a nonlinear solve.  This is
 *  accomplished by computing the time derivative of the state, x_dot,
 *  (through Lambda functions), supplying the current time, and calling
 *  the application transient ModelEvaluator, \f$ f(\dot{x},x,t) \f$.
 *
 *  This class breaks the primary design principle for ModelEvaluators;
 *  it is not stateless!
 */
template <typename Scalar>
class SecondOrderResidualModelEvaluator
 : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:
  typedef Thyra::VectorBase<Scalar>  Vector;

  /// Constructor
  SecondOrderResidualModelEvaluator(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& transientModel)
    : transientModel_(transientModel),
      out_(Teuchos::VerboseObjectBase::getDefaultOStream())
  {
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
  }

  /// Set the underlying transient ModelEvaluator
  void setTransientModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me)
  { 
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
    transientModel_ = me; 
  }

  /// Get the underlying transient model 'f'
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getTransientModel() const
  { 
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
    return transientModel_; 
  }

  /// Set values needed in evalModelImpl 
  void initialize(Teuchos::RCP<const Vector> a, Teuchos::RCP<Vector> v_pred, 
                  Teuchos::RCP<Vector> d_pred, Scalar delta_t, 
                   Scalar t, Scalar beta, Scalar gamma) 
  {
    *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
    a_ = a; v_pred_ = v_pred; d_pred_ = d_pred; 
    delta_t_ = delta_t; t_ = t; beta_ = beta; gamma_ = gamma; 
  } 

  /// \name Overridden from Thyra::StateFuncModelEvaluatorBase
  //@{
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const
    { 
     *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
      return transientModel_->create_W_op(); 
    }

    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
    get_W_factory() const 
    { 
      *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
      return transientModel_->get_W_factory(); 
    }

    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const
    { 
      *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
      return transientModel_->get_f_space(); 
    }

    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int p) const
    { 
      *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
      return transientModel_->get_p_space(p); 
    };

    Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int p) const
    { 
      *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
      return transientModel_->get_p_names(p); 
    }

    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const
    { 
      *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
      return transientModel_->get_x_space(); 
    }

    Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const
    { 
      *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
      Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues = transientModel_->getNominalValues(); 
      nominalValues.set_x(a_); 
      return nominalValues; 
    }

    Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

    void evalModelImpl(
              const Thyra::ModelEvaluatorBase::InArgs<Scalar>  &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;
  //@}
  
private:

  /// Default constructor - not allowed
  SecondOrderResidualModelEvaluator();

private:

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > transientModel_;
  Scalar t_;
  Scalar gamma_;
  Scalar beta_;
  Scalar delta_t_;
  Teuchos::RCP<const Vector> a_; 
  Teuchos::RCP<Vector> d_pred_; 
  Teuchos::RCP<Vector> v_pred_; 
  Teuchos::RCP<Teuchos::FancyOStream> out_;

};

} // namespace Tempus

#endif // Tempus_SecondOrderResidualModelEvaluator_hpp
