#ifndef PIKE_MODEL_EVALUATOR_COMPOSITE_HPP
#define PIKE_MODEL_EVALUATOR_COMPOSITE_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

namespace pike {
  
  /** \brief Forms a Product or Composite ModelEvaluator from multiple model evaluators.

   */ 
  template<typename Scalar>
  class CompositeModelEvaluator : Thyra::StateFuncModelEvaluatorBase<Scalar> {

  public:

    CompositeModelEvaluator();

    void addModel(const Teuchos::RCP<const Thyra::ModelEvaluator>& me);

    //! Build the composite objects.  Once called, no more ModelEvaluators can be added via the addModel method.
    void setup();

    // From Thyra::ModelEvaluator
    int Np() const;
    int Ng() const;
    RCP<const VectorSpaceBase<Scalar> > get_x_space() const;    
    RCP<const VectorSpaceBase<Scalar> > get_f_space() const;
    RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
    RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
    RCP<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
    ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
    RCP<LinearOpWithSolveBase<Scalar> > create_W() const;
    RCP<LinearOpBase<Scalar> > create_W_op() const;
    RCP<PreconditionerBase<Scalar> > create_W_prec() const;
    RCP<const LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
    ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
    ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
    void evalModel(const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
		   const ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  private:
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models_;
    bool setupCalled_;
    
    int Np_;
    int Ng_;

  };

}

#endif
