#ifndef PIKE_MODEL_EVALUATOR_COMPOSITE_HPP
#define PIKE_MODEL_EVALUATOR_COMPOSITE_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"

namespace pike {
  
  /** \brief Forms a Product or Composite ModelEvaluator from multiple model evaluators.

   */ 
  template<typename Scalar>
  class CompositeModelEvaluator : Thyra::StateFuncModelEvaluatorBase<Scalar> {

  public:

    CompositeModelEvaluator();

    //! Register the model evaluators and build the product objects.
    void setModels(const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > >& me);

    // From Thyra::ModelEvaluator
    int Np() const;
    int Ng() const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;    
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
    Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
    Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
    Teuchos::RCP<Thyra::PreconditionerBase<Scalar> > create_W_prec() const;
    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
    Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
    void evalModel(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
		   const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  private:
    std::vector<const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > models_;

    Teuchos::RCP<Thyra::DefaultProductVectorSpace<Scalar> > x_space_;
    Teuchos::RCP<Thyra::DefaultProductVectorSpace<Scalar> > f_space_;
    std::vector<Teuchos::RCP<Thyra::VectorSpace<Scalar> > > p_spaces_;
    std::vector<Teuchos::RCP<Thyra::VectorSpace<Scalar> > > g_spaces_;

    //! Maps a model parameter to a submodel evaluator and index for that me.
    std::vector<std::pair<int,int> > p_map_;
    std::vector<std::pair<int,int> > g_map_;
    
    int Np_;
    int Ng_;

  };

}

#endif
