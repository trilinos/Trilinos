#include "Pike_ModelEvaluator_Composite.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"

namespace pike {

  template<typename Scalar>
  CompositeModelEvaluator::CompositeModelEvaluator()
  {
    
  }

  template<typename Scalar>
  void CompositeModelEvaluator::
  setModels(const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > >& me)
  {
    models_ = me;

    std::vector<Teuchos::RCP<Thyra::VectorSpace<Scalar> > > x_spaces;
    std::vector<Teuchos::RCP<Thyra::VectorSpace<Scalar> > > f_spaces;

    for (int model = 0; model < me.size(); ++model) {

      x_spaces.push_back(me->get_x_space());
      f_spaces.push_back(me->get_f_space());

    }

    x_space_ = Teuchos::rcp(new Thyra::DefaultProductVectorSpace<Scalar>(x_spaces));
    f_space_ = Teuchos::rcp(new Thyra::DefaultProductVectorSpace<Scalar>(f_spaces));
  }

    // From Thyra::ModelEvaluator
  template<typename Scalar>
  int CompositeModelEvaluator::Np() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Np() not implemented yet!");
    return 0;
  }
  
  int CompositeModelEvaluator::Ng() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Ng() not implemented yet!");
    return 0;
  }

    RCP<const VectorSpaceBase<Scalar> > get_x_space() const
    {
      return x_space_;
    }

    RCP<const VectorSpaceBase<Scalar> > get_f_space() const
    {
      return f_space_;
    }

    RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const
    {
      return p_space_;
    }

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

}

#endif
