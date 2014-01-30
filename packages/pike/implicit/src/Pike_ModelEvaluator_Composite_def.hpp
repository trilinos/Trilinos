#include "Pike_ModelEvaluator_Composite.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"

namespace pike {

  CompositeModelEvaluator() :
    setupCalled_(false)
    {

    }

    void addModel(const Teuchos::RCP<const Thyra::ModelEvaluator>& me)
    {
      models_.push_back(me);
    }

    //! Build the composite objects.  Once called, no more ModelEvaluators can be added via the addModel method.
    void setup()
    {
      // Create the x_space

      // Create the f_space

      // Create the composite p structures

      // Create the composite g structure

 
      setupCalled_ = true;
    }

    // From Thyra::ModelEvaluator
    int Np() const
    {
      
    }

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

}

#endif
