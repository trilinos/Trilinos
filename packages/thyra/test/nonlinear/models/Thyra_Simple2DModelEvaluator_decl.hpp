
#ifndef THYRA_SIMPLE_2D_MODEL_EVALUATOR_DECL_HPP
#define THYRA_SIMPLE_2D_MODEL_EVALUATOR_DECL_HPP


#include "Thyra_StateFuncModelEvaluatorBase.hpp"


namespace Thyra {


template<class Scalar> class Simple2DModelEvaluator;


/** \brief Nonmember constuctor.
 *
 * \relates Simple2DModelEvaluator
 */
template<class Scalar>
Teuchos::RCP<Simple2DModelEvaluator<Scalar> >
simple2DModelEvaluator();


/** \brief Simple 2d simulation only ModelEvaluator for f(x) = 0.
 *
 * The equations modeled are:

 \verbatim

    f[0] =       x[0]      + x[1]*x[1] - p[0];
    f[1] = d * ( x[0]*x[0] - x[1]      - p[1] );

 \endverbatim

 * The Matrix <tt>W = d(f)/d(x)</tt> is implemented as a
 * <tt>Thyra::MultiVectorBase</tt> object and the class
 * <tt>Thyra::DefaultSerialDenseLinearOpWithSolveFactory</tt> is used to
 * create the linear solver.
 *
 * This is really more of a mock test driver model for Thyra than an example
 * of implementing a real simulation-constrained ModelEvaluator subclass.
 * That is really what the EpetraExt::ModelEvaluator is for.
 */
template<class Scalar>
class Simple2DModelEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Initializers/Accessors */
  //@{

  /** \brief . */
  void set_d(const Scalar &d);

  /** \brief . */
  void set_p(const Teuchos::ArrayView<const Scalar> &p);

  /** \brief . */
  void set_x0(const Teuchos::ArrayView<const Scalar> &x0);

  /** \brief . */
  void setShowGetInvalidArgs(bool showGetInvalidArg);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

#ifndef TEMPLATE_FRIENDS_NOT_SUPPORTED

  friend Teuchos::RCP<Simple2DModelEvaluator<Scalar> > simple2DModelEvaluator<>();
  
private:

#endif // TEMPLATE_FRIENDS_NOT_SUPPORTED

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private: // data members

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Scalar d_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x0_;
  Teuchos::Array<Scalar> p_;
  bool showGetInvalidArg_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

#ifdef TEMPLATE_FRIENDS_NOT_SUPPORTED
public: // member functions
#else
private: // member functions
#endif

  Simple2DModelEvaluator();

};


} // namespace Thyra


#endif // THYRA_SIMPLE_2D_MODEL_EVALUATOR_DECL_HPP
