#ifndef NOX_THYRA_MODEL_EVALUATOR_2DSIM_DECL_HPP
#define NOX_THYRA_MODEL_EVALUATOR_2DSIM_DECL_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

template<class Scalar> class ModelEvaluator2DSim;

/** \brief Nonmember constuctor.
 *
 * \relates ModelEvaluator2DSim
 */
template<class Scalar>
Teuchos::RCP<ModelEvaluator2DSim<Scalar> >
modelEvaluator2DSim(const Teuchos::RCP<const Epetra_Comm>& comm,
		    const Scalar d = 10.0,
		    const Scalar p0 = 2.0,
		    const Scalar p1 = 0.0,
		    const Scalar x0 = 0.0,
		    const Scalar x1 = 1.0);


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
 */
template<class Scalar>
class ModelEvaluator2DSim
  : public ::Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  ModelEvaluator2DSim(const Teuchos::RCP<const Epetra_Comm>& comm,
		      const Scalar d = 10.0,
		      const Scalar p0 = 2.0,
		      const Scalar p1 = 0.0,
		      const Scalar x0 = 0.0,
		      const Scalar x1 = 1.0);

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

  void set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >& W_factory);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RCP< ::Thyra::LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief . */
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const ::Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private: // data members

  Teuchos::RCP<const Epetra_Comm>  epetra_comm_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Epetra_Map>   x_epetra_map_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Epetra_Map>   f_epetra_map_;

  Teuchos::RCP<Epetra_CrsGraph>  W_graph_;
  
  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory_;

  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Scalar d_;
  Teuchos::RCP< ::Thyra::VectorBase<Scalar> > x0_;
  Teuchos::Array<Scalar> p_;
  bool showGetInvalidArg_;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

};

#include "ModelEvaluator2DSim_def.hpp"

#endif // NOX_THYRA_MODEL_EVALUATOR_2DSIM_DECL_HPP
