// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_THYRA_MODEL_EVALUATOR_HEQ_DECL_HPP
#define NOX_THYRA_MODEL_EVALUATOR_HEQ_DECL_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include <vector>

template<class Scalar> class ModelEvaluatorHeq;

/** \brief Nonmember constuctor.
 *
 * \relates ModelEvaluatorHeq
 */
template<class Scalar>
Teuchos::RCP<ModelEvaluatorHeq<Scalar> >
modelEvaluatorHeq(const Teuchos::RCP<const Epetra_Comm>& comm,
            const int NumGlobalElements,
            const Scalar paramC);


/** \brief Chandrasekhar H-Equation
 *
 * The equation modeled is:

  \f$ F(x)_i = x_i - (1 - \frac{c}{2N}\sum_{j=1}^N \frac{\mu_ix_j}{\mu_i + \mu_j})^{-1} = 0, 1 <= i <= N\f$

 * In this, \f$[0,1]\f$ is discretized into \f$N\f$ uniformly sized
 * cells and \f$\{\mu_j\}\f$ are the midpoints of these cells.
 */
template<class Scalar>
class ModelEvaluatorHeq
  : public ::Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  ModelEvaluatorHeq(const Teuchos::RCP<const Epetra_Comm>& comm,
              const int num_global_elements,
              const Scalar paramC);

  /** \name Initializers/Accessors */
  //@{

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
  /** \brief . */
  Teuchos::RCP< ::Thyra::PreconditionerBase< Scalar > > create_W_prec() const;
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

  const Teuchos::RCP<const Epetra_Comm>  comm_;
  const int num_global_elements_;
  const Scalar paramC_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Epetra_Map>   x_map_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Epetra_Map>   f_map_;

  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory_;

  Teuchos::RCP<Epetra_CrsMatrix> A_heq_;

  mutable Teuchos::RCP<const Epetra_Vector> x_ptr_;

  std::vector<Scalar> node_coordinates_;

  Teuchos::RCP<Epetra_Vector> ones_;

  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP< ::Thyra::VectorBase<Scalar> > x0_;
  bool showGetInvalidArg_;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

};

//==================================================================
#include "ModelEvaluatorHeq_def.hpp"
//==================================================================

#endif
