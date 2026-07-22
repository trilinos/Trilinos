// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SIMPLE_2D_TPETRA_MODEL_EVALUATOR_DECL_HPP
#define SIMPLE_2D_TPETRA_MODEL_EVALUATOR_DECL_HPP


#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Vector.hpp"


/** \brief Simple 2d simulation only ModelEvaluator for f(x) = 0 using Tpetra
 * objects.
 *
 * The equations modeled are:

 \verbatim

    f[0] =       x[0]      + x[1]*x[1] - p[0];
    f[1] = d * ( x[0]*x[0] - x[1]      - p[1] );

 \endverbatim

 * The Matrix <tt>W_op = d(f)/d(x)</tt> is implemented as a
 * <tt>Thyra::TpetraLinearOp</tt> object and all of the other objects are
 * Thyra wrappers for Tpetra objects.
 */
template<class Scalar>
class Simple2DTpetraModelEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/Initializers/Accessors */
  //@{

  /** \brief . */
  Simple2DTpetraModelEvaluator();

  /** \brief . */
  void set_d(const Scalar &d);

  /** \brief . */
  void set_p(const Teuchos::ArrayView<const Scalar> &p);

  /** \brief . */
  void set_x0(const Teuchos::ArrayView<const Scalar> &x0);

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
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

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
  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Scalar d_;
  Teuchos::RCP<Tpetra::Vector<Scalar> > x0_;
  Teuchos::Array<Scalar> p_;
  Teuchos::RCP<Tpetra::CrsGraph<> > W_op_graph_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

};


/** \brief Non-member constructor.
 *
 * \relates Simple2DTpetraModelEvaluator
 */
template<class Scalar>
Teuchos::RCP<Simple2DTpetraModelEvaluator<Scalar> >
simple2DTpetraModelEvaluator()
{
  return Teuchos::rcp(new Simple2DTpetraModelEvaluator<Scalar>);
}


#endif // SIMPLE_2D_TPETRA_MODEL_EVALUATOR_DECL_HPP
