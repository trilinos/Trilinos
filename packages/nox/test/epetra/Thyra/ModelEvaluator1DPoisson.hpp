// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_THYRA_MODEL_EVALUATOR_1DFEM_DECL_HPP
#define NOX_THYRA_MODEL_EVALUATOR_1DFEM_DECL_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

template<class Scalar> class ModelEvaluator1DPoisson;

/** \brief Nonmember constuctor.
 *
 * \relates ModelEvaluator1DPoisson
 */
template<class Scalar>
Teuchos::RCP<ModelEvaluator1DPoisson<Scalar> >
modelEvaluator1DPoisson(const Teuchos::RCP<const Epetra_Comm>& comm,
            const int NumGlobalElements,
            const Scalar z_min,
            const Scalar z_max);


/** \brief 1D Finite Element model for Poisson's equation
 *
 * The equation modeled is:

 \verbatim

   d2u
   --- = 0
   dx2

   subject to:
      u = a  @ x = 0
      u = b @ x = N

 \endverbatim

 * The Matrix <tt>W = d(f)/d(x)</tt> is implemented as a
 * <tt>Thyra::MultiVectorBase</tt> object
 */
template<class Scalar>
class ModelEvaluator1DPoisson
  : public ::Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  ModelEvaluator1DPoisson(const Teuchos::RCP<const Epetra_Comm>& comm,
              const int num_global_elements,
              const Scalar a,
              const Scalar b);

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
  //@}

private:

  /** Allocates and returns the Jacobian matrix graph */
  virtual Teuchos::RCP<Epetra_CrsGraph> createGraph();

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
  const Scalar a_;
  const Scalar b_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > x_space_;
  Teuchos::RCP<const Epetra_Map>   x_owned_map_;
  Teuchos::RCP<const Epetra_Map>   x_ghosted_map_;
  Teuchos::RCP<const Epetra_Import> importer_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<Scalar> > f_space_;
  Teuchos::RCP<const Epetra_Map>   f_owned_map_;

  Teuchos::RCP<Epetra_CrsGraph>  W_graph_;

  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory_;

  Teuchos::RCP<Epetra_Vector> node_coordinates_;
  Teuchos::RCP<Epetra_Vector> ghosted_node_coordinates_;

  mutable Teuchos::RCP<Epetra_Vector> u_ptr;
  mutable Teuchos::RCP<Epetra_Vector> x_ptr;

  mutable Teuchos::RCP<Epetra_Vector> J_diagonal_;

  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP< ::Thyra::VectorBase<Scalar> > x0_;
  Teuchos::Array<Scalar> p_;
  bool showGetInvalidArg_;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

};

//==================================================================
#include "ModelEvaluator1DPoisson_def.hpp"
//==================================================================

#endif
