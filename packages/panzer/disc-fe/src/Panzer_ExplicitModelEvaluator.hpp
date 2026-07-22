// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EXPLICIT_MODEL_EVALUATOR_HPP
#define PANZER_EXPLICIT_MODEL_EVALUATOR_HPP

#include "PanzerDiscFE_config.hpp"

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"

#include "Panzer_ModelEvaluator.hpp"

#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_ModelEvaluator_Epetra.hpp"
#endif

#include "Panzer_MassMatrixModelEvaluator.hpp"

namespace panzer {

/** This is a model evaluator decorator that will take an implicit model evaluator
  * and make it explicit. If the model evaluator is not a panzer model evaluator
  * then there may be problems with constructing the dirichlet conditions in the
  * mass matrix. However, for a pzner model evaluator this has been taken care of.
  */
template<typename Scalar>
class ExplicitModelEvaluator
  : public Thyra::ModelEvaluatorDelegatorBase<Scalar>,
    public panzer::MassMatrixModelEvaluator<Scalar> {
public:

  /** \name Constructors/Initializers/Accessors */
  //@{

  /** Take in a Thyra model evaluator and then turn it into an explicit model evaluator. Assume the
    * mass matrix is constant unless the user specifies otherwise.
    */
  ExplicitModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > & model,
                         bool constantMassMatrix,
                         bool useLumpedMass,
                         bool applyMassInverse=true);

  //@}

  //! Build the nominal values, modifies the underlying models in args slightly
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;

  //! Build the in args, modifies the underlying models in args slightly
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //! Build the out args, modifies the underlying models in args slightly
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;

  //! Get the underlying panzer::ModelEvaluator
  Teuchos::RCP<panzer::ModelEvaluator<Scalar> > getPanzerUnderlyingModel();

  virtual void applyInverseMassMatrix(const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > input, const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > output) const
  {
    Thyra::apply(*invMassMatrix_,Thyra::NOTRANS,*input,output.ptr());
  }

  virtual void applyMassMatrix(const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > input, const Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > output) const
  {
    Thyra::apply(*mass_,Thyra::NOTRANS,*input,output.ptr());
  }

private: // data members

  /** Apply the dirichlet boundary conditions to the vector "f" using the
    * "x" values as the current solution.
    */
  void applyDirichletBCs(const Teuchos::RCP<Thyra::VectorBase<Scalar> > & x,
                         const Teuchos::RCP<Thyra::VectorBase<Scalar> > & f) const
  {
    if(panzerModel_!=Teuchos::null)       { panzerModel_->applyDirichletBCs(x,f); return; }
#ifdef PANZER_HAVE_EPETRA_STACK
    if(panzerEpetraModel_!=Teuchos::null) { panzerEpetraModel_->applyDirichletBCs(x,f); return; }
#endif

    TEUCHOS_ASSERT(false);
  }

  /** This method builds the inverse mass matrix from the underlying model evaluator.
    * Not that this is constant method that modifies a mutable member.
    */
  void buildInverseMassMatrix(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs) const;

  //! Build prototype in/out args
  void buildArgsPrototypes();

  /** \name Private functions overridden from ModelEvaulatorDelegatorBase. */
  //@{

  //! Evaluation of model, applies the mass matrix to the RHS
  void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
                     const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  //@}

  /** Set one time dirichlet beta using the underlying model evaluator.
    * If this model evaluator is a MEDelegator then this method is called
    * recursively, until the call succeeds (finding a panzer::ME) or fails
    * because a delegator is no longer used. If it fails an exception is thrown.
    * Note: The me used in this recursion is constant. This is consistent with the
    * one time dirichlet beta call in the model evaluators.
    */
  void setOneTimeDirichletBeta(double beta,const Thyra::ModelEvaluator<Scalar> & me) const;

  //! Is the mass matrix constant
  bool constantMassMatrix_;

  //! Use mass lumping, or a full solve
  bool massLumping_;

  //! Access to the panzer model evaluator pointer (thyra version)
  Teuchos::RCP<const panzer::ModelEvaluator<Scalar> > panzerModel_;

#ifdef PANZER_HAVE_EPETRA_STACK
  //! Access to the epetra panzer model evaluator pointer
  Teuchos::RCP<const panzer::ModelEvaluator_Epetra> panzerEpetraModel_;
#endif

  mutable Teuchos::RCP<Thyra::LinearOpBase<Scalar> > mass_;
  mutable Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > invMassMatrix_;
  mutable Teuchos::RCP<Thyra::VectorBase<Scalar> > scrap_f_;
  mutable Teuchos::RCP<Thyra::VectorBase<Scalar> > zero_;

  // build prototype in/out args
  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  // hide the default constructor
  ExplicitModelEvaluator();

};

} // end namespace panzer

#endif
