// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_TpetraInverseFactoryOperator_hpp__
#define __Teko_TpetraInverseFactoryOperator_hpp__

#include "Teuchos_ConstNonconstObjectContainer.hpp"

#include "Teko_ConfigDefs.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_TpetraInverseOpWrapper.hpp"

namespace Teko {
namespace TpetraHelpers {

/** \brief A single Epetra wrapper for all operators constructed
 *        from an inverse operator.
 *
 * This class uses the Teko inverse factories to
 * build an Epetra_Operator that behaves like the inverse operatotr.
 * This is done by using the InverseFactory, and letting
 * it build whatever operator is neccessary. Thus the Epetra
 * "layer" is just a single class that handles any generic
 * InverseFactory.
 */
class InverseFactoryOperator : public TpetraInverseOpWrapper {
 public:
  /** \brief Constructor that takes the InverseFactory that will
   *        build the operator.
   *
   * Constructor that takes the InverseFactory that will
   * build the operator.
   */
  InverseFactoryOperator(const Teuchos::RCP<const InverseFactory> &bfp);

  /** \brief Build the underlying data structure for the inverse operator.
   *
   * Build the underlying data structure for the inverse operator. This
   * permits the manipulation of the state object for an inverse operator.
   *
   * \param[in] clearOld If true any previously constructed
   *                     operator will be wiped out and
   *                     a new one created. If false, anoperator
   *                     will be created only if the current one is
   *                     empty (i.e. <code>initPreconditioner</code>
   *                     had not been called).
   */
  virtual void initInverse(bool clearOld = false);

  /** \brief Build this inverse operator from an Epetra_Operator
   * passed in to this object.
   *
   * Build this inverse opeerator from an Epetra_Operator
   * passed in to this object. If this Epetra_Operator
   * is an EpetraOperatorWrapper object then the block Thyra components
   * are extracted.
   *
   * \param[in] A The Epetra source operator.
   * \param[in] clear If true, than any previous state saved by the operator
   *                  is discarded.
   */
  virtual void buildInverseOperator(const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A,
                                    bool clear = true);

  /** \brief Build this inverse operator from an Epetra_Operator
   * passed in to this object.
   *
   * Build this inverse opeerator from an Epetra_Operator
   * passed in to this object. If this Epetra_Operator
   * is an EpetraOperatorWrapper object then the block Thyra components
   * are extracted.
   *
   * \param[in] A The Epetra source operator.
   * \param[in] clear If true, than any previous state saved by the operator
   *                  is discarded.
   */
  virtual void buildInverseOperator(const Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> > &A,
                                    bool clear = true);

  /** \brief Rebuild this inverse from an Epetra_Operator passed
   * in this to object.
   *
   * Rebuild this inverse from an Epetra_Operator passed
   * in this to object.  If <code>buildInverseOperator</code> has not been called
   * the inverse operator will be built instead. Otherwise efforts are taken
   * to only rebuild what is neccessary. Also, that this Epetra_Operator
   * may be an EpetraOperatorWrapper object, so the block Thyra components
   * can be extracted.
   *
   * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
   */
  virtual void rebuildInverseOperator(
      const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A);

  /** \brief Rebuild this inverse from an Epetra_Operator passed
   * in this to object.
   *
   * Rebuild this inverse from an Epetra_Operator passed
   * in this to object.  If <code>buildInverseOperator</code> has not been called
   * the inverse operator will be built instead. Otherwise efforts are taken
   * to only rebuild what is neccessary. Also, that this Epetra_Operator
   * may be an EpetraOperatorWrapper object, so the block Thyra components
   * can be extracted.
   *
   * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
   */
  virtual void rebuildInverseOperator(const Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> > &A);

  /** Extract the forward op used by <code>buildInverseOperator</code>
   * or <code>rebuildInverseOperator</code>.
   */
  Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > getForwardOp() const {
    return fwdOp_.getConstObj();
  }

  /** Extract the forward op used by <code>buildInverseOperator</code>
   * or <code>rebuildInverseOperator</code>.
   */
  Teuchos::RCP<Tpetra::Operator<ST, LO, GO, NT> > getNonconstForwardOp() const {
    return fwdOp_.getNonconstObj();
  }

 protected:
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > extractLinearOp(
      const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A) const;
  Teuchos::RCP<const MappingStrategy> extractMappingStrategy(
      const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A) const;

  InverseFactoryOperator();
  InverseFactoryOperator(const InverseFactoryOperator &);

  Teuchos::RCP<const Teko::InverseFactory> inverseFactory_;
  Teko::ModifiableLinearOp invOperator_;
  bool firstBuildComplete_;

  Teuchos::ConstNonconstObjectContainer<Tpetra::Operator<ST, LO, GO, NT> > fwdOp_;
  bool setConstFwdOp_;
};

}  // end namespace TpetraHelpers
}  // end namespace Teko

#endif
