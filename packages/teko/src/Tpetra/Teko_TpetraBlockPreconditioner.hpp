// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_TpetraBlockPreconditioner_hpp__
#define __Teko_TpetraBlockPreconditioner_hpp__

// #include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_TpetraInverseOpWrapper.hpp"
#include "Teko_ConfigDefs.hpp"

namespace Teko {
namespace TpetraHelpers {

/** \brief A single Epetra wrapper for all the BlockPreconditioners.
 *
 * This class uses the Thyra based preconditioner factories to
 * build an Epetra_Operator that behaves like a preconditioner.
 * This is done by using the BlockPreconditionerFactory, and letting
 * it build whatever preconditioner is neccessary. Thus the Epetra
 * "layer" is just a single class that handles any generic
 * BlockPreconditionerFactory.
 */
class TpetraBlockPreconditioner : public TpetraInverseOpWrapper {
 public:
  /** \brief Constructor that takes the BlockPreconditionerFactory that will
   *        build the preconditioner.
   *
   * Constructor that takes the BlockPreconditionerFactory that will
   * build the preconditioner.
   */
  TpetraBlockPreconditioner(const Teuchos::RCP<const PreconditionerFactory> &bfp);

  /** \brief Build the underlying data structure for the preconditioner.
   *
   * Build the underlying data structure for the preconditioner. This
   * permits the manipulation of the state object for a preconditioner.
   * and is useful in that case some extra data needs to fill the
   * preconditioner state.
   *
   * \param[in] clearOld If true any previously constructed
   *                     preconditioner will be wiped out and
   *                     a new one created. If false, a preconditioner
   *                     will be create only if the current one is
   *                     empty (i.e. <code>initPreconditioner</code>
   *                     had not been called).
   */
  virtual void initPreconditioner(bool clearOld = false);

  /** \brief Build this preconditioner from an Epetra_Operator
   * passed in to this object.
   *
   * Build this preconditioner from an Epetra_Operator
   * passed in to this object. It is assumed that this Epetra_Operator
   * will be a EpetraOperatorWrapper object, so the block Thyra components
   * can be easily extracted.
   *
   * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
   * \param[in] clear If true, than any previous state saved by the preconditioner
   *                  is discarded.
   */
  virtual void buildPreconditioner(const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A,
                                   bool clear = true);

  /** \brief Build this preconditioner from an Epetra_Operator
   * passed in to this object. It is assumed that this Epetra_Operator
   *
   * Build this preconditioner from an Epetra_Operator
   * passed in to this object. It is assumed that this Epetra_Operator
   * will be a EpetraOperatorWrapper object, so the block Thyra components
   * can be easily extracted.
   *
   * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
   * \param[in] mv A vector that was used to build the source operator.
   * \param[in] clear If true, than any previous state saved by the preconditioner
   *                  is discarded.
   */
  virtual void buildPreconditioner(const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A,
                                   const Tpetra::MultiVector<ST, LO, GO, NT> &mv,
                                   bool clear = true);

  /** \brief Rebuild this preconditioner from an Epetra_Operator passed
   * in this to object.
   *
   * Rebuild this preconditioner from an Epetra_Operator passed
   * in this to object.  If <code>buildPreconditioner</code> has not been called
   * the preconditioner will be built instead. Otherwise efforts are taken
   * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
   * will be an EpetraOperatorWrapper object, so the block Thyra components
   * can be easily extracted.
   *
   * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
   */
  virtual void rebuildPreconditioner(
      const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A);

  /** \brief Rebuild this preconditioner from an Epetra_Operator passed
   * in this to object.
   *
   * Rebuild this preconditioner from an Epetra_Operator passed
   * in this to object.  If <code>buildPreconditioner</code> has not been called
   * the preconditioner will be built instead. Otherwise efforts are taken
   * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
   * will be an EpetraOperatorWrapper object, so the block Thyra components
   * can be easily extracted.
   *
   * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
   * \param[in] mv A vector that was used to build the source operator.
   */
  virtual void rebuildPreconditioner(const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A,
                                     const Tpetra::MultiVector<ST, LO, GO, NT> &mv);

  /** Try to get a <code>Teko::PreconditionerState</code> object. This method
   * attempts to cast its internal representation of a preconditioner
   * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a
   * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
   *
   * \returns Get the state object associated with this preconditioner.
   *          If it doesn't exist for this type of preconditioner factory
   *          this method returns null.
   */
  virtual Teuchos::RCP<PreconditionerState> getPreconditionerState();

  /** Try to get a <code>Teko::PreconditionerState</code> object. This method
   * attempts to cast its internal representation of a preconditioner
   * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a
   * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
   *
   * \returns Get the state object associated with this preconditioner.
   *          If it doesn't exist for this type of preconditioner factory
   *          this method returns null.
   */
  virtual Teuchos::RCP<const PreconditionerState> getPreconditionerState() const;

 protected:
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > extractLinearOp(
      const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A) const;
  Teuchos::RCP<const MappingStrategy> extractMappingStrategy(
      const Teuchos::RCP<const Tpetra::Operator<ST, LO, GO, NT> > &A) const;

  TpetraBlockPreconditioner();
  TpetraBlockPreconditioner(const TpetraBlockPreconditioner &);

  // Teuchos::RCP<const BlockPreconditionerFactory> preconFactory_;
  Teuchos::RCP<const PreconditionerFactory> preconFactory_;
  Teuchos::RCP<Thyra::PreconditionerBase<ST> > preconObj_;
  bool firstBuildComplete_;
};

}  // end namespace TpetraHelpers
}  // end namespace Teko

#endif
