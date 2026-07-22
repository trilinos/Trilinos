// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_BlockedEpetraOperator_hpp__
#define __Teko_BlockedEpetraOperator_hpp__

// Epetra includes
#include "Epetra_Operator.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

// Teko includes
#include "Teko_BlockedReordering.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"
#include "Teko_BlockedMappingStrategy.hpp"

namespace Teko {
namespace Epetra {

/** \brief Tear about a user specified Epetra_Operator (CrsMatrix)
 *        using a vector of vectors of GIDs for each block.
 */
class BlockedEpetraOperator : public EpetraOperatorWrapper {
 public:
  /** Build a blocked operator based on a vector of vector
   * of global IDs.
   *
   * \param[in] vars Vector of vectors of global ids specifying
   *                 how the operator is to be blocked.
   * \param[in] content Operator to be blocked
   * \param[in] label Label for name the operator
   */
  BlockedEpetraOperator(const std::vector<std::vector<int> > &vars,
                        const Teuchos::RCP<const Epetra_Operator> &content,
                        const std::string &label = "<ANYM>");

  /** Build a blocked operator based on a vector of vector
   * of global IDs. This function basically sets up the mapping
   * strategy used by this operator.
   *
   * \param[in] vars Vector of vectors of global ids specifying
   *                 how the operator is to be blocked.
   * \param[in] content Operator to be blocked
   */
  virtual void SetContent(const std::vector<std::vector<int> > &vars,
                          const Teuchos::RCP<const Epetra_Operator> &content);

  /** Force a rebuild of the blocked operator from the stored
   * content operator.
   */
  virtual void RebuildOps() { BuildBlockedOperator(); }

  virtual const Teuchos::RCP<const Epetra_Operator> GetContent() const { return fullContent_; }

  virtual const Teuchos::RCP<const Epetra_Operator> GetContent() { return fullContent_; }

  const Teuchos::RCP<const Epetra_Operator> GetBlock(int i, int j) const;

  /** Use a reorder manager to block this operator as desired.
   * Multiple calls to the function reorder only the underlying object.
   */
  void Reorder(const BlockReorderManager &brm);

  //! Remove any reordering on this object
  void RemoveReording();

  /** Write out this operator to matrix market files
   */
  virtual void WriteBlocks(const std::string &prefix) const;

  // functions overloading Epetra_Operator
  ////////////////////////////////////////////////

  // destructor
  virtual ~BlockedEpetraOperator() {}

  // attribute set methods

  // don't use transpose...ever!
  virtual int SetUseTranspose(bool /* useTranspose */) { return -1; }

  virtual int ApplyInverse(const Epetra_MultiVector & /* X */, Epetra_MultiVector & /* Y */) const {
    TEUCHOS_ASSERT(false);
    return -1;
  }

  virtual double NormInf() const {
    TEUCHOS_ASSERT(false);
    return 0.0;
  }

  // attribute access functions
  virtual bool UseTranspose() const { return false; }
  virtual bool HasNormInf() const { return false; }
  virtual const Epetra_Comm &Comm() const { return fullContent_->Comm(); }

  //! Helps perform sanity checks
  bool testAgainstFullOperator(int count, double tol) const;

 protected:
  // gooey center of this shell
  Teuchos::RCP<const Epetra_Operator> fullContent_;
  Teuchos::RCP<BlockedMappingStrategy> blockedMapping_;
  Teuchos::RCP<Thyra::LinearOpBase<double> > blockedOperator_;
  Teuchos::RCP<const BlockReorderManager> reorderManager_;

  std::string label_;

  void BuildBlockedOperator();
};

}  // end namespace Epetra
}  // end namespace Teko

#endif
