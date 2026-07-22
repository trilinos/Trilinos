// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_StridedEpetraOperator_hpp__
#define __Teko_StridedEpetraOperator_hpp__

// Epetra includes
#include "Epetra_Operator.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

// Teko includes
#include "Teko_BlockedReordering.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"
#include "Teko_StridedMappingStrategy.hpp"

namespace Teko {
namespace Epetra {

class StridedEpetraOperator : public EpetraOperatorWrapper {
 public:
  enum eNormType { Inf, One, Frobenius };

  StridedEpetraOperator(int numVars, const Teuchos::RCP<const Epetra_Operator> &content,
                        const std::string &label = "<ANYM>");
  StridedEpetraOperator(const std::vector<int> &vars,
                        const Teuchos::RCP<const Epetra_Operator> &content,
                        const std::string &label = "<ANYM>");

  virtual void SetContent(const std::vector<int> &vars,
                          const Teuchos::RCP<const Epetra_Operator> &content);

  virtual void RebuildOps() { BuildBlockedOperator(); }

  virtual const Teuchos::RCP<const Epetra_Operator> GetContent() const { return fullContent_; }

  // virtual const Teuchos::RCP<Epetra_Operator> GetContent()
  // { return fullContent_; }

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

  /** Print the Norm of the sub matrces. The type of norm
   * is specified by the argument.
   *
   * \param[in] nrmType Type of norm to use
   * \param[in] newline Character to use when a new line
   *                    is needed.
   */
  virtual std::string PrintNorm(const eNormType &nrmType = Frobenius, const char newline = '\n');

  // functions overloading Epetra_Operator
  ////////////////////////////////////////////////

  // destructor
  virtual ~StridedEpetraOperator() {}

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
  Teuchos::RCP<StridedMappingStrategy> stridedMapping_;
  Teuchos::RCP<Thyra::LinearOpBase<double> > stridedOperator_;
  Teuchos::RCP<const BlockReorderManager> reorderManager_;

  std::string label_;

  void BuildBlockedOperator();
};

}  // end namespace Epetra
}  // end namespace Teko

#endif
