// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_MLLinearOp_hpp__
#define __Teko_MLLinearOp_hpp__

#include "Teko_Utilities.hpp"
#include "Teko_BlockImplicitLinearOp.hpp"

// forward declarations!
namespace ML_Epetra {
class MultiLevelPreconditioner;
}

namespace Teko {

namespace Epetra {
class MappingStrategy;
class EpetraOperatorWrapper;
}  // namespace Epetra

class MLLinearOp : public BlockImplicitLinearOp {
 public:
  MLLinearOp(const Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> &mlPrecOp);

  /** @brief Range space of this operator */
  virtual VectorSpace range() const { return productRange_; }

  /** @brief Domain space of this operator */
  virtual VectorSpace domain() const { return productDomain_; }

  /** @brief Perform a matrix vector multiply with this operator.
   *
   * The <code>apply</code> function takes one vector as input
   * and applies the inverse \f$ LDU \f$ decomposition. The result
   * is returned in \f$y\f$. If this operator is reprsented as \f$M\f$ then
   * \f$ y = \alpha M x + \beta y \f$ (ignoring conjugation!).
   *
   * @param[in]     x
   * @param[in,out] y
   * @param[in]     alpha (default=1)
   * @param[in]     beta  (default=0)
   */
  virtual void implicitApply(const BlockedMultiVector &x, BlockedMultiVector &y,
                             const double alpha = 1.0, const double beta = 0.0) const;

  virtual void describe(Teuchos::FancyOStream &out_arg,
                        const Teuchos::EVerbosityLevel verbLevel) const;

  Teuchos::RCP<const ML_Epetra::MultiLevelPreconditioner> getMLPreconditioner() const;
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> getMLPreconditioner();

 protected:
  void extractConversionInformation(ML_Epetra::MultiLevelPreconditioner &mlPrec);

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productRange_;  ///< Range vector space.
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productDomain_;  ///< Domain vector space.

  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> mlPrecOp_;
  Teuchos::RCP<Epetra::EpetraOperatorWrapper> Amat_;
  Teuchos::RCP<const Epetra::MappingStrategy>
      mappingStrategy_;  /// Convert Epetra_Vectors to blocked vectors

  mutable Teuchos::RCP<Epetra_MultiVector> eX_, eY_;  // storage to avoid repeated reallocation

 private:
  // hide me!
  MLLinearOp();
  MLLinearOp(const MLLinearOp &);
};

Teuchos::RCP<const ML_Epetra::MultiLevelPreconditioner> getMLPreconditioner(
    const Teko::LinearOp &lo);

}  // namespace Teko

#endif
