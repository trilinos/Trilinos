// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file Teko_LU2x2InverseOp.hpp
 *
 * File that implements the inverse of a block 2x2 LU decomposition.
 */

#ifndef __Teko_LU2x2InverseOp_hpp__
#define __Teko_LU2x2InverseOp_hpp__

#include "Teko_Utilities.hpp"
#include "Teko_BlockImplicitLinearOp.hpp"

namespace Teko {

/** \brief This linear operator approximates the inverse
 *         of a block \f$ 2\times 2 \f$ operator using a
 *         block \f$ LDU \f$ decomposition.
 *
 * For a matrix that is blocked like
 *
 * \f$ A = \left[\begin{array}{cc}
 *           A_{00} & A_{01} \\
 *           A_{10} & A_{11}
 *           \end{array}\right] \f$
 *
 * this class evaluates the \f$A^{-1}\f$ given \f$A_{00}^{-1}\f$ and the inverse of
 * the Schur complement. The \f$ LDU \f$ factorization is defined as
 *
 * \f$
 * A = \left[ \begin{array}{cc}
 * I & 0  \\
 * A_{10} A_{00}^{-1} & I
 * \end{array} \right]
 * \left[ \begin{array}{cc}
 * A_{00} & 0  \\
 * 0 & -S
 * \end{array} \right]
 * \left[ \begin{array}{cc}
 * I &  A_{00}^{-1} A_{01} \\
 * 0 & I
 * \end{array} \right]\f$
 *
 * where the Schur complement is \f$ S=-A_{11}+A_{10} A_{00}^{-1} A_{01} \f$ .
 * In order to do this 2 evaluations of \f$ A_{00}^{-1} \f$ and a single
 * evalution of \f$ S^{-1} \f$ are needed. For increased flexibility both
 * evaluations of \f$A_{00}^{-1}\f$ can be specified independently.
 * For righthand side vector \f$[f, g]^T\f$ and solution vector \f$[u,v]^T\f$
 * the two inverses (\f$A\f$-hat and \f$A\f$-tilde) are needed to evaluate
 *
 * \f$\hat{A}_{00} u^* = f\f$,
 *
 * \f$\tilde{A}_{00} v = A_{01} v\f$
 *
 * where \f$u^*\f$ is an intermediate step.
 */
class LU2x2InverseOp : public BlockImplicitLinearOp {
 public:
  /** \brief This constructor explicitly takes the parts of \f$ A \f$ required to
   *        build the inverse operator.
   *
   * This constructor explicitly takes the parts of \f$ A \f$ required to build
   * the inverse operator.
   *
   * \param[in] A The block \f$ 2 \times 2 \f$ \f$A\f$ operator.
   * \param[in] invA00  An approximate inverse of \f$ A_{00} \f$, used for both \f$\hat{A}_{00}\f$
   * and \f$\tilde{A}_{00}\f$ \param[in] invS  An approximate inverse of \f$ S = -A_{11} + A_{10}
   * A_{00}^{-1} A_{01} \f$.
   */
  LU2x2InverseOp(const BlockedLinearOp &A, const LinearOp &invA00, const LinearOp &invS);

  /** \brief This constructor explicitly takes the parts of \f$ A \f$ required to
   *        build the inverse operator.
   *
   * This constructor explicitly takes the parts of \f$ A \f$ required to build
   * the inverse operator.
   *
   * \param[in] A The block \f$ 2 \times 2 \f$ \f$A\f$ operator.
   * \param[in] hatInvA00  An approximate inverse of \f$ \hat{A}_{00} \f$
   * \param[in] tildeInvA00  An approximate inverse of \f$ \tilde{A}_{00} \f$
   * \param[in] invS  An approximate inverse of \f$ S = -A_{11} + A_{10} A_{00}^{-1} A_{01} \f$.
   */
  LU2x2InverseOp(const BlockedLinearOp &A, const LinearOp &hatInvA00, const LinearOp &tildeInvA00,
                 const LinearOp &invS);

  //! \name Inherited methods from Thyra::LinearOpBase
  //@{

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
  //@}

  virtual void describe(Teuchos::FancyOStream &out_arg,
                        const Teuchos::EVerbosityLevel verbLevel) const;

 protected:
  using BlockImplicitLinearOp::implicitApply;

  // fundamental operators to use
  const BlockedLinearOp A_;     ///< operator \f$ A \f$
  const LinearOp hatInvA00_;    ///< inverse of \f$ A_{00} \f$
  const LinearOp tildeInvA00_;  ///< inverse of \f$ A_{00} \f$
  const LinearOp invS_;         ///< inverse of \f$ S \f$

  // some blocks of A
  const LinearOp A10_;  ///< operator \f$ A_{10} \f$
  const LinearOp A01_;  ///< operator \f$ A_{01} \f$

  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productRange_;  ///< Range vector space.
  Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >
      productDomain_;  ///< Domain vector space.

 private:
  // hide me!
  LU2x2InverseOp();
  LU2x2InverseOp(const LU2x2InverseOp &);
};

/** \brief Constructor method for building <code>LU2x2InverseOp</code>.
 *
 * Constructor method for building <code>LU2x2InverseOp</code>.
 *
 * \param[in] A      2x2 Operator to be decomposed
 * \param[in] invA00 Approximate inverse of the operators \f$(0,0)\f$ block.
 * \param[in] invS   Approximate inverse of the Schur complement
 *
 * \returns A linear operator that behaves like the inverse of the
 *          LU decomposition.
 *
 * \relates LU2x2InverseOp
 */
inline LinearOp createLU2x2InverseOp(BlockedLinearOp &A, LinearOp &invA00, LinearOp &invS) {
  return Teuchos::rcp(new LU2x2InverseOp(A, invA00, invS));
}

/** \brief Constructor method for building <code>LU2x2InverseOp</code>.
 *
 * Constructor method for building <code>LU2x2InverseOp</code>.
 *
 * \param[in] A      2x2 Operator to be decomposed
 * \param[in] invA00 Approximate inverse of the operators \f$(0,0)\f$ block.
 * \param[in] invS   Approximate inverse of the Schur complement
 * \param[in] str    String to label the operator
 *
 * \returns A linear operator that behaves like the inverse of the
 *          LU decomposition.
 *
 * \relates LU2x2InverseOp
 */
inline LinearOp createLU2x2InverseOp(BlockedLinearOp &A, LinearOp &invA00, LinearOp &invS,
                                     const std::string &str) {
  Teuchos::RCP<Thyra::LinearOpBase<double> > result =
      Teuchos::rcp(new LU2x2InverseOp(A, invA00, invS));
  result->setObjectLabel(str);

  return result;
}

/** \brief Constructor method for building <code>LU2x2InverseOp</code>.
 *
 * Constructor method for building <code>LU2x2InverseOp</code>.
 *
 * \param[in] A           2x2 Operator to be decomposed
 * \param[in] hatInvA00   First approximate inverse of the operators \f$(0,0)\f$ block.
 * \param[in] tildeInvA00 Second approximate inverse of the operators \f$(0,0)\f$ block.
 * \param[in] invS        Approximate inverse of the Schur complement
 *
 * \returns A linear operator that behaves like the inverse of the
 *          LU decomposition.
 *
 * \relates LU2x2InverseOp
 */
inline LinearOp createLU2x2InverseOp(BlockedLinearOp &A, LinearOp &hatInvA00, LinearOp &tildeInvA00,
                                     LinearOp &invS) {
  return Teuchos::rcp(new LU2x2InverseOp(A, hatInvA00, tildeInvA00, invS));
}

/** \brief Constructor method for building <code>LU2x2InverseOp</code>.
 *
 * Constructor method for building <code>LU2x2InverseOp</code>.
 *
 * \param[in] A           2x2 Operator to be decomposed
 * \param[in] hatInvA00   First approximate inverse of the operators \f$(0,0)\f$ block.
 * \param[in] tildeInvA00 Second approximate inverse of the operators \f$(0,0)\f$ block.
 * \param[in] invS        Approximate inverse of the Schur complement
 * \param[in] str         String to label the operator
 *
 * \returns A linear operator that behaves like the inverse of the
 *          LU decomposition.
 *
 * \relates LU2x2InverseOp
 */
inline LinearOp createLU2x2InverseOp(BlockedLinearOp &A, LinearOp &hatInvA00, LinearOp &tildeInvA00,
                                     LinearOp &invS, const std::string &str) {
  Teuchos::RCP<Thyra::LinearOpBase<double> > result =
      Teuchos::rcp(new LU2x2InverseOp(A, hatInvA00, tildeInvA00, invS));
  result->setObjectLabel(str);

  return result;
}

}  // end namespace Teko

#endif
