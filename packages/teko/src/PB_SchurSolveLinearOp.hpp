#ifndef __PB_SchurSolveLinearOp_hpp__
#define __PB_SchurSolveLinearOp_hpp__

#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"

namespace PB {

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
 * evalution of \f$ S^{-1} \f$ are needed.
 */
class SchurSolveLinearOp : public Thyra::LinearOpBase<double> {
public:
   /** @brief This constructor explicitly takes the parts of \f$ A \f$ required to
     *        build the inverse operator.
     *
     * This constructor explicitly takes the parts of \f$ A \f$ required to build
     * the inverse operator. 
     *
     * @param[in] A The block \f$ 2 \times 2 \f$ \f$A\f$ operator.
     * @param[in] invA00  An approximate inverse of \f$ A_{00} \f$.
     * @param[in] invS  An approximate inverse of \f$ S = -A_{11} + A_{10} A_{00}^{-1} A_{01} \f$.
     */
   SchurSolveLinearOp(const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > & A,
                      const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invA00,
                      const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invS);

   /** @name Inherited methods from Thyra::LinearOpBase */
   //@{

   /** @brief Range space of this operator */
   Teuchos::RCP< const Thyra::VectorSpaceBase<double> > range() const { return productRange_; }

   /** @brief Domain space of this operator */
   Teuchos::RCP< const Thyra::VectorSpaceBase<double> > domain() const { return productDomain_; }

   /** @brief Perform a matrix vector multiply with this operator. 
     *
     * The <code>apply</code> function takes one vector as input 
     * and applies the inverse \f$ LDU \f$ decomposition. The result
     * is returned in \f$y\f$. If this operator is reprsented as \f$M\f$ then
     * \f$ y = \alpha M x + \beta y \f$ (ignoring conjugation!).
     *
     * @param[in]     conj (Must be set to Thyra::NONCONJ_ELE)
     * @param[in]     x 
     * @param[in,out] y 
     * @param[in]     alpha (default=1)
     * @param[in]     beta  (default=0)
     */
   void apply(const Thyra::EConj conj, const Thyra::MultiVectorBase<double> & x, Thyra::MultiVectorBase<double> * y,
      const double alpha = Teuchos::ScalarTraits<double>::one(),
      const double beta = Teuchos::ScalarTraits<double>::zero()) const;
   //@}

protected:
   // fundamental operators to use
   const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > A_;  ///< operator \f$ A \f$
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invA00_;    ///< inverse of \f$ A_{00} \f$
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invS_;      ///< inverse of \f$ S \f$

   // some blocks of A
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > A10_;       ///< operator \f$ A_{10} \f$
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > A01_;       ///< operator \f$ A_{01} \f$

   Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > productRange_; ///< Range vector space.
   Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > productDomain_; ///< Domain vector space.

private:
   // hide me!
   SchurSolveLinearOp();
   SchurSolveLinearOp(const SchurSolveLinearOp &);
};

} // end namespace PB

#endif	
