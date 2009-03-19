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
 *           F & U \\
 *           L & D
 *           \end{array}\right] \f$
 *
 * this class evaluates the inverse given inv(F) and the inverse of
 * the Schur complement. The \f$ LDU \f$ is defined as
 *
 * \f$
 * A = \left[ \begin{array}{cc}
 * I & 0  \\
 * L F^{-1} & I
 * \end{array} \right]
 * \left[ \begin{array}{cc}
 * F & 0  \\
 * 0 & -S
 * \end{array} \right]
 * \left[ \begin{array}{cc}
 * I &  F^{-1} U \\
 * 0 & I
 * \end{array} \right]\f$
 *
 * where the Schur complement is \f$ S=-D+L F^{-1} U \f$ .
 * In order to do this 2 evaluations of \f$ F^{-1} \f$ and a single
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
     * @param[in] fwdOp The block \f$ 2 \times 2 \f$ \f$A\f$ operator.
     * @param[in] invF  An approximate inverse of \f$ F \f$.
     * @param[in] invS  An approximate inverse of \f$ S = -D + L F^{-1} U \f$.
     */
   SchurSolveLinearOp(const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > & fwdOp,
                      const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
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
     * is returned in Y. If this operator is reprsented as \f$M\f$ then
     * \f$ Y = \alpha M X + \beta Y \f$ (ignoring conjugation!).
     *
     * @param[in]     conj
     * @param[in]     X 
     * @param[in,out] Y 
     * @param[in]     alpha (default=1)
     * @param[in]     beta  (default=0)
     */
   void apply(const Thyra::EConj conj, const Thyra::MultiVectorBase<double> &X, Thyra::MultiVectorBase<double> *Y,
      const double alpha = Teuchos::ScalarTraits<double>::one(),
      const double beta = Teuchos::ScalarTraits<double>::zero()) const;
   //@}

protected:
   // fundamental operators to use
   const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > fwdOp_; ///< operator \f$ A \f$
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invF_;         ///< inverse of \f$ F \f$
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invS_;         ///< inverse of \f$ S \f$

   // some blocks of A
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > fwdL_;         ///< operator \f$ L = A(1,0) \f$
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > fwdU_;         ///< operator \f$ U = A(0,1) \f$

   Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > productRange_; ///< Range vector space.
   Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > productDomain_; ///< Domain vector space.

private:
   // hide me!
   SchurSolveLinearOp();
   SchurSolveLinearOp(const SchurSolveLinearOp &);
};

} // end namespace PB

#endif	
