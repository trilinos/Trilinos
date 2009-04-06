#ifndef __PB_LU2x2InverseOp_hpp__
#define __PB_LU2x2InverseOp_hpp__

#include "PB_Utilities.hpp"
#include "PB_BlockImplicitLinearOp.hpp"

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
class LU2x2InverseOp : public BlockImplicitLinearOp {
public:
   /** \brief This constructor explicitly takes the parts of \f$ A \f$ required to
     *        build the inverse operator.
     *
     * This constructor explicitly takes the parts of \f$ A \f$ required to build
     * the inverse operator. 
     *
     * \param[in] A The block \f$ 2 \times 2 \f$ \f$A\f$ operator.
     * \param[in] invA00  An approximate inverse of \f$ A_{00} \f$.
     * \param[in] invS  An approximate inverse of \f$ S = -A_{11} + A_{10} A_{00}^{-1} A_{01} \f$.
     */
   LU2x2InverseOp(const BlockedLinearOp & A,
                       const LinearOp & invA00,
                       const LinearOp & invS);

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
   virtual void apply(const BlockedMultiVector & x, BlockedMultiVector & y,
              const double alpha = 1.0, const double beta = 0.0) const;
   //@}

protected:
   // fundamental operators to use
   const BlockedLinearOp A_;  ///< operator \f$ A \f$
   const LinearOp invA00_;    ///< inverse of \f$ A_{00} \f$
   const LinearOp invS_;      ///< inverse of \f$ S \f$

   // some blocks of A
   const LinearOp A10_;       ///< operator \f$ A_{10} \f$
   const LinearOp A01_;       ///< operator \f$ A_{01} \f$

   Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > productRange_; ///< Range vector space.
   Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > productDomain_; ///< Domain vector space.

private:
   // hide me!
   LU2x2InverseOp();
   LU2x2InverseOp(const LU2x2InverseOp &);
};

inline LinearOp createNewLU2x2InverseOp(BlockedLinearOp & A,LinearOp & invA00,LinearOp & invS)
{
   return Teuchos::rcp(new LU2x2InverseOp(A,invA00,invS));
}

} // end namespace PB

#endif	
