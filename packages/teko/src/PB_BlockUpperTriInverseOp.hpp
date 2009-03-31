#ifndef __PB_BlockUpperTriInverseOp_hpp__
#define __PB_BlockUpperTriInverseOp_hpp__

#include "PB_Utilities.hpp"
#include "PB_BlockImplicitLinearOp.hpp"

namespace PB {

/** \brief This linear operator computes the inverse
  *        of a upper triangular matrix.
  *
  * This linear operator computes the inverse
  * of an upper triangular matrix. This requires,
  * the upper triangular blocks, as well as the 
  * inverse of the operators on the diagonal.
  */
class BlockUpperTriInverseOp : public BlockImplicitLinearOp {
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
   BlockUpperTriInverseOp(BlockedLinearOp & U,const std::vector<LinearOp> & invDiag);

   //! @name Inherited methods from Thyra::LinearOpBase
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

   virtual void describe(Teuchos::FancyOStream &out_arg,
                         const Teuchos::EVerbosityLevel verbLevel) const;

protected:
   // fundamental operators to use
   const BlockedLinearOp U_;  ///< operator \f$ U \f$
   std::vector<LinearOp> invDiag_; ///< (Approximate) Inverses of the diagonal operators

   Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > productRange_; ///< Range vector space.
   Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > productDomain_; ///< Domain vector space.

private:
   // hide me!
   BlockUpperTriInverseOp();
   BlockUpperTriInverseOp(const BlockUpperTriInverseOp &);
};

inline LinearOp createNewBlockUpperTriInverseOp(BlockedLinearOp & U,const std::vector<LinearOp> & invDiag)
{
   return Teuchos::rcp(new BlockUpperTriInverseOp(U,invDiag));
}

} // end namespace PB

#endif	
