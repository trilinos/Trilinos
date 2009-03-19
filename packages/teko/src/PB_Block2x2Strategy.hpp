#ifndef __PB_Block2x2Strategy_hpp__
#define __PB_Block2x2Strategy_hpp__

#include "Teuchos_RCP.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace PB {

/** @brief Abstract strategy for computing inv(F) and inv(S) in the
 *         Block2x2PreconditionerFactory.
 *
 * This should be paired with a Block2x2PreconditionerFactory,
 * it build the Inv F and InvSchur opreators. Building an approximate
 * inverse of this system
 * 
 * \f$
 * A = \left[ 
 * \begin{array}{cc}
 * F & U \\
 * L & D
 * \end{array}
 * \right]
 * \f$
 *
 * using a 2x2 block LDU decomposition gives
 *
 * \f$
 * A = \left[ 
 * \begin{array}{cc}
 * I & 0  \\
 * L F^{-1} & I
 * \end{array}
 * \right]
 * \left[ 
 * \begin{array}{cc}
 * F & 0  \\
 * 0 & -S
 * \end{array}
 * \right]
 * \left[ 
 * \begin{array}{cc}
 * I &  F^{-1} U \\
 * 0 & I
 * \end{array}
 * \right]
 * \f$
 *
 * where the Shur complement is \f$ S = -D + L F^{-1} U \f$
 * To invert \f$ A \f$, \f$ F^{-1} \f$ and \f$ S^{-1} \f$ are required.  The idea
 * of this strategy is to give those operators.
 */
class Block2x2Strategy {
public:
   /** returns an (approximate) inverse of F */
   virtual const Teuchos::RCP<const Thyra::LinearOpBase<double> >
   getInvF(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;

   /** returns an (approximate) inverse of S = -D + L*inv(F)*U */
   virtual const Teuchos::RCP<const Thyra::LinearOpBase<double> >
   getInvSchur(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;
};

/** @brief A simple strategy for use with Block2x2PreconditionerFactory, that
 *         offers static objects for inv(F) and inv(S)
 *
 * This is a simple startegy for a Block2x2PreconditionerFactory
 * it simply returns statically set RCP pointers to the passed in
 * inv(F) and inv(Schur) operators. Note this will not permit 
 * efficient implementations when the preconditioner has to be rebuilt
 */
class StaticBlock2x2Strategy : public Block2x2Strategy {
public:
   /** @brief Constructor that takes the static inv(F) and inv(S) objects
     *
     * Constructor that takes the static inv(F) and inv(S) objects.
     *
     * @param[in] invF Inverse of the (1,1) block in the source matrix.
     * @param[in] invS Inverse of the Schur complement of the source matrix.
     */
   StaticBlock2x2Strategy(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                          const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invS)
      : invF_(invF), invSchur_(invS)
   {}
 
   /** @name Methods inherited from Block2x2Strategy. */
   //@{

   /** returns a static (approximate) inverse of F */
   virtual const Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvF(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invF_; }

   /** returns a static (approximate) inverse of S = -D + L*inv(F)*U */
   virtual const Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvSchur(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invSchur_; }

   //@}

protected:
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invF_;  /**< Stored value of inv(F) */
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invSchur_; /**< Stored value of inv(S) */

private:
   // hide me!
   StaticBlock2x2Strategy() {}
};

} // end namespace PB

#endif
