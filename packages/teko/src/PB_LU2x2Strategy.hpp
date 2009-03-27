#ifndef __PB_LU2x2Strategy_hpp__
#define __PB_LU2x2Strategy_hpp__

#include "Teuchos_RCP.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace PB {

/** @brief Abstract strategy for computing inv(F) and inv(S) in the
 *         LU2x2PreconditionerFactory.
 *
 * This should be paired with a LU2x2PreconditionerFactory,
 * it build the \f$A_{00}^{-1}\f$ and \f$S\f$ opreators. Building an approximate
 * inverse of this system
 * 
 * \f$
 * A = \left[ 
 * \begin{array}{cc}
 * A_{00} & A_{01} \\
 * A_{10} & A_{11}
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
 * A_{10} A_{00}^{-1} & I
 * \end{array}
 * \right]
 * \left[ 
 * \begin{array}{cc}
 * A_{00} & 0  \\
 * 0 & -S
 * \end{array}
 * \right]
 * \left[ 
 * \begin{array}{cc}
 * I &  A_{00}^{-1} A_{01} \\
 * 0 & I
 * \end{array}
 * \right]
 * \f$
 *
 * where the Shur complement is \f$ S = -A_{11} + A_{10} A_{00}^{-1} A_{01} \f$
 * To invert \f$ A \f$, \f$ A_{00}^{-1} \f$ and \f$ S^{-1} \f$ are required.  The idea
 * of this strategy is to give those operators.
 */
class LU2x2Strategy {
public:
   /** returns an (approximate) inverse of \f$A_{00}\f$ */
   virtual const Teuchos::RCP<const Thyra::LinearOpBase<double> >
   getInvA00(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;

   /** returns an (approximate) inverse of \f$S = -A_{11} + A_{10} A_{00}^{-1} A_{01} */
   virtual const Teuchos::RCP<const Thyra::LinearOpBase<double> >
   getInvS(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;
};

/** @brief A simple strategy for use with LU2x2PreconditionerFactory, that
 *         offers static objects for inv(F) and inv(S)
 *
 * This is a simple startegy for a LU2x2PreconditionerFactory
 * it simply returns statically set RCP pointers to the passed in
 * inv(F) and inv(Schur) operators. Note this will not permit 
 * efficient implementations when the preconditioner has to be rebuilt
 */
class StaticLU2x2Strategy : public LU2x2Strategy {
public:
   /** @brief Constructor that takes the static \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ objects
     *
     * Constructor that takes the static \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ objects.
     *
     * @param[in] invA00 Inverse of \f$A_{00}\f$ in the source matrix.
     * @param[in] invS Inverse of the Schur complement of the source matrix.
     */
   StaticLU2x2Strategy(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invA00,
                          const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invS)
      : invA00_(invA00), invS_(invS)
   {}
 
   /** @name Methods inherited from LU2x2Strategy. */
   //@{

   /** returns a static (approximate) inverse of F */
   virtual const Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvA00(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invA00_; }

   /** returns a static (approximate) inverse of S = -D + L*inv(F)*U */
   virtual const Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvS(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invS_; }

   //@}

protected:
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invA00_;  /**< Stored value of \f$A_{00}^{-1}\f$ */
   const Teuchos::RCP<const Thyra::LinearOpBase<double> > invS_; /**< Stored value of \f$S^{-1}\f$ */

private:
   // hide me!
   StaticLU2x2Strategy() {}
};

} // end namespace PB

#endif
