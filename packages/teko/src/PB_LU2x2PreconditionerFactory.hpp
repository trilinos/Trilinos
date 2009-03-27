#ifndef __PB_LU2x2PreconditionerFactory_hpp__
#define __PB_LU2x2PreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_LU2x2Strategy.hpp"

namespace PB {

/** \brief Construct a preconditioner using a LDU dcomposition of a block
 *  2x2 matrix.
 *
 * This produces a preconditioner using the block-LDU decomposition of
 * the matrix. The general assumption made is that the matrix is 2x2 
 * and the block factorization can be constructed (i.e. assumptions about
 * the invertability of some blocks). The pattern used, and the one you
 * should follow if you want to use this software is
 *
 * \f$
 * A = \left[ 
 * \begin{array}{cc}
 * A_{00} & A_{01} \\
 * A_{10} & A_{11}
 * \end{array}
 * \right]
 * = \left[ 
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
 * where the Schur complement is \f$ S=-A_{11}+A_{10} A_{00}^{-1} A_{01} \f$ .
 *
 * In order to facilate using this class in a nonlinear solve (or for a 
 * time-dependent problem) the additional abstraction of a ``Strategy''
 * has been added. This strategy, abstractly represented as the LU2x2Strategy,
 * provides the \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ operators. Typical usage for this class
 * is to build a LU2x2Strategy and pass it into the primary constructor. 
 * Additional constructors are provided both for convenience and to ease
 * adoption. Underneath the hood all these constructors do is invoke the
 * corresponding strategy object.
 *
 * For example, assume that you have the particularly nice case that
 * your approximations of \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ are indedent of the source
 * operator. Then, one way to instantiate a LU2x2PreconditionerFactory
 * is

   <code>
      RCP<LinearOpBase<double> > invA00 = buildInvA00(...);\n
      RCP<LinearOpBase<double> > invS   = buildInvS(...);\n
      RCP<LU2x2PreconditionerFactory> precFact = rcp(new LU2x2PreconditionerFactory(invA00,invS));
   </code>

 * Now using the strategy constructor, an entirely equivalent factory
 * object can be constructed by

   <code>
      RCP<LinearOpBase<double> > invA00 = buildInvA00(...);\n
      RCP<LinearOpBase<double> > invS   = buildInvS(...);\n
      RCP<LU2x2Strateghy> precStrat = rcp(new StaticLU2x2Strategy(invA00,invS));\n
      RCP<LU2x2PreconditionerFactory> precFact = rcp(new LU2x2PreconditionerFactory(precStrat));
   </code>
 
 * Notice that the StaticLU2x2Strategy takes the same objects 
 * as the original constructor, it acts as an intermediary to tell the 
 * LU2x2PreconditionerFactory what those operators are.
 **/
class LU2x2PreconditionerFactory : public BlockPreconditionerFactory {
   public:
      //! @name Constructors.
      //@{
 
      /** @brief Build a simple static LU2x2 preconditioner */
      LU2x2PreconditionerFactory(LinearOp & invA00,LinearOp & invS);

      /** @brief Constructor that permits the most generality in building \f$A_{00}^{-1}\f$ and
        *        \f$S^{-1}\f$.
        *
        * Constructor that permits the most generality in building \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$.
        *
        * @param[in] strategy  Strategy object that takes a 2x2 block matrix and
        *                      and constructs the \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ objects.
        */
      LU2x2PreconditionerFactory(const Teuchos::RCP<const LU2x2Strategy> & strategy);

      //@}

      /** \brief Create the LU 2x2 preconditioner operator.
        *
        * This method breaks apart the BlockLinearOp and builds a block
        * LU preconditioner. This will require two applications of the inverse
        * of the (0,0) block and one application of the inverse Schur complement.
        */
      LinearOp buildPreconditionerOperator(BlockedLinearOp & blo) const;
 
   protected: 
      //! some members
      Teuchos::RCP<const LU2x2Strategy> invOpsStrategy_;
};

} // end namespace PB

#endif
