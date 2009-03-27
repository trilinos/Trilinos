#ifndef __PB_JacobiPreconditionerFactory_hpp__
#define __PB_JacobiPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_BlockInvDiagonalStrategy.hpp"

namespace PB {

class JacobiPreconditionerFactory : public BlockPreconditionerFactory {
   public:
      //! @name Constructors.
      //@{

      /*! Construct a PreconditionerFactory assuming a specific block
          \f$2\times2\f$ matrix. This case is a common one.
      */ 
      JacobiPreconditionerFactory(const LinearOp & invD0,const LinearOp & invD1);

      /*! The most flexible JacobiPreconditionerFactory constructor.
          Pass in a generally defined BlockInvDiagonalStrategy to use the
          full generality of this class.
      */
      JacobiPreconditionerFactory(const RCP<const BlockInvDiagonalStrategy> & strategy);

      //@}

      /** \brief Create the Jacobi preconditioner operator.
        *
        * This method breaks apart the BlockLinearOp and builds a block
        * diagonal preconditioner. The inverse of the diagonals are specified
        * by the BlockInvDiagonalStrategy object.
        */
      LinearOp buildPreconditionerOperator(BlockedLinearOp & blo) const;
 
   protected: 
      //! some members
      Teuchos::RCP<const BlockInvDiagonalStrategy> invOpsStrategy_;
};

} // end namespace PB

#endif
