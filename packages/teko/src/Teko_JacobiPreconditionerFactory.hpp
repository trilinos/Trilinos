#ifndef __Teko_JacobiPreconditionerFactory_hpp__
#define __Teko_JacobiPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_BlockInvDiagonalStrategy.hpp"

namespace Teko {

class JacobiPreconditionerFactory : public BlockPreconditionerFactory {
   public:
      //! @name Constructors.
      //@{

      /** Construct a PreconditionerFactory assuming a specific block
        * \f$2\times2\f$ matrix. This case is a simple one.
        */ 
      JacobiPreconditionerFactory(const LinearOp & invD0,const LinearOp & invD1);

      /** The most flexible JacobiPreconditionerFactory constructor.
        * Pass in a generally defined BlockInvDiagonalStrategy to use the
        * full generality of this class.
        */
      JacobiPreconditionerFactory(const RCP<const BlockInvDiagonalStrategy> & strategy);

      /** Build an empty Jacobi preconditioner factory.
        */
      JacobiPreconditionerFactory();

      //@}

      /** \brief Create the Jacobi preconditioner operator.
        *
        * This method breaks apart the BlockLinearOp and builds a block
        * diagonal preconditioner. The inverse of the diagonals are specified
        * by the BlockInvDiagonalStrategy object.
        */
      LinearOp buildPreconditionerOperator(BlockedLinearOp & blo,BlockPreconditionerState & state) const;
 
   protected: 
      //! some members
      Teuchos::RCP<const BlockInvDiagonalStrategy> invOpsStrategy_;

      //! Initialize from a parameter list
      virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);
};

} // end namespace Teko

#endif
