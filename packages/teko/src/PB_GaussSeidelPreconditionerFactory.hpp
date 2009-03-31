#ifndef __PB_GaussSeidelPreconditionerFactory_hpp__
#define __PB_GaussSeidelPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_BlockInvDiagonalStrategy.hpp"
#include "PB_Utilities.hpp"

namespace PB {

typedef enum {GS_UseLowerTriangle,GS_UseUpperTriangle} TriSolveType;

class GaussSeidelPreconditionerFactory : public BlockPreconditionerFactory {
   public:
  
      //! @name Constructors.
      //@{

      /*! Construct a PreconditionerFactory assuming a specific block
          \f$2\times2\f$ matrix. This case is a simple one.
      */ 
      GaussSeidelPreconditionerFactory(TriSolveType solveType,const LinearOp & invD0,const LinearOp & invD1);

      /*! The most flexible JacobiPreconditionerFactory constructor.
          Pass in a generally defined BlockInvDiagonalStrategy to use the
          full generality of this class.
      */
      GaussSeidelPreconditionerFactory(TriSolveType solveType,const RCP<const BlockInvDiagonalStrategy> & strategy);

      //@}

      /** \brief Create the Gauss-Seidel preconditioner operator.
        *
        * This method breaks apart the BlockLinearOp and builds a block
        * diagonal preconditioner. The inverse of the diagonals are specified
        * by the BlockInvDiagonalStrategy object.
        */
      LinearOp buildPreconditionerOperator(BlockedLinearOp & blo) const;
 
   protected: 
      //! some members
      Teuchos::RCP<const BlockInvDiagonalStrategy> invOpsStrategy_;
      TriSolveType solveType_;
};

} // end namespace PB

#endif
