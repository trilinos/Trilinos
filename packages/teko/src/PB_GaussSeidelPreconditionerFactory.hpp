#ifndef __PB_GaussSeidelPreconditionerFactory_hpp__
#define __PB_GaussSeidelPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_BlockInvDiagonalStrategy.hpp"
#include "PB_Utilities.hpp"

namespace PB {

typedef enum {GS_UseLowerTriangle,GS_UseUpperTriangle} TriSolveType;

/** \brief A factory that creates a block Gauss Seidel preconditioner.
  *        The user must specify the solvers (or preconditioners) to use
  *        to approximately invert the diagonal operators.
  *
  * A factory that creates a block Gauss Seidel preconditioner.
  * The user must specify the solvers (or preconditioners) to use
  * to approximately invert the diagonal operators.
  * 
  * To invoke this preconditioner using the XML file a diagonal inverse
  * needs to be specified. For example the following XML code creates
  * a Gauss-Seidel preconditioner called "GS-Outer" using Amesos 
  * (a direct solver) to invert the diagonal blocks.
  *
    \verbatim
    <ParameterList name="GS-Outer">
       <Parameter name="Type" type="string" value="Block Gauss-Seidel"/>
       <Parameter name="Inverse Type" type="string" value="Amesos"/>
    </ParameterList>
    \endverbatim
  */
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

      /** Build an empty Gauss-Seidel preconditioner factory
        */
      GaussSeidelPreconditionerFactory();

      //@}

      /** \brief Create the Gauss-Seidel preconditioner operator.
        *
        * This method breaks apart the BlockLinearOp and builds a block
        * diagonal preconditioner. The inverse of the diagonals are specified
        * by the BlockInvDiagonalStrategy object.
        */
      LinearOp buildPreconditionerOperator(BlockedLinearOp & blo,BlockPreconditionerState & state) const;

   protected: 
      //! some members
      Teuchos::RCP<const BlockInvDiagonalStrategy> invOpsStrategy_;
      TriSolveType solveType_;

      //! Initialize from a parameter list
      virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);
};

} // end namespace PB

#endif
