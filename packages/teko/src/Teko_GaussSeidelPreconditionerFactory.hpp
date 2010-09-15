/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// Export of this program may require a license from the United States
// Government.
//  
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//  
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//  
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//  
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission. 
//  
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#ifndef __Teko_GaussSeidelPreconditionerFactory_hpp__
#define __Teko_GaussSeidelPreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_BlockInvDiagonalStrategy.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {

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

} // end namespace Teko

#endif
