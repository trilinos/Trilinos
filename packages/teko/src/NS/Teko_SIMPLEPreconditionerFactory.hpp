/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software. 
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

#ifndef __Teko_SIMPLEPreconditionerFactory_hpp__
#define __Teko_SIMPLEPreconditionerFactory_hpp__

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_InverseFactory.hpp"

namespace Teko {
namespace NS {

// Declaration of the preconditioner factory
/** The basic XML parameter list for SIMPLE looks like.
   
   \code
    <ParameterList name="SIMPLEC">
      <Parameter name="Type" type="string" value="NS SIMPLE"/>

      <!-- Inverse operations to use -->
      <Parameter name="Inverse Velocity Type" type="string" value="ML_GS-Velocity"/>
      <Parameter name="Inverse Pressure Type" type="string" value="ML_GS-Pressure"/>

      <Parameter name="Explicit Velocity Inverse Type" type="string" value="AbsRowSum"/>

      <Parameter name="Alpha" type="double" value="0.9"/>
    </ParameterList>
  \endcode
  */
class SIMPLEPreconditionerFactory : public BlockPreconditionerFactory {
public:
   // Constructor
   SIMPLEPreconditionerFactory(const Teuchos::RCP<InverseFactory> & inverse,
                               double alpha);

   // Constructor
   SIMPLEPreconditionerFactory(const Teuchos::RCP<InverseFactory> & invVelFactory,
                               const Teuchos::RCP<InverseFactory> & invPrsFactory,
                               double alpha);

   //! Default constructor
   SIMPLEPreconditionerFactory();

   // Function inherited from BlockPreconditionerFactory
   LinearOp buildPreconditionerOperator(BlockedLinearOp & blo,
                                        BlockPreconditionerState & state) const;

   //! Set the mass matrix for this factory
   virtual void setMassMatrix(Teko::LinearOp & mass)
   { massMatrix_ = mass; }

   //! For assisting in construction of the preconditioner
   virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const;

   //! For assisting in construction of the preconditioner
   virtual bool updateRequestedParameters(const Teuchos::ParameterList & pl);

   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);
    
protected:
   // class members
   Teuchos::RCP<InverseFactory> customHFactory_;
   Teuchos::RCP<InverseFactory> invVelFactory_;
   Teuchos::RCP<InverseFactory> invPrsFactory_;
   double alpha_;
   DiagonalType fInverseType_;
   // enum FInverseType {Diagonal,Lumped,AbsRowSum,Custom} fInverseType_;

   bool useMass_;
   Teko::LinearOp massMatrix_;
   
   // Info for the block-diagonal approximation to H if used.
   mutable Teuchos::ParameterList BlkDiagList_;
};
 
} // end namespace NS
} // end namespace Teko

#endif
