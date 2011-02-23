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

#ifndef __Teko_DiagonallyScaledPreconditionerFactory_hpp__
#define __Teko_DiagonallyScaledPreconditionerFactory_hpp__

#include "Teuchos_Time.hpp"

// Teko includes
#include "Teko_PreconditionerFactory.hpp"

namespace Teko {

/** \brief Using an absolute row sum approximation of the 
  *        matrix this factory creates an inverse using
  *        the explicity scaled matrix. The inverse of the
  *        scaling operation is automatically applied
  *
  * The current behavior for this factory on an operator
  * \f$A\f$ is to compute the absolute row sum
  * of \f$A\f$ and define the scaling matrix \f$D\f$. Then
  * explicitly right muliply by \f$D^{-1}\f$ (\f$A\f$ is not
  * modified, and invert the linear system \f$ A D^{-1}\f$
  * using the user specified inverse factory. Then the scaling
  * is removed by second (implicit) multiplication by \f$D^{-1}\f$.
  * Thus the operator returned by the factory is 
  * \f$D^{-1} \, \mbox{userInv} (A D^{-1})^{-1})\f$.
  *
  * For construction purposes this class can be initialized
  * using a parameter list. Most often used in conjuncition with
  * an InverseLibrary object. In particular the relevant parameters are
  *
  \code
     <Parameter name="Type" type="string" value="Diagonal Scaling"/>
     <Parameter name="Inverse Factory" type="string" value="<Some Inverse Factory>"/>
     <Parameter name="Scaling Type" type="string" value="<Row,Column>"/> <!-- Row is default -->
     <Parameter name="Diagonal Type" type="string" value="<AbsRowSum,Diagonal,Lumped>"/> <!-- AbsRowSum is default -->
  \endcode
  */
class DiagonallyScaledPreconditionerFactory 
   : public virtual Teko::PreconditionerFactory {
public:
   typedef enum {COLUMN_SCALING, ROW_SCALING} ScalingType;

   //! Default constructor, for use with the AutoClone class.
   DiagonallyScaledPreconditionerFactory();

   /** Construct a preconditioner factory that scales the operator.
     *
     * \param[in] invFactory Factory to perform the inverse
     */
   DiagonallyScaledPreconditionerFactory(const Teuchos::RCP<Teko::InverseFactory> & invFactory,ScalingType columnScaling=COLUMN_SCALING,DiagonalType diagonalType=AbsRowSum);

   //! default destructor: prints out diagnostic string
   virtual ~DiagonallyScaledPreconditionerFactory();

   /** \brief Function that is called to build the preconditioner
     *        for the linear operator that is passed in.
     *
     * This function builds a preconditioner based on the passed
     * in LinearOp. 
     *
     * \param[in] lo    Source linear operator that is to be preconditioned.
     * \param[in] state An object associated with this operator to store
     *                  the preconditioner state.
     * 
     * \returns The preconditioner as a linear operator (i.e. to perform
    *           a matrix-vector operation simply call "apply").
     */
   virtual LinearOp buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const;

   //! @name Methods for construction from a parameter list entry
   //@{

   /** \brief This function builds the internals of the preconditioner factory
     *        from a parameter list.
     *        
     * This function builds the internals of the preconditioner factory
     * from a parameter list. Furthermore, it allows a preconditioner factory
     * developer to easily add a factory to the build system. This function
     * is required for building a preconditioner from a parameter list.
     *
     * \param[in] settings Parmaeter list to use as the internal settings
     */
   virtual void initializeFromParameterList(const Teuchos::ParameterList & settings);
   
   //@}

private:
   Teuchos::RCP<Teko::InverseFactory> invFactory_;
   ScalingType scalingType_;
   DiagonalType diagonalType_;
};

} // end namespace Teko

#endif
