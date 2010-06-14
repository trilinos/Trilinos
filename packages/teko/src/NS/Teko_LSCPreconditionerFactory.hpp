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

#ifndef __Teko_LSCPreconditionerFactory_hpp__
#define __Teko_LSCPreconditionerFactory_hpp__

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_LSCStrategy.hpp"

namespace Teko {
namespace NS { // Navier-Stokes specialization

/** \brief Preconditioner state for the LSC factory. 
  *
  * Preconditioner state for the LSC factory. This is based
  * on the notation and concepts found in
  *
  * Elman, Howle, Shadid, Silvester, and Tuminaro, "Least Squares Preconditioners
  * for Stabilized Discretizations of the Navier-Stokes Euqations," SISC-2007.
  */
class LSCPrecondState : public BlockPreconditionerState {
public:
   LSCPrecondState() {}

   //! Inverse mass operator (\f$Q_u^{-1}\f$)
   LinearOp invMass_;

   /** \f$B Q_u^{-1} B^T\f$
     */
   ModifiableLinearOp BQBt_;

   /** \f$B H B^T\f$
     */
   ModifiableLinearOp BHBt_;

   /** \f$B Q_u^{-1} B^T-\gamma C\f$
     */
   LinearOp BQBtmC_;
   InverseLinearOp invBQBtmC_;

   /** \f$B H B^T-\gamma C\f$
     */
   LinearOp BHBtmC_;
   InverseLinearOp invBHBtmC_;

   //! \f$\alpha D^{-1}\f$ where
   LinearOp aiD_;

   //! \f$\gamma = \rho(Q_u^{-1} F / 3)\f$
   double gamma_;

   /** \f$\alpha = 1/\rho(B \; diag(F)^{-1} B^T D^{-1})\f$ where
     * \f$D = diag(B \; diag(F)^{-1} B^T + C)\f$.
     */
   double alpha_;
};

class LSCPreconditionerFactory 
   : public BlockPreconditionerFactory {
public:
   //! \name Constructors
   //@{

   //! Staiblized constructor
   LSCPreconditionerFactory(const LinearOp & invF,const LinearOp & invBQBtmC,
                            const LinearOp & invD,const LinearOp & invMass);
 
   //! Stable constructor
   LSCPreconditionerFactory(const LinearOp & invF,
                            const LinearOp & invBQBtmC,
                            const LinearOp & invMass);

   //! fully generic constructor
   LSCPreconditionerFactory(const Teuchos::RCP<LSCStrategy> & strategy);

   //! Default constructor
   LSCPreconditionerFactory();
   //@}

   //! for PreconditionerFactoryBase
   virtual LinearOp buildPreconditionerOperator(BlockedLinearOp & blo,BlockPreconditionerState & state) const;

   //! Build the LSCPrecondState object
   virtual RCP<PreconditionerState> buildPreconditionerState() const
   { return rcp(new LSCPrecondState()); }

   //! For assiting in construction of the preconditioner
   virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const;

   //! For assiting in construction of the preconditioner
   virtual bool updateRequestedParameters(const Teuchos::ParameterList & pl);

protected:
   // Gimmie object
   Teuchos::RCP<LSCStrategy> invOpsStrategy_;
   bool isSymmetric_;

   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);

public:
   /** \brief Builder function for creating strategies.
     *
     * Builder function for creating strategies.
     * 
     * \param[in] name     String name of strategy to build
     * \param[in] settings Parameter list describing the parameters for the
     *                     strategy to build
     * \param[in] invLib   Inverse library for the strategy to use.
     *
     * \returns If the name is associated with a strategy
     *          a pointer is returned, otherwise Teuchos::null is returned.
     */
   static RCP<LSCStrategy> 
   buildStrategy(const std::string & name, 
                 const Teuchos::ParameterList & settings,
                 const RCP<const InverseLibrary> & invLib,
                 const RCP<RequestHandler> & rh);

   /** \brief Add a strategy to the builder. This is done using the
     *        clone pattern. 
     *
     * Add a strategy to the builder. This is done using the
     * clone pattern. If your class does not support the Cloneable interface then
     * you can use the AutoClone class to construct your object.
     *
     * \note If this method is called twice with the same string, the latter clone pointer
     *       will be used.
     *
     * \param[in] name String to associate with this object
     * \param[in] clone Pointer to Cloneable object
     */
   static void addStrategy(const std::string & name,const RCP<Cloneable> & clone);

private:
   //! for creating the strategy objects
   static CloneFactory<LSCStrategy> strategyBuilder_;

   //! This is where the default objects are put into the strategyBuilder_
   static void initializeStrategyBuilder();
};

} // end namespace NS
} // end namespace Teko

#endif
