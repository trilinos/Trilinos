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

#ifndef __Teko_PCDStrategy_hpp__
#define __Teko_PCDStrategy_hpp__

#include "Teko_LU2x2Strategy.hpp"

// Teuchos includes
#include "Teuchos_Time.hpp"

namespace Teko {
namespace NS {

/** @brief Strategy for computing implementation of the 
 *         Pressure convection diffusion preconditioner.
 *
 * This requires the user to add a "PCD Operator", 
 * "Pressure Laplace Operator" and "Press Mass Operators"
 * into the preconditioner state. The user is notified by
 * boolean fields returned from the requested paramters list.
 */
class PCDStrategy : public LU2x2Strategy {
public:
   //! default Constructor 
   PCDStrategy();

   //! Constructor to set the inverse factories.
   PCDStrategy(const Teuchos::RCP<InverseFactory> & invFA,
                         const Teuchos::RCP<InverseFactory> & invS);

   //! Destructor (does nothing)
   virtual ~PCDStrategy() {}

   /** returns the first (approximate) inverse of \f$A_{00}\f$ */
   virtual const Teko::LinearOp
   getHatInvA00(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const;

   /** returns the second (approximate) inverse of \f$A_{00}\f$ */
   virtual const Teko::LinearOp
   getTildeInvA00(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const;

   /** returns an (approximate) inverse of \f$S = -A_{11} + A_{10} \mbox{diag}(A_{00})^{-1} A_{01}\f$ */
   virtual const Teko::LinearOp
   getInvS(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const;

   /** \brief This function builds the internals of the state from a parameter list.
     *        
     * \param[in] settings Parameter list to use as the internal settings
     * \param[in] invLib Inverse library to use for building inverse factory objects
     *
     * \note The default implementation does nothing.
     */
   virtual void initializeFromParameterList(const Teuchos::ParameterList & settings,
                                            const InverseLibrary & invLib);

   /** \brief Request the additional parameters this preconditioner factory
     *        needs.
     *
     * Request the additonal parameters needed by this preconditioner factory.
     * The parameter list will have a set of fields that can be filled with
     * the requested values. These fields include all requirements, even those
     * of the sub-solvers if there are any.  Once correctly filled the object
     * can be updated by calling the updateRequestedParameters with the filled
     * parameter list.
     *
     * For the PCD strategy the following fields are required to be set to
     * true, they are passed to the user as false. The user acknowledges that
     * the operators are required by updating the parameters to true.
     * <ul>
     * <li><ParameterList name="PCD Operator" type="bool" value="false"/></li>
     * <li><ParameterList name="Pressure Laplace Operator" type="bool" value="false"/></li>
     * <li><ParameterList name="Pressure Mass Operator" type="bool" value="false"/></li>
     * </ul>
     *
     * \returns A parameter list with the requested parameters.
     *
     * \note The default implementation returns Teuchos::null.
     */
   virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const;

   /** \brief Update this object with the fields from a parameter list.
     *
     * Update the requested fields using a parameter list. This method is
     * expected to pair with the getRequestedParameters method (i.e. the fields
     * requested are going to be update using this method).
     *
     * For the PCD strategy the following fields are required to be set to
     * true. Essentially, the user is acknowledging that the operators are
     * required.
     * <ul>
     * <li><ParameterList name="PCD Operator" type="bool" value="true"/></li>
     * <li><ParameterList name="Pressure Laplace Operator" type="bool" value="true"/></li>
     * <li><ParameterList name="Pressure Mass Operator" type="bool" value="true"/></li>
     * </ul>
     *
     * \param[in] pl Parameter list containing the requested parameters.
     *
     * \returns If the method succeeded (found all its required parameters) this
     *          method returns true, otherwise it returns false.
     *
     * \note The default implementation returns true (it does nothing!).
     */
   virtual bool updateRequestedParameters(const Teuchos::ParameterList & pl);

protected:
   /** Build timers for this type of object.
     */
   static void buildTimers();

   /** Initialize the operator's state. This builds the Schur complement and the inverse
     * operators. If the state has already been initialized this method does nothing.
     *
     * \param[in] A Operator to intialize with.
     * \param[in] state Storage object for this operator.
     */
   void initializeState(const Teko::BlockedLinearOp & A,BlockPreconditionerState & state) const;

   // how to invert the matrices
   Teuchos::RCP<InverseFactory> invFactoryF_; // for \tilde{A_00}\f$
   Teuchos::RCP<InverseFactory> invFactoryS_;

   DiagonalType massInverseType_;

   //! Passed to application for construction of laplace operator
   Teuchos::RCP<Teuchos::ParameterList> lapParams_; 

   //! Passed to application for construction of PCD operator
   Teuchos::RCP<Teuchos::ParameterList> pcdParams_; 

   bool schurCompOrdering_;

   static Teuchos::RCP<Teuchos::Time> initTimer_;
   static Teuchos::RCP<Teuchos::Time> invSTimer_;
   static Teuchos::RCP<Teuchos::Time> invFTimer_;
   static Teuchos::RCP<Teuchos::Time> opsTimer_;

public:
   // some static functions for determining strings

   static std::string getPCDString() { return "PCD Operator"; }
   static std::string getPressureLaplaceString() { return "Pressure Laplace Operator"; }
   static std::string getPressureMassString() { return "Pressure Mass Matrix"; }
};

} // end namespace NS
} // end namespace Teko

#endif
