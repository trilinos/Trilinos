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

#ifndef __Teko_InverseFactory_hpp__
#define __Teko_InverseFactory_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Teko_Config.h"
#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerState.hpp"
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestHandlerContainer.hpp"

namespace Teko {

/** \brief Abstract class for building an inverse operator
  *
  * Abstract class for building an inverse operator. It pairs
  * with a linear operator and gives you a new operator that
  * behaves like its inverse.
  */
class InverseFactory : public RequestHandlerContainer {
public:
   virtual ~InverseFactory() {}

   /** \brief Build an inverse operator
     *
     * Build the inverse operator using this factory.
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp) const = 0;

   /** \brief Build a preconditioned inverse operator
     * 
     * Build the inverse operator using this factory and a user specified
     * preconditioning operator. The default behavior is to call buildInverse
     * ignoring the preconditioner. 
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     * \param[in] precOp Preconditioning operator
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp,const LinearOp & precOp) const
   { return buildInverse(linearOp); }

   #if 0
   /** \brief Build an inverse operator and make sure it aware of some parents state
     *        This functionality is only useful for Teko::PreconditionerFactory inverses.
     *
     * Build an inverse operator and make sure it aware of some parents state
     * This functionality is only useful for Teko::PreconditionerFactory inverses.
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     * \param[in] parentState Current state object to be used. Only useful for preconditioners.
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp, const PreconditionerState & parentState) const
   { return buildInverse(linearOp); }

   /** \brief Build an inverse operator and make sure it aware of some parents state
     *        This functionality is only useful for Teko::PreconditionerFactory inverses.
     *
     * Build an inverse operator using a preconditioning operator and make sure it
     * is aware of some parents state. The default behavior is to call buildInverse
     * ignoring the preconditioner. This functionality is only useful for
     * Teko::PreconditionerFactory inverses.
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     * \param[in] precOp Preconditioning operator
     * \param[in] parentState Current state object to be used. Only useful for preconditioners.
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp, const LinearOp & precOp, const PreconditionerState & parentState) const
   { return buildInverse(linearOp,precOp); }
   #endif

   /** \brief Pass in an already constructed inverse operator. Update
     *        the inverse operator based on the new source operator.
     *
     * Pass in an already constructed inverse operator. Update
     * the inverse operator based on the new source operator.
     *
     * \param[in]     source Source operator to be inverted.
     * \param[in,out] dest   Pre constructed inverse operator to be
     *                        rebuilt using the <code>source</code>
     *                        object.
     */
   virtual void rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const = 0;

   /** \brief Pass in an already constructed inverse operator. Update
     *        the inverse operator based on the new source operator.
     *
     * Pass in an already constructed inverse operator. Update
     * the inverse operator based on the new source operator.
     *
     * \param[in]     source Source operator to be inverted.
     * \param[in]     precOp Preconditioning operator
     * \param[in,out] dest   Pre constructed inverse operator to be
     *                        rebuilt using the <code>source</code>
     *                        object.
     */
   virtual void rebuildInverse(const LinearOp & source,const LinearOp & precOp,InverseLinearOp & dest) const
   { rebuildInverse(source,dest); }

   /** \brief A function that permits inspection of the parameters used to create
     *        this object.
     *
     * A function that permits inspection of the parameters used to create this
     * object. Useful for determining defaults and settings used.
     *
     * \returns A list used to parameterize this object.
     */
   virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const = 0;

   /** Return a string that describes this factory */
   virtual std::string toString() const = 0;

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
     * \returns A parameter list with the requested parameters.
     *
     * \note The default implementation returns Teuchos::null.
     */
   virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const
   { return Teuchos::null; }
   
   /** \brief Update this object with the fields from a parameter list.
     *
     * Update the requested fields using a parameter list. This method is
     * expected to pair with the getRequestedParameters method (i.e. the fields
     * requested are going to be update using this method).
     *
     * \param[in] pl Parameter list containing the requested parameters.
     *
     * \returns If the method succeeded (found all its required parameters) this
     *          method returns true, otherwise it returns false.
     *
     * \note The default implementation returns true (it does nothing!).
     */
   virtual bool updateRequestedParameters(const Teuchos::ParameterList & pl)
   { return true; }

   //! Set the request handler with pointers to the appropriate callbacks
   void setRequestHandler(const Teuchos::RCP<RequestHandler> & rh)
   { callbackHandler_ = rh; }

   //! Get the request handler with pointers to the appropriate callbacks
   Teuchos::RCP<RequestHandler> getRequestHandler() const 
   { return callbackHandler_; }

protected:
   //! For handling requests and send requests back to the user
   Teuchos::RCP<RequestHandler> callbackHandler_;
};

class StaticOpInverseFactory : public InverseFactory {
public:
   //! \name Constructors
   //@{ 
   
   /** \brief Constructor that takes a linear operator and
     *        uses it as a static inverse
     *
     * Constructor that takes a linear operator and
     * uses it as a static inverse
     * 
     * \param[in] inv Linear operator to use as the inverse.
     */
   StaticOpInverseFactory(const LinearOp inv) 
      : inverse_(inv) {}

   //! Copy constructor
   StaticOpInverseFactory(const StaticOpInverseFactory & saFactory)
      : inverse_(saFactory.inverse_) {}
   //@}

   virtual ~StaticOpInverseFactory() {}

   /** \brief Build an inverse operator
     *
     * Build the inverse operator using this factory. This also tacks
     * on extra data to the RCP called "prec". This is the
     * PreconditionerBase object, and it is need for <code>rebuildInverse</code>
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp) const
   { return Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(inverse_); }

   /** \brief Pass in an already constructed inverse operator. Update
     *        the inverse operator based on the new source operator.
     *
     * Pass in an already constructed inverse operator. Update
     * the inverse operator based on the new source operator. This
     * method assumes the <code>dest</code> object also contains
     * the associated PreconditionerBase object as "prec" as extra
     * data in the RCP.
     *
     * \param[in]     source Source operator to be inverted.
     * \param[in,out] dest   Pre constructed inverse operator to be
     *                        rebuilt using the <code>source</code>
     *                        object.
     */
   virtual void rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const
   { }

   /** \brief A function that permits inspection of the parameters used to create
     *        this object.
     *
     * A function that permits inspection of the parameters used to create this
     * object. Useful for determining defaults and settings used.
     *
     * \returns A list used to parameterize this object.
     */
   virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const
   { return Teuchos::null; }

   /** Return a string that describes this factory */
   virtual std::string toString() const { return inverse_->description(); }

protected:
   Teko::LinearOp inverse_;

private:
   // hide me!
   StaticOpInverseFactory();
};

//! @name Functions for constructing and initializing solvers
//@{

/** Build an inverse operator using a factory and a linear operator
  *
  * \param[in] factory The inverse factory used to construct the inverse
  *                    operator
  * \param[in] A       Linear operator whose inverse is required
  *
  * \returns An (approximate) inverse operator is returned for the operator <code>A</code>.
  *
  * \relates InverseFactory
  */
InverseLinearOp buildInverse(const InverseFactory & factory,const LinearOp & A);

/** Build an inverse operator using a factory and a linear operator
  *
  * \param[in] factory The inverse factory used to construct the inverse
  *                    operator
  * \param[in] precOp  Preconditioning operator
  * \param[in] A       Linear operator whose inverse is required
  *
  * \returns An (approximate) inverse operator is returned for the operator <code>A</code>.
  *
  * \relates InverseFactory
  */
InverseLinearOp buildInverse(const InverseFactory & factory,const LinearOp & A,const LinearOp & precOp);

#if 0
/** Build an inverse operator using a factory and a linear operator
  * This functionality is only useful for Teko::PreconditionerFactory inverses.
  *
  * \param[in] factory The inverse factory used to construct the inverse
  *                    operator
  * \param[in] A       Linear operator whose inverse is required
  * \param[in] parentState Current state object to be used. Only useful for preconditioners.
  *
  * \returns An (approximate) inverse operator is returned for the operator <code>A</code>.
  *
  * \relates PreconditionerInverseFactory InverseFactory
  */
InverseLinearOp buildInverse(const InverseFactory & factory,const LinearOp & linearOp,
                             const PreconditionerState & parentState);

/** Build an inverse operator using a factory and a linear operator
  * This functionality is only useful for Teko::PreconditionerFactory inverses.
  *
  * \param[in] factory The inverse factory used to construct the inverse
  *                    operator
  * \param[in] A       Linear operator whose inverse is required
  * \param[in] precOp  Preconditioning operator
  * \param[in] parentState Current state object to be used. Only useful for preconditioners.
  *
  * \returns An (approximate) inverse operator is returned for the operator <code>A</code>.
  *
  * \relates PreconditionerInverseFactory InverseFactory
  */
InverseLinearOp buildInverse(const InverseFactory & factory,const LinearOp & linearOp,
                             const LinearOp & precOp,
                             const PreconditionerState & parentState);
#endif

/** Using a prebuilt linear operator, use factory to build an inverse operator
  * given a new forward operator.
  *
  * \note This function sometimes fails depending on the underlying type
  *       of the inverse factory.  Use with caution.
  *
  * \param[in] factory The inverse factory used to construct the inverse
  *                    operator
  * \param[in] A       Linear operator whose inverse is required
  * \param[in] invA    The inverse operator that is to be rebuilt using
  *                    the <code>A</code> operator.
  *
  * \relates InverseFactory
  */
void rebuildInverse(const InverseFactory & factory, const LinearOp & A, InverseLinearOp & invA);

/** Using a prebuilt linear operator, use factory to build an inverse operator
  * given a new forward operator.
  *
  * \note This function sometimes fails depending on the underlying type
  *       of the inverse factory.  Use with caution.
  *
  * \param[in] factory The inverse factory used to construct the inverse
  *                    operator
  * \param[in] A       Linear operator whose inverse is required
  * \param[in] precOp  Preconditioning operator
  * \param[in] invA    The inverse operator that is to be rebuilt using
  *                    the <code>A</code> operator.
  *
  * \relates InverseFactory
  */
void rebuildInverse(const InverseFactory & factory, const LinearOp & A,const LinearOp & precOp, InverseLinearOp & invA);

/** \brief Build an InverseFactory object from a ParameterList, as specified in Stratimikos.
  *
  * Build an InverseFactory object from a ParameterList, as specified in Stratimikos.
  * The specific inverse routine (either solver or preconditioner) to be chosen is specified
  * by a string.
  *
  * \note It is preferred that the <code>InverseLibrary</code> is used to construct an
  *       <code>InverseFactory</code> instead.
  *
  * \param[in] list ParameterList that describes the available solvers/preconditioners.
  * \param[in] type String saying which solver/preconditioner to use.
  *
  * \returns An inverse factory using the specified inverse operation.
  *
  * \relates InverseFactory
  */
Teuchos::RCP<InverseFactory> invFactoryFromParamList(const Teuchos::ParameterList & list,const std::string & type);

/** \brief Get a valid parameter list for the inverse factory class.
  *
  * Get a valid parameter list for the inverse factory class. This will
  * specify the set of parameters for each possible "inverse".
  *
  * \note It is preferred that the <code>InverseLibrary</code> is used 
  *       to get paramter lists for <code>InverseFactory</code> construction.
  *
  * \returns A parameter list is returned that is suitable to be passed
  *          to <code>invFactoryFromParamList</code>.
  *
  * \relates InverseFactory
  */
Teuchos::RCP<const Teuchos::ParameterList> invFactoryValidParameters();

//@}

} // end namespace Teko

#endif
