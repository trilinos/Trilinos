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

#ifndef __Teko_DiagnosticPreconditionerFactory_hpp__
#define __Teko_DiagnosticPreconditionerFactory_hpp__

#include "Teuchos_Time.hpp"

// Teko includes
#include "Teko_PreconditionerFactory.hpp"

namespace Teko {

/** \brief A class which builds a diagnostic operator
  *        to wrap the application of the inverse operator.
  *        Also times the construction of the inverse operator.
  *
  * For construction purposes this class can be initialized
  * using a parameter list. Most often used in conjuncition with
  * an InverseLibrary object. In particular the relevant parameters are
  *
  \code
       <Parameter name="Type" type="string" value="Diagnostic Inverse"/>
       <Parameter name="Inverse Factory" type="string" value="<Some Inverse Factory>"/>
       <Parameter name="Descriptive Label" type="string" value="<Some Label>"/>
       <Parameter name="Print Residual" type="bool" value="false"/>
  \endcode
  */
class DiagnosticPreconditionerFactory 
   : public virtual Teko::PreconditionerFactory {
public:
   //! Default constructor, for use with the AutoClone class.
   DiagnosticPreconditionerFactory();

   /** Construct a preconditioner factory that prints diagnostics about 
     * a particualar inverse operator.
     *
     * \param[in] invFactory Factory and operator to use diagnostics
     * \param[in] label      Label to give to factory and operator
     */
   DiagnosticPreconditionerFactory(const Teuchos::RCP<Teko::InverseFactory> & invFactory,const std::string & label,
                                   const Teuchos::RCP<std::ostream> & os=Teuchos::null,bool printResidual=false);

   //! default destructor: prints out diagnostic string
   virtual ~DiagnosticPreconditionerFactory();

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
     */
   virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const;
   
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
     */
   virtual bool updateRequestedParameters(const Teuchos::ParameterList & pl);
   
   //@}

   int numInitialBuilds() const { return buildTimer_->numCalls(); }
   double totalInitialBuildTime() const { return buildTimer_->totalElapsedTime(); }

   int numRebuilds() const { return rebuildTimer_->numCalls(); }
   double totalRebuildTime() const { return rebuildTimer_->totalElapsedTime(); }

private:
   void initTimers(const std::string & label);

   Teuchos::RCP<std::ostream> outputStream_;
   Teuchos::RCP<Teko::InverseFactory> invFactory_;
   std::string diagString_;
   bool printResidual_;

   mutable Teuchos::RCP<Teuchos::Time> buildTimer_;   // only first pass construction time (no rebuild)
   mutable Teuchos::RCP<Teuchos::Time> rebuildTimer_; // rebuild-construction timer (no build)
};

} // end namespace Teko

#endif
