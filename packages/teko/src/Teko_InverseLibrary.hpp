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

#ifndef __Teko_InverseLibrary_hpp__
#define __Teko_InverseLibrary_hpp__

#include "Teko_InverseFactory.hpp"

// Teuchos includes
#include "Teuchos_ParameterList.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Teko includes
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestHandlerContainer.hpp"

namespace Teko {

class InverseLibrary : public RequestHandlerContainer {
public:
   InverseLibrary();

   InverseLibrary(const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> & strat);

   //! add an unspecified inverse to the library
   void addInverse(const std::string & label,const Teuchos::ParameterList & pl);

   //! Add a Stratimikos solver with a label to the library
   void addStratSolver(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl);

   //! Add a Stratimikos preconditioner with a label to the library
   void addStratPrecond(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl);

   //! Add a Teko preconditioner to the library with a label
   void addBlockPrecond(const std::string & label,const std::string & type,const Teuchos::ParameterList & pl);

   /** Get the fully constructed parameter list for a particular label
     *
     * \param[in] label Name used for the desired solver.
     *
     * \returns If the label is found in the library the corresponding parameter list
     *          is returned, otherwise <code>Teuchos::null</code> is returned.
     */
   Teuchos::RCP<const Teuchos::ParameterList> getParameterList(const std::string & label) const;

   //! Get the inverse factory associated with a particular label
   Teuchos::RCP<InverseFactory> getInverseFactory(const std::string & label) const;

   //! Print the inverses and parameter lists available for use
   void PrintAvailableInverses(std::ostream & os) const;

   //! Set the request handler with pointers to the appropriate callbacks
   void setRequestHandler(const Teuchos::RCP<RequestHandler> & rh)
   { callbackHandler_ = rh; }

   //! Get the request handler with pointers to the appropriate callbacks
   Teuchos::RCP<RequestHandler> getRequestHandler() const 
   { return callbackHandler_; }

protected:

   // stratimikos type Inverse objects: mapping the label to a parameter list
   std::map<std::string,Teuchos::RCP<const Teuchos::ParameterList> > stratSolver_;
   std::map<std::string,Teuchos::RCP<const Teuchos::ParameterList> > stratPrecond_;
   std::map<std::string,Teuchos::RCP<const Teuchos::ParameterList> > blockPrecond_;

   // vectors showing which string types are in Stratimikos
   std::vector<std::string> stratValidSolver_;
   std::vector<std::string> stratValidPrecond_;
   std::vector<std::string> blockValidPrecond_;

   //! For handling requests and send requests back to the user
   Teuchos::RCP<RequestHandler> callbackHandler_;

   //! This is the default builder used by stratimikos
   Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> defaultBuilder_;
    
public:

   /** \brief Build an inverse library from a parameter list.
     * 
     * Build an inverse library from a parameter list. This will
     * contain all the labeled inverses specified.
     *
     * \param[in] pl Parameter list to build the library from
     * \param[in] useStratDefaults Also add the default parameters from Stratimikos
     *
     * \returns A pointer to the inverse library created.
     */
   static Teuchos::RCP<InverseLibrary> buildFromParameterList(const Teuchos::ParameterList & pl,bool useStratDefaults=true);

   /** \brief Build an inverse library from a parameter list.
     * 
     * Build an inverse library from a parameter list. This will
     * contain all the labeled inverses specified.
     *
     * \param[in] pl Parameter list to build the library from
     * \param[in] strat Stratimikos object to use
     *
     * \returns A pointer to the inverse library created.
     */
   static Teuchos::RCP<InverseLibrary> buildFromParameterList(const Teuchos::ParameterList & pl,
                                                              const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> & strat);

   /** \brief Build an inverse library from Stratimikos
     * 
     * Build an inverse library from Stratimkos. The labels
     * will just be the names in Stratimikos.
     *
     * \param[in] strat Stratimikos object to use
     *
     * \returns A pointer to the inverse library created.
     */
   static Teuchos::RCP<InverseLibrary> buildFromStratimikos(
         const Stratimikos::DefaultLinearSolverBuilder & strat=Stratimikos::DefaultLinearSolverBuilder());

   /** \brief Build an inverse library from Stratimikos
     * 
     * Build an inverse library from Stratimkos. The labels
     * will just be the names in Stratimikos.
     *
     * \param[in] strat Stratimikos pointer to use
     *
     * \returns A pointer to the inverse library created.
     */
   static Teuchos::RCP<InverseLibrary> buildFromStratimikos(
         const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> & strat);
};

} // end namespace Teko

#endif
