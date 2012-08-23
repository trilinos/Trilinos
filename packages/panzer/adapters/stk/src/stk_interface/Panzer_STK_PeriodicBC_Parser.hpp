// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_STK_PeriodicBC_Parser_hpp__
#define __Panzer_STK_PeriodicBC_Parser_hpp__

#include "Panzer_STK_PeriodicBC_Matcher.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_RCP.hpp"

#include <vector>
#include <string>

namespace panzer_stk {

/** Read a parameter list to describe the periodic
  * boundary conditions. This object then provides
  * a vector of the PeriodicBC_Matcher objects.
  */
class PeriodicBC_Parser : Teuchos::ParameterListAcceptor {
public:
   PeriodicBC_Parser();

   /** Return a vector containing all the periodic boundary conditions.
     */
   const std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > &
   getMatchers() const;

   // parameterlistacceptor required functions
   /////////////////////////////////////
   virtual void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);

   Teuchos::RCP<Teuchos::ParameterList>  getNonconstParameterList();

   Teuchos::RCP<Teuchos::ParameterList>  unsetParameterList();

   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   // specific to this class
   ///////////////////////////

   /** Get valid parameters given a count parameter for total number
     * of boundary conditions.
     *
     * \param[in] count Number of periodic boundary conditions
     */
   Teuchos::RCP<Teuchos::ParameterList> getValidParameters(int count) const;

   /** Build a periodic matcher object given a string.
     *
     * \param[in] buildStr String specifying the matcher to build. 
     *                     Format: "MatchCondition bndry1;bndry2"
     */
   Teuchos::RCP<const PeriodicBC_MatcherBase> buildMatcher(const std::string & buildStr) const;

   /** Parse a string describing the periodic boundary condition
     * Format: "MatchCondition bndry1;bndry2"
     */
   void buildMatcher_Tokenize(const std::string & buildStr,
                             std::string & matcher,
                             std::string & bndry1,
                             std::string & bndry2) const;

   /** Parse a string describing the periodic boundary condition
     * Format: "MatchCondition paramA, paramB, paraC, ... : bndry1;bndry2"
     *
     * \returns False if no parameters are found (defaults to old style of
     *          periodic BC input)
     */
   bool buildMatcher_Tokenize_withParams(const std::string & buildStr,
                                         std::string & matcher,
                                         std::vector<std::string> & params,
                                         std::string & bndry1,
                                         std::string & bndry2) const;
private:
   //! stored parameter list
   Teuchos::RCP<Teuchos::ParameterList> storedPL_;

   //! matchers constructed by "setParameterList"
   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > matchers_;

   // stored string values
   const std::string countStr_;
   const std::string condPrefix_;
};

}

#endif
