// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

   /** Get the flag indicating if the bounding box search is used
   * when matching periodic Ids. */
   const bool & useBoundingBoxSearch() const;

   // parameterlistacceptor required functions
   /////////////////////////////////////
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & pl);

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

   /** Build a periodic matcher object given a string
     *
     * \param[in] buildStr String specifying the matcher to build. 
     *                     Format: "MatchCondition bndry1;bndry2"
     */
   Teuchos::RCP<const PeriodicBC_MatcherBase>
   buildMatcher(const std::string & buildStr) const;

   /** Return the string for the type of matcher: coord, edge, face, or all
     *    and the dimension: 2 or 3
     *
     * */
   std::pair<std::string, unsigned int> getMatcherTypeAndDim(const std::string & buildStr) const; 

   /** Replace "all" with specific type in matcher string:
     *    coord, edge, or face
     * */
   std::string replaceMatcherType(const std::string & buildStr, const std::string & matcherType) const; 

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
   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > edgeMatchers_;
   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > faceMatchers_;

   // stored string values
   const std::string countStr_;
   const std::string condPrefix_;
   const std::string searchStr_;
   
   // stored flag indicating if bounding box search is used for periodic DOF matching
   bool useBBoxSearch_;
};

}

#endif
