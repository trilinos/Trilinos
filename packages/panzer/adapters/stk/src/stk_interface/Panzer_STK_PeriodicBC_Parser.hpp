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
