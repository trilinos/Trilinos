// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Util.cpp
  
  \brief Non-templated helper functions.

*/

#include <Zoltan2_Util.hpp>

namespace Zoltan2 {

/*! \brief Search a parameter list for an ostream type of parameter.

    \param pl   the ParameterList to search
    \param key  the name of the parameter to look for
    \param os   a pointer to an ostream for output
    \param defaultValue the value to use if the key is not found in the pl

    \result os is written with the pointer to the ostream in the parameter list,
               or with the defaultValue if the key was not found.

    Note: Parameter lists store pointers to ostreams, not ostreams themselves.
 */

void getOutputStreamFromParameterList( Teuchos::ParameterList &pl, 
  std::string key, std::ostream *&os, std::ostream &defaultValue)
{
  std::ostream **ostreamPtr=NULL;
  std::ofstream **ofstreamPtr=NULL;
  std::ostringstream **ostringstreamPtr=NULL;

  Teuchos::ParameterEntry *entry = pl.getEntryPtr(key);

  if (entry){
    if (entry->isType<std::ostream *>())
      os = entry->getValue(ostreamPtr);
  
    else if (entry->isType<std::ofstream *>())
      os = entry->getValue(ofstreamPtr);
  
    else if (entry->isType<std::ostringstream *>())
      os = entry->getValue(ostringstreamPtr);
  
    else{
      os = &defaultValue;
    }
  }
  else{
    os = &defaultValue;
  }
}

}  //namespace Zoltan2
