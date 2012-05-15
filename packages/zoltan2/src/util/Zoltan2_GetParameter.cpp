// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GetParameter.cpp
 *  \brief Convenience methods for working with the parameter list.
 */

#include <Zoltan2_GetParameter.hpp>

namespace Zoltan2{

/*! \brief Return a sublist of the given parameter list.
 *   \param pl  a Teuchos::ParameterList.
 *   \param listName  the name of a parameter list that may be a sublist
 *                          of pl.
 *   \return the requested parameter list if it exists, an empty list
 *                    otherwise.
 *
 *  If the sublist does not exist, an empty parameter list named "emptyList"
 *  is returned.  If the input list \c pl is such an empty list, then
 *  the empty list is returned.  In this way getList() can be nested when
 *  it is not known if the intermediate lists exist.  For example:
 *
   \code
         getList(getList(getParameters(), "partitioning"), "geometric")
   \endcode
 *
 * will work (by returning an empty list) even if there is
 * no "partitioning" list.
 */

const Teuchos::ParameterList & getParameterList(
  const Teuchos::ParameterList &superList, const char *listName)
{
  static Teuchos::ParameterList emptyList("emptyList");

  if (superList.name() == std::string("emptyList"))
    return superList;

  const Teuchos::ParameterEntry *sublist = superList.getEntryPtr(listName);

  if (!sublist || !sublist->isList()){
    return emptyList;
  }

  return superList.sublist(listName);
}

} // namespace Zoltan2

