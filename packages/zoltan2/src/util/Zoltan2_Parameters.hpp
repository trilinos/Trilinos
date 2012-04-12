// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Parameters.hpp
 *  \brief Convenience methods for working with the parameter list.
 */

#ifndef ZOLTAN2_PARAMETERS_HPP
#define ZOLTAN2_PARAMETERS_HPP

#include <Zoltan2_Standards.hpp>
#include <Teuchos_ParameterList.hpp>

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
  const Teuchos::ParameterList &superList, const char *listName) const
{
  static Teuchos::ParameterList emptyList("emptyList");

  if (superList.name() == std::string("emptyList"))
    return superList;

  const Teuchos::ParameterEntry *sublist = superList.getEntryPtr(listName);

  if (!sublist || !sublist->isList()){
    return &emptyList;
  }

  return superList.sublist(listName);
}

/*! \brief Find the value of the named parameter in the list.
 *  \param pl A parameter list that may contain the parameter.
 *  \param name  The name of the parameter entry.
 *  \param set  On return, true if the parameter is set and false if
 *               it is not set (does not appear in the parameter list).
 *  \param value On return, if the entry was found, this will be set
 *                    to the value of the entry.  Otherwise it is
 *                    undefined.
 */

template <typename T>
  void getParameterValue(const Teuchos::ParameterList &pl,
    const char *name, bool &isSet, T &value) const
{
  isSet = false;
  value = T(0);

  if (pl.name() == std::string("emptyList"))
    return;

  const Teuchos::ParameterEntry *pe = pl.getEntryPtr(name);

  if (!pe)
    return;

  isSet = true;
  value = pe->getValue<T>(&value);
}

/*! \brief Find the value of a second level parameter.
 *  \param name1  The sublist containing the parameter.
 *  \param name2  The name of the parameter entry.
 *  \param set  On return, true if the parameter is set and false if
 *               it is not set (does not appear in the parameter list).
 *  \param value On return, if the entry was found, this will be set
 *                    to the value of the entry.  Otherwise it is
 *                    undefined.
 */

template <typename T>
  void getParameterValue(const Teuchos::ParameterList &pl,
    const char *name1, const char *name2, 
    bool &isSet, T &value) const
{
  getParameterValue(getParameterList(pl, name1), name2, isSet, value);
}

/*! \brief Find the value of a third level parameter.
 *  \param name1  The top level sublist.
 *  \param name2  The sublist in \c name1 that contains the parameter.
 *  \param name3  The name of the parameter entry.
 *  \param set  On return, true if the parameter is set and false if
 *               it is not set (does not appear in the parameter list).
 *  \param value On return, if the entry was found, this will be set
 *                    to the value of the entry.  Otherwise it is
 *                    undefined.
 */

template <typename T>
  void getParameterValue(const Teuchos::ParameterList &pl,
    const char *name1, const char *name2, 
    bool &isSet, T &value) const
{
  getParameterValue(
    getParameterList(getParameterList(pl, name1), name2), 
    name3, isSet, value);
}

} // namespace Zoltan2

#endif
