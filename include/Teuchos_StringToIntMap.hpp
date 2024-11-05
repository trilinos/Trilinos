// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STRING_TO_INT_MAP_HPP
#define TEUCHOS_STRING_TO_INT_MAP_HPP

#include "Teuchos_Assert.hpp"
#include <map>
#include <string>

namespace Teuchos {

/** \brief Map a std::string to an enumeration.
 *
 * The purpose of this class is to simplify mapping a standard std::string
 * to an integer which can be interpreted as an enumeration.
 *
 * Here is an example of its use.
 \verbatim

  const int n_opt = 3;
  enum MyOptEnum {
    OPT_ONE
    ,OPT_TWO
    ,OPT_THREE
  };  // NOTE: Must be 0, 1,..., n_opt - 1
  const char* MyOptStrings[n_opt] = {
    "OPT_ONE
    ,"OPT_TWO"
    ,"OPT_THREE"
  }; // NOTE: parallels enums in MyOptEnum
  StringToIntMap my_enum_map( "opt_map", n_opt, NyOptStrings );
  ...
  switch( my_enum_map.get<MyOptEnum>("OPT_ONE") ) {
    case OPT_ONE:
      // do stuff
    case OPT_TWO:
      // do stuff
    case OPT_THREE:
      // do stuff
    default:
      // ???
  }

 \endverbatim
 *
 * The number of strings passed to the constructor must equal the number of
 * options in the enumeration.  If there are duplicate strings
 * (capitalization concidered) then the std::exception <tt>AlreadyExists</tt> is
 * throw.  If a std::string that was not passed in the constructor if given to
 * <tt>operator()( const std::string& str )</tt> then the std::exception
 * <tt>DoesNotExist</tt> is thrown.
 *
 * In the constructor, <tt>defaultGroupName</tt> is used in error messages in
 * the exceptions thrown to help make since out of the message.
 *
 * The default constructor is not defined and not to be called.
 */
class StringToIntMap {
public:

  /** \brief . */
  class AlreadyExists : public std::logic_error
  {public: AlreadyExists(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  class DoesNotExist : public std::logic_error
  {public: DoesNotExist(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /** \brief . */
  StringToIntMap( const std::string& defaultGroupName, int n, const char* strings[] );

  /** \brief . */
  int get( const std::string& option, const std::string& groupName = "" ) const;

  /** \brief . */
  template<class EnumType>
  EnumType get( const std::string& option, const std::string& groupName = "" ) const;

  /** \brief . */
  const std::string& defaultGroupName() const;

private:

  typedef std::map< std::string, int > map_t;  // all share implementation.
  std::string defaultGroupName_;
  map_t map_;

  std::string validSelections() const;

  // not defined and not to be called.
  StringToIntMap();

};  // end class StringToIntMap

/** \brief Nonmember get function.
 * \relates StringToIntMap
 */
template<class EnumType>
inline
EnumType get(
  StringToIntMap const& theMap
  ,std::string const& option
  ,std::string const& groupName = ""
  )
{
  return static_cast<EnumType>(theMap.get(option,groupName));
}

// ////////////////////////////////////////////
// Inline declarations

template<class EnumType>
inline
EnumType StringToIntMap::get( const std::string& option, const std::string& groupName ) const
{
  return static_cast<EnumType>(get(option,groupName));
}

inline
const std::string& StringToIntMap::defaultGroupName() const
{
  return defaultGroupName_;
}

}  // end namespace Teuchos

#endif  // TEUCHOS_STRING_TO_INT_MAP_HPP
