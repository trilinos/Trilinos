// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Parameters.hpp
    \brief Defines Parameter related enumerators, methods, and validators.
*/

#ifndef _ZOLTAN2_PARAMETERS_HPP_
#define _ZOLTAN2_PARAMETERS_HPP_

#include <Zoltan2_config.h>

#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cctype>

#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Comm.hpp>

// Had to redefine this type from Teuchos_ParameterEntryValidator.hpp.
// Compiler stumbled on it.
typedef Teuchos::RCP<const Teuchos::Array<std::string> > ValidStringsList;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////
// Parameter-related namespace methods 

/*! \brief Create a ParameterList suitable for validating another list.
 *
 *   \param  plIn   a ParameterList that has not been validated yet.
 *   \param  plOut  on return, plOut has exactly the same parameters are
 *                      plIn, but it also has a validator for each parameter.
 *
 *   The validating ParameterList is required to verify that the
 *   input ParameterList is valid.
 */
void createValidatorList(
   const Teuchos::ParameterList &plIn, Teuchos::ParameterList &plOut);

/*! \brief Print out the ParameterList's built in documenation.
 *
 *   \param  pl      a ParameterList
 *   \param  os      the output stream to which the documentation should go
 *   \param listNames  leave this unset - it is set in recursive calls
 */
void printListDocumentation( const Teuchos::ParameterList &pl, std::ostream &os,
  std::string listNames=std::string(""));

////////////////////////////////////////////////////////////////////
// Parameter-related enumerated types.
//
//  If you change these enumerators, change their documentation
//  in Zoltan2_Parameters.cpp.
//

/*! \brief Level of error checking or assertions desired.
 *
 *  Each assertion in the code must have a level. Tests for
 *  logic errors should always be level DEBUG_MODE_ASSERTION.
 *  Quick tests are BASIC, longer tests for common errors are
 *  COMPLEX, and tests for unlikely errors are only done in DEBUG_MODE.
 *  The user sets the assertion level with the parameter \c error_check_level.
 */

enum AssertionLevel {
  BASIC_ASSERTION,    /*!< \brief checks that should always be done (user input) */
  COMPLEX_ASSERTION,  /*!< \brief checks that take extra time (validate a graph) */
  DEBUG_MODE_ASSERTION, /*!< \brief done when checking everything incl logic errors */
  NUM_ASSERTION_LEVELS};

/*! \brief The amount of debugging or status output to print.
 *
 *  Each debug/status message must have an output level.  The
 *  user specfies the level desired with the \c debug_level parameter.
 *
 *  If Zoltan2 is compiled with \b Z2_OMIT_ALL_STATUS_MESSAGES, no
 *  messages will be displayed, \c debug_level is ignored,
 *  and status message code is ifdef'd out.
 */
 
enum MessageOutputLevel {
  NO_STATUS,                 /*!< \brief don't display status/debug messages */
  BASIC_STATUS,              /*!< \brief the status at each high level step */
  DETAILED_STATUS, /*!< \brief include sub-steps, plus each method's entry and exit */
  VERBOSE_DETAILED_STATUS,   /*!< \brief include more detail about sub-steps */
  NUM_STATUS_OUTPUT_LEVELS};

/*! \brief Whether profiling information should be local or should include
 *             global reductions.
 *
 *  This is unused.
 * \todo It would be good for timing and memory profiling to include an
 *          option to print out local values (which would be faster) or
 *          includes global min, max, average, total, etc.
 */
enum MessageSummaryLevel{
  LOCAL_SUMMARY,              /*!< \brief messages should display local info only */
  GLOBAL_SUMMARY,             /*!< \brief include global min, max, avg, etc. */
  NUM_STATUS_SUMMARY_LEVELS};

////////////////////////////////////////////////////////////////////
// A validator for integer range lists.


/*! \brief Codes that indicate how to interpret the Array<int> representing
 *            the user's integer range list
 */
enum RangeType {
  RANGE_INCLUDES_ALL,    /*!< \brief all values were listed */
  RANGE_IS_EMPTY,        /*!< \brief no values were listed */
  RANGE_IS_LISTED,       /*!< \brief the listed values are in the Array<int> */
  NUM_RANGE_TYPES};

/*! \class IntegerRangeListValidator
 *  \brief A ParameterList validator for integer range lists
 *
 * An integer range list is a concise way to provide a list of
 * identifiers.  It is set as a string.  Valid values are:
 *
 * The template parameter is the data type of the values in the list.
 *
 * \li an integer
 * \li a range of integers given as two integers separated by a dash
 * \li the word "all"
 * \li any comma separated list of the above three
 *
 * Examples:
 *    \li "1,5,12,30-39,101"
 *    \li "all"
 *
 * Redundant specifiers are and'ed:  "1,5,all" is just "all".
 *
 * Typical use cases for an integer range list are:
 *   \li the list of processes that are to output debugging information
 *   \li the list of fixed vertex IDs in a partitioning operation
 *
 * At the call to validateAndModify(), the integer range list parameter 
 * value is changed from a string to an Array<Integral> which encodes
 * the meaning of the string.
 *
 * Helper functions for interpreting the integer range list after it has
 * been validated are:
 *
 *\li validIntegralRangeList(const Teuchos::Array<Integral> &vals)
 *
 *\li allValuesAreInRangeList(const Teuchos::Array<Integral> &vals)
 *
 *\li allValuesAreInRangeList(const Teuchos::ParameterEntry &e)
 *
 *\li noValuesAreInRangeList(const Teuchos::Array<Integral> &vals)
 *
 *\li noValuesAreInRangeList(const Teuchos::ParameterEntry &e)
 *
 *\li IsInRangeList(const Integral val, const Teuchos::Array<Integral> &valList)
 *
 *\li IsInRangeList(const Integral val, const Teuchos::ParameterEntry &e)
 *
 *\li printIntegralRangeList(std::ostream &os, Teuchos::Array<Integral> &irl)
 */

template <typename Integral>
  class IntegerRangeListValidator : public Teuchos::ParameterEntryValidator
{
private:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  Integral min_;
  Integral max_;

  static const std::string listDelim_;
  static const std::string rangeDelim_;
  static const std::string allText_;

  static void checkValid(char c); 
  static bool listSaysAll(std::string &l);
  static int breakRange(std::string &range, std::string &from, std::string &to);

#endif

public:
  /*! \brief Constructor: any Integral is valid
   */
  IntegerRangeListValidator();

  /*! \brief Constructor: only Integrals in the specified range are valid.
   *   \param validMin  all values implied by the integer range list
   *                          must be bounded by this minimum.
   *   \param validMax  all values implied by the integer range list
   *                          must be bounded by this maximum.
   */
  IntegerRangeListValidator(Integral validMin, Integral validMax); 

  // Implementation of ParameterEntryValidator interface

  const std::string getXMLTypeName() const; 

  void printDoc(std::string const& docString, std::ostream &out) const;

  ValidStringsList validStringValues() const ;

  void validate( Teuchos::ParameterEntry  const& entry,
    std::string const& paramName, std::string const& sublistName
    ) const;

  void validateAndModify( std::string const& paramName,
    std::string const& sublistName, Teuchos::ParameterEntry * entry
    ) const;
}; // end class

/*! \brief A helper function that indicates whether an array is a valid
 *           integer range list.
 *
 *    \param vals   An array that may encode an integer range list.
 *    \return        true if the array encodes such a list, false otherwise.
 */
template <typename Integral>
  bool validIntegralRangeList(const Teuchos::Array<Integral> &vals);

/*! \brief  A helper function that determines if all values are in the list.
 *
 *    \param vals   An array encoding an integer range list.
 *    \return        true if the array encoding implies all values,
 *                  false otherwise.
 *
 *  If the array is not a valid encoding of an integer range list,
 *  a std::runtime_error will be thrown.
 *
 *   If the user's parameter value was "all", then after validation
 *   the parameter value will be an array of size one containing
 *   a code that indicates all values are included, and this
 *   function will return \c true.  
 */
template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::Array<Integral> &vals);

/*! \brief  A helper function that determines if all values are in the list.
 *
 *    \param e    A parameter entry
 *    \return     true if the entry value is an array encoding all values,
 *                  false otherwise.
 *
 *  If the entry value is not a valid encoding of an integer range list,
 *  a std::runtime_error will be thrown.
 *
 *   If the user's parameter value was "all", then after validation
 *   the parameter value will be an array of size one containing
 *   a code that indicates all values are included, and this
 *   function will return \c true.  
 */
template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::ParameterEntry &e);

/*! \brief  A helper function that determines if no values are in the list.
 *
 *    \param vals   An array encoding an integer range list.
 *    \return        true if the array encoding implies no values,
 *                  false otherwise.
 *
 *  If the array is not a valid encoding of an integer range list,
 *  a std::runtime_error will be thrown.
 *
 *   If the user's parameter value was empty, then after validation
 *   the parameter value will be an array of size one containing
 *   a code that indicates no values are included, and this
 *   function will return \c true.  
 */
template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::Array<Integral> &vals);

/*! \brief  A helper function that determines if no values are in the list.
 *
 *    \param e    A parameter entry
 *    \return     true if the entry value is an array encoding no values,
 *                  false otherwise.
 *
 *  If the entry value is not a valid encoding of an integer range list,
 *  a std::runtime_error will be thrown.
 *
 *   If the user's parameter value was empty, then after validation
 *   the parameter value will be an array of size one containing
 *   a code that indicates no values are included, and this
 *   function will return \c true.  
 */
template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::ParameterEntry &e);

/*! \brief  A helper function that determines if a value is in the list.
 *
 *    \param val  A value that could be in the list.
 *    \param e    A parameter entry
 *    \return     true if the entry value implies \c val is in the list,
 *                  false otherwise.
 *
 *  If the entry value is not a valid encoding of an integer range list,
 *  a std::runtime_error will be thrown.
 */
template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::ParameterEntry &e);

/*! \brief  A helper function that determines if a value is in the list.
 *
 *    \param val  A value that could be in the list.
 *    \param vals   An array encoding an integer range list.
 *    \return     true if the encoding of \c vals implies \c val is in the list,
 *                  false otherwise.
 *
 *  If \c vals is not a valid encoding of an integer range list,
 *  a std::runtime_error will be thrown.
 */
template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::Array<Integral> &vals);

/*! \brief  A helper function that prints the meaning of an encoded
 *             integer range list
 *
 *    \param os   An output stream to which the list will be printed.
 *    \param val  A value of an integer range list parameter after it has
 *                     been validated.
 *
 */
template <typename Integral>
  void printIntegralRangeList(std::ostream &os, Teuchos::Array<Integral> &irl);

////////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////////

template <typename Integral>
const std::string IntegerRangeListValidator<Integral>::listDelim_(",");

template <typename Integral>
const std::string IntegerRangeListValidator<Integral>::rangeDelim_("-");

template <typename Integral>
const std::string IntegerRangeListValidator<Integral>::allText_("all");

template <typename Integral>
  void IntegerRangeListValidator<Integral>::checkValid(char c)
{
  if (std::isspace(c) || std::isdigit(c) || (c == ',') || (c == '-'))
    return;
  else
    throw std::runtime_error("invalid integer range list");
}

template <typename Integral>
  bool IntegerRangeListValidator<Integral>::listSaysAll(std::string &l)
{
  std::transform(l.begin(), l.end(), l.begin(), (int(*)(int)) tolower);
  if (l.find(allText_) != std::string::npos)
    return true;  // "all" is in the string
  else
    return false;
}

template <typename Integral>
  int IntegerRangeListValidator<Integral>::breakRange(
    std::string &range, std::string &from, std::string &to)
{
  from.clear();
  to.clear();
  std::string::size_type loc = range.find(rangeDelim_);
  if (loc == std::string::npos){
    from = range;
  }
  else{
    from = range.substr(0, loc);
    to = range.substr(loc+1, range.size());
  }
  long a, b;
  std::istringstream iss1(from);
  iss1 >> a;
  b = a;
  if (to.size() > 0){
    std::istringstream iss2(to);
    iss2 >> b;
    if (b < a)
      std::swap(from,to);
  }
  return (b != a) ? 2 : 1;
}


template <typename Integral>
  IntegerRangeListValidator<Integral>::IntegerRangeListValidator(): 
    min_(1), max_(0)
{
}

template <typename Integral>
  IntegerRangeListValidator<Integral>::IntegerRangeListValidator(
    Integral validMin, Integral validMax) :
      min_(validMin), max_(validMax)
{
  if (min_ < max_) std::swap(min_,max_);
}

  // Implementation of ParameterEntryValidator interface

template <typename Integral>
  const std::string 
    IntegerRangeListValidator<Integral>::getXMLTypeName() const 
{
  std::string className("IntegerRangeListValidator");
  std::string classType("("+Teuchos::TypeNameTraits<Integral>::name()+")");
  return std::string(className + classType);
}

template <typename Integral>
  void IntegerRangeListValidator<Integral>::printDoc(
    std::string const& docString, std::ostream &out) const  
{
  Teuchos::StrUtils::printLines(out,"# ",docString);
  out << "#\tAn integer range list is a string which can contain:\n";
  out << "#\t\tthe text \"all\", which indicates all values\n";
  out << "#\t\ta list of integer ranges separated by commas.\n";
  out << "#\tA range is one value, or two values separated by a dash.\n";
  out << "#\tExample: \"all\" or \"1-10\" or \"3, 10-12\" or \"25\"\n";
  if (max_ >= min_){
    out << "#\tThe range of valid integers is [";
    out << min_ << "," << max_ << "]\n";
  }
}

template <typename Integral>
  ValidStringsList IntegerRangeListValidator<Integral>::validStringValues() const 
{ 
  return Teuchos::null; 
}

template <typename Integral>
  void IntegerRangeListValidator<Integral>::validate( 
    Teuchos::ParameterEntry  const& entry,
    std::string const& paramName, std::string const& sublistName) const
{
  if (!entry.isType<std::string>()){
    return;  // already converted to an an array
  }
  std::string *sptr=NULL;
  std::string &rangeList = entry.getValue(sptr);
  std::string inValue(rangeList);

  if (listSaysAll(inValue))
    return;  // "all" is in the string

  // throw error if invalid integer range list
  std::for_each(inValue.begin(), inValue.end(), checkValid);

  if (max_ >= min_){
    std::string::const_iterator rangeBegin = inValue.begin();
    std::string::const_iterator valueEnd = inValue.end();

    while (rangeBegin != valueEnd){
      std::string::const_iterator rangeEnd = std::search(
        rangeBegin, valueEnd, listDelim_.begin(), listDelim_.end());
      std::string range(rangeBegin, rangeEnd);
      std::string aHalf, bHalf;
      int count = breakRange(range, aHalf, bHalf);

      Integral a, b;
      std::istringstream iss1(aHalf);
      iss1 >> a;
      if (count > 1){
        std::istringstream iss2(bHalf);
        iss2 >> b;
      }
      else
        b = a;

      if ((a < min_) || (b > max_)){
        std::ostringstream oss;
        oss << "input range [" << a << "," << b << "] ";
        oss << "exceeds valid range [" << min_ << "," << max_ << "] ";
        throw std::runtime_error(oss.str());
      }
      if (rangeEnd == valueEnd)
        rangeBegin = rangeEnd;
      else
        rangeBegin = ++rangeEnd;
    }
  }
}

template <typename Integral>
  void IntegerRangeListValidator<Integral>::validateAndModify( 
    std::string const& paramName, std::string const& sublistName, 
    Teuchos::ParameterEntry * entry) const
{
  typedef typename Teuchos::Array<Integral>::size_type arraySize_t;
  if (!entry->isType<std::string>()){
    return;
  }
  
  std::string *sptr=NULL;
  std::string &rangeList = entry->getValue(sptr);
  Teuchos::Array<Integral> valueList;

  std::string inValue(rangeList);

  if (listSaysAll(inValue)){
    valueList.push_back(RANGE_INCLUDES_ALL);
  }
  else{
    // throw error if invalid integer range list
    std::for_each(inValue.begin(), inValue.end(), checkValid);

    std::string::const_iterator rangeBegin = inValue.begin();
    std::string::const_iterator valueEnd = inValue.end();

    while (rangeBegin != valueEnd){
      std::string::const_iterator rangeEnd = std::search(rangeBegin,
        valueEnd, listDelim_.begin(), listDelim_.end());
      std::string range(rangeBegin, rangeEnd);
      std::string aHalf, bHalf;
      int count = breakRange(range, aHalf, bHalf);

      Integral a, b;
      std::istringstream iss1(aHalf);
      iss1 >> a;
      if (count > 1){
        std::istringstream iss2(bHalf);
        iss2 >> b;
      }
      else
        b = a;

      if ((max_ >= min_) && ((a < min_) || (b > max_))){
        std::ostringstream oss;
        oss << "input range [" << a << "," << b << "] ";
        oss << "exceeds valid range [" << min_ << "," << max_ << "] ";
        throw std::runtime_error(oss.str());
      }
      for (Integral i=a; i <=b; i++)
        valueList.push_back(i);
      if (rangeEnd == valueEnd)
       rangeBegin = rangeEnd;
      else
        rangeBegin = ++rangeEnd;
    }
    if (valueList.size() > 1){  // sort & remove duplicates
      std::sort(valueList.begin(), valueList.end());
      arraySize_t listEnd = 0;
      arraySize_t length = valueList.size();
      for (arraySize_t i=1; i < length; i++){
        if (valueList[i] > valueList[listEnd]){
          listEnd++;
          if (listEnd != i)
            valueList[listEnd] = valueList[i];
        }
      }
      if (++listEnd < length)
        valueList.resize(listEnd);
    }

    Integral flag = RANGE_IS_LISTED;
    if (valueList.size() == 0){
      flag = RANGE_IS_EMPTY;
    }
    else if (max_ >= min_){
      Integral allSize = max_ - min_ + 1;
      if (valueList.size() == allSize){
        flag = RANGE_INCLUDES_ALL;
        valueList.clear();
      }
    }
    valueList.push_back(flag);
  }
  entry->setValue(valueList);
}

// Helpers for IntegralRangeList parameter type

template <typename Integral>
  bool validIntegralRangeList(const Teuchos::Array<Integral> &vals)
{
  typedef typename Teuchos::Array<Integral>::size_type arraySize_t;
  arraySize_t len = vals.size();
  if (len==0)
    return false;

  Integral flag = vals[len-1];
  if ((flag != RANGE_INCLUDES_ALL) && (flag != RANGE_IS_EMPTY) &&
      (flag != RANGE_IS_LISTED))
    return false;

  return true;
}

template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::Array<Integral> &vals)
{
  if (!validIntegralRangeList(vals))
    throw std::runtime_error("list is not a valid range list");
  return vals[vals.size()-1] == RANGE_INCLUDES_ALL;
}

template <typename Integral>
  bool allValuesAreInRangeList(const Teuchos::ParameterEntry &e)
{
  if (!e.isType<Teuchos::Array<Integral> >())
    throw std::runtime_error("Should not call until modified");

  Teuchos::Array<Integral> *valPtr=NULL;
  Teuchos::Array<Integral> &vals = e.getValue(valPtr);
  return allValuesAreInRangeList(vals);
}

template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::Array<Integral> &vals)
{
  if (!validIntegralRangeList(vals))
    throw std::runtime_error("list is not a valid range list");
  return vals[vals.size()-1] == RANGE_IS_EMPTY;
}

template <typename Integral>
  bool noValuesAreInRangeList(const Teuchos::ParameterEntry &e)
{
  if (!e.isType<Teuchos::Array<Integral> >())
    throw std::runtime_error("Should not call until modified");

  Teuchos::Array<Integral> *valPtr=NULL;
  Teuchos::Array<Integral> &vals = e.getValue(valPtr);
  return noValuesAreInRangeList(vals);
}

// TODO :
// inList(std::vector<Integral> &val, std::vector<bool> &result)

template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::Array<Integral> &valList)
{
  if (allValuesAreInRangeList(valList))
    return true;
  else if (noValuesAreInRangeList(valList))
    return false;

  typename Teuchos::Array<Integral>::const_iterator flag = valList.end();
  --flag;
  if (std::binary_search(valList.begin(), flag, val))
    return true;
  else
    return false;
}

template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::ParameterEntry &e)
{
  if (!e.isType<Teuchos::Array<Integral> >())
    throw std::runtime_error("Should not call until modified");

  Teuchos::Array<Integral> *valPtr=NULL;
  Teuchos::Array<Integral> &valList = e.getValue(valPtr);

  return IsInRangeList(val, valList);
} 

template <typename Integral>
  void printIntegralRangeList(std::ostream &os, Teuchos::Array<Integral> &irl)
{
  if (Zoltan2::allValuesAreInRangeList(irl))
    os << "all";
  else if (Zoltan2::noValuesAreInRangeList(irl))
    os << "empty";
  else{
    Teuchos::ArrayView<const Integral> view = 
      irl.view(0, irl.size()-1); // skip last value, it's a flag
    os << view;
  }    
}

}  // end of namespace Zoltan2

#endif
