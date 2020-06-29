// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_IntegerRangeList.hpp
    \brief Define IntegerRangeList validator.
*/

#ifndef _ZOLTAN2_INTEGERRANGELIST_HPP_
#define _ZOLTAN2_INTEGERRANGELIST_HPP_

#include <Zoltan2_config.h>

#include <Teuchos_ParameterEntryValidator.hpp>
#include <Teuchos_ValidatorXMLConverter.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_ValidatorMaps.hpp>
#include <Teuchos_DummyObjectGetter.hpp>
#include <Teuchos_StrUtils.hpp>

#ifdef _MSC_VER
// for isspace(int), else one gets isspace(_Elem _Ch, const locale& _Loc) from <locale>
#include <cctype>
#endif

// Had to redefine this type from Teuchos_ParameterEntryValidator.hpp.
// Compiler stumbled on it.
typedef Teuchos::RCP<const Teuchos::Array<std::string> > ValidStringsList;

namespace Zoltan2{

/*! \brief Codes that indicate how to interpret the Array<int> representing
 *            the user's integer range list
 */
enum RangeType {
  RANGE_INCLUDES_ALL,    /*!< \brief all values were listed */
  RANGE_IS_EMPTY,        /*!< \brief no values were listed */
  RANGE_IS_LISTED,       /*!< \brief the listed values are in the Array<int> */
  NUM_RANGE_TYPES
};

template <typename T> class IntegerRangeListValidator;

////////////////////////////////////////////////////////////////////
// Helper functions for interpreting integer range list arrays.
////////////////////////////////////////////////////////////////////

/*! \brief A helper function that indicates whether an array is a valid
 *           integer range list.
 *
 *    \param vals   An array that may encode an integer range list.
 *    \return        true if the array encodes such a list, false otherwise.
 */
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
  bool allValuesAreInRangeList(const Teuchos::Array<Integral> &vals)
{
  if (!validIntegralRangeList(vals))
    throw std::runtime_error("list is not a valid range list");
  return vals[vals.size()-1] == RANGE_INCLUDES_ALL;
}

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
  bool allValuesAreInRangeList(const Teuchos::ParameterEntry &e)
{
  if (!e.isType<Teuchos::Array<Integral> >())
    throw std::runtime_error("Should not call until modified");

  Teuchos::Array<Integral> *valPtr=NULL;
  Teuchos::Array<Integral> &vals = e.getValue(valPtr);
  return allValuesAreInRangeList(vals);
}

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
  bool noValuesAreInRangeList(const Teuchos::Array<Integral> &vals)
{
  if (!validIntegralRangeList(vals))
    throw std::runtime_error("list is not a valid range list");
  return vals[vals.size()-1] == RANGE_IS_EMPTY;
}

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

/*! \brief  A helper function that determines if a value is in the list.
 *
 *    \param val  A value that could be in the list.
 *    \param vals   An array encoding an integer range list.
 *    \param sorted  Set to false if the integer range list is not sorted
 *    \return     true if the encoding of \c vals implies \c val is in the list,
 *                  false otherwise.
 *
 *  If \c vals is not a valid encoding of an integer range list,
 *  a std::runtime_error will be thrown.
 */
template <typename Integral>
  bool IsInRangeList(const Integral val, const Teuchos::Array<Integral> &valList, bool sorted=true)
{
  if (allValuesAreInRangeList(valList))
    return true;
  else if (noValuesAreInRangeList(valList))
    return false;

  if (sorted){
    typename Teuchos::Array<Integral>::const_iterator flag = valList.end();
    --flag;
    if (std::binary_search(valList.begin(), flag, val))
      return true;
    else
      return false;
  }
  else{
    for (typename Teuchos::Array<Integral>::size_type i=0; i < valList.size()-1; i++){
      if (valList[i] == val)
        return true;
    }
    return false;
  }
}

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
  bool IsInRangeList(const Integral val, const Teuchos::ParameterEntry &e)
{
  if (!e.isType<Teuchos::Array<Integral> >())
    throw std::runtime_error("Should not call until modified");

  Teuchos::Array<Integral> *valPtr=NULL;
  Teuchos::Array<Integral> &valList = e.getValue(valPtr);

  typedef IntegerRangeListValidator<Integral> irl_t;
  RCP<const irl_t> irl;
  bool fail = false;
 
  RCP<const Teuchos::ParameterEntryValidator> valtor = e.validator();
  if (!valtor.is_null()){
    try{
      irl = rcp_dynamic_cast<const irl_t>(valtor);
    }
    catch (...){
      fail = true;
    }
  }
  else{
    fail = true;
  }

  if (fail)
    throw std::runtime_error("wrong type of parameter entry");

  bool sorted = irl->inputListWillBeSorted();

  return IsInRangeList(val, valList, sorted);
} 

/*! \brief  A helper function that returns a view of the list.
 *
 *    \param irl The value of an integer range list parameter after it has
 *                 been validated
 *    \return If the user's parameter value did not imply "all" or "none",
 *        then a view of the array of integers implied by the value is
 *        returned.  Otherwise an empty array is returned.
 */   
template <typename Integral>
  Teuchos::ArrayView<Integral> getList(Teuchos::Array<Integral> &irl)
{
  Teuchos::ArrayView<Integral> av;  // av.size() == 0

  if (!Zoltan2::allValuesAreInRangeList(irl) &&
      !Zoltan2::noValuesAreInRangeList(irl))
    av = irl.view(0, irl.size()-1); // skip last value, it's a flag

  return av;
}

/*! \brief  A helper function that prints the meaning of an encoded
 *             integer range list
 *
 *    \param os   An output stream to which the list will be printed.
 *    \param val  A value of an integer range list parameter after it has
 *                     been validated.
 *
 */
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

////////////////////////////////////////////////////////////////////
// Validator class
////////////////////////////////////////////////////////////////////

/*! \class IntegerRangeListValidator
 *  \brief A ParameterList validator for integer range lists
 *
 * An integer range list is a concise way to provide a list of
 * integers.  It is set as a string.  Valid values are:
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
 * A constructor flag determines how the list is processed at
 * validateAndModify() time.  Either the list is sorted and duplicates
 * are removed, or the list remains as it was entered by the user.
 * For example, if the list is to be modified:

        - "1,5,2,1" becomes "1,2,5"
        - "1-10,9-15" becomes 1-15"

 * Typical use cases for an integer range list are:

 *   - the list of processes that print status information
 *   - the list of fixed vertex IDs in a partitioning operation
 *
 * At the call to validateAndModify(), the integer range list parameter 
 * value is changed from a string to an Array<Integral> which encodes
 * the meaning of the string. (The last value in the array specifies 
 * whether the values are listed, whether the user requested "all",
 * or whether the first and last value of a range is provided.)
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
 *\li IsInRangeList(const Integral val, const Teuchos::Array<Integral> &valList,
 *      bool sorted=true)
 *
 *\li IsInRangeList(const Integral val, const Teuchos::ParameterEntry &e)
 *
 *\li Teuchos::ArrayView<Integral> getList(Teuchos::Array<Integral> &irl)
 *
 *\li printIntegralRangeList(std::ostream &os, Teuchos::Array<Integral> &irl)
 *
 * The template parameter is the data type of the values in the list.
 */

template <typename Integral>
  class IntegerRangeListValidator : public Teuchos::ParameterEntryValidator
{
private:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  Integral min_;
  Integral max_;
  bool unsorted_;

  static const std::string listDelim_;
  static const std::string rangeDelim_;
  static const std::string allText_;

  static void checkValid(char c); 
  static bool listSaysAll(std::string &l);
  static int breakRange(std::string &range, std::string &from, std::string &to);

#endif

public:
  /*! \brief Constructor: any Integral is valid
   *   \param unsorted normally the input integers will be sorted and
   *       duplicates will be removed.  If
   *        this is not the desired behavior, then set \c unsorted to true.
   */
  IntegerRangeListValidator(bool unsorted=false);

  /*! \brief Constructor: only Integrals in the specified range are valid.
   *   \param validMin  all values implied by the integer range list
   *                          must be bounded by this minimum.
   *   \param validMax  all values implied by the integer range list
   *                          must be bounded by this maximum.
   *   \param unsorted normally the input integers will be sorted and
   *       duplicates will be removed..  If
   *        this is not the desired behavior, then set \c unsorted to true.
   */
  IntegerRangeListValidator(Integral validMin, Integral validMax,
    bool unsorted=false); 

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

  // IntegerRangeListValidator methods

  /*! \brief Return the minimum value permitted in the list.
   *
   *  If getAllowedMinimum() > getAllowedMaximum(), then there are
   *  no limits on the integer values in the list.
   */
  Integral getAllowedMinimum() const { return min_;}

  /*! \brief Return the maximum value permitted in the list.
   *
   *  If getAllowedMinimum() > getAllowedMaximum(), then there are
   *  no limits on the integer values in the list.
   */
  Integral getAllowedMaximum() const { return max_;}

  /*! \brief Return whether the list is sorted or not.
   *
   *   By default, when the parameter value (a list of integers and integer ranges),
   *   is processed it is sorted and duplicates are removed.  A constructor argument
   *   can be set so that the list is not sorted and duplicates are not removed.
   */  
  bool inputListWillBeSorted() const { return !unsorted_;}

}; // end class

////////////////////////////////////////////////////////////////////
// Class definitions
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
  IntegerRangeListValidator<Integral>::IntegerRangeListValidator(
    bool unsorted): min_(1), max_(0), unsorted_(unsorted)
{
}

template <typename Integral>
  IntegerRangeListValidator<Integral>::IntegerRangeListValidator(
    Integral validMin, Integral validMax, bool unsorted) :
      min_(validMin), max_(validMax), unsorted_(unsorted)
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
    std::string const& /* paramName */, std::string const& /* sublistName */) const
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
    std::string const& /* paramName */, std::string const& /* sublistName */, 
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
    if (!unsorted_ && valueList.size() > 1){  // sort & remove duplicates
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

////////////////////////////////////////////////////////////////////
// Parameter entry validator <-> XML conversion
////////////////////////////////////////////////////////////////////

/*! \brief XML conversion code for IntegerRangeListValidator
 *
 * The valid XML representation of an IntegerRangeListValidator is:
  \code
  <Validator type="IntegerRangeListValidator(template-arg)"
   min="optional minimum value"
   max="optional maximum value"
   unsorted="true if input list should remain unsorted"
   validatorId="Validator Id"
  />
  \endcode
 */
template <typename Integral>
class IntegerRangeListValidatorXMLConverter : public Teuchos::ValidatorXMLConverter
{

public:

  RCP<Teuchos::ParameterEntryValidator> convertXML(
    const Teuchos::XMLObject& xmlObj,
    const Teuchos::IDtoValidatorMap& validatorIDsMap) const;

  void convertValidator(
    const RCP<const Teuchos::ParameterEntryValidator> validator,
    Teuchos::XMLObject& xmlObj,
    const Teuchos::ValidatortoIDMap& validatorIDsMap) const;

#ifdef HAVE_TEUCHOS_DEBUG
  /** \brief . */
  RCP<const Teuchos::ParameterEntryValidator> getDummyValidator() const{
    return Teuchos::DummyObjectGetter<IntegerRangeListValidator<Integral> >::getDummyObject();
  }
#endif
};

template <typename Integral>
  RCP<Teuchos::ParameterEntryValidator>
   IntegerRangeListValidatorXMLConverter<Integral>::convertXML(
     const Teuchos::XMLObject& xmlObj, 
     const Teuchos::IDtoValidatorMap& /*validatorIDsMap*/) const
{
  Integral minValue=0, maxValue=0;
  bool unsorted=false, hasMin=false, hasMax=false;

  if (xmlObj.hasAttribute(std::string("min"))) {
    minValue = xmlObj.getRequired<Integral>(std::string("min"));
    hasMin = true;
  }

  if (xmlObj.hasAttribute(std::string("max"))) {
    maxValue = xmlObj.getRequired<Integral>(std::string("max"));
    hasMax = true;
  }

  if (xmlObj.hasAttribute(std::string("unsorted"))) 
    unsorted = xmlObj.getRequired<bool>(std::string("unsorted"));

  RCP<Teuchos::ParameterEntryValidator> toReturn;

  if (hasMin && hasMax)
    toReturn = rcp(new IntegerRangeListValidator<Integral>(minValue, maxValue, unsorted));
  else if (!hasMin && !hasMax)
    toReturn = rcp(new IntegerRangeListValidator<Integral>(unsorted));
  else
    throw std::runtime_error("invalid XML representation");

  return toReturn;
}

template<typename Integral>
void IntegerRangeListValidatorXMLConverter<Integral>::convertValidator(
  const RCP<const Teuchos::ParameterEntryValidator > validator,
  Teuchos::XMLObject& xmlObj,
  const Teuchos::ValidatortoIDMap& /*validatorIDsMap*/) const
{
  RCP<const IntegerRangeListValidator<Integral> > castedValidator =
    rcp_dynamic_cast<const IntegerRangeListValidator<Integral> >(
      validator, true);

  Integral minValue = castedValidator->getAllowedMinimum();
  Integral maxValue = castedValidator->getAllowedMaximum();
  bool unsorted = castedValidator->inputListWillBeSorted();

  if (minValue < maxValue){
    xmlObj.addAttribute<Integral>(std::string("min"), minValue);
    xmlObj.addAttribute<Integral>(std::string("max"), maxValue);
  }

  xmlObj.addAttribute<bool>(std::string("unsorted"), unsorted);
}

}
 // end of namespace Zoltan2

#endif
