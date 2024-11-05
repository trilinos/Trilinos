// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PARAMETER_ENTRY_VALIDATOR_H
#define TEUCHOS_PARAMETER_ENTRY_VALIDATOR_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Describable.hpp"

namespace Teuchos {


#ifndef DOXYGEN_SHOULD_SKIP_THIS
class ParameterEntry;
#endif

/** \brief Abstract interface for an object that can validate a
 *  ParameterEntry's value.
 *
 * Not only can a validator validate and entry but it can also help to set
 * and/or adjust the default value.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterEntryValidator : public Describable
{
public:

  /** \name Public types */
  //@{

  /** \brief . */
  typedef unsigned int ValidatorID;

  /** \brief . */
  typedef RCP<const Array<std::string> > ValidStringsList;

  //@}

  /** \brief Default Constructor */
  ParameterEntryValidator() {}

  /** \brief Get a string that should be used as a value of the type attribute
   * when serializing it to XML.
   *
   * \return a string that should be used as a tag for this validator
   * when serializing it to XML.
   */
  virtual const std::string getXMLTypeName() const=0;

  /** \brief Print documentation for this parameter.
   *
   * \param docString [in] (Multi-line) documentation std::string.
   *
   * \param out [out] The std::ostream used for the output
   *
   * The purpose of this function is to augment what is
   * in <tt>docString</tt>
   * with some description of what valid values this parameter
   * validator will accept.
   */
  virtual void printDoc(
    std::string const& docString,
    std::ostream &out
    ) const = 0;

  /** \brief Return an array of strings of valid values if applicable.
   *
   * If there is no such array of std::string values that makes since, just return
   * <tt>return.get()==NULL</tt>.
   *
   * The returned strings must not contain any newlines (i.e. no <tt>'\n'</tt>
   * characters) and must be short enough to fit on one line and be readable.
   */
  virtual ValidStringsList validStringValues() const = 0;

  /** \brief Validate a parameter entry value and throw std::exception (with a
   * great error message) if validation fails.
   *
   * \param  entry
   *            [in] The ParameterEntry who's type and value is being validated
   * \param  paramName
   *            [in] The name of the ParameterEntry that is used to build error messages.
   * \param  sublistName
   *            [in] The name of the ParameterList that <tt>paramName</tt> exists in
   *            that is used to build error messages.
   */
  virtual void validate(
    ParameterEntry  const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const = 0;

  /** \brief Validate and perhaps modify a parameter entry's value.
   *
   * \param paramName [in] The name of the ParameterEntry that is used to
   * build error messages.
   *
   * \param sublistName [in] The name of the ParameterList that
   * <tt>paramName</tt> exists in that is used to build error messages.
   *
   * \param entry [in/out] The ParameterEntry who's type and value is being
   * validated and perhaps even changed as a result of calling this function.
   *
   * The default implementation simply calls <tt>this->validate()</tt>.
   */
  virtual void validateAndModify(
    std::string const& paramName,
    std::string const& sublistName,
    ParameterEntry * entry
    ) const
    {
      TEUCHOS_TEST_FOR_EXCEPT(0==entry);
      this->validate(*entry,paramName,sublistName);
    }

  double convertStringToDouble(std::string str) const
  {
    #ifdef HAVE_TEUCHOSCORE_CXX11
      size_t idx = 0;
      double result = std::stod(str, &idx); // can throw std::invalid_argument
      if(idx != str.length()) { // check for extra bad format characters
        throw std::invalid_argument( "String: '" + str + "' had bad formatting for converting to a double." );
      }
      return result;
    #else
      return std::atof(str.c_str());
    #endif
  }

  int convertStringToInt(std::string str) const
  {
    #ifdef HAVE_TEUCHOSCORE_CXX11
      size_t idx = 0;
      int result = std::stoi(str, &idx); // can throw std::invalid_argument
      if(idx != str.length()) { // check for extra bad format characters
        throw std::invalid_argument( "String: '" + str + "' had bad formatting for converting to an int." );
      }
      return result;
    #else
      return std::atoi(str.c_str());
    #endif
  }

  int convertStringToLongLong(std::string str) const
  {
    size_t idx = 0;
    long long result = std::stoll(str, &idx); // can throw std::invalid_argument
    if(idx != str.length()) { // check for extra bad format characters
      throw std::invalid_argument( "String: '" + str + "' had bad formatting for converting to a long long." );
    }
    return result;
  }

};


} // namespace Teuchos


#endif // TEUCHOS_PARAMETER_ENTRY_VALIDATOR_H
