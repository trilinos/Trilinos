// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
};


} // namespace Teuchos


#endif // TEUCHOS_PARAMETER_ENTRY_VALIDATOR_H
