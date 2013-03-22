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

#ifndef TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
#define TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP

/*! \file Teuchos_StandardParameterEntryXMLConverters.hpp
 * \brief A collection of standard ParameterEntryXMLConverters.
*/


#include "Teuchos_ParameterEntryXMLConverter.hpp"


namespace Teuchos {


/** \brief A last resort converter for when no others will do.
 *
 * Writes out a raw string representation to xml and sets ParameterEntryValues
 * as strings when they are read back in.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT AnyParameterEntryConverter : public ParameterEntryXMLConverter{

public:

  /** \name Overridden from ParameterEntryXMLConverter */
  //@{

  /** \brief . */
  const std::string getTypeAttributeValue() const;

  /** \brief . */
  const std::string getValueAttributeValue(
    RCP<const ParameterEntry> entry) const;

  /** \brief . */
  any getAny(const XMLObject& xmlObj) const; 
  //@}

};


/** \brief A standard ParameterEntryXMLConverter for most data types.
 *
 * This converter is appropriate for most data types.
 */
template<class T>
class StandardTemplatedParameterConverter : 
  public ParameterEntryXMLConverter
{

public:

  /** \name Overridden from ParameterEntryXMLConverter */
  //@{

  /** \brief . */
  virtual const std::string getTypeAttributeValue() const{
    return TypeNameTraits<T>::name();
  }

  /** \brief . */
  virtual const std::string getValueAttributeValue(
    RCP<const ParameterEntry> entry) const {
    return toString(any_cast<T>(entry->getAny(false)));
  }

  /** \brief . */
  any getAny(const XMLObject& xmlObj) const{
    return any(xmlObj.getRequired<T>(getValueAttributeName()));
  }
  
  //@}

};



} // namespace Teuchos


#endif // TEUCHOS_STANDARDPARAMETERENTRYXMLCONVERTERS_HPP
