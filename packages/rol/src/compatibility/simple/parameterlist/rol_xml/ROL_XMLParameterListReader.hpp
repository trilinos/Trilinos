// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_XMLPARAMETERLISTREADER_H
#define ROL_XMLPARAMETERLISTREADER_H

/*! \file ROL_XMLParameterListReader.hpp
    \brief Writes an XML object to a parameter list
*/

#include "Teuchos_ParameterList.hpp"
#include "ROL_XMLObject.hpp"


namespace ROL {


/** \brief Writes an XML object to a parameter list.
 *
 * \ingroup XML
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT XMLParameterListReader{

public:

  /** \name Public Types */
  //@{

  /** \brief . */
  static const std::string& getParameterListTagName(){
    static const std::string parameterListTagName = "ParameterList";
    return parameterListTagName;
  }

  /** \brief . */
  static const std::string& getNameAttributeName(){
    static const std::string nameAttributeName = "name";
    return nameAttributeName;
  }

  //@}

  //! @name Constructors
  //@{
  /** \brief . */
  XMLParameterListReader();
  //@}

  /** Write the given XML object to a parameter list */
  ParameterList toParameterList(const XMLObject& xml) const;

  /** \brief Set policy regarding duplicated sublists
    *
    * The default behavior of this class is to allow duplicated sublists,
    * although the resulting
    * ParameterList is undefined for the duplicated sublists (in most
    * cases, they will be merged in the order they are encountered in the
    * XML character stream).
    *
    * If set \c false, then duplicated sublists in the XML tree
    * will result in the Teuchos::DuplicateParameterSublist
    * exception being thrown.
    *
    * If set \c true, the default behavior is restored.
    */
  void setAllowsDuplicateSublists(bool policy);

  /** \brief Specifies the current policy regarding duplicated sublists.
      See setAllowsDuplicateSublists() for more details.
  */
  bool getAllowsDuplicateSublists() const;

private:

  bool _allowDuplicateSublists;

  /** \brief Write the given XML object to a parameter list 
   */
  void convertParameterList(const XMLObject& xml,
    RCP<ParameterList> parentList) const;
};



} // namespace ROL


#endif
