// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Teuchos_STRINGINPUTSOURCE_H
#define Teuchos_STRINGINPUTSOURCE_H

/*! \file Teuchos_StringInputSource.hpp
    \brief Definition of XMLInputSource derived class for reading XML from
	a std::string
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLInputSource.hpp"


namespace Teuchos
{

  using std::string;

  /** \ingroup XML
   * \brief Instantiation of XMLInputSource class for reading XML from a std::string
   */
  class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT StringInputSource : public XMLInputSource
    {
    public:
      /** \brief Constructor */
      StringInputSource(const std::string& text);

      /** \brief Destructor */
      virtual ~StringInputSource(){;}

      /** \brief Create a StringInputStream */
      virtual RCP<XMLInputStream> stream() const;

    private:
      std::string text_;
    };


}
#endif

