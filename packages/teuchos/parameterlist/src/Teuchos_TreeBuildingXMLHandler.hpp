// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Teuchos_TREEBUILDINGXMLHANDLER_H
#define Teuchos_TREEBUILDINGXMLHANDLER_H

/*! \file Teuchos_TreeBuildingXMLHandler.hpp
    \brief Defines a class for assembling an XMLObject from XML input
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLObject.hpp"
#include <list>
#include <stack>


namespace Teuchos
{

  /** \ingroup XML
   * \brief TreeBuildingXMLHandler assembles a XMLObject from your XML input.
   */
  class TreeBuildingXMLHandler
    {
      typedef Teuchos::map<std::string, std::string> Map;
    public:
      /** \brief Empty constructor */
      TreeBuildingXMLHandler();

      /** \brief Retrieve the entire XML tree */
      const XMLObject& getObject() const
        {
          // valid XML requires a root object; we will allow empty XML files
          TEUCHOS_TEST_FOR_EXCEPTION(
            root_.isEmpty(), EmptyXMLError
            ,"Error, An empty XML std::string or file was specified!"
            "  The XML specification requires at minimum the presence"
            " of a root element!" );
          return root_.getChild(0);
        }

      /** \brief Process character data */
      void characters(const std::string& chars);

      /** \brief Receive notification of the end of an element */
      int endElement(const std::string& tag);

      /** \brief Receive notification of the start of an element */
      void startElement(const std::string& tag,
                        const Map& attributes);

    private:
      XMLObject root_;
      XMLObject current_;
      std::stack<XMLObject> path_;
    };
}


#endif

