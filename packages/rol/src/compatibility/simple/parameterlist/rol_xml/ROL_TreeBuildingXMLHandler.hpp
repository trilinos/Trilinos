// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TREEBUILDINGXMLHANDLER_H
#define ROL_TREEBUILDINGXMLHANDLER_H

/*! \file ROL_TreeBuildingXMLHandler.hpp
    \brief Defines a class for assembling an XMLObject from XML input
*/

#include "Teuchos_ConfigDefs.hpp"
#include "ROL_XMLObject.hpp"
#include <list>
#include <stack>

namespace ROL
{

  /** \ingroup XML
   * \brief TreeBuildingXMLHandler assembles a XMLObject from your XML input.
   */
  class TreeBuildingXMLHandler
    {
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
                        const std::map<std::string, std::string>& attributes);

    private:
      XMLObject root_;
      XMLObject current_;
      std::stack<XMLObject> path_;
    };
}

#endif

