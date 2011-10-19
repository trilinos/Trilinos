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
          TEST_FOR_EXCEPTION(
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

