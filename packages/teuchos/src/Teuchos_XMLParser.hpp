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

#ifndef TEUCHOS_XMLPARSER_H
#define TEUCHOS_XMLPARSER_H

/*! \file Teuchos_XMLParser.hpp
    \brief A class providing a simple XML parser. Methods can be overloaded 
           to exploit external XML parsing libraries.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLInputStream.hpp"

namespace Teuchos
{
  /** 
   * \brief XMLParser consumes characters from an XMLInputStream object,
   * parsing the XML and using a TreeBuildingXMLHandler to construct an
   * XMLObject.
   */
  class TEUCHOS_LIB_DLL_EXPORT XMLParser
    {
    public:
     
      /** \brief Constructor */
      XMLParser(RCP<XMLInputStream> is) : _is(is) {;}
      
      /** \brief Destructor */
      ~XMLParser(){;}
      
      /** \brief Consume the XMLInputStream to build an XMLObject. */
      XMLObject parse();
    private:
//use pragmas to disable some false-positive warnings for windows sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
      RCP<XMLInputStream> _is;
      Teuchos::map<std::string,string> _entities;
#ifdef _MSC_VER
#pragma warning(pop)
#endif
      
      /** \brief Determine whether \c c matches the <tt>Letter</tt> production according to the XML specification.*/
      inline static bool isLetter(unsigned char c);
      /** \brief Determine whether \c c matches the <tt>NameChar</tt> production according to the XML specification.*/
      inline static bool isNameChar(unsigned char c);
      /** \brief Determine whether \c c matches the <tt>Char</tt> production according to the XML specification.*/
      inline static bool isChar(unsigned char c);
      /** \brief Determine whether \c c matches the <tt>Space</tt> production according to the XML specification.*/
      inline static bool isSpace(unsigned char c);

      /** \brief Consume a <tt>ETag</tt> production according to the XML specification.
       *  <tt>getETag</tt> throws an std::exception if the input does not match the production rule.
       *  
       *  @param tag
       *         [out] On output, will be set to the tag name of the closing tag.
       */
      void getETag(std::string &tag);

      /** \brief Consume a <tt>STag</tt> production according to the XML specification.
       *  <tt>getSTag</tt> throws an std::exception if the input does not match the production rule.
       *  
       *  @param lookahead
       *         [in] Contains the first character of the tag name.
       * 
       *  @param tag
       *         [out] On output, will be set to the tag name of the opening tag.
       * 
       *  @param attrs
       *         [out] On output, contains the attributes of the tag.
       *
       *  @param emptytag
       *         [out] On output, specifies if this was an empty element tag.
       *
       */
      void getSTag(unsigned char lookahead, std::string &tag, Teuchos::map<std::string,string> &attrs, bool &emptytag);

      /** \brief Consume a <tt>Comment</tt> production according to the XML specification.
       *  <tt>getComment</tt> throws an std::exception if the input does not match the production rule.
       */
      void getComment();

      /** \brief Consumes a <tt>Space</tt> (block of whitepace) production according to the XML specification.
       *
       *  @param lookahead
       *         [out] On output, specifies the first character after the whitespace.
       *
       *  @return Returns non-zero if the input stream was exhausted while reading whitespace.
       */
      int getSpace(unsigned char &lookahead);

      /** \brief Consumes a <tt>Reference</tt> production according to the XML specification.
       *
       *  @param refstr
       *         [out] On output, specifies the decoded reference.
       *
       */
      void getReference(std::string &refstr);

      /** \brief Determines if the next character on the stream 
       *
       *  @param cexp
       *         [in] The expected character.
       *
       *  @return Returns non-zero if the next character on the stream is not \c cexp.
       * 
       */
      int assertChar(unsigned char cexp);
    };
}

#endif
