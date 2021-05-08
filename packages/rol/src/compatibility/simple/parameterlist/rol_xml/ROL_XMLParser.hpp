#ifndef ROL_XMLPARSER_H
#define ROL_XMLPARSER_H

/*! \file ROL_XMLParser.hpp
    \brief A class providing a simple XML parser. Methods can be overloaded
           to exploit external XML parsing libraries.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "ROL_FileInputStream.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL
{
  /**
   * \brief XMLParser consumes characters from an XMLInputStream object,
   * parsing the XML and using a TreeBuildingXMLHandler to construct an
   * XMLObject.

    Note: per the XML standard, certain characters must be represented with escape codes within fields: 

    @verbatim
         character            description             code
             <                 less-than              &lt;
             &                 ampersand              &amp;
    @endverbatim

   */
  class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT XMLParser
    {
    public:

      /** \brief Constructor */
      XMLParser(FileInputStream& is) : _is(is), _lineNo(1) {;}

      /** \brief Destructor */
      ~XMLParser(){;}

      /** \brief Consume the XMLInputStream to build an parameter list. */
      Teuchos::ParameterList parse();
    private:
//use pragmas to disable some false-positive warnings for windows sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
      FileInputStream& _is;
      std::map<std::string,std::string> _entities;
      long _lineNo;
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
      void getSTag(unsigned char lookahead, std::string &tag, std::map<std::string,std::string> &attrs, bool &emptytag);

      /** \brief Consume a <tt>Comment</tt> production according to the XML specification.
       *  <tt>getComment</tt> throws an std::exception if the input does not match the production rule.
       */
      void getComment(long startLine);

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

      /** \brief Ignore the rest of an XML declaration tag.
       */
      void ignoreXMLDeclaration();

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
