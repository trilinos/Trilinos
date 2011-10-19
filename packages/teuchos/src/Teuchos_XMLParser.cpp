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

// BUGS: There is a bug in Teuchos_XMLObjectImplem.cpp, line 82
// when printing attribute values, one must check if the value contains quote
// or apost; 
// a quot'd attval cannot contain literal quot
// a apos'd attval cannot contain literal apos
// either they have to be matched appropriately or (easier) all quot and apos must
// be replaced by &quot; and &apos;

#include "Teuchos_XMLParser.hpp"
#include "Teuchos_TreeBuildingXMLHandler.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;

// this parser currently does not support:
// * XML declaration
// * processing instructions
// * XML schemas
// * CDATA sections...see http://www.w3.org/TR/2004/REC-xml-20040204/#dt-cdsection
// * full Unicode support (we read unsigned bytes, so we get only 0x00 through 0xFF)

// it currently does support:
// * comments
// * empty element tags, e.g.   <hello />
// * entity references: &amp; &lt; &gt; &apos; &quot;
// * numeric character references: &#32;
// * std::exception/error handling on parse errors


/* From the W3C XML 1.0 Third Edition
   http://www.w3.org/TR/2004/REC-xml-20040204/
  
   The following productions specify well-formed XML documents.
   These have been reduced to the support anticipated for support by this parser.
        
     element      ::=  EmptyElemTag
                       | STag content ETag 
     STag         ::=  '<' Name (S Attribute)* S? '>' 
     Attribute    ::=  Name Eq AttValue 
     ETag         ::=  '</' Name S? '>'
     content      ::=  CharData? ((element | Reference | CDSect | Comment) CharData?)*
     EmptyElemTag ::=  '<' Name (S Attribute)* S? '/>'
     
     AttValue     ::=  '"' ([^<&"] | Reference)* '"'
                       | "'" ([^<&'] | Reference)* "'"
     
     CharRef   ::= '&#' [0-9]+ ';'
     EntityRef ::= '&' Name ';'
     Reference ::= EntityRef | CharRef
     
     #x20 (space)
     #x9  (horizontal tab)
     #xD  (carriage return)
     #xA  (new line, new line line feed)
     
     S        ::=  (#x20 | #x9 | #xD | #xA)+
     Eq       ::=   S? '=' S?
     NameChar ::=  Letter | Digit | '.' | '-' | '_' | ':' | #x00B7
     Name     ::=  (Letter | '_' | ':') (NameChar)*
     
     Letter   ::= [#x0041-#x005A] | [#x0061-#x007A] 
                  | [#x00C0-#x00D6] | [#x00D8-#x00F6] 
                  | [#x00F8-#x00FF]
     Digit    ::= [#x0030-#x0039]
     
     Char      ::=  #x9 | #xA | #xD | [#x20-#xFF]   
     CharData  ::= [^<&]* - ([^<&]* ']]>' [^<&]*)
                   that is, some std::string of characters not containing '<' or '&' or ']]>'
     Comment   ::= '<!--' ((Char - '-') | ('-' (Char - '-')))* '-->'
                   that is, '<!--' txt '-->', where txt does not contain '--' 
     
     CDSect    ::= CDStart CData CDEnd
     CDStart   ::= '<![CDATA['
     CData     ::= (Char* - (Char* ']]>' Char*))
     CDEnd     ::= ']]>'
     
     document  ::=   prolog element Misc*
     prolog    ::=   XMLDecl? Misc*
     XMLDecl   ::=   '<?xml' VersionInfo EncodingDecl? SDDecl? S? '?>'
     Misc      ::=   Comment | S
        
*/

XMLObject XMLParser::parse() 
{
  
  RCP<TreeBuildingXMLHandler> handler = rcp(new TreeBuildingXMLHandler());
  
  _entities.clear();
  _entities["apos"] = "'";
  _entities["quot"] = "\"";
  _entities["lt"]   = "<";
  _entities["gt"]   = ">";
  _entities["amp"]  = "&";
  
  bool done = false;
  int curopen = 0;  // number of currently open tags, or "do we process character data?"
  bool gotRoot = false;

  while (!done) {
    
    std::string tag, cdata;
    unsigned char c1, c2;
    Teuchos::map<std::string,string> attrs;
    
    // Consume any whitespace
    if (curopen == 0) {
      // this will leave a lookahead in c1
      if ( getSpace(c1) ) {
        done = true;
        break;
      }
    }
    else {
      // need to manually lookahead
      if (_is->readBytes(&c1,1) < 1) {
        done = true;
        break;
      }
    }

    if (c1 == '<') {
      // determine if it is a STag/EmptyElemTag or ETag or Comment
      // get lookahead
      TEST_FOR_EXCEPTION( _is->readBytes(&c2,1) < 1 , std::runtime_error, "XMLParser::parse(): stream ended in tag begin/end");

      if (c2 == '/') {
        // we have: </
        // try to get an ETag
        getETag(tag);
        TEST_FOR_EXCEPTION( handler->endElement(tag)!=0, std::runtime_error,
          "XMLParser::getETag(): document not well-formed: end element"
          " tag = '"<<tag<<"' did not match start element");
        curopen--;
      }
      else if (isLetter(c2) || c2==':' || c2=='_') {
        // it looks like a STag or an EmptyElemTag
        bool emptytag;
        getSTag(c2, tag, attrs, emptytag);
        handler->startElement(tag,attrs);
        if (curopen == 0) {
          TEST_FOR_EXCEPTION(gotRoot == true, std::runtime_error,
            "XMLParser::getETag(): document not well-formed: more than one root element specified");
          gotRoot = true;
        }
        curopen++;
        if (emptytag) {
          TEST_FOR_EXCEPTION( handler->endElement(tag)!=0, std::runtime_error,
            "XMLParser::getETag(): document not well-formed: end element tag did not match start element");
          curopen--;
        }
      }
      else if (c2 == '!') {
        // it is starting to look like a comment; we need '--'
        // if we don't get this, it means
        // * the document is not well-formed
        // * the document employs a feature not supported by this parser, 
        //   e.g. <!ELEMENT...  <!ATTLIST...  <!DOCTYPE...  <![CDATA[...
        TEST_FOR_EXCEPTION( assertChar('-')!=0, std::runtime_error,
            "XMLParser::parse(): element not well-formed or exploits unsupported feature" );
        TEST_FOR_EXCEPTION( assertChar('-')!=0 , std::runtime_error,
            "XMLParser::parse(): element not well-formed or exploits unsupported feature" );
        getComment();
      }
      else {
        TEST_FOR_EXCEPTION(1, std::runtime_error, "XMLParser::parse(): element not well-formed or exploits unsupported feature" );
      }
    }
    else if ( (curopen > 0) && (c1 == '&') ) {
      std::string chars = "";
      getReference(chars);
      handler->characters(chars);
    }
    else if ( (curopen > 0) ) {
      std::string chars = "";
      chars.push_back(c1);
      handler->characters(chars);
    }
    else {
      TEST_FOR_EXCEPTION(1,std::runtime_error,"XMLParser::parse(): document not well-formed");
    }
  }

  TEST_FOR_EXCEPTION( curopen != 0 , std::runtime_error, "XMLParser::parse(): document not well-formed: elements not matched" );

  return handler->getObject();

}


void XMLParser::getETag(std::string &tag)
{
  /* Recall from the specification:
        ETag  ::=  '</' Name S? '>'
        Name  ::=  (Letter | '_' | ':') (NameChar)*
    
     We have already consumed: </
  */
  
  bool tagover = false;
  unsigned char c;
  // clear tag
  tag = "";
  TEST_FOR_EXCEPTION( _is->readBytes(&c,1) < 1 , std::runtime_error , "XMLParser::getETag(): EOF before end element was terminated");
  TEST_FOR_EXCEPTION( !isLetter(c) && c!='_' && c!=':' , std::runtime_error , "XMLParser::getETag(): tag not well-formed");
  tag.push_back(c);
  while (1) {
    TEST_FOR_EXCEPTION( _is->readBytes(&c,1) < 1 , std::runtime_error , "XMLParser::getETag(): EOF before end element was terminated");
    if ( isNameChar(c) ) {
      if (tagover) {
        TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getETag(): end element not well-formed: expected '>'");
      }
      tag.push_back(c);
    }
    else if (isSpace(c)) {
      // mark the end of the tag and consume the whitespace
      tagover = true;
    }
    else if (c == '>') {
      break; 
    }
    else {
      TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getETag(): end element not well-formed");
    }
  }
}


void XMLParser::getSTag(unsigned char lookahead, std::string &tag, Teuchos::map<std::string,string> &attrs, bool &emptytag) 
{
  
  /* Recall from the specification:
        
        STag         ::=  '<' Name (S Attribute)* S? '>' 
        EmptyElemTag ::=  '<' Name (S Attribute)* S? '/>'
        Name         ::=  (Letter | '_' | ':') (NameChar)*
        NameChar     ::=  Letter | Digit | '.' | '-' | '_' | ':' | #x00B7
        
        S            ::=  (#x20 | #x9 | #xD | #xA)+
        Attribute    ::=  Name Eq AttValue 
        Eq           ::=   S? '=' S?
        AttValue     ::=  '"' ([^<&"] | Reference)* '"'
                          | "'" ([^<&'] | Reference)* "'"
        Reference ::= EntityRef | CharRef
        CharRef   ::= '&#' [0-9]+ ';'
        EntityRef ::= '&' Name ';'
        
     We have already consumed: <lookahead
  */
  
  unsigned char c;
  attrs.clear();
  
  tag = lookahead;
  // get the rest of the tag: (NameChar)*
  while (1) {
    TEST_FOR_EXCEPTION( _is->readBytes(&c,1) < 1 , std::runtime_error , "XMLParser::getSTag(): EOF before start element was terminated");
    if (isNameChar(c)) {
      tag.push_back(c);
    }
    else {
      break; 
    }
  }
  
  // after the name: should be one of the following
  // (S Attribute) | S? '>' | S? '/>'
  do {
    
    bool hadspace = false;
    
    // if space, consume the whitespace
    if ( isSpace(c) ) {
      hadspace = true;
      TEST_FOR_EXCEPTION( getSpace(c)!=0, std::runtime_error,
        "XMLParser::getSTag(): EOF before start element was terminated");
    }
    
    // now, either Attribute | '>' | '/>'
    if ( (isLetter(c) || c=='_' || c==':') && hadspace ) {
      
      // Attribute
      // get attribute name, starting with contents of c
      std::string attname, attval;
      attname = c;
      do {
        TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getSTag(): EOF before start element was terminated");
        if ( isNameChar(c) ) {
          attname.push_back(c);
        }
        else if ( isSpace(c) || c=='=' ) {
          break; 
        }
        else {
          TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getSTag(): attribute not well-formed: expected whitespace or '='");
        }
      } while (1);
      
      // if whitespace, consume it
      if (isSpace(c)) {
        getSpace(c);  
      }
      // should be on '='
      if (c != '=') {
        TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getSTag(): attribute not well-formed: expected '='");
      }
      
      // get any whitespace following the '='
      TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getSTag(): EOF before start element was terminated");
      if (isSpace(c)) {
        getSpace(c);
      }
      
      // now get the quoted attribute value
      bool apost;
      attval = "";
      if (c == '\'') {
        apost = true;
      }
      else if (c == '\"') {
        apost = false;
      }
      else {
        TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getSTag(): attribute value must be quoted with either ''' or '\"'");
      }
      do {
        TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getSTag(): EOF before start element was terminated");
        if (apost && c=='\'') {
          // end of attval
          break;
        }
        else if (!apost && c=='\"') {
          // end of attval
          break;
        }
        else if ( c == '&' ) {
          // finish: need to add support for Reference
          std::string refstr;
          getReference(refstr);
          attval += refstr;
        }
        else if ( c!='<' ) {
          // valid character for attval
          attval.push_back(c);
        }
        else {
          TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getSTag(): invalid character in attribute value");
        }
      } while(1);
      
      // add attribute to list
      TEST_FOR_EXCEPTION( attrs.find(attname) != attrs.end() , std::runtime_error , "XMLParser::getSTag(): cannot have two attributes with the same name");
      attrs[attname] = attval;
    }
    else if (c == '>') {
      emptytag = false;
      break;
    }
    else if (c == '/') {
      TEST_FOR_EXCEPTION(assertChar('>')!=0, std::runtime_error,
        "XMLParser::getSTag(): empty element tag not well-formed: expected '>'");
      emptytag = true;
      break;
    }
    else {
      TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getSTag(): start element not well-formed: invalid character");
    }
  
    // get next char
    TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getSTag(): EOF before start element was terminated");
  
  } while(1);
}


void XMLParser::getComment() 
{
  /* Recall from the specification:
        Comment   ::= '<!--' ((Char - '-') | ('-' (Char - '-')))* '-->'
                      that is, '<!--' txt '-->', where txt does not contain '--' 
     We have already consumed: <!--
     
     Be wary here of the fact that c=='-' implies isChar(c)
  */
  unsigned char c;
  while (1) {
    TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getComment(): EOF before comment was terminated");
    // if we have a -
    if (c=='-') {
      // then it must be the end of the comment or be a Char
      TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getComment(): EOF before comment was terminated");
      if (c=='-') {
        // this had better be leading to the end of the comment
        TEST_FOR_EXCEPTION( assertChar('>')!=0, std::runtime_error,
            "XMLParser::getComment(): comment not well-formed: expected '>'");
        break;
      }
      else if (!isChar(c)) {
        TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getComment(): comment not well-formed: invalid character");
      }
    }
    else if (!isChar(c)) {
      TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getComment(): comment not well-formed: invalid character");
    }
  } 
}


void XMLParser::getReference(std::string &refstr) {
  // finish: does CharRef support only dec, or hex as well?
  unsigned char c;
  unsigned int num, base;
  refstr = "";
  TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getReference(): EOF before reference was terminated");
  if (c == '#') {
    // get a CharRef
    // CharRef   ::= '&#' [0-9]+ ';'
    //               | '&#x' [0-9]+ ';'
    // get first number
    TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getReference(): EOF before reference was terminated");
    if (c == 'x') {
      base = 16;
      num = 0;
    }
    else if ('0' <= c && c <= '9') {
      base = 10;
      num = c - '0';
    }
    else {
      TEST_FOR_EXCEPTION(1, std::runtime_error, "XMLParser::getReference(): invalid character in character reference: expected 'x' or [0-9]");
    }

    do {
      TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getReference(): EOF before reference was terminated");
      TEST_FOR_EXCEPTION( c != ';' && !('0' <= c && c <= '9') , std::runtime_error , "XMLParser::getReference(): invalid character in character reference: expected [0-9] or ';'");
      if (c == ';') {
        break;
      }
      num = num*base + (c-'0');
    } while (1);
    TEST_FOR_EXCEPTION(num > 0xFF, std::runtime_error , "XMLParser::getReference(): character reference value out of range");
    refstr.push_back( (unsigned char)num );
  }
  else if (isLetter(c) || c=='_' || c==':') {
    // get an EntityRef
    // EntityRef ::= '&' Name ';'
    std::string entname = "";
    entname.push_back(c);
    do {
      TEST_FOR_EXCEPTION(_is->readBytes(&c,1) < 1, std::runtime_error , "XMLParser::getReference(): EOF before reference was terminated");
      if (c==';') {
        break;
      }
      else if ( isLetter(c) || ('0' <= c && c <= '9')
                || c=='.' || c=='-' || c=='_' || c==':' 
                || c==0xB7 ) {
        entname.push_back(c);
      }
      else {
        TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getReference(): entity reference not well-formed: invalid character");
      }
    } while (1);
    TEST_FOR_EXCEPTION( _entities.find(entname) == _entities.end(), std::runtime_error , "XMLParser::getReference(): entity reference not well-formed: undefined entity");
    refstr = _entities[entname];  
  }
  else {
    TEST_FOR_EXCEPTION(1, std::runtime_error , "XMLParser::getReference(): reference not well-formed: expected name or '#'");
  }
}


int XMLParser::getSpace(unsigned char &lookahead) {
  // if space, consume the whitespace
  do {
    if (_is->readBytes(&lookahead,1) < 1) {
      return 1; // inform caller that we reached the end
    }
  }
  while (isSpace(lookahead));
  return 0;
}


bool XMLParser::isLetter(unsigned char c) {
  if ( (0x41 <= c && c <= 0x5A) || (0x61 <= c && c <= 0x7A) ||
       (0xC0 <= c && c <= 0xD6) || (0xD8 <= c && c <= 0xF6) ||
       (0xF8 <= c) /* unsigned char must be <= 0xFF */         )
  {
    return true;
  }
  return false;
}


bool XMLParser::isNameChar(unsigned char c) {
  if ( isLetter(c) || ('0' <= c && c <= '9') ||
       c=='.' || c=='-' || c=='_' || c==':' || c==0xB7 ) 
  {
    return true;
  }
  return false;
}


bool XMLParser::isSpace(unsigned char c) {
  if ( c==0x20 || c==0x9 || c==0xD || c==0xA )
  {
    return true;
  }
  return false;
}


bool XMLParser::isChar(unsigned char c) {
  if ( c==0x9 || c==0xA || c==0xD || 0x20 <= c) {  // unsigned char must be <= 0xFF
    return true;
  }
  return false;
}


int XMLParser::assertChar(unsigned char cexp) 
{
  // pull the next character off the stream and verify that it is what is expected
  // if not, return an error to the caller
  unsigned char c;
  if (_is->readBytes(&c,1) < 1) {
    return 1;
  }
  if (c != cexp) {
    return 2;
  }
  return 0; 
}

