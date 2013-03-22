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
#include "Teuchos_Assert.hpp"
#include <stack>

using namespace Teuchos;

// this parser currently does not support:
// * processing instructions
// * XML schemas
// * CDATA sections...see http://www.w3.org/TR/2004/REC-xml-20040204/#dt-cdsection
// * full Unicode support (we read unsigned bytes, so we get only 0x00 through 0xFF)
// 
// it tolerates (read: ignores) xml declarations, at any point in the file where a tag would be valid
//
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
     
     document     ::=   prolog element Misc*
     prolog       ::=   XMLDecl? Misc*
     XMLDecl      ::=   '<?xml' VersionInfo EncodingDecl? SDDecl? S? '?>'
     Misc         ::=   Comment | S

     VersionInfo  ::=     S 'version' Eq ("'" VersionNum "'" | '"' VersionNum '"')
     Eq           ::=     S? '=' S?
     VersionNum   ::=    '1.' [0-9]+
     Misc         ::=     Comment | S


        
*/

#define XMLPARSER_TFE( T , S ) \
  TEUCHOS_TEST_FOR_EXCEPTION( T, std::runtime_error, "XML parse error at line " << _lineNo << ": " << S )

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
  std::stack<long> tagLineStarts;
  std::stack<string> tags;

  while (!done) {
    
    std::string tag, cdata;
    unsigned char c1, c2;
    Teuchos::map<std::string,string> attrs;
    
    // Consume any whitespace
    if (curopen == 0) {
      // this will leave a lookahead in c1
      c1 = '\0';
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
      if (c1 == '\n') ++_lineNo; // a newline while processing character data; not an error
    }

    if (c1 == '<') {
      // determine if it is a STag/EmptyElemTag or ETag or Comment
      // get lookahead
      XMLPARSER_TFE( _is->readBytes(&c2,1) < 1 , "stream ended in tag begin/end");

      if (c2 == '/') {
        // we have: </
        // try to get an ETag
        getETag(tag);
        // have to check whether we have an enclosing, otherwise tags and tagLineStarts have no top()
        XMLPARSER_TFE( curopen == 0,  "document not well-formed: encountered end element '" << tag << "' while not enclosed." );
        XMLPARSER_TFE( handler->endElement(tag)!=0, "document not well-formed: end element tag = '" << tag << "'"
                                                    << " did not match start element '" << tags.top() 
                                                    << "' from line " << tagLineStarts.top() );
        curopen--;
        tagLineStarts.pop();
        tags.pop();
      }
      else if (isLetter(c2) || c2==':' || c2=='_') {
        // it looks like a STag or an EmptyElemTag
        bool emptytag;
        tagLineStarts.push(_lineNo);
        getSTag(c2, tag, attrs, emptytag);
        tags.push(tag);
        handler->startElement(tag,attrs);
        if (curopen == 0) {
          XMLPARSER_TFE(gotRoot == true, "document not well-formed: more than one root element specified" );
          gotRoot = true;
        }
        curopen++;
        if (emptytag) {
          // we just open this tag, so we should have any trouble closing it
          XMLPARSER_TFE( handler->endElement(tag)!=0, "unknown failure from handler while processing tag '" << tag << "'" );
          curopen--;
          tagLineStarts.pop();
          tags.pop();
        }
      }
      else if (c2 == '?') {
        // it is starting to look like an xml declaration
        XMLPARSER_TFE( assertChar('x') != 0 , "was expecting an XML declaration; element not well-formed or exploits unsupported feature" );
        XMLPARSER_TFE( assertChar('m') != 0 , "was expecting an XML declaration; element not well-formed or exploits unsupported feature" );
        XMLPARSER_TFE( assertChar('l') != 0 , "was expecting an XML declaration; element not well-formed or exploits unsupported feature" );
        ignoreXMLDeclaration();
      }
      else if (c2 == '!') {
        // it is starting to look like a comment; we need '--'
        // if we don't get this, it means
        // * the document is not well-formed
        // * the document employs a feature not supported by this parser, 
        //   e.g. <!ELEMENT...  <!ATTLIST...  <!DOCTYPE...  <![CDATA[...
        XMLPARSER_TFE( assertChar('-') != 0 , "element not well-formed or exploits unsupported feature" );
        XMLPARSER_TFE( assertChar('-') != 0 , "element not well-formed or exploits unsupported feature" );
        getComment(_lineNo);
      }
      else {
        XMLPARSER_TFE(true,  "element not well-formed or exploits unsupported feature" );
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
      XMLPARSER_TFE(1 , "document not well-formed: character data outside of an enclosing tag");
    }
  }

  XMLPARSER_TFE( curopen != 0 ,  "file ended before closing element '" << tags.top() << "' from line " << tagLineStarts.top() );

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
  XMLPARSER_TFE( _is->readBytes(&c,1) < 1 ,  "EOF before end element was terminated");
  XMLPARSER_TFE( !isLetter(c) && c!='_' && c!=':' ,  "tag not well-formed");
  tag.push_back(c);
  while (1) {
    XMLPARSER_TFE( _is->readBytes(&c,1) < 1 ,  "EOF before end element was terminated");
    if ( isNameChar(c) ) {
      if (tagover) {
        XMLPARSER_TFE(1,  "end element not well-formed: expected '>'");
      }
      tag.push_back(c);
    }
    else if (isSpace(c)) {
      // mark the end of the tag and consume the whitespace
      // if it is ia newline, it isn't an error
      if (c == '\n') ++_lineNo;
      tagover = true;
    }
    else if (c == '>') {
      break; 
    }
    else {
      XMLPARSER_TFE(1,  "end element not well-formed");
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
    XMLPARSER_TFE( _is->readBytes(&c,1) < 1 ,  "EOF before start element was terminated");
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
      XMLPARSER_TFE( getSpace(c)!=0, "EOF before start element was terminated");
    }
    
    // now, either Attribute | '>' | '/>'
    if ( (isLetter(c) || c=='_' || c==':') && hadspace ) {
      
      // Attribute
      // get attribute name, starting with contents of c
      std::string attname, attval;
      attname = c;
      do {
        XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before start element was terminated");
        if ( isNameChar(c) ) {
          attname.push_back(c);
        }
        else if ( isSpace(c) || c=='=' ) {
          break; 
        }
        else {
          XMLPARSER_TFE(1,  "attribute not well-formed: expected whitespace or '='");
        }
      } while (1);
      
      // if whitespace, consume it
      if (isSpace(c)) {
        getSpace(c);  
      }
      // should be on '='
      if (c != '=') {
        XMLPARSER_TFE(1,  "attribute not well-formed: expected '='");
      }
      
      // get any whitespace following the '='
      XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before start element was terminated");
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
        XMLPARSER_TFE(1,  "attribute value must be quoted with either ''' or '\"'");
      }
      do {
        XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before start element was terminated");
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
          XMLPARSER_TFE(1,  "invalid character in attribute value");
        }
      } while(1);
      
      // add attribute to list
      XMLPARSER_TFE( attrs.find(attname) != attrs.end() ,  "cannot have two attributes with the same name");
      attrs[attname] = attval;
    }
    else if (c == '>') {
      emptytag = false;
      break;
    }
    else if (c == '/') {
      XMLPARSER_TFE(assertChar('>')!=0, "empty element tag not well-formed: expected '>'");
      emptytag = true;
      break;
    }
    else {
      XMLPARSER_TFE(1,  "start element not well-formed: invalid character");
    }
  
    // get next char
    XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before start element was terminated");
  
  } while(1);
}


void XMLParser::getComment(long startLine) 
{
  /* Recall from the specification:
        Comment   ::= '<!--' ((Char - '-') | ('-' (Char - '-')))* '-->'
                      that is, '<!--' txt '-->', where txt does not contain '--' 
     We have already consumed: <!--
     
     Be wary here of the fact that c=='-' implies isChar(c)
  */
  unsigned char c;
  while (1) {
    XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before terminating comment begun at line " << _lineNo );
    if (c == '\n') ++_lineNo;
    // if we have a -
    if (c=='-') {
      // then it must be the end of the comment or be a Char
      XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before terminating comment begun at line " << _lineNo );
      if (c == '\n') ++_lineNo;
      if (c=='-') {
        // this had better be leading to the end of the comment
        XMLPARSER_TFE( assertChar('>')!=0, "comment not well-formed: missing expected '>' at line " << _lineNo );
        break;
      }
      else if (!isChar(c)) {
        XMLPARSER_TFE(1,  "comment not well-formed: invalid character at line " << _lineNo );
      }
    }
    else if (!isChar(c)) {
      XMLPARSER_TFE(1,  "comment not well-formed: invalid character at line " << _lineNo );
    }
  } 
}


void XMLParser::getReference(std::string &refstr) {
  // finish: does CharRef support only dec, or hex as well?
  unsigned char c;
  unsigned int num, base;
  refstr = "";
  // none of these bytes read are allowed to be a newline, so don't do any incrementing of _lineNo
  XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before reference was terminated");
  if (c == '#') {
    // get a CharRef
    // CharRef   ::= '&#' [0-9]+ ';'
    //               | '&#x' [0-9]+ ';'
    // get first number
    XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before reference was terminated");
    if (c == 'x') {
      base = 16;
      num = 0;
    }
    else if ('0' <= c && c <= '9') {
      base = 10;
      num = c - '0';
    }
    else {
      XMLPARSER_TFE(1,  "invalid character in character reference: expected 'x' or [0-9]");
    }

    do {
      XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before reference was terminated");
      XMLPARSER_TFE( c != ';' && !('0' <= c && c <= '9') ,  "invalid character in character reference: expected [0-9] or ';'");
      if (c == ';') {
        break;
      }
      num = num*base + (c-'0');
    } while (1);
    XMLPARSER_TFE(num > 0xFF,  "character reference value out of range");
    refstr.push_back( (unsigned char)num );
  }
  else if (isLetter(c) || c=='_' || c==':') {
    // get an EntityRef
    // EntityRef ::= '&' Name ';'
    std::string entname = "";
    entname.push_back(c);
    do {
      XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before reference was terminated");
      if (c==';') {
        break;
      }
      else if ( isLetter(c) || ('0' <= c && c <= '9')
                || c=='.' || c=='-' || c=='_' || c==':' 
                || c==0xB7 ) {
        entname.push_back(c);
      }
      else {
        XMLPARSER_TFE(1,  "entity reference not well-formed: invalid character");
      }
    } while (1);
    XMLPARSER_TFE( _entities.find(entname) == _entities.end(),  "entity reference not well-formed: undefined entity");
    refstr = _entities[entname];  
  }
  else {
    XMLPARSER_TFE(1,  "reference not well-formed: expected name or '#'");
  }
}


int XMLParser::getSpace(unsigned char &lookahead) {
  // if space, consume the whitespace
  do {
    if (lookahead == '\n') ++_lineNo;
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
  // don't worry about newlines; assertChar is always wrapped in TEST_FOR_EXCEPTION, so we don't want to advance the line counter
  if (_is->readBytes(&c,1) < 1) {
    return 1;
  }
  if (c != cexp) {
    return 2;
  }
  return 0; 
}

void XMLParser::ignoreXMLDeclaration() 
{
  /* Be a little lax on the spec here; read until we get to '?', then assert '>'
     We have already consumed: <xml
  */
  unsigned char c;
  while (1) {
    XMLPARSER_TFE(_is->readBytes(&c,1) < 1,  "EOF before terminating XML declaration begun at line " << _lineNo );
    if (c == '\n') ++_lineNo;
    // if we have a -
    if (c=='?') {
      // this had better be leading to the end of the declaration
      XMLPARSER_TFE( assertChar('>')!=0, "XML declaration not well-formed: missing expected '>' at line " << _lineNo );
      break;
    }
  } 
}
