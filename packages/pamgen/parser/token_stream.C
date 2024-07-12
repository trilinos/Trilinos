// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "token_stream.h"
#include "InputBlock.h"

#include <string>
#include <iostream>
#include <sstream>
#include <cstring>

using namespace std;
namespace PAMGEN_NEVADA{
  /*****************************************************************************/
  Token_Stream::Token_Stream(istream &in, ostream &out,
      InputBlock* ib, int depth, int verb, bool ia)
    /*****************************************************************************/
    // Create a Token_Stream whose text is taken from in and echoed to out.
    // Only processor 0 does this echoing.
    //
    // depth indicates the indentation depth to use in the echo file, which is
    // useful for grammars in which an input file can refer to another input
    // file.  The included file can be more highly indented to set it off from
    // the including file.
    //
    // verb indicates the verbosity level for diagnostic messages.  A value
    // of zero, the default, causes only the basic diagnostic message to print.
    // A higher value enables printing of supplemental information, such as a
    // list of recognized keywords or recommendations for correcting the
    // diagnosed problem.
    //
    // ia flags whether the builtin include keyword should be active.  If it
    // is active, then whenever the scanner scans the identifier INCLUDE,
    // the scanner scans a filename, opens this file as a new broadcast_ifstream,
    // and takes its next token from this new stream.  When the end of this
    // stream is reached, input from the previous stream resumes.  If the
    // INCLUDE identifier is not followed by a string, it is a parse error.
    // If the stream cannot be opened, it is a semantics error.
    //
    // Upon successful construction, the Token_Stream points to the first token
    // in the input stream.  In other words, Lookahead() returns the first token
    // in the input stream, and Shift() returns this token and advances to the
    // second token in the input stream.
    //
    // in must already be opened and in a good state for reading.  out must
    // already be opened and in a good state for writing on processor zero.
    // On other processors, out is not used and need not be open or in a
    // good state.
    : output(out),
    all_input(0),
    indentation_depth(depth),
    verbosity(verb),
    include_active(ia),
    token_needed(true),
    position(0),
    newline(true),
    line_number(0),
    error_count(0),
    recovery_flag(false)
  {
    assert (in);

    input = &in;
    if(ib)blockstack.push(ib);

    input_buffer.reserve(80);        // Most likely maximum token length

    assert(Line_Number()==0);
    assert(Indentation_Depth()==depth);
    assert(Verbosity()==verb);
    assert(Error_Count()==0);
    assert(Recovery_Flag()==false);
  }

  /*****************************************************************************/
  Token Token_Stream::Shift()
    /*****************************************************************************/
    // Return the current token in the stream and advance to the next token.
  {
    Token tok;

    if (token_needed){
      tok = read_token();
    }
    else
    {
      token_needed = true;
      tok = lookahead;
    }

    InputBlock* iblk = blockstack.top();

    if (tok.Type() == TK_IDENTIFIER) {
      iblk->add( tok.As_String() );
    }
    else if (tok.Type() == TK_STRING) {
      iblk->add( tok.As_Stripped_String().c_str() );
    }
    else if (tok.Type() == TK_INTEGER) {
      ostringstream oss; oss << tok.As_Int();
      iblk->add( oss.str().c_str() );
    }
    else if (tok.Type() == TK_REAL)
    {
      ostringstream oss; oss.setf( ios_base::scientific ); oss.precision(12);
      oss << tok.As_Real();
      iblk->add( oss.str().c_str() );
    }
    else if (tok.Type() == TK_LP)       iblk->add("(");
    else if (tok.Type() == TK_RP)       iblk->add(")");
    else if (tok.Type() == TK_PER)      iblk->add("/");
    else if (tok.Type() == TK_CARET)    iblk->add("^");
    else if (tok.Type() == TK_LT)       iblk->add("<");
    else if (tok.Type() == TK_GT)       iblk->add(">");
    else if (tok.Type() == TK_PLUS)     iblk->add("+");
    else if (tok.Type() == TK_NONE) ;
    else if (tok.Type() == TK_END)  ;
    else
      cout << "magic: type " << tok.Type() << endl;

    return tok;
  }

  /*****************************************************************************/
  Token Token_Stream::BlockShift()
    /*****************************************************************************/
    // Return the current token in the stream and advance to the next token.
  {
    if (token_needed){
      return read_token();
    }else{
      token_needed = true;
      return lookahead;
    }
  }

  /*****************************************************************************/
  Token Token_Stream::Lookahead(void)
    /*****************************************************************************/
    // Return the current token on the stream, but do not advance the stream.
  {
    if (token_needed){
      lookahead = read_token();
      token_needed = false;
    }
    return lookahead;
  }



  /*****************************************************************************/
  bool Token_Stream::At_Integer()
    /*****************************************************************************/
  {
    Token token = Lookahead();
    return token.Type() == TK_INTEGER;
  }

  /*****************************************************************************/
  int Token_Stream::Parse_Integer()
    /*****************************************************************************/
    // Parse an integer quantity.  Unlike real quantities,
  {
    Token token = Shift();
    if (token.Type() != TK_INTEGER){
      Parse_Error("Integer expected instead of keyword before this message");
    }
    return token.As_Int();
  }

  /*****************************************************************************/
  bool Token_Stream::At_Real()
    /*****************************************************************************/
  {
    Token token = Lookahead();
    return token.Type() == TK_REAL;
  }

  /*****************************************************************************/
  double Token_Stream::Parse_Real()
    /*****************************************************************************/
  {
    Token token = Shift();
    if (token.Type() != TK_REAL && token.Type() != TK_INTEGER)
      Parse_Error("number expected instead of the keyword before this message");
    double rtn = token.As_Real();
    return rtn;
  }

  /**************************************************************************/
  void Token_Stream::Set_Recovery_Context(jmp_buf buf)
    /**************************************************************************/
    // This is part of the mechanism whereby we implement error recovery in
    // a parser without having to use exceptions.  (Exceptions are not yet
    // adequately supported on all platforms we must support.)
  {
    memcpy(recovery_context, buf, sizeof(recovery_context));
  }

  /*************************************************************************/
  bool Token_Stream::At_String()
    /*************************************************************************/
  {
    Token token = Lookahead();
    return token.Type()==TK_STRING;
  }

  /*************************************************************************/
  string Token_Stream::Parse_String()
    /*************************************************************************/
  {
    Token token = Shift();
    if (token.Type() != TK_STRING)
      Parse_Error("Quoted string expected instead of keyword before this message");
    return token.As_String();
  }

  /*************************************************************************/
  Token Get_Integer_Token(Token_Stream *token_stream, int)
    /*************************************************************************/
  {
    Token token = token_stream->Shift();
    if (token.Type() != TK_INTEGER)
      token_stream->Parse_Error("integer expected instead of keyword before this message");
    return token;
  }

  /*************************************************************************/
  Token Get_Real_Token(Token_Stream *token_stream, int)
    /*************************************************************************/
  {
    Token token = token_stream->Shift();
    if (token.Type() != TK_REAL && token.Type() != TK_INTEGER)
      token_stream->Parse_Error("number expected instead of keyword before this message");
    return token;
  }

  /*************************************************************************/
  Token Get_String_Token(Token_Stream *token_stream, int)
    /*************************************************************************/
  {
    Token token = token_stream->Shift();
    if (token.Type() != TK_STRING)
      token_stream->Parse_Error("quoted string expected instead of keyword before this message");
    return token;
  }

  /*************************************************************************/
  Token Get_Identifier_Token(Token_Stream *token_stream, int)
    /*************************************************************************/
  {
    Token token = token_stream->Shift();
    if (token.Type() != TK_IDENTIFIER)
      token_stream->Parse_Error("keyword expected instead of keyword before this message");
    return token;
  }

  /*************************************************************************/
  Token Get_No_Token(Token_Stream *, int)
    /*************************************************************************/
    // Used in tables of parsing functions as a nop
  {
    return Token(TK_NONE, Token::no_value);
  }

  /*************************************************************************/
  Token_Stream::~Token_Stream()
    /*************************************************************************/
    // Destroy any streams remaining in the input stack.
  {
    while (input_stack.size()){
      istream *stream = input_stack.top();
      delete stream;
      input_stack.pop();
    }
  }

  void Token_Stream::pushNewInputBlock( const char* name )
  {
    InputBlock* ib = new InputBlock;
    ib->set( name, line_number );
    blockstack.top()->add( ib );
    blockstack.push(ib);
  }

  void Token_Stream::popInputBlock()
  {
    assert( blockstack.size() > 0 );
    blockstack.pop();
  }

  InputBlock* Token_Stream::getInputBlock()
  {
    assert( blockstack.size() > 0 );
    return blockstack.top();
  }
}//end namespace
