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

#include <fstream>
#include <string>

#include <ctype.h>
#include <stdlib.h>
#include <unistd.h>
#include <search.h>
#include <string.h>

using namespace std;

namespace PAMGEN_NEVADA{

static void makeUpperCase( string& s )
{
  for (unsigned int i = 0; i < s.size(); i = i + 1)
    s[i] = toupper(s[i]);
}

/*****************************************************************************/
char Token_Stream::get_next_char()
/*****************************************************************************/
{
  see_next_char();
  return input_buffer[position++];
}

/*****************************************************************************/
char Token_Stream::see_next_char()
/*****************************************************************************/
{
  if (position>=input_buffer.size()){
    if (input_buffer.size()==input_buffer.capacity())
      // Extend the input buffer
      input_buffer.reserve(input_buffer.capacity()+80);
    // Read another character
    if (newline){
      output.width(5);
      output << ++line_number << ": ";
      output.flush();
      newline = false;
    }
    char c = input->get();
    if (c=='\r') c = input->get(); // allow reading DOS formatted input files
    if (c=='\n') newline = true;
    if (all_input) (*all_input).append(1,c);
    
    if (input->eof() || input->fail()){ 
      c = 0;
    }
    else{
      output << c;
    }
    
    input_buffer += c;
  }
  return input_buffer[position];
}

/*****************************************************************************/
Token Token_Stream::read_token()
/*****************************************************************************/
  // Routine read_token checks for comments, ints, floats, and identifiers, 
  // Both routines check for whitespace and eof
{
  
  Token_Value yylval;

  // Begin scanning the next token.

  // Discard initial whitespace
  char c = see_next_char();
  while (iswhite(c))
  {
    get_next_char();
    c = see_next_char();
  }
  if (position){
    input_buffer = string(input_buffer, position,
                          input_buffer.size() - position);
    position = 0;
  }

  // Discard comment: *, $, # to next newline.
  if (c=='*' || c=='$' || c=='#'){
    while (c!='\n' && c!=0){
      get_next_char();
      c = see_next_char();
    }
    if (c=='\n'){
      get_next_char();
    }
    input_buffer = string(input_buffer, position,
                          input_buffer.size() - position);
    position = 0;
    return read_token();
  }

  //check for eof
  if (c==0){
    if (input->eof()){
      input_buffer = string(input_buffer, position,
                            input_buffer.size() - position);
      position = 0;
      if (!input_stack.size()){
        return Token(TK_EXIT, Token::no_value);
      } else {
//         pop_include();
        return read_token();
      }
    }else{
      Parse_Error("Token_Stream::read_token: unable to read input file");
      return Token(TK_EXIT, Token::no_value);
    }  
  }else if (c=='+' || c=='-' || c=='.' || isdigit(c)){
    if (c=='-' || c=='+'){
      assert(position==0);
      char c1 = get_next_char();
      c = see_next_char();
      if (!isdigit(c) && c!='.'){
        input_buffer = string(input_buffer, position,
                              input_buffer.size() - position);
        position = 0;
        return Token((Token_Type)c1, Token::no_value);
      }
      if (c=='.'){
        // Must be a float if a number at all
        get_next_char();
        c = see_next_char();
        if (!isdigit(c)){
          input_buffer = string(input_buffer, 1, input_buffer.size() - 1);
          position = 0;
          return Token((Token_Type)c1, Token::no_value);
        }else{
          position = 1;
          c = see_next_char();
        }
      }
    }
    while (isdigit(c)){
      get_next_char();
      c = see_next_char();
    }
    if (c=='.' || c=='e' || c=='E'){
      if (c=='.'){
        // it's a float; look for tail, if any
        get_next_char();
        c = see_next_char();
        while (isdigit(c)){
          get_next_char();
          c = see_next_char();
        }
      }
      if (c=='e' || c=='E'){
        int wpos = position;
        get_next_char();
        c = see_next_char();
        if (c=='+' || c=='-'){
          get_next_char();
          c = see_next_char();
        }
        if (!isdigit(c)){
          position = wpos;
        }else{
          while (isdigit(c)){
            get_next_char();
            c = see_next_char();
          }
        }
      }
      string text = string(input_buffer, 0, position);
      input_buffer = string(input_buffer, position,
                            input_buffer.size() - position);
      position = 0;
      yylval.fval = atof(text.c_str());
      return Token(TK_REAL, yylval);
    }else{
      // nope, it's an integer
      string text = string(input_buffer, 0, position);
      input_buffer = string(input_buffer, position,
                          input_buffer.size() - position);
      position = 0;
      yylval.ival = atoi(text.c_str());
      return Token(TK_INTEGER, yylval);
    }
  //check for indication of include file
  }else if (c=='\"' || c=='<' || c=='\''){
    char cterm = (c=='<'? '>' : c);
    get_next_char();
    c = see_next_char();
    while (c!=cterm && c!=0){
      get_next_char();
      c = see_next_char();
    }
    if (c==0){
      Parse_Error("unexpected end of file; did you forget a closing quote?");
    }
    get_next_char();
    string text = string(input_buffer, 0, position);
    input_buffer = string(input_buffer, position,
                          input_buffer.size() - position);
    position = 0;
    yylval.sval = (char*)text.c_str();
    return Token(TK_STRING, yylval);
  //check for character string (underscores allowed) 
  }else if (isalpha(c) || c=='_'){
    while (1){
      while (isalnum(c) || c=='_'){
        get_next_char();
        c = see_next_char();
      }
      int wpos = position;
      while (c==' ' || c=='\t'){
        get_next_char();
        c = see_next_char();
      }
      if (!isalpha(c) && c!='_'){
        string token_text = string(input_buffer, 0, wpos);
	makeUpperCase(token_text);
        input_buffer = string(input_buffer, wpos,
                              input_buffer.size() - wpos);
        position = 0;
        if (token_text=="END"){
          return Token(TK_END, Token::no_value);
        }else if (token_text=="EXIT"){
          if (!input_stack.size()){
            return Token(TK_EXIT, Token::no_value);
          } else {
//             pop_include();
            return read_token();
          }
//         }else if (include_active && token_text=="INCLUDE"){
//           push_include();
//           return read_token();
        }else{
          yylval.sval = (char*)token_text.c_str();
          return Token(TK_IDENTIFIER, yylval);
        }
      }
    }
  }else if (c=='/'){
    get_next_char();
    if (see_next_char()=='*'){
      // a C-style comment
      c = get_next_char();
      while (c!=0 && see_next_char()!=0 && (c!='*' || see_next_char()!='/'))
        c = get_next_char();
      if (c=='*' && see_next_char()=='/'){
        get_next_char();
        return read_token();
      }else{
	Parse_Error("Token_Stream::read_token: End of file encountered within C-style comment");
        return Token(TK_EXIT, Token::no_value);
      }
    }else{
      return Token((Token_Type)'/', Token::no_value);
    }
  }else{
    input_buffer = string(input_buffer, position + 1,
                          input_buffer.size() - position - 1);
    position = 0;
    return Token((Token_Type)c, Token::no_value);
  }
}

/*****************************************************************************/
string Token_Stream::Get_Line()
/*****************************************************************************/
  // Read the rest of the current line from the input stream.
{
  string line;
  line.reserve(80);

  char c = get_next_char();
  line += c;
  while (c!='\n' && c!=0) {
    c = get_next_char();
    line += c;
  }
 
  /* skip leading characters which are whitespace */
  const char *cptr=line.c_str();
  while(iswhite(*cptr) && *cptr != '\n') {
    cptr++;
  }   

  blockstack.top()->add(cptr);

  input_buffer = string(input_buffer, position,
                        input_buffer.size() - position);
  position = 0;
  return string(cptr);
}
}//end namespace PAMGEN_NEVADA
