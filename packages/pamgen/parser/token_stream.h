// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef token_streamH
#define token_streamH

#include "token.h"

#include <sstream>
#include <stack>
#include <string>
#include <stdio.h>
#include <setjmp.h>


namespace PAMGEN_NEVADA{
class InputBlock;
/*****************************************************************************/
class Token_Stream
/*****************************************************************************/
// Represents a stream of tokens.  These are extracted from the associated
// istream and echoed to the associated ostream.
{
  public:
  Token_Stream(std::istream &in, std::ostream &out,
                 InputBlock* inputstructure,
                 int depth, int verb, bool include_active);

  private:
    // Not defined.  A Token_Stream is not copyable or assignable.
    Token_Stream(const Token_Stream&);
    Token_Stream& operator=(const Token_Stream&);

  public:
    ~Token_Stream();

// Access functions

    int Error_Count() const {return error_count;}
    int Indentation_Depth() const {return indentation_depth;}
    int Line_Number() const {return line_number;}
    int Recovery_Flag() const {return recovery_flag;}
    int Verbosity() const {return verbosity;}

    std::ostream &Output_Stream() const {return output;}

// Definitions

    void Set_Recovery_Flag(bool flag){recovery_flag=flag;}
    // Normally used by a client to tell the Token_Stream that recovery
    // from the last error is complete.  May also be used to manually suppress
    // error messages by turning the recover flag on, but normally the recovery
    // flag is turned on by the Parse_Error function.

    void Set_Verbosity(int v){verbosity=v;}

    const std::string& RunID() const { return runid; }
    const std::string& CurDir() const { return curdir; }

// General methods

    // Examine the next token in the stream without extracting it.
    Token Lookahead();

    // Extract the next token in the stream.
    Token Shift();

    // Extract the next token and treat as a sub-block keyword.
    Token BlockShift();

    // Save a context to which the Parse_Error function should execute
    // a longjmp upon detecting an error.  Note that the client must
    // call setjmp himself and then pass the resulting jmp_buf into
    // the Token_Stream.  It would be nice to wrap this call, but
    // unfortunately the setjmp/longjmp mechanism works only if the
    // function calling setjmp has not yet returned when longjmp is
    // called.  (See man -s 3c setjmp for more info.)
    void Set_Recovery_Context(jmp_buf);

    // Print a nicely formatted message to the output stream associated
    // with this Token_Stream.  The message indicates a grammatical error
    // in the input stream, such as an unrecognized identifier, and it will
    // be properly synchronized with other output to the output stream.
    //
    // Upon printing the synchronized message, the function uses a longjmp
    // to return to the last saved context.  This allows error recovery.
    //
    // The first argument should contain a very brief description of the
    // error, such as "integer expected".  The second argument, if non-empty,
    // should be a more verbose explanation of the error or a hint on what
    // the correct grammar might be.  The second argument is inserted into
    // the output stream only if verbosity is nonzero.
    void Parse_Error(const std::string &s, const std::string &v = "");

    // Similar to Parse_Error, except that parsing continues without an
    // intervening call to longjmp.
    void Semantics_Error(const std::string &s, const std::string &v = "");

    // Read the rest of the current line in the associated input stream.
    // If this consists only of whitespace, read the next line.
    std::string Get_Line();

// Specific elementary parsing methods

    bool At_Integer();
    int  Parse_Integer();

    bool At_Real();
    double Parse_Real();

    bool At_String();
    std::string Parse_String();


// Define the set of whitespace characters
    static bool iswhite(char c){
      return c==' ' || c=='\t' || c==',' || c=='=' || c==';' || c==':' ||
         c=='\n';
    }

    void setRunID( const std::string& s ) { runid = s; }
    void setCurDir( const std::string& d ) { curdir = d; }

    void setCollectInput( std::string* is ) { all_input = is; }

    void pushNewInputBlock( const char* name );
    void popInputBlock();
    InputBlock* getInputBlock();  // returns the current block

  private:
    std::istream *input;             // not owned
    std::ostream &output;
    std::string* all_input;    // collects the entire input stream
    std::stack< InputBlock* > blockstack;
    int indentation_depth;
    int verbosity;
    bool include_active;

    std::string runid;
    std::string curdir;

    bool token_needed;
    Token lookahead;

    std::string input_buffer;   // Characters read so far from input streams
    unsigned  position;         // Current position in input buffer
    bool newline;

    int line_number;
    int error_count;
    bool recovery_flag;
    jmp_buf recovery_context;

    // Include support
    std::stack<std::istream*> input_stack;
    std::stack<int> line_number_stack;

    Token read_token();
    char get_next_char();
    char see_next_char();

    // Parse a file name from the current input stream, push the current
    // input stream, and open a new input stream to the named file.  Increment
    // the indentation depth.  Input now croms from the new input stream. When
    // the new stream reaches EOF, close it, decrement the indentation depth,
    // and pop the old stream.  This sequence of actions can be nested.
    void push_include();
    void pop_include();
};

// The following functions have the signature required for a keyword function
// (see keyword.h) but tie directly into the corresponding token stream member
// function.

Token Get_Real_Token(Token_Stream *, int);
Token Get_Integer_Token(Token_Stream *, int);
Token Get_String_Token(Token_Stream *, int);
Token Get_Identifier_Token(Token_Stream *, int);

// The following is the keyword function version of a no-op.
Token Get_No_Token(Token_Stream *, int);
}//end namespace PAMGEN_NEVADA
#endif
