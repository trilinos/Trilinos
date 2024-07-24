// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "keyword.h"
#include "parser.h"
#include "parse_routines.h"

#include <sstream>
#include <ctype.h>
#include <stdio.h>
#include <search.h>
#include <stdarg.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>

using namespace std;

namespace PAMGEN_NEVADA {

/*****************************************************************************/
int Parse( Token_Stream *token_stream,
           const Parse_Table & parse_table,
           Token_Type end_token )
/*****************************************************************************/
  // Given a keyword table, parse the input token stream.  Although this
  // function stands alone, it ought not to be used this way; use of the
  // Parse_Table class is preferred, since it provides a cleaner interface,
  // many more services, and a means to resolve ambiguities in the keyword
  // table.  Standalone ::Parse is legacy code.
{
  const Keyword * key_table = parse_table.begin();
  unsigned key_table_size = parse_table.size();

  assert( token_stream != 0 );
  assert( key_table != 0 );
  assert( key_table_size > 0 );
  assert( PAMGEN_NEVADA::Check_Keyword_Table(key_table, key_table_size) == 0 );

  while(1){
    jmp_buf recovery_context;
    setjmp(recovery_context);
    token_stream->Set_Recovery_Context(recovery_context);

    Token token = token_stream->BlockShift();
    if (token.Type() == end_token) break;
    if (token.Type() == TK_EXIT) {
      token_stream->Semantics_Error("End of file reached;"
                                    "did you forget a delimiter"
                                    "(such as END?)");
      break;
    }

    if (token.Type() != TK_IDENTIFIER)
      token_stream->Parse_Error("expected a keyword",
                  PAMGEN_NEVADA::Concatenate_Legal_Commands(key_table, key_table_size));

    const char *name = token.As_String();
    const Keyword *match = (Keyword*)bsearch( name,
                                              (char*)key_table,
                                              key_table_size,
                                              sizeof(Keyword),
                                              PAMGEN_Cstring_Keyword_Compare
                                            );

    if (!match){
      token_stream->
        Parse_Error(string("this keyword is unrecognized: ") + name,
                    PAMGEN_NEVADA::Concatenate_Legal_Commands(key_table,
                                               key_table_size));
      continue;
    }

    PAMGEN_NEVADA::Check_for_Ambiguity(token_stream, name, match, key_table, key_table_size);

    token_stream->Set_Recovery_Flag(false);
    assert(match->func != 0);
    token_stream->pushNewInputBlock( match->name );
    Allow_New_Mesh_Specification();
    match->func(token_stream, match->argument);
    token_stream->popInputBlock();
  }

  return token_stream->Error_Count();
}

/*****************************************************************************/
int Parse( Token_Stream *token_stream,
           const Parse_Table * parse_table,
           Token_Type end_token )
/*****************************************************************************/
{
  return Parse( token_stream, *parse_table, end_token );
}


/*****************************************************************************/
string Concatenate_Legal_Commands(const Keyword *table, int N)
/*****************************************************************************/
  // Called whenever a keyword is not recognized.  Generates a list of
  // all the keywords in the keyword table so the user can see what the
  // possibilities are.  Since this can be a very long list, the ::Parse
  // function passes it as the second message to Token_Stream::Parse_Error,
  // which ensures that it is not actually printed unless the user has
  // specified verbose mode.
{
  string list = "Select from:\n";
  for (int i=0; i<N; i++)
    list += string("  ") + string(table[i].name) + "\n";
  return list;
}

/*****************************************************************************/
void Check_for_Ambiguity(Token_Stream *token_stream,
                         const char *name,
                         const Keyword *&match,
                         const Keyword *table,
                         int N)
/*****************************************************************************/
// Check to see that name matches match unambiguously.  If more than one
// keyword matches name, but one matches exactly, replace match with the
// exact matching keyword.  Otherwise, if more than one keyword matches,
// report an error.
{
  assert(name!=NULL);
  assert(match!=NULL);
  assert(table!=NULL);
  assert(N>0);
  assert(!PAMGEN_NEVADA::Check_Keyword_Table(table, N));
  // The keyword table must be sorted if this function is to work properly.

  int error_flag = 0;
  string matches = " " + string(match->name);

  // Check if match is an exact match.  There can only be one exact match
  // unless the keyword table itself is ambiguous, which is ruled out
  // by the last precondition (above.)

  if(!strcmp(name, match->name)) {
    return;
  }

  // Check for ambiguous or exact matches behind input match.

  const Keyword *p = match-1;
  while (p >= table && !Token::Token_Match(name, p->name)){
    if(!strcmp(name, p->name)) {
      match = p;
      return;
    }
    error_flag++;
    matches = " " + string(p->name) + " " + matches;
    --p;
  }

  // Check for ambiguous or exact matches ahead of input match.

  p = match+1;
  while (p - table < N && !Token::Token_Match(name, p->name)){
    if(!strcmp(name, p->name)) {
      match = p;
      return;
    }
    error_flag++;
    matches += " " + string(p->name);
    p++;
  }

  if (error_flag)
    token_stream->Parse_Error(string("this abbreviation is ambiguous: ") + name +
                              " possible matches: " + matches);
}


/*****************************************************************************/
Keyword* Search_For_Keyword(const char* name, Keyword* table, int length)
/*****************************************************************************/
  // At present (010521), this function is used only for handcrafted parsers,
  // such as the ALEGRA adaptivity parser.  Such use is discouraged.  We
  // would like this function to go away eventually.
{
  for (int i = 0; i < length; ++i)
    if (Token::Token_Match(name, table[i].name) == 0)
      return &table[i];
  return 0;
}



string MainKeyword(const char* name)
{
  string s = string("$") + string("\n");
  s += string("$ Example Syntax for ") + string(name) + string("\n");
  s += string("$") + string("\n");
  s += string(name) + string("\n");

  return s;
}

string MainKeyword(const char* name, const char* line)
{
  string s = string("$") + string("\n");
  s += string("$ Example Syntax for ") + string(name) + string("\n");
  s += string("$") + string("\n");
  s += string(name) + string(": ") + string(line) + string("\n");

  return s;
}

string MainKeyword(string name, string line)
{
  string s = string("$") + string("\n");
  s += string("$ Example Syntax for ") + name + string("\n");
  s += string("$") + string("\n");
  s += name + string(": ") + line + string("\n");

  return s;
}

string MainKeywordId(const char* name) {
  string s = string("$") + string("\n");
  s += string("$ Example Syntax for ") + string(name) + string("\n");
  s += string("$") + string("\n");
  s += string(name);
  s += string(" integer $integer-id for main keyword") + string("\n");

  return s;
}

string MainKeywordOptId(const char* name) {
  string s = string("$") + string("\n");
  s += string("$ Example Syntax for ") + string(name) + string("\n");
  s += string("$") + string("\n");
  s += string(name);
  s += string(" [integer] $optional integer-id for main keyword") + string("\n");

  return s;
}

string SubKeyword(const char* name) {
  string s = string("  ") + string(name);
  s += string("\n");
  return s;
}

string SubKeyword(const char* name, const char* opts) {
  string s = string("  ") + string(name);
  s += string(" options = [") + string(opts);
  s += string("]") + string("\n");
  return s;
}

string SubKeyword(const char* name, const char* opts, const char* def) {
  string s = string("  ") + string(name);
  s += string(" options = [") + string(opts);
  s += string("]   $default = ") + string(def);
  s += string("\n");
  return s;
}

string SubKeyword(const char* name, double d) {
  ostringstream oss;
  oss << "  " << name << " real" << "      $default = " << d << "\n";
  return oss.str();
}

string SubKeyword(const char* name, int d) {
  ostringstream oss;
  oss << "  " << name << " integer   $default = " << d << "\n";
  return oss.str();
}

string EndKeyword() {
  string s = string("END") + string("\n");
  s += string("*****************************************") + string("\n");
  return s;
}

string NoEndKeyword() {
  string s = string("*****************************************") + string("\n");
  return s;
}

string SubSubKeyword(const char* name) {
  string s = string("  ") + SubKeyword(name);
  return s;
}

string SubSubKeyword(const char* name, const char* opts) {
  string s = string("  ") + SubKeyword(name,opts);
  return s;
}

string SubSubKeyword(const char* name, const char* opts, const char* def) {
  string s = string("  ") + SubKeyword(name,opts,def);
  return s;
}

string SubSubKeyword(const char* name, double d) {
  string s = string("  ") + SubKeyword(name,d);
  return s;
}

string SubSubKeyword(const char* name, int d) {
  string s = string("  ") + SubKeyword(name,d);
  return s;
}

string SubEndKeyword() {
  string s = string("  END") + string("\n");
  return s;
}
}  // end namespace PAMGEN_NEVADA


extern "C" {

// these two functions taken out of line to remove compiler warnings (rrd)

/*****************************************************************************/
int PAMGEN_Cstring_Keyword_Compare(const void *token, const void *key)
/*****************************************************************************/
  // This is a bsearch-compatible comparison function.  The first pointer
  // must point to a valid const char* and the second must point to a valid
  // const Keyword struct.
{
  const char *t = (const char*)token;
  const char *k = ((const PAMGEN_NEVADA::Keyword*)key)->name;
  return PAMGEN_NEVADA::Token::Token_Match(t, k);
}

/*****************************************************************************/
int PAMGEN_Keyword_Compare(const void *data1, const void *data2)
/*****************************************************************************/
  // This is a qsort-compatible comparison function for keywords.
  // The two arguments must point to valid const Keyword structs.
{
  const char *string1 = ((const PAMGEN_NEVADA::Keyword*)data1)->name;
  const char *string2 = ((const PAMGEN_NEVADA::Keyword*)data2)->name;
  return PAMGEN_NEVADA::Token::Token_Match(string1,string2);
}

}  /* end extern C */
