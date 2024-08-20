// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef parserH
#define parserH

#include "keyword.h"
#include "token_stream.h"
#include "parse_table.h"
#include <string>

// Some helpful routines for parsing from raw keyword tables.

// Check that a keyword table has no duplicate keywords and is properly sorted.
namespace PAMGEN_NEVADA {
bool Check_Keyword_Table(const Keyword *table, int table_length);

// Create a nicely formatted std::string summarizing the keywords in a table.
std::string Concatenate_Legal_Commands(const Keyword *table, int table_length);

// See if a C string matches a particular keyword in a table unambiguously.
// If the original match is exact, or if there is no other match, return
// normally.  Otherwise, if there is an exact match, substitute this for
// the original match and return normally.  Otherwise, call 
// token_stream->Parse_Error.

void Check_for_Ambiguity(Token_Stream *token_stream, 
                         const char *cstring,
                         const Keyword *&original_match,
                         const Keyword *keyword_table,
                         int keyword_table_size); 

//! Parse a Token_Stream using a keyword table.  The end_token is the token
//! to use as a terminator--usually TK_END. The value returned is the number
//! of errors detected.  They are embedded in the NEVADA namespace because
//! Parse is such a common name.
//@{
int Parse( Token_Stream *token_stream, 
           const Parse_Table & parse_table,
           Token_Type end_token );
int Parse( Token_Stream *token_stream, 
           const Parse_Table * parse_table,
           Token_Type end_token );

//@}

Keyword* Search_For_Keyword(const char* name, Keyword* table, int length);

// Some useful functions for describing correct syntax of keywords to the
// user.
std::string MainKeyword(const char* name);
std::string MainKeyword(const char* name, const char* line);
std::string MainKeyword(std::string name, std::string line);
std::string MainKeywordId(const char* name);
std::string MainKeywordOptId(const char* name);
std::string SubKeyword(const char* name);
std::string SubKeyword(const char* name, const char* opts);
std::string SubKeyword(const char* name, const char* opts, const char* def);
std::string SubKeyword(const char* name, double def);
std::string SubKeyword(const char* name, int def);
std::string EndKeyword();
std::string NoEndKeyword();
std::string SubSubKeyword(const char* name);
std::string SubSubKeyword(const char* name, const char* opts);
std::string SubSubKeyword(const char* name, const char* opts, const char* def);
std::string SubSubKeyword(const char* name, double def);
std::string SubSubKeyword(const char* name, int def);
std::string SubEndKeyword();
}//end namespace PAMGEN_NEVADA
#endif
