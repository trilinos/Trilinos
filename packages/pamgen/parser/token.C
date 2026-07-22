// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <ctype.h>

#include "token.h"
#include <algorithm>
#include <string>
#include <cstring>

using namespace std;
namespace PAMGEN_NEVADA {
/*****************************************************************************/
Token::Token() : type(TK_NONE) 
/*****************************************************************************/
{
}

/*****************************************************************************/
Token::Token(Token_Type t, Token_Value& sv)
/*****************************************************************************/
:
  type(t)
{
  switch (type){
    case TK_IDENTIFIER:
      s.sval = new char[strlen(sv.sval)+1];
      {
        char *p = sv.sval;
        char *q = s.sval;
        while (*p){
          if (islower(*p)) *q = toupper(*p);
          else             *q = *p;
          ++q; ++p;
        }
        *q = 0;
      }
      break;
    case TK_STRING:
      s.sval = new char[strlen(sv.sval)+1];
      {
        char *p = sv.sval;
        char *q = s.sval;
        while (*p){
          //We don't really want to change file names to upper case
          //if (islower(*p)) *q = toupper(*p);
          //else             *q = *p;
          *q = *p;
          ++q; ++p;
        }
        *q = 0;
      }
      break;
    default:
      s = sv;
      break;
  }
}

/*****************************************************************************/
Token::Token(Token_Type t, No_Value)
/*****************************************************************************/
:
  type(t)
{
  assert(t>=TK_EXIT && t<=TK_ERROR);
  assert(t!=TK_IDENTIFIER);
  assert(t!=TK_INTEGER);
  assert(t!=TK_REAL);
  assert(t!=TK_STRING);

  assert(Type()==t);
}

/*****************************************************************************/
Token::Token(const Token& t)
/*****************************************************************************/
:
  type(t.type)
{
  switch (type){
    case TK_IDENTIFIER:
    case TK_STRING:
      s.sval = new char[strlen(t.s.sval)+1];
      strcpy(s.sval, t.s.sval);
      break;
    default:
      s = t.s;
      break;
  }
}

/*****************************************************************************/
Token& Token::operator=(const Token& t)
/*****************************************************************************/
{
  switch (type){
    case TK_IDENTIFIER:
    case TK_STRING:
      if(s.sval) delete[] s.sval;
    break;
    default: break;
  }

  switch (t.type){
    case TK_IDENTIFIER:
    case TK_STRING:
      s.sval = new char[strlen(t.s.sval)+1];
      strcpy(s.sval, t.s.sval);
      break;
    default:
      s = t.s;
      break;
  }
  type = t.type;
  return *this;
}

/*****************************************************************************/
Token::~Token(void)
/*****************************************************************************/
{
  switch (type){
    case TK_IDENTIFIER:
    case TK_STRING:
      if(s.sval) delete[] s.sval;
      break;
    default:
      break;
  }
}

/*****************************************************************************/
string Token::As_Stripped_String() const
/*****************************************************************************/
  // Strip off the first and last characters, which are assumed to be
  // quotation marks "
{
  assert(Type()==TK_STRING);
  assert(strlen(As_String())>=2);
  assert((As_String()[0]=='\"' && As_String()[strlen(As_String())-1]=='\"') || 
               (As_String()[0]=='\'' && As_String()[strlen(As_String())-1]=='\''));

  string src = As_String();
  string Result(src,1,src.size()-2);

  assert(Result.size()==strlen(As_String())-2);
  assert(!strncmp(Result.c_str(), As_String()+1,Result.size()));
  return Result;
}

/*****************************************************************************/
int Token::Token_Match(const char *t, const char *k)
/*****************************************************************************/
    // Matches two character strings according to the following rules:
    // Two character strings t and k match if each word of t matches
    // each word of k.
    // Two words match if word from t matches the first part of word 
    // from k.
{
  assert(t!=NULL);
  assert(k!=NULL);

// kgb (990831):  new version is not limited in character or word count
  int Result;
  const char *tb = t;
  const char *kb = k;
  while (isspace(*tb)) tb++;
  while (isspace(*kb)) kb++;
  while (*tb && *kb){
    const char *te = tb;
    const char *ke = kb;
    while (*te && !isspace(*te)) te++;
    while (*ke && !isspace(*ke)) ke++;
    if (te-tb<3 || ke-kb<3) {
      // short keywords or tokens must match exactly
      Result = strncmp(tb, kb, std::max(te-tb,ke-kb));
    }else{
      Result = strncmp(tb, kb, std::min(te-tb, ke-kb));
    }
    if (Result) return Result;
    tb = te;
    kb = ke;
    while (isspace(*tb)) tb++;
    while (isspace(*kb)) kb++;
  }
  if (*tb) return 1;
  if (*kb) return -1;
  return 0;
}

/*****************************************************************************/
bool operator!=(const Token &a, const Token &b)
/*****************************************************************************/
{
  if (a.type!=b.type) return true;
  switch (a.type){
    case TK_IDENTIFIER:
    case TK_STRING:
      return strcmp(a.s.sval, b.s.sval);
    case TK_INTEGER:
      return a.s.ival != b.s.ival;
    case TK_REAL:
      return a.s.fval != b.s.fval;
    default:
      return false;
  }
}
/*****************************************************************************/
bool operator==(const Token &a, const Token &b)
/*****************************************************************************/
{
  return !operator!=(a, b);
}
}
