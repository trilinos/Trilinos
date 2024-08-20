// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_commonRTC.hh"

#include <string>
#include <cmath>
#include <typeinfo>
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;
using namespace PG_RuntimeCompiler;

////////////////HELPER FUNCTIONS////////////////////

/*****************************************************************************/
bool PG_RuntimeCompiler::isAssignable(ObjectType type)
/*****************************************************************************/
{
  return (type == ScalarVarOT || type == ArrayIndexOT);
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isValue(ObjectType type)
/*****************************************************************************/
{
  return (type != OperatorOT);
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isVariable(ObjectType type)
/*****************************************************************************/
{
  return (type == ScalarVarOT || type == ArrayVarOT);
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isInt(const string& str)
/*****************************************************************************/
{
  for(unsigned int i = 0; i < str.length(); ++i) {
    char c = str[i];
    if ( !(c >= '0' && c <= '9'))
      return false;
  }
  return true;
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isDouble(const string& str)
/*****************************************************************************/
{
  bool afterDecimal = false, afterExponent = false;

  for(unsigned int i = 0; i < str.length(); ++i) {
    char c = str[i];
    if ( !((c >= '0' && c <= '9') || c == '.' || c == 'e' || c == 'E' ||
           c == '-' || c == '+'))
      return false;
    if (c == '.' && !afterDecimal)
      afterDecimal = true;
    else if (isExp(c) && !afterExponent)
      afterExponent = true;
    else if (c == '.' && afterDecimal)
      return false;
    else if (isExp(c) && afterExponent)
      return false;
    else if (c == '.' && afterExponent)
      return false;
  }

  return afterDecimal || afterExponent;
}

/*****************************************************************************/
char PG_RuntimeCompiler::getNextRealChar(const string& str, unsigned int i)
/*****************************************************************************/
{
  while (i < str.size() && (isWhitespace(str[i]) || str[i] == '\\'))
    ++i;
  if (i < str.size())
    return str[i];
  else
    return 0;
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isString(const string& str)
/*****************************************************************************/
{
  return (str[0] == '\"' && str[str.size()-1] == '\"');
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isChar(const string& str)
/*****************************************************************************/
{
  return (str.length() == 3 && str[0] == '\'' && str[2] == '\'');
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isLetter(char c)
/*****************************************************************************/
{
  return ( (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'));
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isValidVariableChar(char c)
/*****************************************************************************/
{
  return ((c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') ||
          (c >= 'a' && c <= 'z') || (c == '_') || (c == '.'));
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isValidNumericChar(char c)
/*****************************************************************************/
{
  return ((c >= '0' && c <= '9') || c == '.');
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isValidOperatorChar(char c)
/*****************************************************************************/
{
  return ( c=='+' || c=='-' || c=='|' || c=='&' || c=='=' || c=='!' || c=='<'||
           c=='>' || c=='^' || c=='/' || c=='*' || c=='%' || c=='(' || c==')');
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isWhitespace(char o)
/*****************************************************************************/
{
  return (o == 9 || o == 10 || o == 32);
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isExp(char c)
/*****************************************************************************/
{
  return (c == 'E' || c == 'e');
}

/*****************************************************************************/
bool PG_RuntimeCompiler::isLengthTwoOp(const string& op)
/*****************************************************************************/
{
  return (op=="||" || op=="&&" || op=="==" || op=="!=" || op=="<=" ||op==">=");
}

/*****************************************************************************/
string PG_RuntimeCompiler::intToString(int i)
/*****************************************************************************/
{
  static char buf[32];

  sprintf(buf, "%i", i);
  string rv(buf);
  return rv;
}

/*****************************************************************************/
string PG_RuntimeCompiler::tokenTypeToString(TokenType tt)
/*****************************************************************************/
{
    if      (tt == CONSTANT)    return "CONSTANT";
    else if (tt == VAR)         return "VAR";
    else if (tt == OPERATOR)    return "OPERATOR";
    else if (tt == OPENINDEX)   return "OPENINDEX";
    else if (tt == CLOSEINDEX)  return "CLOSEINDEX";
    else if (tt == COMMA)       return "COMMA";
    else if (tt == DECL)        return "DECL";
    else if (tt == BLOCKOPENER) return "BLOCKOPENER";
    else if (tt == FUNC)        return "FUNC";
    else if (tt == WHITE)       return "WHITE";
    else if (tt == OPENBRACE)   return "OPENBRACE";
    else if (tt == CLOSEBRACE)  return "CLOSEBRACE";
    else                        return "SEMICOLON";
}

/*****************************************************************************/
string PG_RuntimeCompiler::typeToString(Type t)
/*****************************************************************************/
{
  if (t == CharT)
    return "char";
  else if (t == FloatT)
    return "float";
  else if (t == DoubleT)
    return "double";
  else if (t == IntT)
    return "int";
  else if (t == LongT)
    return "long";
  else
    return "unknown type";
}


#include "RTC_ExecutableRTC.hh"

using namespace std;

ostream& PG_RuntimeCompiler::operator<<(ostream& os, const Executable& obj)
{
  return obj.operator<<(os);
}
#include "RTC_ObjectRTC.hh"

using namespace std;

ostream& PG_RuntimeCompiler::operator<<(ostream& os, const Object& obj)
{
  return obj.operator<<(os);
}
