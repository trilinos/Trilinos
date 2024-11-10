// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _COMMONRTC_H
#define _COMMONRTC_H

#include <string>
#include <typeinfo>
#include <cmath>

namespace PG_RuntimeCompiler {

//errs must be defined
#define CHECKERR(c,m) \
    if(c) { errs += "Error at line " + intToString(tokens.lineNum()) + ": " + m; return; }

#define CHECKARGERR(c,m) \
    if(c) { errs += m; return; }

/**
 * TokenType is an enum that tell us what type a token is
 */
enum TokenType {CONSTANT,    /**< constant value - ex: 4.5, 3, 'a', "str"  */
                VAR,         /**< variable name  - ex: i, yTy32            */
                OPERATOR,    /**< operator - ex: +, -, ^                   */
                OPENINDEX,   /**< open square brace:   [                   */
                CLOSEINDEX,  /**< closed square brace: ]                   */
                COMMA,       /**< a comma: ,                               */
                DECL,        /**< a variable declaration keyword - ex: int */
                BLOCKOPENER, /**< a blockopening keyword - ex: if, while   */
                FUNC,        /**< a function call - ex: sin(8)             */
                WHITE,       /**< whitespace                               */
                OPENBRACE,   /**< open curly brace:   {                    */
                CLOSEBRACE,  /**< closed curly brace: }                    */
                SEMICOLON    /**< a semicolon: ;                           */
};

std::string tokenTypeToString(TokenType tt);

/**
 * A Token is a single item in our program
 */
class Token {
 private:
  TokenType _type;  //!< The type of the token
  std::string _value; //!< The value of the token
  int _lineLoc; //!< The line number location of this token

 public:

  /**
   * Constructor -> Trivial
   *
   * @param tt    - The token's type
   * @param value - The token's value
   * @param loc   - The line number location of this token
   */
  Token(TokenType Ttt, const std::string& Tvalue, int Tloc) {
    _value   = Tvalue;
    _type    = Ttt;
    _lineLoc = Tloc;
  }

  Token() { }

  /**
   * get_type -> Returns the type of the token
   */
  TokenType type() const {return _type;}

  bool operator!=(const Token& rhs) {
    return (_value != rhs._value || _type != rhs._type);
  }

  bool operator==(const Token& rhs) {
    return (_value == rhs._value && _type == rhs._type);
  }

  /**
   * value -> Returns the value of the token
   */
  const std::string& value() const {return _value;}

  int lineLoc() const {return _lineLoc;}

  std::string toString() const {
    return (tokenTypeToString(_type) + "::" + _value);
  }
};

/**
 * Type
 * An enum to specify the type of a value
 */
enum Type {CharT, IntT, LongT, FloatT, DoubleT};

/**
 * Some meta functions to help with Type
 */
template <typename T>
struct TypeToTypeT;

template<> struct TypeToTypeT<char>   { static const Type value = CharT; };
template<> struct TypeToTypeT<int>    { static const Type value = IntT; };
template<> struct TypeToTypeT<long>   { static const Type value = LongT; };
template<> struct TypeToTypeT<float>  { static const Type value = FloatT; };
template<> struct TypeToTypeT<double> { static const Type value = DoubleT; };

/**
 * ObjectType
 * An enum to specify what an object really is
 */
enum ObjectType {OperatorOT, ScalarNumberOT, ArrayNumberOT, ScalarVarOT,
                 ArrayIndexOT, ArrayVarOT, FunctionOT};


////////////////HELPER FUNCTIONS////////////////////

/**
 * isAssignable -> Tells us if a certain object type can be assigned to
 *
 * @param type - The type we are inquiring about
 */
bool isAssignable(ObjectType type);

/**
 * isValue -> Tells us if a certain object type can be interpreted as a value
 *
 * @param type - The type we are inquiring about
 */
bool isValue(ObjectType type);

/**
 * isVariable -> Tells us if a certain object type can be considered to be a
 *               variable
 *
 * @param type - The type we are inquiring about
 */
bool isVariable(ObjectType type);

/**
 * isValidVariableChar -> This method returns true if character c can be part
 *                        of a variable name.
 *
 * @param c - The character we are looking at
 */
bool isValidVariableChar(char c);

/**
 * isValidOperatorChar -> This method returns true if character c can be part
 *                        of an operator.
 *
 * @param c - The character we are looking at
 */
bool isValidOperatorChar(char c);

/**
 * isValidNumericChar -> This method returns true if character c can be part
 *                       of a number.
 *
 * @param c - The character we are looking at
 */
bool isValidNumericChar(char c);

/**
 * isWhitespace -> This method returns true if character c is whitespace
 *
 * @param c - The character we are looking at
 */
bool isWhitespace(char c);

/**
 * isLengthTwoOp -> This method returns true if op is a 2-char operator
 *
 * @param op - The string representation of the operator we are looking at
 */
bool isLengthTwoOp(const std::string& op);

/**
 * isInt -> This method returns true if s contains a valid int
 *
 * @param s - The string we are looking at
 */
bool isInt(const std::string& s);

/**
 * isDouble -> This method returns true if s contains a valid double
 *
 * @param s - The string we are looking at
 */
bool isDouble(const std::string& s);

/**
 * isChar -> This method returns true if s contains a valid char
 *
 * @param s - The string we are looking at
 */
bool isChar(const std::string& s);

/**
 *
 */
bool isString(const std::string& s);

/**
 * intToString -> This method returns int i in string form. 46 -> "46"
 *
 * @param i - The integer we are converting into a string
 */
std::string intToString(int i);

std::string typeToString(Type t);

char getNextRealChar(const std::string& str, unsigned int i);

bool isLetter(char c);

bool isExp(char c);
}
#endif
