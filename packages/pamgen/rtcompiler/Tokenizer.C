// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "RTC_TokenizerRTC.hh"
#include "RTC_commonRTC.hh"

#include <string>
#include <vector>
#include <iostream>
#include <stack>
#include <cassert>

using namespace std;
using namespace PG_RuntimeCompiler;

/*****************************************************************************/
Tokenizer::Tokenizer(vector<Token> tokens)
/*****************************************************************************/
{
  _tokenizedLines.push_back(tokens);

  _currLine  = _tokenizedLines.begin();
  _currToken = _currLine->begin();
}

/*****************************************************************************/
Tokenizer::Tokenizer(const string& body, string& errors)
/*****************************************************************************/
{
  assert(body != "");
  int line = 1;
  unsigned int index = 0;
  bool insideForStatement = false;
  bool blockOpeningStatement = false;

  while(index < body.size()) {
    vector<Token> tokens;
    string currStr = "";
    string lastStr = "";
    TokenType previous = WHITE;

    while(index < body.size()) {

      ////////////////////////////////////////////////////////////////////////
      //process what appears to be a number
      ////////////////////////////////////////////////////////////////////////
      if      (isValidNumericChar(body[index])) {
        bool lastCharWasExp = false;
        for ( ; index < body.size() &&
                (isValidNumericChar(body[index]) ||
                 isExp(body[index]) ||
                 ((body[index]=='-' || body[index]=='+') && lastCharWasExp));
              ++index) {
          currStr += body[index];
          lastCharWasExp = isExp(body[index]);
        }
        //number format will be checked later
        if (index < body.size() && isLetter(body[index])) {
          errors += "Error at line: " + intToString(line) +
                    " letters embedded in a number.";
          return;
        }
        if (previous != OPERATOR && previous != OPENINDEX &&
            previous != COMMA && previous != WHITE) {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede: " + currStr;
          return;
        }
        tokens.push_back(Token(CONSTANT, currStr, line));
        previous = CONSTANT;
        lastStr  = currStr;
      }

      ////////////////////////////////////////////////////////////////////////
      //process a name of some sort, could be variable, keyword, or function
      ////////////////////////////////////////////////////////////////////////
      else if (isValidVariableChar(body[index])) {
        for (;index < body.size() && isValidVariableChar(body[index]); ++index)
          currStr += body[index];

        if (currStr == "int" || currStr == "long" || currStr == "float" ||
            currStr == "double" || currStr == "char") {
          if (lastStr != "(" && previous != WHITE && previous != COMMA) {
            errors += "Error at line: " + intToString(line) + " :" + lastStr +
              " cannot precede: " + currStr;
            return;
          }
          tokens.push_back(Token(DECL, currStr, line));
          previous = DECL;
        }
        //Note: Block opener checks must come before function check
        else if (currStr == "for" || currStr == "while" || currStr == "if" ||
                 currStr == "else") {
          if (previous != WHITE && !(currStr == "if" && lastStr == "else")) {
            errors += "Error at line: " + intToString(line) + " :" + lastStr +
              " cannot precede: " + currStr;
            return;
          }
          if (currStr != "else" && getNextRealChar(body, index) != '(') {
            errors += "Error at line: " + intToString(line) + " block opening statements (if, while etc) must be followed by a (";
            return;
          }
          tokens.push_back(Token(BLOCKOPENER, currStr, line));
          if (currStr == "for")
            insideForStatement = true;
	  //Important: else if should be a single token
	  if (currStr == "if" && lastStr == "else" && previous == BLOCKOPENER){
	    tokens.pop_back();
	    tokens.pop_back();
	    tokens.push_back(Token(BLOCKOPENER, "else if", line));
	  }
          previous = BLOCKOPENER;
          blockOpeningStatement = true;
        }
        else if (getNextRealChar(body, index) == '(') {
          if (previous != OPERATOR && previous != WHITE &&
              previous != COMMA && previous != OPENINDEX) {
            errors += "Error at line: " + intToString(line) + " :" + lastStr +
              " cannot precede: " + currStr;
            return;
          }
          tokens.push_back(Token(FUNC, currStr, line));
          previous = FUNC;
        }
        else {
          if (previous != DECL && previous != OPERATOR && previous != WHITE &&
              previous != COMMA && previous != OPENINDEX) {
            errors += "Error at line: " + intToString(line) + " :" + lastStr +
              " cannot precede: " + currStr;
            return;
          }
          tokens.push_back(Token(VAR, currStr, line));
          previous = VAR;
        }
        lastStr = currStr;
      }

      ////////////////////////////////////////////////////////////////////////
      //process operator
      ////////////////////////////////////////////////////////////////////////
      else if (isValidOperatorChar(body[index])) {
        currStr = body[index++];
        if (index < body.size() && isLengthTwoOp(currStr+body[index]))
          currStr += body[index++];
        else if (index < body.size() && ((currStr+body[index]) == "/*")) {
          //we have a comment;
          ++index;
          while (index+1 < body.size() &&
                 (body[index] != '*' || body[index+1] != '/')) {
            if (body[index] == '\n')
              ++line;
            ++index; //skip over comment
          }
          if (index+1 == body.size()) {
            errors += "Error at line: " + intToString(line) +
              ": Comment was never terminated.";
            return;
          }
          else {
            assert(body[index] == '*');
            ++index;
            assert(body[index] == '/');
            ++index;
            continue;
          }
        }

        if (previous == DECL) {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede: " + currStr;
          return;
        }

        //Important: The difference between negation and minus is detected
        //here. I am assuming that if a minus symbol is preceded by an operator
        //other than a closing parentheses, nothing, an open index, or a comma,
        //then it means negation.
        if ( currStr == "-" &&
             ( (previous == OPERATOR && lastStr != ")") ||
               (previous == OPENINDEX) || (previous == COMMA) ||
               (previous == WHITE)))
          currStr = "_"; //different symbol for negation
        //Important: for-loops violate many of the conventions followed in the
        //rest of the code. To compensate, the first ( after the for statement
        //will be ignored as will the last ) in the for statement.
        bool skip = ( (currStr == "(" && tokens.size() == 1 && previous == BLOCKOPENER && lastStr == "for") || (currStr == ")" && insideForStatement && getNextRealChar(body,index+1) == '{'));
        if (!skip) {
          if (currStr == "(" || currStr == ")") {
            if (!checkStack(currStr)) {
              errors += "Error at line: " + intToString(line) +
                " parentheses error.";
              return;
            }
          }
          tokens.push_back(Token(OPERATOR, currStr, line));
        }
        previous = OPERATOR;
        lastStr = currStr;
      }

      ////////////////////////////////////////////////////////////////////////
      //process whitespace -- skips it
      ////////////////////////////////////////////////////////////////////////
      else if (isWhitespace(body[index])) {
        for(; index < body.size() && isWhitespace(body[index]); ++index) {
          if (body[index] == '\n')
            ++line;
        }
      }

      ////////////////////////////////////////////////////////////////////////
      //process the character opener
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == '\'') {
        if (previous != OPERATOR && previous != OPENINDEX &&
            previous != COMMA && previous != WHITE) {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede \'";
          return;
        }
        currStr = body.substr(index, 3); //safe, substr wont overrun array
        if (isChar(currStr))
          tokens.push_back(Token(CONSTANT, currStr, line));
        else {
          errors += "Error at line: " + intToString(line) +
            " incorrect use of character tick.";
          return;
        }
        index += 3;
        previous = CONSTANT;
        lastStr = currStr;
      }

      ////////////////////////////////////////////////////////////////////////
      //process the string opener
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == '\"') {
        if (lastStr != "(" && lastStr != ",") {
          errors += "Error at line: " + intToString(line) +
            ": strings can only be used as function arguments.";
          return;
        }
        currStr += body[index++];
        while (index < body.size() && body[index] != '\"') {
          currStr += body[index++];
        }
        currStr += body[index];
        if (index == body.size()) {
          errors += "Error at line: " + intToString(line) +
            " string not terminated.";
          return;
        }
        tokens.push_back(Token(CONSTANT, currStr, line));
        assert(body[index] == '\"');
        index++;
        previous = CONSTANT;
        lastStr = currStr;
      }

      ////////////////////////////////////////////////////////////////////////
      //process the backslash - the prefix of an aprepro-safe brace
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == '\\') {
        //I believe we can simply ignore it
        ++index;
      }

      ////////////////////////////////////////////////////////////////////////
      //process opening brace
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == '{') {
        if (lastStr != ")" && lastStr != "else") {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede: {";
          return;
        }
        checkStack("{");
        if (!blockOpeningStatement) {
          errors += "Error at line: " + intToString(line) +
            " only block opening statements(if,while etc) can end with {.";
          return;
        }
        insideForStatement = false;
        blockOpeningStatement = false;
        tokens.push_back(Token(OPENBRACE, "{", line));
        ++index;
        break; // { finishes the line
      }

      ////////////////////////////////////////////////////////////////////////
      //process closing brace
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == '}') {
        if (previous != WHITE) {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede: }";
          return;
        }
        if (!checkStack("}")) {
          errors += "Error at line: " + intToString(line) +
            " mismatched }.";
          return;
        }
        tokens.push_back(Token(CLOSEBRACE, "}", line));
        ++index;
        break; // } finishes the line
      }

      ////////////////////////////////////////////////////////////////////////
      //process opening index
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == '[') {
        if (previous != VAR) {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede: [";
          return;
        }
        checkStack("[");
        tokens.push_back(Token(OPENINDEX, "[", line));
        ++index;
        previous = OPENINDEX;
        lastStr = "[";
      }

      ////////////////////////////////////////////////////////////////////////
      //process closing index
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == ']') {
        if (previous != CONSTANT && previous != VAR && previous != OPERATOR) {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede: ]";
          return;
        }
        if (!checkStack("]")) {
          errors += "Error at line: " + intToString(line) +
            " mismatched ].";
          return;
        }
        tokens.push_back(Token(CLOSEINDEX, "]", line));
        ++index;
        previous = CLOSEINDEX;
        lastStr = "]";
      }

      ////////////////////////////////////////////////////////////////////////
      //process comma
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == ',') {
        if (previous != CONSTANT && previous != VAR && previous != OPERATOR &&
            previous != CLOSEINDEX) {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede: ,";
          return;
        }
        tokens.push_back(Token(COMMA, ",", line));
        ++index;
        previous = COMMA;
        lastStr = ",";
      }

      ////////////////////////////////////////////////////////////////////////
      //process semicolon
      ////////////////////////////////////////////////////////////////////////
      else if (body[index] == ';') {
        if (previous == OPENINDEX || previous == COMMA ||
            previous == BLOCKOPENER) {
          errors += "Error at line: " + intToString(line) + " :" + lastStr +
            " cannot precede: ;";
          return;
        }
        if (blockOpeningStatement && !insideForStatement) {
          errors += "Error at line: " + intToString(line) +
            " block opening statements(if,while etc) must end with {.";
          return;
        }
        tokens.push_back(Token(SEMICOLON, ";", line));
        ++index;
        break; // ; finishes the line
      }

      else {
        errors += "Error at line: " + intToString(line) +
          " unknown symbol: " + body[index];
        return;
      }

      currStr = "";
    }

    //check stack -- should only contain braces
    for(unsigned int i = 0; i < _stack.size(); ++i) {
      if (_stack[i] != "}" && _stack[i] != "{") {
        errors += "Error at line: " + intToString(line) +
          " mismatched: " + _stack[i];
        return;
      }
    }

    if (tokens.size() > 0)
      _tokenizedLines.push_back(tokens);
  }

  if (!_stack.empty()) {
    errors += "Error mismatched {} in your function";
    return;
  }

  _currLine  = _tokenizedLines.begin();
  _currToken = _currLine->begin();
}

/*****************************************************************************/
bool Tokenizer::checkStack(const string& op)
/*****************************************************************************/
{
  if (_stack.empty() && (op == ")" || op == "]" || op == "}"))
    return false;

  if (op == ")") {
    if (_stack.back() != "(")
      return false;
    _stack.pop_back();
  }
  else if (op == "}") {
    if (_stack.back() != "{")
      return false;
    _stack.pop_back();
  }
  else if (op == "]") {
    if (_stack.back() != "[")
      return false;
    _stack.pop_back();
  }
  else if (op == "(" || op == "{" || op == "[") {
    _stack.push_back(op);
  }
  else {
    cout << "Unknown stack symbol: " << op << endl;
    return false;
  }
  return true;
}

/*****************************************************************************/
bool Tokenizer::check(const string& testName, const string& errs,
                      unsigned int expectedNumLines,
                      unsigned int* expectedNumTokens,
                      Token** expectedTokens)
/*****************************************************************************/
{
  cout << testName << " begin." << endl;
  if (errs != "") {
    cerr << "The test had errors: " << errs << endl;
    return false;
  }
  if (_tokenizedLines.size() != expectedNumLines) {
    cerr << "Expected numLines = " << expectedNumLines << " actual = "
         << _tokenizedLines.size() << endl;
    return false;
  }
  for (unsigned int i = 0; i < expectedNumLines; ++i) {
    if (_tokenizedLines[i].size() != expectedNumTokens[i]) {
      cerr << "Expected numTokens for line" << i << ": = "
           << expectedNumTokens[i] << " actual = " << _tokenizedLines[i].size()
           << endl;
      return false;
    }
    for (unsigned int j = 0; j < expectedNumTokens[i]; ++j) {
      if (_tokenizedLines[i][j] != expectedTokens[i][j]) {
        cerr << "Token mismatch at line: " << i << " and token: " << j << endl;
        cerr << "Expected: " << expectedTokens[i][j].toString()
             << " actual: " << _tokenizedLines[i][j].toString() << endl;
      }
    }
  }
  cout << testName << " passed." << endl;
  return true;
}

/*****************************************************************************/
void Tokenizer::print()
/*****************************************************************************/
{
  for(unsigned int i = 0; i < _tokenizedLines.size(); ++i) {
    for (unsigned int j = 0; j < _tokenizedLines[i].size(); ++j) {
      cout << _tokenizedLines[i][j].toString() << " ";
    }
    cout << endl;
  }
}

/*****************************************************************************/
void Tokenizer::printCurrLine()
/*****************************************************************************/
{
  assert(_currLine != _tokenizedLines.end());
  for (unsigned int i = 0; i < _currLine->size(); ++i)
    cout << (*_currLine)[i].toString() << " ";
  cout << endl;
}

/*****************************************************************************/
void Tokenizer::nextLine()
/*****************************************************************************/
{
  assert(_currLine != _tokenizedLines.end());

  _currLine++;
  if (_currLine != _tokenizedLines.end())
    _currToken = _currLine->begin();
}

/*****************************************************************************/
void Tokenizer::nextToken()
/*****************************************************************************/
{
  assert(_currLine != _tokenizedLines.end());
  assert(_currToken != _currLine->end());

  _currToken++;
}

/*****************************************************************************/
void Tokenizer::previousToken()
/*****************************************************************************/
{
  assert(_currLine != _tokenizedLines.end());
  assert(_currToken != _currLine->begin());

  _currToken--;
}

/*****************************************************************************/
int Tokenizer::lineNum()
/*****************************************************************************/
{
  assert(_currLine != _tokenizedLines.end());
  if (_currToken == _currLine->end()) {
    _currToken--;
    int rv =  _currToken->lineLoc();
    _currToken++;
    return rv;
  }
  else
    return _currToken->lineLoc();
}


/*****************************************************************************/
void Tokenizer::test()
/*****************************************************************************/
{
  string errs = "";
  Token** lines;

  string body = "     i = 5.4 * varname;";
  Tokenizer test1(body, errs);
  unsigned int expectedNumTokens = 6;
  lines = new Token*[1];
  Token tokens1[6] = {Token(VAR, "i",0), Token(OPERATOR, "=",0),
                      Token(CONSTANT, "5.4",0), Token(OPERATOR, "*",0),
                      Token(VAR, "varname",0), Token(SEMICOLON, ";",0)};
  lines[0] = tokens1;
  if (!test1.check("test1", errs, 1, &expectedNumTokens, lines))
    return;


  body = "var = 2341.12113e-123;";
  Tokenizer test2(body, errs);
  expectedNumTokens = 4;
  Token tokens2[4] = {Token(VAR, "var",0), Token(OPERATOR, "=",0),
                      Token(CONSTANT, "2341.12113e-123",0),Token(SEMICOLON,";",0)};
  lines[0] = tokens2;
  if (!test2.check("test2", errs, 1, &expectedNumTokens, lines))
    return;


  body = "int i_3[5*6] = func(var/((err-2) * 7), phi);";
  Tokenizer test3(body, errs);
  expectedNumTokens = 25;
  Token tokens3[25] = {Token(DECL, "int",0), Token(VAR, "i_3",0), Token(OPENINDEX, "[",0), Token(CONSTANT, "5",0), Token(OPERATOR, "*",0), Token(CONSTANT, "6",0), Token(CLOSEINDEX, "]",0), Token(OPERATOR, "=",0), Token(FUNC, "func",0), Token(OPERATOR, "(",0), Token(VAR, "var",0), Token(OPERATOR, "/",0), Token(OPERATOR, "(",0), Token(OPERATOR, "(",0), Token(VAR, "err",0), Token(OPERATOR, "-",0), Token(CONSTANT, "2",0), Token(OPERATOR, ")",0), Token(OPERATOR, "*",0), Token(CONSTANT, "7",0), Token(OPERATOR, ")",0), Token(COMMA, ",",0), Token(VAR, "phi",0), Token(OPERATOR, ")",0), Token(SEMICOLON, ";",0)};
  lines[0] = tokens3;
  if (!test3.check("test3", errs, 1, &expectedNumTokens, lines))
    return;


  body = "for (int i = 'e'; i <= 5; i = i + 1)   { }";
  Tokenizer test4(body, errs);
  delete[] lines;
  lines = new Token*[4];
  unsigned int expectedTokens[4] = {6, 4, 6, 1};
  Token tokens4[6] = {Token(BLOCKOPENER, "for",0), Token(DECL, "int",0),
                      Token(VAR, "i",0), Token(OPERATOR, "=",0),
                      Token(CONSTANT, "'e'",0), Token(SEMICOLON, ";",0)};
  Token tokens5[4] = {Token(VAR, "i",0), Token(OPERATOR, "<=",0),
                      Token(CONSTANT, "5",0), Token(SEMICOLON, ";",0)};
  Token tokens6[6] = {Token(VAR, "i",0), Token(OPERATOR, "=",0), Token(VAR, "i",0),
                      Token(OPERATOR, "+",0), Token(CONSTANT, "1",0),
                      Token(OPENBRACE, "{",0)};
  Token tokens7[1] = {Token(CLOSEBRACE, "}",0)};
  lines[0] = tokens4;
  lines[1] = tokens5;
  lines[2] = tokens6;
  lines[3] = tokens7;
  if (!test4.check("test4", errs, 4, expectedTokens, lines))
    return;


  body = "-x = -func(-(e-y), -array[-r*-8]);";
  Tokenizer test5(body, errs);
  expectedNumTokens = 24;
  Token tokens8[24] = {Token(OPERATOR, "_",0), Token(VAR, "x",0), Token(OPERATOR, "=",0), Token(OPERATOR, "_",0), Token(FUNC, "func",0), Token(OPERATOR, "(",0), Token(OPERATOR, "_",0), Token(OPERATOR, "(",0), Token(VAR, "e",0), Token(OPERATOR, "-",0), Token(VAR, "y",0), Token(OPERATOR, ")",0), Token(COMMA, ",",0), Token(OPERATOR, "_",0), Token(VAR, "array",0), Token(OPENINDEX, "[",0), Token(OPERATOR, "_",0), Token(VAR, "r",0), Token(OPERATOR, "*",0), Token(OPERATOR, "_",0), Token(CONSTANT, "8",0), Token(CLOSEINDEX, "]",0), Token(OPERATOR, ")",0), Token(SEMICOLON, ";",0) };
  lines[0] = tokens8;
  if (!test5.check("test5", errs, 1, &expectedNumTokens, lines))
    return;


  body = "else if (condition) \\{ \\}";
  Tokenizer test6(body, errs);
  unsigned int expectedToks[2] = {5, 1};
  Token tokens9[5] = {Token(BLOCKOPENER, "else if",0), Token(OPERATOR, "(",0),
		      Token(VAR, "condition",0), Token(OPERATOR, ")",0),
		      Token(OPENBRACE, "{",0)};
  Token tokens10[1] = {Token(CLOSEBRACE, "}",0)};
  lines[0] = tokens9;
  lines[1] = tokens10;
  if(!test6.check("test6", errs, 2, expectedToks, lines))
    return;
  delete[] lines;


  //check to see if errors are generated when expected
  body = "1231asda;";
  Tokenizer test7(body, errs);
  assert(errs != "");
  errs = "";

  body = "char c = 'asda';";
  Tokenizer test8(body, errs);
  assert(errs != "");
  errs = "";

  body = "int i = array[123(]);";
  Tokenizer test9(body, errs);
  assert(errs != "");
  errs = "";
}
