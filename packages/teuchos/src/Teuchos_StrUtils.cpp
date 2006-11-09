// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TestForException.hpp"


using namespace Teuchos;


Array<string> StrUtils::readFile(istream& is, char comment)
{
  string line;
  Array<string> rtn(0);

  while (readLine(is, line))
    {
      if (line.length() > 0) rtn.append(before(line, comment));
      line="";
    }
	
  return rtn;
}

Array<string> StrUtils::splitIntoLines(const string& input)
{
  int begin = 0;
  Array<string> rtn;

  for (unsigned int p=0; p<input.length(); ++p) {
    const bool isEnd = p==input.length()-1;
    if( input[p]=='\n' || input[p]=='\0' || input[p]=='\r' || isEnd )
    {
      if (p-begin > 1) rtn.append(subString(input, begin, p+(isEnd?1:0)));
      begin = p+1;
    }
  }
  return rtn;
}

Array<Array<string> > StrUtils::tokenizeFile(istream& is, char comment)
{
  string line;
  Array<Array<string> > rtn(0);
  Array<string> lines = readFile(is, comment);
  rtn.reserve(lines.length());
	
  int count = 0;
  for (int i=0; i<lines.length(); i++)
    {
      if (lines[i].length() == 0) continue;
      Array<string> tokens = stringTokenizer(lines[i]);
      if (tokens.length() == 0) continue;
      rtn.append(tokens);
      count++;
    }
	
  return rtn;
}

bool StrUtils::readLine(istream& is, string& line)
{
  char c[500];
  if (line.length() > 0) line[0] = '\0';
	
  if (is.eof()) return false;
  if (is.getline(c, 499))
    {
      line = string(c);
    }
	
  return true;
}

	

Array<string> StrUtils::getTokensPlusWhitespace(const string& str){
  Array<string> rtn(0);
  unsigned int start = 0;
	
  while(start < str.length())
    {
      unsigned int wordStart =  findNextNonWhitespace(str, start);
			/* add any preceding whitespace */
			if (wordStart > start)
				{
					rtn.append(subString(str, start, wordStart));
				}
			start = wordStart;
			/* add the next word */
      int stop = findNextWhitespace(str, start);
      if (start-stop == 0) return rtn;
      string sub = subString(str, start, stop);
      rtn.append(sub);
			start = stop;// findNextNonWhitespace(str, stop);
    }
  return rtn;
}

Array<string> StrUtils::stringTokenizer(const string& str){
  Array<string> rtn(0);
  unsigned int start = 0;
	
  while(start < str.length())
    {
      start =  findNextNonWhitespace(str, start);
      int stop = findNextWhitespace(str, start);
      if (start-stop == 0) return rtn;
      string sub = subString(str, start, stop);
      rtn.append(sub);
      start =  findNextNonWhitespace(str, stop);
    }
  return rtn;
}

string StrUtils::reassembleFromTokens(const Array<string>& tokens, int iStart)
{
  string rtn;

  for (int i=iStart; i<tokens.length(); i++) 
    {
      rtn += tokens[i];
      if (i < (tokens.length()-1)) rtn += " ";
    }
  return rtn;
}

void StrUtils::splitList(const string& big, Array<string>& list) 
{
  if (subString(big, 0,1)!="[") 
    {
      list.resize(1);
      list[0] = big;
      return;
    }
	
  int parenDepth = 0;
  int localCount = 0;
  string tmp(big);
  list.resize(0);

  // start at 1 to ignore '[';
	
  for (unsigned int i=1; i<big.length(); i++)
    {
      if (big[i]=='(') parenDepth++;
      if (big[i]==')') parenDepth--;
      if (big[i]==']') 
	{
	  tmp[localCount]='\0'; 
	  list.append(tmp);
	  break;
	}
      if (big[i]==',' && parenDepth==0)
	{
	  tmp[localCount]='\0';
	  list.append(tmp);
	  tmp = big;
	  localCount = 0;
	  continue;
	}
      tmp[localCount] = big[i];
      localCount++;
    }
}
							

// return the position of the next whitespace in a string. 
// If no whitespace, return -1;

int StrUtils::findNextWhitespace(const string& str, int offset)
{
  for (unsigned int i=0; i<(str.length()-offset); i++)
    {
      if (str[i+offset]==' ' || str[i+offset]=='\t' || str[i+offset]=='\n')
	{
	  return i+offset;
	}
    }
  return str.length();
}

int StrUtils::findNextNonWhitespace(const string& str, int offset)
{
  for (unsigned int i=0; i<(str.length()-offset); i++)
    {
      if (!(str[i+offset]==' ' || str[i+offset]=='\t' || str[i+offset]=='\n'))
	{
	  return i+offset;
	}
    }
  return str.length();
}


string StrUtils::varTableSubstitute(const string& rawLine,
				    const Array<string>& varNames,
				    const Array<string>& varValues)
{
  TEST_FOR_EXCEPTION(varNames.length() != varValues.length(),
                     runtime_error,
                     "mismatched variable tables in varTableSubstitute");
                     
  string line = rawLine;
  for (int i=0; i<varNames.length(); i++)
    {
      line = varSubstitute(line, varNames[i], varValues[i]);
    }
  return line;
}




string StrUtils::varSubstitute(const string& rawLine, 
			       const string& varName, 
			       const string& varValue)
{
  string line = rawLine;
  
  // iterate because there might be more than one occurance on this line
  while (find(line, varName) >= 0)
    {
      string b = before(line, varName);
      string a = after(line, varName);
      line = b + varValue + a;
    }
  return line;
}


string StrUtils::before(const string& str, char sub)
{
  char c[2];
  c[0] = sub;
  c[1] = 0;
  return before(str, c);
}

string StrUtils::before(const string& str, const string& sub)
{
  TEST_FOR_EXCEPTION(sub.c_str()==0,
                     runtime_error, "String::before: arg is null pointer");

  char* p = strstr((char*) str.c_str(), (char*) sub.c_str());
  if (p==0) return str;
  int subLen = p-str.c_str();
  string rtn(str.c_str(), subLen);
  return rtn;
}

string StrUtils::after(const string& str, const string& sub)
{
  TEST_FOR_EXCEPTION(sub.c_str()==0,
                     runtime_error, "String::after: arg is null pointer");

  // find beginning of substring
  char* p = strstr((char*) str.c_str(), (char*) sub.c_str()) ;
  // if substring not found, return empty string
  if (p==0) return string();
  // offset to end of substring
  p+= strlen(sub.c_str());
  return string(p);
}

int StrUtils::find(const string& str, const string& sub)
{
  char* p = strstr((char*) str.c_str(), (char*) sub.c_str());
  if (p==0) return -1;
  return p-str.c_str();
}

bool StrUtils::isWhite(const string& str)
{
  for (unsigned int i=0; i<str.length(); i++)
    {
      unsigned char c = str[i];
      if (c >= 33 && c <= 126)
        {
          return false;
        }
    }
  return true;
}

string StrUtils::fixUnprintableCharacters(const string& str)
{
  string rtn = str;
  for (unsigned int i=0; i<rtn.length(); i++)
    {
      unsigned char c = rtn[i];
      if (c < 33 || c > 126) 
        {
          if (c != '\t' && c != '\n'&& c != '\r' && c != '\f' && c != ' ')
            {
              rtn[i] = ' ';
            }
        }
    }
  return rtn;
}

string StrUtils::between(const string& str, const string& begin,
			 const string& end, string& front,
			 string& back)
{
  front = before(str, begin);
  string middle = before(after(str, begin), end);
  back = after(str, end);
  return middle;
}


string StrUtils::subString(const string& str, int begin, int end)
{
return string(str.c_str()+begin, end-begin);
}

string StrUtils::readFromStream(istream& is)
{
  TEST_FOR_EXCEPTION(true, logic_error, 
                     "StrUtils::readFromStream isn't implemented yet");

	return "";
}

string StrUtils::allCaps(const string& s)
{
  string rtn = s;
  for (unsigned int i=0; i<rtn.length(); i++)
    {
      rtn[i] = toupper(rtn[i]);
    }
  return rtn;
}

double StrUtils::atof(const string& s)
{
	return ::atof(s.c_str());
}

int StrUtils::atoi(const string& s)
{
	return ::atoi(s.c_str());
}

std::ostream& StrUtils::printLines(
  std::ostream             &os
  ,const std::string       &linePrefix
  ,const std::string       &lines
  )
{
  typedef Teuchos::Array<string> array_t;
  array_t linesArray = splitIntoLines(lines);
  for( int i = 0; i < static_cast<int>(linesArray.size()); ++i )
  {
    os << linePrefix << linesArray[i] << "\n";
  }
  return os;
}
  

	
