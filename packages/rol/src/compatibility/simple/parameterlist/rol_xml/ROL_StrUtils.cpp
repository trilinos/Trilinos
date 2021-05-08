#include "ROL_StrUtils.hpp"
#include "Teuchos_Assert.hpp"

namespace ROL {


std::vector<std::string> StrUtils::readFile(std::istream& is, char comment)
{
  std::string line;
  std::vector<std::string> rtn(0);

  while (readLine(is, line))
  {
    if (line.length() > 0) rtn.push_back(before(line, comment));
    line="";
  }

  return rtn;
}


std::vector<std::string> StrUtils::splitIntoLines(const std::string& input)
{
  int begin = 0;
  std::vector<std::string> rtn;
  const unsigned int len = input.length();
  for (unsigned int p=0; p<len; ++p) {
    const bool isEnd = p==len-1;
    if( input[p]=='\n' || input[p]=='\0' || input[p]=='\r' || isEnd )
    {
      if (p-begin > 1)
        rtn.push_back(
          subString( input, begin, p+(isEnd?(input[len-1]=='\n'?0:1):0) )
          );
      begin = p+1;
    }
  }
  return rtn;
}


std::vector<std::vector<std::string> > StrUtils::tokenizeFile(std::istream& is, char comment)
{
  std::string line;
  std::vector<std::vector<std::string> > rtn(0);
  std::vector<std::string> lines = readFile(is, comment);
  rtn.reserve(lines.size());

  int count = 0;
  for (unsigned int i=0; i<lines.size(); i++)
  {
    if (lines[i].length() == 0) continue;
    std::vector<std::string> tokens = stringTokenizer(lines[i]);
    if (tokens.size() == 0) continue;
    rtn.push_back(tokens);
    count++;
  }

  return rtn;
}


bool StrUtils::readLine(std::istream& is, std::string& line)
{
  char c[500];
  if (line.length() > 0) line[0] = '\0';

  if (is.eof()) return false;
  if (is.getline(c, 499))
  {
    line = std::string(c);
  }

  return true;
}


std::vector<std::string> StrUtils::getTokensPlusWhitespace(const std::string& str){
  std::vector<std::string> rtn(0);
  unsigned int start = 0;

  while(start < str.length())
  {
    unsigned int wordStart =  findNextNonWhitespace(str, start);
    /* add any preceding whitespace */
    if (wordStart > start)
    {
      rtn.push_back(subString(str, start, wordStart));
    }
    start = wordStart;
    /* add the next word */
    int stop = findNextWhitespace(str, start);
    if (start-stop == 0) return rtn;
    std::string sub = subString(str, start, stop);
    rtn.push_back(sub);
    start = stop;// findNextNonWhitespace(str, stop);
  }
  return rtn;
}


std::vector<std::string> StrUtils::stringTokenizer(const std::string& str){
  std::vector<std::string> rtn(0);
  unsigned int start = 0;

  while(start < str.length())
  {
    start =  findNextNonWhitespace(str, start);
    int stop = findNextWhitespace(str, start);
    if (start-stop == 0) return rtn;
    std::string sub = subString(str, start, stop);
    rtn.push_back(sub);
    start =  findNextNonWhitespace(str, stop);
  }
  return rtn;
}


std::string StrUtils::reassembleFromTokens(const std::vector<std::string>& tokens,
  int iStart)
{
  std::string rtn;

  for (unsigned int i=iStart; i<tokens.size(); i++)
  {
    rtn += tokens[i];
    if (i < (tokens.size()-1)) rtn += " ";
  }
  return rtn;
}


void StrUtils::splitList(const std::string& big, std::vector<std::string>& list)
{
  if (subString(big, 0,1)!="[")
  {
    list.resize(1);
    list[0] = big;
    return;
  }

  int parenDepth = 0;
  int localCount = 0;
  std::string tmp(big);
  list.resize(0);

  // start at 1 to ignore '[';

  for (unsigned int i=1; i<big.length(); i++)
  {
    if (big[i]=='(') parenDepth++;
    if (big[i]==')') parenDepth--;
    if (big[i]==']')
    {
      tmp[localCount]='\0';
      list.push_back(tmp);
      break;
    }
    if (big[i]==',' && parenDepth==0)
    {
      tmp[localCount]='\0';
      list.push_back(tmp);
      tmp = big;
      localCount = 0;
      continue;
    }
    tmp[localCount] = big[i];
    localCount++;
  }
}


// return the position of the next whitespace in a std::string.
// If no whitespace, return -1;

int StrUtils::findNextWhitespace(const std::string& str, int offset)
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


int StrUtils::findNextNonWhitespace(const std::string& str, int offset)
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


std::string StrUtils::varTableSubstitute(const std::string& rawLine,
  const std::vector<std::string>& varNames,
  const std::vector<std::string>& varValues)
{
  TEUCHOS_TEST_FOR_EXCEPTION(varNames.size() != varValues.size(),
    std::runtime_error,
    "mismatched variable tables in varTableSubstitute");

  std::string line = rawLine;
  for (unsigned int i=0; i<varNames.size(); i++)
  {
    line = varSubstitute(line, varNames[i], varValues[i]);
  }
  return line;
}


std::string StrUtils::varSubstitute(const std::string& rawLine,
  const std::string& varName,
  const std::string& varValue)
{
  std::string line = rawLine;

  // iterate because there might be more than one occurance on this line
  while (find(line, varName) >= 0)
  {
    std::string b = before(line, varName);
    std::string a = after(line, varName);
    line = b + varValue + a;
  }
  return line;
}


std::string StrUtils::before(const std::string& str, char sub)
{
  char c[2];
  c[0] = sub;
  c[1] = 0;
  return before(str, c);
}


std::string StrUtils::before(const std::string& str, const std::string& sub)
{
  TEUCHOS_TEST_FOR_EXCEPTION(sub.c_str()==0,
    std::runtime_error, "String::before: arg is null pointer");

  char* p = std::strstr((char*) str.c_str(), (char*) sub.c_str());
  if (p==0) return str;
  int subLen = p-str.c_str();
  std::string rtn(str.c_str(), subLen);
  return rtn;
}


std::string StrUtils::after(const std::string& str, const std::string& sub)
{
  TEUCHOS_TEST_FOR_EXCEPTION(sub.c_str()==0,
    std::runtime_error, "String::after: arg is null pointer");

  // find beginning of substring
  char* p = std::strstr((char*) str.c_str(), (char*) sub.c_str()) ;
  // if substring not found, return empty std::string
  if (p==0) return std::string();
  // offset to end of substring
  p+= std::strlen(sub.c_str());
  return std::string(p);
}


int StrUtils::find(const std::string& str, const std::string& sub)
{
  char* p = std::strstr((char*) str.c_str(), (char*) sub.c_str());
  if (p==0) return -1;
  return p-str.c_str();
}


bool StrUtils::isWhite(const std::string& str)
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


std::string StrUtils::fixUnprintableCharacters(const std::string& str)
{
  std::string rtn = str;
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


std::string StrUtils::between(const std::string& str, const std::string& begin,
  const std::string& end, std::string& front,
  std::string& back)
{
  front = before(str, begin);
  std::string middle = before(after(str, begin), end);
  back = after(str, end);
  return middle;
}


std::string StrUtils::subString(const std::string& str, int begin, int end)
{
  return std::string(str.c_str()+begin, end-begin);
}


std::string StrUtils::readFromStream(std::istream& is)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "StrUtils::readFromStream isn't implemented yet");

  // NOTE (mfh 15 Sep 2014): Most compilers have figured out that the
  // return statement below is unreachable.  Some older compilers
  // might not realize this.  That's why the return statement was put
  // there, so that those compilers don't warn that this function
  // doesn't return a value.  If it's a choice between one warning and
  // another, I would prefer the choice that produces less code and
  // doesn't have unreachable code (which never gets tested).

  //return "";
}


std::string StrUtils::allCaps(const std::string& s)
{
  std::string rtn = s;
  for (unsigned int i=0; i<rtn.length(); i++)
  {
    rtn[i] = toupper(rtn[i]);
  }
  return rtn;
}


double StrUtils::atof(const std::string& s)
{
        return std::atof(s.c_str());
}


int StrUtils::atoi(const std::string& s)
{
        return std::atoi(s.c_str());
}


std::ostream& StrUtils::printLines(
  std::ostream             &os
  ,const std::string       &linePrefix
  ,const std::string       &lines
  )
{
  std::vector<std::string> linesArray = splitIntoLines(lines);
  for(unsigned int i = 0; i < linesArray.size(); ++i )
  {
    os << linePrefix << linesArray[i] << "\n";
  }
  return os;
}


std::string StrUtils::removeAllSpaces(std::string stringToClean)
{
  std::string::size_type pos=0;
  bool spacesLeft = true;

  while(spacesLeft){
    pos = stringToClean.find(" ");
    if(pos != std::string::npos){
      stringToClean.erase(pos,1);
    }
    else{
      spacesLeft = false;
    }
  }
  return stringToClean;
}


} // namespace ROL
