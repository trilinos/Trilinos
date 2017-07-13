#include "CatalystParserInterface.h"

int parseCatalystFile(const std::string& filepath,
                      CatalystParserInterface::parse_info& pinfo);

int parseCatalystString(const std::string& s,
                        CatalystParserInterface::parse_info& pinfo);

int CatalystParserInterface::parseFile(const std::string& filepath,
                                       parse_info& pinfo)
{
  return parseCatalystFile(filepath, 
                           pinfo);
}

int CatalystParserInterface::parseString(const std::string& s,
                                         parse_info& pinfo)
{
  return parseCatalystString(s, 
                             pinfo);
}
