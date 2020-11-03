// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystParserInterface.h"

int parseCatalystFile(const std::string &filepath, CatalystParserInterface::parse_info &pinfo);

int parseCatalystString(const std::string &s, CatalystParserInterface::parse_info &pinfo);

int CatalystParserInterface::parseFile(const std::string &filepath, parse_info &pinfo)
{
  return parseCatalystFile(filepath, pinfo);
}

int CatalystParserInterface::parseString(const std::string &s, parse_info &pinfo)
{
  return parseCatalystString(s, pinfo);
}
