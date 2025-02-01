// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __PhactoriParserInterface_h
#define __PhactoriParserInterface_h

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Iovs {

  class PhactoriParserInterface
  {

  public:
    enum VarType {
      SCALAR,
      VECTOR,
      TENSOR,
      ALLTYPE,
      SCALAR_USED,
      VECTOR_USED,
      TENSOR_USED,
      ALLTYPE_USED
    };

    typedef std::map<std::string, VarType>        varMap;
    typedef std::pair<unsigned int, unsigned int> idRange;

    struct ParseInfo
    {
      varMap      nodeVars;
      varMap      elementVars;
      varMap      globalVars;
      bool        checkNodeIDRange = false;
      idRange     nodeIDRange;
      bool        checkElementIDRange = false;
      idRange     elementIDRange;
      std::string separator       = "_";
      std::string jsonParseResult = "";
      bool        parseFailed     = true;
    };

    static void parseFile(const std::string &filepath, ParseInfo &pinfo);

    static void parseString(const std::string &s, ParseInfo &pinfo);
  };

} // namespace Iovs

#endif
