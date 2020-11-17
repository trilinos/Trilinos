/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef __CatalystParserInterface_h
#define __CatalystParserInterface_h

#include <map>
#include <string>
#include <utility>
#include <vector>

class CatalystParserInterface
{
public:
  enum var_type {
    SCALAR,
    VECTOR,
    TENSOR,
    ALLTYPE,
    SCALAR_USED,
    VECTOR_USED,
    TENSOR_USED,
    ALLTYPE_USED
  };

  typedef std::map<std::string, var_type>       var_map;
  typedef std::pair<unsigned int, unsigned int> id_range;

  class parse_info
  {
  public:
    parse_info()
    {
      this->node_vars    = 0;
      this->element_vars = 0;
      this->global_vars  = 0;
      this->nodeIDs      = 0;
      this->elementIDs   = 0;
      this->separator    = "_";
    };

    std::string json_result;
    var_map *   node_vars;
    var_map *   element_vars;
    var_map *   global_vars;
    id_range *  nodeIDs;
    id_range *  elementIDs;
    std::string separator;
  };

  static int parseFile(const std::string &filepath, parse_info &pinfo);

  static int parseString(const std::string &s, parse_info &pinfo);
};

#endif
