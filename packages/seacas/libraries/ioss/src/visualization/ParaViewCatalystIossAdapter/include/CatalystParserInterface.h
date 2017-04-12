
#ifndef __CatalystParserInterface_h
#define __CatalystParserInterface_h

#include <string>
#include <vector>
#include <map>
#include <utility>

class CatalystParserInterface
{
public:

  enum var_type { SCALAR, VECTOR, TENSOR, ALLTYPE, 
                  SCALAR_USED, VECTOR_USED, TENSOR_USED, 
                  ALLTYPE_USED };

  typedef std::map<std::string, var_type> var_map;
  typedef std::pair<unsigned int, unsigned int> id_range;

  class parse_info
  {
  public:

    parse_info()
      {
      this->node_vars = 0;
      this->element_vars = 0;
      this->global_vars = 0;
      this->nodeIDs = 0;
      this->elementIDs = 0;
      this->separator = "_";
      };

    std::string json_result;
    var_map* node_vars;
    var_map* element_vars;
    var_map* global_vars;
    id_range* nodeIDs;
    id_range* elementIDs;
    std::string separator;
  };

  static int parseFile(const std::string& filepath,
                       parse_info& pinfo);

  static int parseString(const std::string& s,
                         parse_info& pinfo);
};

#endif
