#include <string>
#include <stk_unit_test_utils/getOption.h>

namespace stk
{
namespace unit_test_util
{

inline
std::string get_mesh_spec(const std::string &optionName)
{
    std::string meshSpec("generated:");
    std::string dim = unitTestUtils::getOption(optionName, "20");
    meshSpec += dim+"x"+dim+"x"+dim;
    return meshSpec;
}

}
}
