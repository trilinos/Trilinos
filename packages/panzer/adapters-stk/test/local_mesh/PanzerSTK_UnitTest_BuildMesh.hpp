#ifndef PANZER_STK_UTILITIES_BUILD_MESH_HPP_
#define PANZER_STK_UTILITIES_BUILD_MESH_HPP_

#include "Teuchos_RCP.hpp"

#include <vector>

namespace Teuchos
{
class ParameterList;
}

namespace panzer_stk
{

class STK_Interface;

Teuchos::RCP<panzer_stk::STK_Interface>
buildMesh(const std::vector<int> & N,
          const std::vector<int> & B,
          const std::vector<double> & L);

}

#endif
