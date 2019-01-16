#ifndef PANZER_STK_STK_INTERFACE_GENERATOR_HPP_
#define PANZER_STK_STK_INTERFACE_GENERATOR_HPP_

#include "Teuchos_RCP.hpp"

namespace Teuchos
{
class ParameterList;
}

namespace panzer_stk
{

class STK_Interface;

/**\brief Generalized mesh interface
 *
 * Generates a mesh
 *
 * \param[in] parameter_list The mesh description should be in this
 *
 * \return A mesh
 */
Teuchos::RCP<panzer_stk::STK_Interface>
generateMesh(const Teuchos::ParameterList & parameter_list);

}

#endif
