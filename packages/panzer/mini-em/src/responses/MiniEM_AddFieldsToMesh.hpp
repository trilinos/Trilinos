#ifndef _MiniEM_AddFieldsToMesh_hpp_
#define _MiniEM_AddFieldsToMesh_hpp_

#include "Teuchos_ParameterList.hpp"
#include "Panzer_STK_Interface.hpp"

namespace mini_em {

  void addFieldsToMesh(panzer_stk::STK_Interface & mesh,
                       const Teuchos::ParameterList & output_list);

}

#endif
