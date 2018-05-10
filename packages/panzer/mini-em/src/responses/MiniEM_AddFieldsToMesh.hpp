#ifndef _MiniEM_AddFieldsToMesh_hpp_
#define _MiniEM_AddFieldsToMesh_hpp_

namespace panzer_stk {
  class STK_Interface;
}
namespace Teuchos {
  class ParameterList;
}

namespace mini_em {

  void addFieldsToMesh(panzer_stk::STK_Interface & mesh,
                       const Teuchos::ParameterList & output_list);

}

#endif
