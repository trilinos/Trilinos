// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
