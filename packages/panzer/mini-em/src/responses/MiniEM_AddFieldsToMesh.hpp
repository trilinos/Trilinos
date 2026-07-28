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

  /** \brief Registers the cell/nodal fields requested by an output
    * ParameterList with the mesh so they can be written to Exodus.
    *
    * \param[in,out] mesh the mesh to register fields with.
    * \param[in] output_list a ParameterList with (a subset of) the sublists "Cell Average Quantities", "Cell Average Vectors", "Cell Quantities", "Nodal Quantities", and "Allocate Nodal Quantities", each mapping an element block ID to a comma-separated list of field names.
    */
  void addFieldsToMesh(panzer_stk::STK_Interface & mesh,
                       const Teuchos::ParameterList & output_list);

}

#endif
