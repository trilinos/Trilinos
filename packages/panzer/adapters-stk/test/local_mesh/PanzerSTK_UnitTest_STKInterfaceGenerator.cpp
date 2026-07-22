// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerSTK_UnitTest_STKInterfaceGenerator.hpp"

#include "Teuchos_Assert.hpp"

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tpetra_Core.hpp"

#include "Panzer_STK_Interface.hpp"

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_LineMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

#include <string>

namespace panzer_stk
{

Teuchos::RCP<panzer_stk::STK_Interface>
generateMesh(const Teuchos::ParameterList & parameter_list_in)
{

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Setup mesh
  Teuchos::RCP<panzer_stk::STK_Interface> mesh;
  Teuchos::RCP<Teuchos::ParameterList> mesh_params;
  Teuchos::ParameterList parameter_list(parameter_list_in);

  // In case we don't have a wrapper around the mesh
  if(parameter_list.isSublist("Mesh")){
    mesh_params.reset(new Teuchos::ParameterList(parameter_list.sublist("Mesh")));
  } else if(parameter_list.name() == "Mesh"){
    // mesh_params does not pass ownership outside of this function so we can wrap parameter_list in an RCP
    mesh_params = Teuchos::rcpFromRef(parameter_list);
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "generate_mesh : Could not find 'Mesh' parameter list in input parameter list");
  }


  const std::string mesh_type = mesh_params->get<std::string>("Mesh Type");

  // Have to remove this before we send it along to the interpreter
  mesh_params->remove("Mesh Type");

  Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm,true);

  MPI_Comm raw_comm = *(mpi_comm->getRawMpiComm());

  if(mesh_type == "Hex"){
    panzer_stk::CubeHexMeshFactory meshFactory;
    meshFactory.setParameterList(mesh_params);
    mesh = meshFactory.buildMesh(raw_comm);
  } else if(mesh_type == "Tet"){
    panzer_stk::CubeTetMeshFactory meshFactory;
    meshFactory.setParameterList(mesh_params);
    mesh = meshFactory.buildMesh(raw_comm);
  } else if(mesh_type == "Quad"){
    panzer_stk::SquareQuadMeshFactory meshFactory;
    meshFactory.setParameterList(mesh_params);
    mesh = meshFactory.buildMesh(raw_comm);
  } else if(mesh_type == "Tri"){
    panzer_stk::SquareTriMeshFactory meshFactory;
    meshFactory.setParameterList(mesh_params);
    mesh = meshFactory.buildMesh(raw_comm);
  } else if(mesh_type == "Line"){
    panzer_stk::LineMeshFactory meshFactory;
    meshFactory.setParameterList(mesh_params);
    mesh = meshFactory.buildMesh(raw_comm);
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "generate_mesh : Mesh Type '"<<mesh_type<<"' not found. Options: Line, Tri, Quad, Tet, and Hex.");
  }

  return mesh;

}

}
