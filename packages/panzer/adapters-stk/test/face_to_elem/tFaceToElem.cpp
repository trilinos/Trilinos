// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ParameterList.hpp>

#include "PanzerCore_config.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"
#include "Panzer_Relations.hpp"

#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopology.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Import.hpp>

#include "SillyEquationSet.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

typedef Traits::Residual ResidualType;
typedef Traits::Jacobian JacobianType;

namespace unit_test {

void setParamList(RCP<Teuchos::ParameterList>  mesh_pl, int nx, int ny, int nz=0) {

  mesh_pl->set("X Blocks", 1);
  mesh_pl->set("Y Blocks", 1);
  mesh_pl->set("X Elements", nx);
  mesh_pl->set("Y Elements", ny);

  mesh_pl->set("X0", 0.0);
  mesh_pl->set("Y0", 0.0);
  mesh_pl->set("Xf", 1.0);
  mesh_pl->set("Yf", 1.0);

  if ( nz > 0 ) {
    mesh_pl->set("Z Blocks", 1);
    mesh_pl->set("Z Elements", nz);
    mesh_pl->set("Z0", 0.0);
    mesh_pl->set("Zf", 1.0);
  }
}
Teuchos::RCP<panzer::PhysicsBlock> createPhysicsBlock(RCP<panzer_stk::STK_Interface> mesh) {
  std::string eblock = "eblock-0_0_0";
  if ( !mesh->validBlockId(eblock)) // Try 2D
    eblock = "eblock-0_0";
  int workset_size=100;
  panzer::CellData cellData(workset_size,mesh->getCellTopology(eblock));
  Teuchos::RCP<MyFactory> eqset_factory = Teuchos::rcp(new MyFactory);

  int default_int_order =1;
  Teuchos::RCP<Teuchos::ParameterList> ipb = Teuchos::parameterList();
  ipb->setName("test physics");
  {
    Teuchos::ParameterList& p = ipb->sublist("a");
    p.set("Type","MeshCoords");
    p.set("Prefix","");
    p.set("Model ID","solid");
    p.set("Basis Type","HGrad");
    p.set("Basis Order",1);
    p.set("Integration Order",1);
  }

  Teuchos::RCP<panzer::GlobalData> gd = panzer::createGlobalData();
  Teuchos::RCP<panzer::PhysicsBlock> physicsBlock =
    Teuchos::rcp(new PhysicsBlock(ipb,eblock,default_int_order,cellData,eqset_factory,gd,false));
  return physicsBlock;
}

TEUCHOS_UNIT_TEST(tFullRelations, test_face_tet)
{
  RCP<Teuchos::ParameterList>  mesh_pl(new Teuchos::ParameterList("Mesh"));
  int nx=4, ny=3, nz=2;
  setParamList(mesh_pl, nx, ny, nz);

  RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::CubeTetMeshFactory);
  mesh_factory->setParameterList(mesh_pl);
  RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  const Teuchos::RCP<panzer::ConnManager>
    conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  FaceToElems rel(conn_manager);

  int num_faces_correct = 4*(nx*ny+nx*nz+ny*nz);
  int num_faces_local = rel.numberBoundaryFaces(), num_faces_global;
  MPI_Allreduce(&num_faces_local, &num_faces_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  TEST_EQUALITY(num_faces_correct, num_faces_global);
}

TEUCHOS_UNIT_TEST(tFullRelations, test_face_hex)
{
  RCP<Teuchos::ParameterList>  mesh_pl(new Teuchos::ParameterList("Mesh"));
  int nx=4, ny=3, nz=2;
  setParamList(mesh_pl, nx, ny, nz);

  RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::CubeHexMeshFactory);
  mesh_factory->setParameterList(mesh_pl);
  RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  const Teuchos::RCP<panzer::ConnManager>
    conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  FaceToElems rel(conn_manager);

  int num_faces_correct = 2*(nx*ny+nx*nz+ny*nz);
  int num_faces_local = rel.numberBoundaryFaces(), num_faces_global;
  MPI_Allreduce(&num_faces_local, &num_faces_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  TEST_EQUALITY(num_faces_correct, num_faces_global);
}

TEUCHOS_UNIT_TEST(tFullRelations, test_face_quad)
{
  RCP<Teuchos::ParameterList>  mesh_pl(new Teuchos::ParameterList("Mesh"));
  int nx=4, ny=3, nz=-1;
  setParamList(mesh_pl, nx, ny, nz);

  RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::SquareQuadMeshFactory);
  mesh_factory->setParameterList(mesh_pl);
  RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  const Teuchos::RCP<panzer::ConnManager>
    conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  FaceToElems rel(conn_manager);

  int num_faces_correct = 2*(nx+ny);
  int num_faces_local = rel.numberBoundaryFaces(), num_faces_global;
  MPI_Allreduce(&num_faces_local, &num_faces_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  TEST_EQUALITY(num_faces_correct, num_faces_global);
}

TEUCHOS_UNIT_TEST(tFullRelations, test_face_tri)
{
  RCP<Teuchos::ParameterList>  mesh_pl(new Teuchos::ParameterList("Mesh"));
  int nx=4, ny=3, nz=-1;
  setParamList(mesh_pl, nx, ny, nz);

  RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::SquareTriMeshFactory);
  mesh_factory->setParameterList(mesh_pl);
  RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);
  const Teuchos::RCP<panzer::ConnManager>
    conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  FaceToElems rel(conn_manager);

  int num_faces_correct = 2*(nx+ny);
  int num_faces_local = rel.numberBoundaryFaces(), num_faces_global;
  MPI_Allreduce(&num_faces_local, &num_faces_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  TEST_EQUALITY(num_faces_correct, num_faces_global);
}

TEUCHOS_UNIT_TEST(tFullRelations, test_build_workset)
{
  int nx=3, ny=2, nz=4;
  RCP<Teuchos::ParameterList>  mesh_pl(new Teuchos::ParameterList("Mesh"));

  setParamList(mesh_pl, nx, ny, nz);
  RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::CubeTetMeshFactory);
  mesh_factory->setParameterList(mesh_pl);
  RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

  const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  FaceToElems rel(conn_manager);

  Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = createPhysicsBlock(mesh);
  Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
                                                                                          physicsBlock->getWorksetNeeds());                    
}

TEUCHOS_UNIT_TEST(tFullRelations, test_tet_normals)
{
  int nx=3, ny=2, nz=4;
  RCP<Teuchos::ParameterList>  mesh_pl(new Teuchos::ParameterList("Mesh"));

  setParamList(mesh_pl, nx, ny, nz);
  RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::CubeTetMeshFactory);
  mesh_factory->setParameterList(mesh_pl);
  RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

  const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  FaceToElems rel(conn_manager);

  Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = createPhysicsBlock(mesh);
  Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
                                                                                          physicsBlock->getWorksetNeeds());                    

  rel.setNormals(work_sets);

  std::vector<double> norm;
  double nxny2=sqrt(nx*nx+ny*ny), nxnz2=sqrt(nx*nx+nz*nz),  nynz2=sqrt(ny*ny+nz*nz);

  // face 0 is in the y/z plane
  rel.getNormal(4, 0, norm);
  TEST_FLOATING_EQUALITY( norm[1], -ny/nynz2, 1e-14);
  TEST_FLOATING_EQUALITY( norm[2], +nz/nynz2, 1e-14);

  // face 1 is in the x/z plane
  rel.getNormal(4, 1, norm);
  TEST_FLOATING_EQUALITY( norm[0], -nx/nxnz2, 1e-14);
  TEST_FLOATING_EQUALITY( norm[2], -nz/nxnz2, 1e-14);

  // face 2 is in the x/y plane
  rel.getNormal(4, 2, norm);
  TEST_FLOATING_EQUALITY( norm[0], -nx/nxny2, 1e-14);
  TEST_FLOATING_EQUALITY( norm[1], +ny/nxny2, 1e-14);
  
  // and the last face is pointing to the right.
  rel.getNormal(4, 3, norm);
  TEST_FLOATING_EQUALITY( norm[0], +1.0, 1e-14);

}
TEUCHOS_UNIT_TEST(tFullRelations, test_tri_normals)
{
  int nx=3, ny=2, nz=-1;
  RCP<Teuchos::ParameterList>  mesh_pl(new Teuchos::ParameterList("Mesh"));

  setParamList(mesh_pl, nx, ny, nz);
  RCP<panzer_stk::STK_MeshFactory> mesh_factory = rcp(new panzer_stk::SquareTriMeshFactory);
  mesh_factory->setParameterList(mesh_pl);
  RCP<panzer_stk::STK_Interface> mesh = mesh_factory->buildUncommitedMesh(MPI_COMM_WORLD);
  mesh_factory->completeMeshConstruction(*mesh,MPI_COMM_WORLD);

  const Teuchos::RCP<panzer::ConnManager>
      conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  FaceToElems rel(conn_manager);

  Teuchos::RCP<panzer::PhysicsBlock> physicsBlock = createPhysicsBlock(mesh);
  Teuchos::RCP<std::vector<panzer::Workset> > work_sets = panzer_stk::buildWorksets(*mesh,physicsBlock->elementBlockID(),
                                                                                          physicsBlock->getWorksetNeeds());                    

  rel.setNormals(work_sets);
  std::vector<double> norm;
  double nxny2=sqrt(nx*nx+ny*ny);

  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &rank);
  if (rank > 0) return;

  // face 0 is pointing up.
  rel.getNormal(4, 0, norm);
  TEST_FLOATING_EQUALITY( norm[1], 1.0, 1e-14);

  // face 1 is pointing to the left
  rel.getNormal(4, 1, norm);
  TEST_FLOATING_EQUALITY( norm[0], -1.0, 1e-14);

  // face 2 is in the x/y plane
  rel.getNormal(4, 2, norm);
  TEST_FLOATING_EQUALITY( norm[0], +nx/nxny2, 1e-14);
  TEST_FLOATING_EQUALITY( norm[1], -ny/nxny2, 1e-14);

}

} // end unit test
} // end panzer
