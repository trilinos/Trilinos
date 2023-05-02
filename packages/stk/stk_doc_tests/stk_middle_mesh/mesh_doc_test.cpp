#include "gtest/gtest.h"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/variable_size_field.hpp"

using namespace stk::middle_mesh::mesh;

//BEGIN_FIELD_ACCESS
void foo(std::shared_ptr<Mesh> mesh, FieldPtr<double> fieldPtr)
{
  auto& field = *fieldPtr;
  for (MeshEntityPtr vert : mesh->get_vertices())
    if (vert)
    {
      double val = field(vert, 0, 0);
      // do something with val
      std::cout << "field value for vert " << vert << " = " << val << std::endl;
    }
}
//END_FIELD_ACCESS


TEST(Mesh, CreateMesh)
{
  //BEGIN_CREATE_MESH
  MPI_Comm comm = MPI_COMM_WORLD;
  std::shared_ptr<Mesh> mesh = make_empty_mesh(comm);
  //END_CREATE_MESH

  //BEGIN_CREATE_VERTS
  MeshEntityPtr v1 = mesh->create_vertex(0, 0, 0);
  MeshEntityPtr v2 = mesh->create_vertex(1, 0, 0);
  MeshEntityPtr v3 = mesh->create_vertex(0, 1, 0);
  //END_CREATE_VERTS


  //BEGIN_CREATE_EDGES
  MeshEntityPtr edge1 = mesh->create_edge(v1, v2);
  MeshEntityPtr edge2 = mesh->create_edge(v2, v3);
  MeshEntityPtr edge3 = mesh->create_edge(v3, v1);
  //END_CREATE_EDGES


  //BEGIN_CREATE_TRI_FROM_EDGES
  MeshEntityPtr tri = mesh->create_triangle(edge1, edge2,  edge3, EntityOrientation::Standard);
  //END_CREATE_TRI_FROM_EDGES

  int dim = 0;
  //BEGIN_MESH_ITERATION_FUNCTIONS
  mesh->get_vertices();
  mesh->get_edges();
  mesh->get_elements();
  mesh->get_mesh_entities(dim);  // returns same as one of the above functions
                                 // depending on dim
  //END_MESH_ITERATION_FUNCTIONS

  //BEGIN_ITERATE_MESH
  for (MeshEntityPtr vert : mesh->get_vertices())
    if (vert)
    {
      // do stuff with vert
    }  
  //END_ITERATE_MESH

  int i = 0;
  //BEGIN_MESH_ENTITY_FUNCTIONS
  MeshEntityType type = edge1->get_type();   // returns enum telling if this entity is a vert, 
                                             // edge, tri, or quad
  int localId         = edge1->get_id();     //  returns the local ID of this entity
  MeshEntityPtr vert  = edge1->get_down(i);  // returns the i'th downward adjacent entity
  int numVerts        = edge1->count_down(); // returns the number of downward adjacent
                                             // entities
  MeshEntityPtr triUp = edge1->get_up(i);    // returns the i'th upward adjacent entity
  int numEls          = edge1->count_up();   // returns the number of upward adjacent entities  
  //END_MESH_ENTITY_FUNCTIONS
  EXPECT_EQ(type, MeshEntityType::Edge);
  EXPECT_EQ(localId, 0);
  EXPECT_EQ(vert->get_id(), 0);
  EXPECT_EQ(numVerts, 2);
  EXPECT_EQ(triUp->get_id(), 0);
  EXPECT_EQ(numEls, 1);

  int dimEl = 2;
  //BEGIN_ADJACENCY_FUNCTIONS
  std::array<MeshEntityPtr, MAX_DOWN> downEntities;
  int numDown = get_downward(tri, dim, downEntities.data());

  std::vector<MeshEntityPtr> upEntities;
  int numUp = get_upward(v1, dimEl, upEntities);
  //END_ADJACENCY_FUNCTIONS
  EXPECT_EQ(numDown, 3);
  EXPECT_EQ(numUp, 1);

  int remoteRank = 0, remoteId = 0;
  //BEGIN_REMOTE_ENTITY_FUNCTIONS
  // registers that the given `RemoteSharedEntity` represents the same mesh entity as this one
  v1->add_remote_shared_entity(RemoteSharedEntity{remoteRank, remoteId}); 

  // returns the i'th `RemoteSharedEntity`  
  const RemoteSharedEntity& remote = v1->get_remote_shared_entity(i);  

  // returns the number of remote shared entities  
  int numRemotes                   = v1->count_remote_shared_entities(); 
  //END_REMOTE_ENTITY_FUNCTIONS
  EXPECT_EQ(remote.remoteRank, 0);
  EXPECT_EQ(remote.remoteRank, 0);
  EXPECT_EQ(numRemotes, 1);

  {
    //BEGIN_REMOTE_ENTITY_FREE_FUNCTIONS
    int ownerRank = get_owner(mesh, v1);
    RemoteSharedEntity remote2 = get_remote_shared_entity(v1, remoteRank);
    //END_REMOTE_ENTITY_FREE_FUNCTIONS
    EXPECT_EQ(ownerRank, 0);
    EXPECT_EQ(remote2.remoteRank, remoteRank);
  }

  int componentsPerNode = 1;
  double init = 0;
  //BEGIN_FIELD_CREATION
  FieldPtr<double> field = create_field<double>(mesh, FieldShape(1, 0, 0), componentsPerNode, init);
  //END_FIELD_CREATION


  foo(mesh, field);


  //BEGIN_VARIABLE_FIELD_CREATION
  VariableSizeFieldPtr<double> variField = create_variable_size_field<double>(mesh, FieldShape(1, 0, 0));
  //END_VARIABLE_FIELD_CREATION


  int node = 0;
  double val = 0;
  {
    //BEGIN_VARIABLE_FIELD_INSERT
    variField->insert(v1, node, val);
    //END_VARIABLE_FIELD_INSERT
  }

  {
    //BEGIN_VARIABLE_FIELD_NUM_COMP
    int numComp = variField->get_num_comp(v1, node);
    //END_VARIABLE_FIELD_NUM_COMP
    EXPECT_EQ(numComp, 1);
  }

  {
    int comp = 0;
    //BEGIN_VARIABLE_FIELD_ACCESS
    auto& variFieldRef = *field;
    double val2 = variFieldRef(v1, node, comp);
    //END_VARIABLE_FIELD_ACCESS
    EXPECT_EQ(val2, 0);
  }

}

TEST(Mesh, TriFromVerts)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  std::shared_ptr<Mesh> mesh = make_empty_mesh(comm);

  MeshEntityPtr v1 = mesh->create_vertex(0, 0, 0);
  MeshEntityPtr v2 = mesh->create_vertex(1, 0, 0);
  MeshEntityPtr v3 = mesh->create_vertex(0, 1, 0);

  //BEGIN_CREATE_TRI_FROM_VERTS
  MeshEntityPtr tri = mesh->create_triangle_from_verts(v1, v2, v3);
  //END_CREATE_TRI_FROM_VERTS
  EXPECT_EQ(tri->get_id(), 0);

  {
    //BEGIN_FIELD2_CREATION
    FieldPtr<double> field = create_field<double>(mesh, FieldShape(0, 0, 3), 5);
    //END_FIELD2_CREATION
  }  
}