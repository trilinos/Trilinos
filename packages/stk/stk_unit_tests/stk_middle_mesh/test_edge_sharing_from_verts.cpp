#include "gtest/gtest.h"
#include "stk_middle_mesh/create_sharing_from_verts.hpp"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/utils.hpp"

using namespace stk::middle_mesh;

TEST(EdgeSharingFromVerts, 2Procs)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 2)
    GTEST_SKIP();

  double delta = 0.5;
  double xoffset = delta * utils::impl::comm_rank(MPI_COMM_WORLD);
  auto mesh = mesh::make_empty_mesh();

  auto v1 = mesh->create_vertex(xoffset,         0);
  auto v2 = mesh->create_vertex(xoffset + delta, 0);
  auto v3 = mesh->create_vertex(xoffset,         delta);
  auto v4 = mesh->create_vertex(xoffset + delta, delta);
  auto v5 = mesh->create_vertex(xoffset,         2*delta);
  auto v6 = mesh->create_vertex(xoffset + delta, 2*delta);

  auto el1 = mesh->create_quad_from_verts(v1, v2, v4, v3);
  auto el2 = mesh->create_quad_from_verts(v3, v4, v6, v5);

  if (utils::impl::comm_rank(MPI_COMM_WORLD) == 0)
  {
    v2->add_remote_shared_entity({1, 0});
    v4->add_remote_shared_entity({1, 2});
    v6->add_remote_shared_entity({1, 4});
  } else
  {
    v1->add_remote_shared_entity({0, 1});
    v3->add_remote_shared_entity({0, 3});
    v5->add_remote_shared_entity({0, 5});
  }

  mesh::impl::CreateSharingFromVert sharingCreator(mesh);
  sharingCreator.create_sharing_from_verts();

  if (utils::impl::comm_rank(MPI_COMM_WORLD) == 0)
  {
    mesh::MeshEntityPtr edge1 = el1->get_down(1);
    mesh::MeshEntityPtr edge2 = el2->get_down(1);

    EXPECT_EQ(edge1->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge1->get_remote_shared_entity(0), mesh::RemoteSharedEntity(1, 3));

    EXPECT_EQ(edge2->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge2->get_remote_shared_entity(0), mesh::RemoteSharedEntity(1, 6));
  } else
  {
    mesh::MeshEntityPtr edge1 = el1->get_down(3);
    mesh::MeshEntityPtr edge2 = el2->get_down(3);

    EXPECT_EQ(edge1->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge1->get_remote_shared_entity(0), mesh::RemoteSharedEntity(0, 1));

    EXPECT_EQ(edge2->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge2->get_remote_shared_entity(0), mesh::RemoteSharedEntity(0, 4));    
  }
}


TEST(EdgeSharingFromVerts, 4Procs)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();


  auto mesh = mesh::make_empty_mesh();

  int myrank = utils::impl::comm_rank(MPI_COMM_WORLD);
  int iproc = myrank % 2, jproc = myrank / 2;
  double delta = 0.5;
  double xoffset = iproc * delta, yoffset = jproc * delta;

  auto v1 = mesh->create_vertex(xoffset,         yoffset);
  auto v2 = mesh->create_vertex(xoffset + delta, yoffset);
  auto v3 = mesh->create_vertex(xoffset,         yoffset + delta);
  auto v4 = mesh->create_vertex(xoffset + delta, yoffset + delta);
  auto el1 = mesh->create_quad_from_verts(v1, v2, v4, v3);

  if (myrank == 0)
  {
    v2->add_remote_shared_entity({1, 0});

    v4->add_remote_shared_entity({1, 2});
    v4->add_remote_shared_entity({2, 1});
    v4->add_remote_shared_entity({3, 0});

    v3->add_remote_shared_entity({2, 0});
  } else if (myrank == 1)
  {
    v1->add_remote_shared_entity({0, 1});

    v3->add_remote_shared_entity({0, 3});
    v3->add_remote_shared_entity({2, 1});
    v3->add_remote_shared_entity({3, 0});

    v4->add_remote_shared_entity({3, 1});
  } else if (myrank == 2)
  {
    v1->add_remote_shared_entity({0, 2});

    v2->add_remote_shared_entity({0, 3});
    v2->add_remote_shared_entity({1, 2});
    v2->add_remote_shared_entity({3, 0});

    v4->add_remote_shared_entity({3, 2});
  } else if (myrank == 3)
  {
    v2->add_remote_shared_entity({1, 3});

    v1->add_remote_shared_entity({0, 3});
    v1->add_remote_shared_entity({1, 2});
    v1->add_remote_shared_entity({2, 1});

    v3->add_remote_shared_entity({2, 3});
  }

  mesh::impl::CreateSharingFromVert sharingCreator(mesh);
  sharingCreator.create_sharing_from_verts();

  if (myrank == 0)
  {
    mesh::MeshEntityPtr edge1 = el1->get_down(1);
    mesh::MeshEntityPtr edge2 = el1->get_down(2);

    EXPECT_EQ(edge1->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge1->get_remote_shared_entity(0), mesh::RemoteSharedEntity(1, 3));

    EXPECT_EQ(edge2->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge2->get_remote_shared_entity(0), mesh::RemoteSharedEntity(2, 0));
  } else if (myrank == 1)
  {
    mesh::MeshEntityPtr edge1 = el1->get_down(3);
    mesh::MeshEntityPtr edge2 = el1->get_down(2);

    EXPECT_EQ(edge1->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge1->get_remote_shared_entity(0), mesh::RemoteSharedEntity(0, 1));

    EXPECT_EQ(edge2->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge2->get_remote_shared_entity(0), mesh::RemoteSharedEntity(3, 0));    
  } else if (myrank == 2)
  {
    mesh::MeshEntityPtr edge1 = el1->get_down(0);
    mesh::MeshEntityPtr edge2 = el1->get_down(1);

    EXPECT_EQ(edge1->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge1->get_remote_shared_entity(0), mesh::RemoteSharedEntity(0, 2));

    EXPECT_EQ(edge2->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge2->get_remote_shared_entity(0), mesh::RemoteSharedEntity(3, 3));    
  } else if (myrank == 3)
  {
    mesh::MeshEntityPtr edge1 = el1->get_down(0);
    mesh::MeshEntityPtr edge2 = el1->get_down(3);

    EXPECT_EQ(edge1->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge1->get_remote_shared_entity(0), mesh::RemoteSharedEntity(1, 2));

    EXPECT_EQ(edge2->count_remote_shared_entities(), 1);
    EXPECT_EQ(edge2->get_remote_shared_entity(0), mesh::RemoteSharedEntity(2, 1));      
  }
}

TEST(ComputingSharingFromVerts, Elements)
{
  int myrank   = utils::impl::comm_rank(MPI_COMM_WORLD);
  int commsize = utils::impl::comm_size(MPI_COMM_WORLD);
  auto mesh    = mesh::make_empty_mesh();

  auto v1 = mesh->create_vertex(0, 0);
  auto v2 = mesh->create_vertex(1, 0);
  auto v3 = mesh->create_vertex(1, 1);
  auto v4 = mesh->create_vertex(0, 0);
  std::array<mesh::MeshEntityPtr, 4> verts = {v1, v2, v3, v4};

  auto el1 = mesh->create_quad_from_verts(v1, v2, v3, v4);

  for (int rank=0; rank < commsize; ++rank)
  {
    if (rank != myrank)
    {
      for (int i=0; i < 4; ++i)
      {
        verts[i]->add_remote_shared_entity({rank, i});
      }
    }
  }

  mesh::impl::CreateSharingFromVert sharingCreator(mesh, 2);
  sharingCreator.create_sharing_from_verts();

  EXPECT_EQ(el1->count_remote_shared_entities(), commsize - 1);
  for (int rank=0; rank < commsize; ++rank)
  {
    if (rank != myrank)
    {
      mesh::RemoteSharedEntity remote = mesh::get_remote_shared_entity(el1, rank);
      EXPECT_EQ(remote, mesh::RemoteSharedEntity(rank, 0));
    }
  }

}