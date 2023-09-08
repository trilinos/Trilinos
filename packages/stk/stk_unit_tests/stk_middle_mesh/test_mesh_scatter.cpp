#include "gtest/gtest.h"
#include "stk_middle_mesh/mesh.hpp"
#include "stk_middle_mesh/mesh_entity.hpp"
#include "stk_middle_mesh/field.hpp"
#include "stk_middle_mesh/mesh_scatter.hpp"
#include "stk_middle_mesh/variable_size_field.hpp"

using namespace stk::middle_mesh;


namespace {

void expect_near(const utils::Point& pt1, const utils::Point& pt2, double tol)
{
  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(pt1[i], pt2[i], tol);
}

mesh::MeshEntityPtr get_closest_entity(std::shared_ptr<mesh::Mesh> mesh, int dim, const utils::Point& centroid, double tol=1e-12)
{
  double minDist = tol;
  mesh::MeshEntityPtr minEntity = nullptr;
  for (auto& entity : mesh->get_mesh_entities(dim))
    if (entity)
    {
      utils::Point disp = mesh::compute_centroid(entity) - centroid;
      double dist = std::sqrt(dot(disp, disp));
      if (dist < minDist)
      {
        minDist = dist;
        minEntity = entity;
      }
    }

  return minEntity;
}

mesh::RemoteSharedEntity get_remote(mesh::VariableSizeFieldPtr<mesh::RemoteSharedEntity> fieldPtr, mesh::MeshEntityPtr entity, int rank)
{
  auto& field = *fieldPtr;
  for (int i=0; i < field.get_num_comp(entity, 0); ++i)
    if (field(entity, 0, i).remoteRank == rank)
      return field(entity, 0, i);

  throw std::runtime_error("unable to find remote");
}


class RemoteChecker
{
  struct RemoteData
  {
    mesh::RemoteSharedEntity remote = {-1, -1};
    utils::Point centroid = {0, 0, 0};
  };

  using Exchanger = stk::DataExchangeKnownPatternNonBlockingBuffer<RemoteData>;
  public:
    void check_remotes(std::shared_ptr<mesh::Mesh> mesh)
    {
      for (int dim=0; dim < 2; ++dim)
      {
        Exchanger exchanger(mesh->get_comm());
        pack_buffers(exchanger, mesh, dim);

        exchanger.start_nonblocking();
        auto f = [](int rank, const std::vector<RemoteData>& buf) {};
        exchanger.complete_receives(f);

        check_remotes(exchanger, mesh, dim);
      }
    }

  private:
    void pack_buffers(Exchanger& exchanger, std::shared_ptr<mesh::Mesh> mesh, int dim)
    {
      for (mesh::MeshEntityPtr entity : mesh->get_mesh_entities(dim))
        if (entity)
        {
          utils::Point centroid = mesh::compute_centroid(entity);
          for (int i=0; i < entity->count_remote_shared_entities(); ++i)
          {
            mesh::RemoteSharedEntity remote = entity->get_remote_shared_entity(i);
            RemoteData data{remote, centroid};
            exchanger.get_send_buf(remote.remoteRank).push_back(data);
            exchanger.get_recv_buf(remote.remoteRank).emplace_back();
          }
        }
    }

    void check_remotes(Exchanger& exchanger, std::shared_ptr<mesh::Mesh> mesh, int dim)
    {
      mesh::FieldShape fshape;
      fshape.count[dim] = 1;
      auto fieldPtr = mesh::create_field<int>(mesh, fshape, 1, 0);
      auto& recvCounts = *fieldPtr;

      for (int rank=0; rank < utils::impl::comm_size(mesh->get_comm()); ++rank)
      {
        std::set<int> seen_entities;
        for (RemoteData& data : exchanger.get_recv_buf(rank))
        {
          EXPECT_EQ(seen_entities.count(data.remote.remoteId), 0u);
          mesh::MeshEntityPtr entity = mesh->get_mesh_entities(dim)[data.remote.remoteId];
          utils::Point entityCentroid = mesh::compute_centroid(entity);

          auto disp = entityCentroid - data.centroid;
          double dist = std::sqrt(dot(disp, disp));
          EXPECT_NEAR(dist, 0, 1e-13);
          recvCounts(entity, 0, 0) += 1;
          seen_entities.insert(entity->get_id());
        }
      }

      for (mesh::MeshEntityPtr entity : mesh->get_mesh_entities(dim))
      {
        EXPECT_EQ(recvCounts(entity, 0, 0), entity->count_remote_shared_entities());
      }
    }
};

void check_remotes_centroids(std::shared_ptr<mesh::Mesh> mesh)
{
  RemoteChecker checker;
  checker.check_remotes(mesh);
}

}

// Tests to write
// 2 x 2 mesh, scatter to 4 procs
//  - test coordinates of vertices of element
//  - check remote info for verts
//  - run consistency checks

TEST(MeshScatter, 2x2FromOneProc)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 5)
    GTEST_SKIP();

  std::shared_ptr<mesh::Mesh> mesh;
  std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, mesh);
  MPI_Comm unionComm = MPI_COMM_WORLD;
  MPI_Comm meshComm;
  int color = utils::impl::comm_rank(unionComm) == 0 ? 0 : 1;
  MPI_Comm_split(unionComm, color, 0, &meshComm);

  if (color == 0)
  {
    mesh = mesh::make_empty_mesh(meshComm);
    auto v1 = mesh->create_vertex(0,   0,   0);
    auto v2 = mesh->create_vertex(0.5, 0,   0);
    auto v3 = mesh->create_vertex(1,   0,   0);
    auto v4 = mesh->create_vertex(0,   0.5, 0);
    auto v5 = mesh->create_vertex(0.5, 0.5, 0);
    auto v6 = mesh->create_vertex(1,   0.5, 0);
    auto v7 = mesh->create_vertex(0,   1,   0);
    auto v8 = mesh->create_vertex(0.5, 1,   0);
    auto v9 = mesh->create_vertex(1,   1,   0);
    
    auto el1 = mesh->create_quad_from_verts(v1, v2, v5, v4);
    auto el2 = mesh->create_quad_from_verts(v2, v3, v6, v5);
    auto el3 = mesh->create_quad_from_verts(v4, v5, v8, v7);  
    auto el4 = mesh->create_quad_from_verts(v5, v6, v9, v8);

    scatterSpec->add_destination(el1, 1);
    scatterSpec->add_destination(el2, 2);
    scatterSpec->add_destination(el3, 3);
    scatterSpec->add_destination(el4, 4);
  }

  mesh::impl::MeshScatter scatter(scatterSpec, mesh, color == 1 ? meshComm : MPI_COMM_NULL, true);
  auto meshScattered  = scatter.scatter();
  auto elementOrigins = scatter.get_element_origins();
  auto entityOrigins  = scatter.get_entity_origins();
  auto entityDests    = scatter.get_entity_destinations();

  if (color == 0)
  {
    EXPECT_EQ(entityOrigins, nullptr);
    EXPECT_NE(entityDests, nullptr);

    auto v1 = get_closest_entity(mesh, 0, {0,   0,   0});
    auto v2 = get_closest_entity(mesh, 0, {0.5, 0,   0});
    auto v3 = get_closest_entity(mesh, 0, {1.0, 0,   0});
    auto v4 = get_closest_entity(mesh, 0, {0,   0.5, 0});
    auto v5 = get_closest_entity(mesh, 0, {0.5, 0.5, 0});
    auto v6 = get_closest_entity(mesh, 0, {1.0, 0.5, 0});
    auto v7 = get_closest_entity(mesh, 0, {0,   1.0, 0});
    auto v8 = get_closest_entity(mesh, 0, {0.5, 1.0, 0});
    auto v9 = get_closest_entity(mesh, 0, {1.0, 1.0, 0});

    auto edge1  = get_closest_entity(mesh, 1, {0.25,   0,   0});
    auto edge2  = get_closest_entity(mesh, 1, {0.50, 0.25,   0});
    auto edge3  = get_closest_entity(mesh, 1, {0.25, 0.50,   0});
    auto edge4  = get_closest_entity(mesh, 1, {0.0,  0.25,   0});

    auto edge5  = get_closest_entity(mesh, 1, {0.75,   0,   0});
    auto edge6  = get_closest_entity(mesh, 1, {1.0,  0.25,   0});
    auto edge7  = get_closest_entity(mesh, 1, {0.75, 0.50,   0});

    auto edge8  = get_closest_entity(mesh, 1, {0.50, 0.75,   0});
    auto edge9  = get_closest_entity(mesh, 1, {0.25, 1.0,   0});
    auto edge10 = get_closest_entity(mesh, 1, {0,    0.75,   0});

    auto edge11 = get_closest_entity(mesh, 1, {1.0,  0.75,   0});
    auto edge12 = get_closest_entity(mesh, 1, {0.75, 1.0,   0});

    auto el1 = get_closest_entity(mesh, 2, {0.25, 0.25, 0});
    auto el2 = get_closest_entity(mesh, 2, {0.75, 0.25, 0});
    auto el3 = get_closest_entity(mesh, 2, {0.25, 0.75, 0});
    auto el4 = get_closest_entity(mesh, 2, {0.75, 0.75, 0});


    EXPECT_EQ(entityDests->get_num_comp(v1, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(v2, 0), 2);
    EXPECT_EQ(entityDests->get_num_comp(v3, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(v4, 0), 2);
    EXPECT_EQ(entityDests->get_num_comp(v5, 0), 4);
    EXPECT_EQ(entityDests->get_num_comp(v6, 0), 2);
    EXPECT_EQ(entityDests->get_num_comp(v7, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(v8, 0), 2);
    EXPECT_EQ(entityDests->get_num_comp(v9, 0), 1);

    EXPECT_EQ(get_remote(entityDests, v1, 1), mesh::RemoteSharedEntity(1, 0));

    EXPECT_EQ(get_remote(entityDests, v2, 1), mesh::RemoteSharedEntity(1, 1));
    EXPECT_EQ(get_remote(entityDests, v2, 2), mesh::RemoteSharedEntity(2, 0));

    EXPECT_EQ(get_remote(entityDests, v3, 2), mesh::RemoteSharedEntity(2, 1));

    EXPECT_EQ(get_remote(entityDests, v4, 1), mesh::RemoteSharedEntity(1, 2));
    EXPECT_EQ(get_remote(entityDests, v4, 3), mesh::RemoteSharedEntity(3, 0));

    EXPECT_EQ(get_remote(entityDests, v5, 1), mesh::RemoteSharedEntity(1, 3));
    EXPECT_EQ(get_remote(entityDests, v5, 2), mesh::RemoteSharedEntity(2, 2));
    EXPECT_EQ(get_remote(entityDests, v5, 3), mesh::RemoteSharedEntity(3, 1));
    EXPECT_EQ(get_remote(entityDests, v5, 4), mesh::RemoteSharedEntity(4, 0));

    EXPECT_EQ(get_remote(entityDests, v6, 2), mesh::RemoteSharedEntity(2, 3));
    EXPECT_EQ(get_remote(entityDests, v6, 4), mesh::RemoteSharedEntity(4, 1));

    EXPECT_EQ(get_remote(entityDests, v7, 3), mesh::RemoteSharedEntity(3, 2));

    EXPECT_EQ(get_remote(entityDests, v8, 3), mesh::RemoteSharedEntity(3, 3));
    EXPECT_EQ(get_remote(entityDests, v8, 4), mesh::RemoteSharedEntity(4, 2));

    EXPECT_EQ(get_remote(entityDests, v9, 4), mesh::RemoteSharedEntity(4, 3));

    EXPECT_EQ(entityDests->get_num_comp(edge1, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(edge2, 0), 2);
    EXPECT_EQ(entityDests->get_num_comp(edge3, 0), 2);
    EXPECT_EQ(entityDests->get_num_comp(edge4, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(edge5, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(edge6, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(edge7, 0), 2);
    EXPECT_EQ(entityDests->get_num_comp(edge8, 0), 2);
    EXPECT_EQ(entityDests->get_num_comp(edge9, 0), 1);    
    EXPECT_EQ(entityDests->get_num_comp(edge10, 0), 1);    
    EXPECT_EQ(entityDests->get_num_comp(edge11, 0), 1);    
    EXPECT_EQ(entityDests->get_num_comp(edge12, 0), 1);

    
    EXPECT_EQ(get_remote(entityDests, edge1, 1), mesh::RemoteSharedEntity(1, 0));

    EXPECT_EQ(get_remote(entityDests, edge2, 1), mesh::RemoteSharedEntity(1, 1));
    EXPECT_EQ(get_remote(entityDests, edge2, 2), mesh::RemoteSharedEntity(2, 0));
    
    EXPECT_EQ(get_remote(entityDests, edge3, 1), mesh::RemoteSharedEntity(1, 2));
    EXPECT_EQ(get_remote(entityDests, edge3, 3), mesh::RemoteSharedEntity(3, 0));

    EXPECT_EQ(get_remote(entityDests, edge4, 1), mesh::RemoteSharedEntity(1, 3));

    EXPECT_EQ(get_remote(entityDests, edge5, 2), mesh::RemoteSharedEntity(2, 1));

    EXPECT_EQ(get_remote(entityDests, edge6, 2), mesh::RemoteSharedEntity(2, 2));

    EXPECT_EQ(get_remote(entityDests, edge7, 2), mesh::RemoteSharedEntity(2, 3));
    EXPECT_EQ(get_remote(entityDests, edge7, 4), mesh::RemoteSharedEntity(4, 0));

    EXPECT_EQ(get_remote(entityDests, edge8, 3), mesh::RemoteSharedEntity(3, 1));
    EXPECT_EQ(get_remote(entityDests, edge8, 4), mesh::RemoteSharedEntity(4, 1));

    EXPECT_EQ(get_remote(entityDests, edge9, 3), mesh::RemoteSharedEntity(3, 2));

    EXPECT_EQ(get_remote(entityDests, edge10, 3), mesh::RemoteSharedEntity(3, 3));

    EXPECT_EQ(get_remote(entityDests, edge11, 4), mesh::RemoteSharedEntity(4, 2));

    EXPECT_EQ(get_remote(entityDests, edge12, 4), mesh::RemoteSharedEntity(4, 3));

    EXPECT_EQ(entityDests->get_num_comp(el1, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(el2, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(el3, 0), 1);
    EXPECT_EQ(entityDests->get_num_comp(el4, 0), 1);

    EXPECT_EQ(get_remote(entityDests, el1, 1), mesh::RemoteSharedEntity(1, 0));
    EXPECT_EQ(get_remote(entityDests, el2, 2), mesh::RemoteSharedEntity(2, 0));
    EXPECT_EQ(get_remote(entityDests, el3, 3), mesh::RemoteSharedEntity(3, 0));
    EXPECT_EQ(get_remote(entityDests, el4, 4), mesh::RemoteSharedEntity(4, 0));

    
  } else if (color == 1)
  {
    EXPECT_EQ(mesh::count_valid(meshScattered->get_vertices()), 4);
    EXPECT_EQ(mesh::count_valid(meshScattered->get_edges()), 4);
    EXPECT_EQ(mesh::count_valid(meshScattered->get_elements()), 1);

    int myRank = utils::impl::comm_rank(meshComm);
    utils::Point lowerLeftCorner;
    switch (myRank)
    {
      case 0: {lowerLeftCorner = {0,   0,   0}; break; }
      case 1: {lowerLeftCorner = {0.5, 0,   0}; break; }
      case 2: {lowerLeftCorner = {0,   0.5, 0}; break; }
      case 3: {lowerLeftCorner = {0.5, 0.5, 0}; break; }
    }

    auto el = meshScattered->get_elements()[0];
    std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
    int nverts = mesh::get_downward(el, 0, verts.data());
    EXPECT_EQ(nverts, 4);
    expect_near(verts[0]->get_point_orig(0), lowerLeftCorner, 1e-13);
    expect_near(verts[1]->get_point_orig(0), lowerLeftCorner + utils::Point(0.5, 0,   0), 1e-13);
    expect_near(verts[2]->get_point_orig(0), lowerLeftCorner + utils::Point(0.5, 0.5, 0), 1e-13);
    expect_near(verts[3]->get_point_orig(0), lowerLeftCorner + utils::Point(0,   0.5, 0), 1e-13);

    if (myRank == 0)
    {
      EXPECT_EQ(verts[0]->count_remote_shared_entities(), 0);
      EXPECT_EQ(verts[1]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[2]->count_remote_shared_entities(), 3);
      EXPECT_EQ(verts[3]->count_remote_shared_entities(), 1);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 3));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 2));      

      EXPECT_EQ(el->get_down(0)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(1)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(2)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(3)->count_remote_shared_entities(), 0);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(1), 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(2), 2));

      EXPECT_EQ(entityOrigins->get_num_comp(verts[0], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[1], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[2], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[3], 0), 1);

      EXPECT_EQ(get_remote(entityOrigins, verts[0], 0), mesh::RemoteSharedEntity(0, 0));
      EXPECT_EQ(get_remote(entityOrigins, verts[1], 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ(get_remote(entityOrigins, verts[2], 0), mesh::RemoteSharedEntity(0, 4));
      EXPECT_EQ(get_remote(entityOrigins, verts[3], 0), mesh::RemoteSharedEntity(0, 3));

      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(0), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(1), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(2), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(3), 0), 1);

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(0), 0), mesh::RemoteSharedEntity(0, 0));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(2), 0), mesh::RemoteSharedEntity(0, 2));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 0), mesh::RemoteSharedEntity(0, 3));

      EXPECT_EQ(entityOrigins->get_num_comp(el, 0), 1);
      EXPECT_EQ(get_remote(entityOrigins, el, 0), mesh::RemoteSharedEntity(0, 0));

    } else if (myRank == 1)
    {
      EXPECT_EQ(verts[0]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[1]->count_remote_shared_entities(), 0);
      EXPECT_EQ(verts[2]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[3]->count_remote_shared_entities(), 3);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 3));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 3));

      EXPECT_EQ(el->get_down(0)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(1)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(2)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(3)->count_remote_shared_entities(), 1);   

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(2), 3));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(3), 0));

      EXPECT_EQ(entityOrigins->get_num_comp(verts[0], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[1], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[2], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[3], 0), 1);

      EXPECT_EQ(get_remote(entityOrigins, verts[0], 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ(get_remote(entityOrigins, verts[1], 0), mesh::RemoteSharedEntity(0, 2));
      EXPECT_EQ(get_remote(entityOrigins, verts[2], 0), mesh::RemoteSharedEntity(0, 5));
      EXPECT_EQ(get_remote(entityOrigins, verts[3], 0), mesh::RemoteSharedEntity(0, 4));

      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(0), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(1), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(2), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(3), 0), 1);

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(0), 0), mesh::RemoteSharedEntity(0, 4));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 0), mesh::RemoteSharedEntity(0, 5));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(2), 0), mesh::RemoteSharedEntity(0, 6));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 0), mesh::RemoteSharedEntity(0, 1));

      EXPECT_EQ(entityOrigins->get_num_comp(el, 0), 1);
      EXPECT_EQ(get_remote(entityOrigins, el, 0), mesh::RemoteSharedEntity(0, 1));      
    } else if (myRank == 2)
    {
      EXPECT_EQ(verts[0]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[1]->count_remote_shared_entities(), 3);
      EXPECT_EQ(verts[2]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[3]->count_remote_shared_entities(), 0);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 3));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 3));

      EXPECT_EQ(el->get_down(0)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(1)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(2)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(3)->count_remote_shared_entities(), 0);   

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(0), 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(1), 3)); 

      EXPECT_EQ(entityOrigins->get_num_comp(verts[0], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[1], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[2], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[3], 0), 1);

      EXPECT_EQ(get_remote(entityOrigins, verts[0], 0), mesh::RemoteSharedEntity(0, 3));
      EXPECT_EQ(get_remote(entityOrigins, verts[1], 0), mesh::RemoteSharedEntity(0, 4));
      EXPECT_EQ(get_remote(entityOrigins, verts[2], 0), mesh::RemoteSharedEntity(0, 7));
      EXPECT_EQ(get_remote(entityOrigins, verts[3], 0), mesh::RemoteSharedEntity(0, 6));

      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(0), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(1), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(2), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(3), 0), 1);

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(0), 0), mesh::RemoteSharedEntity(0, 2));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 0), mesh::RemoteSharedEntity(0, 7));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(2), 0), mesh::RemoteSharedEntity(0, 8));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 0), mesh::RemoteSharedEntity(0, 9));

      EXPECT_EQ(entityOrigins->get_num_comp(el, 0), 1);
      EXPECT_EQ(get_remote(entityOrigins, el, 0), mesh::RemoteSharedEntity(0, 2));                 
    } else if (myRank == 3)
    {
      EXPECT_EQ(verts[0]->count_remote_shared_entities(), 3);
      EXPECT_EQ(verts[1]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[2]->count_remote_shared_entities(), 0);
      EXPECT_EQ(verts[3]->count_remote_shared_entities(), 1);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 2));

      EXPECT_EQ(el->get_down(0)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(1)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(2)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(3)->count_remote_shared_entities(), 1);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(0), 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(3), 2));   

      EXPECT_EQ(entityOrigins->get_num_comp(verts[0], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[1], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[2], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[3], 0), 1);

      EXPECT_EQ(get_remote(entityOrigins, verts[0], 0), mesh::RemoteSharedEntity(0, 4));
      EXPECT_EQ(get_remote(entityOrigins, verts[1], 0), mesh::RemoteSharedEntity(0, 5));
      EXPECT_EQ(get_remote(entityOrigins, verts[2], 0), mesh::RemoteSharedEntity(0, 8));
      EXPECT_EQ(get_remote(entityOrigins, verts[3], 0), mesh::RemoteSharedEntity(0, 7));

      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(0), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(1), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(2), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(3), 0), 1);

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(0), 0), mesh::RemoteSharedEntity(0, 6));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 0), mesh::RemoteSharedEntity(0, 10));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(2), 0), mesh::RemoteSharedEntity(0, 11));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 0), mesh::RemoteSharedEntity(0, 7));

      EXPECT_EQ(entityOrigins->get_num_comp(el, 0), 1);
      EXPECT_EQ(get_remote(entityOrigins, el, 0), mesh::RemoteSharedEntity(0, 3));        
    }

    mesh::check_topology(meshScattered);
    EXPECT_EQ((*elementOrigins)(el, 0, 0), mesh::RemoteSharedEntity(0, myRank));

    check_remotes_centroids(meshScattered);
  }


  MPI_Comm_free(&meshComm);

}


TEST(MeshScatter, 2x2FromTwoProcs)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 6)
    GTEST_SKIP();

  std::shared_ptr<mesh::Mesh> mesh;
  std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, mesh);
  MPI_Comm unionComm = MPI_COMM_WORLD;
  MPI_Comm meshComm;
  int color = utils::impl::comm_rank(unionComm) <= 1 ? 0 : 1;
  MPI_Comm_split(unionComm, color, 0, &meshComm);
  double xStart = utils::impl::comm_rank(meshComm) * 0.5;

  if (color == 0)
  {
    mesh = mesh::make_empty_mesh(meshComm);
    int myRank = utils::impl::comm_rank(meshComm);
    auto v1 = mesh->create_vertex(xStart + 0,   0,   0);
    auto v2 = mesh->create_vertex(xStart + 0.5, 0,   0);
    auto v3 = mesh->create_vertex(xStart + 0,   0.5, 0);
    auto v4 = mesh->create_vertex(xStart + 0.5, 0.5, 0);
    auto v5 = mesh->create_vertex(xStart + 0,   1,   0);
    auto v6 = mesh->create_vertex(xStart + 0.5, 1,   0);
    
    auto el1 = mesh->create_quad_from_verts(v1, v2, v4, v3);
    auto el2 = mesh->create_quad_from_verts(v3, v4, v6, v5);

    if (myRank == 0)
    {
      v2->add_remote_shared_entity({1, 0});
      v4->add_remote_shared_entity({1, 2});
      v6->add_remote_shared_entity({1, 4});

      el1->get_down(1)->add_remote_shared_entity({1, 3});
      el2->get_down(1)->add_remote_shared_entity({1, 6});
    } else
    {
      v1->add_remote_shared_entity({0, 1});
      v3->add_remote_shared_entity({0, 3});
      v5->add_remote_shared_entity({0, 5});

      el1->get_down(3)->add_remote_shared_entity({0, 1});
      el2->get_down(3)->add_remote_shared_entity({0, 4});      
    }

    scatterSpec->add_destination(el1, 2*myRank + 2);
    scatterSpec->add_destination(el2, 2*myRank + 3);
  }

  mesh::impl::MeshScatter scatter(scatterSpec, mesh, color == 1 ? meshComm : MPI_COMM_NULL, true);
  auto meshScattered  = scatter.scatter();
  auto elementOrigins = scatter.get_element_origins();
  auto entityOrigins  = scatter.get_entity_origins();
  auto entityDests    = scatter.get_entity_destinations();

  if (color == 0)
  {
    EXPECT_EQ(entityOrigins, nullptr);
    EXPECT_NE(entityDests, nullptr);

    auto v1 = get_closest_entity(mesh, 0, {xStart,       0,   0});
    auto v2 = get_closest_entity(mesh, 0, {xStart + 0.5, 0,   0});
    auto v3 = get_closest_entity(mesh, 0, {xStart,       0.5,   0});
    auto v4 = get_closest_entity(mesh, 0, {xStart + 0.5, 0.5, 0});
    auto v5 = get_closest_entity(mesh, 0, {xStart,       1.0, 0});
    auto v6 = get_closest_entity(mesh, 0, {xStart + 0.5, 1.0, 0});

    double x1 = xStart, x2 = xStart + 0.5;
    auto edge1 = get_closest_entity(mesh, 1, {(x1 + x2)/2, 0,    0});
    auto edge2 = get_closest_entity(mesh, 1, {x1 + 0.5,    0.25, 0});
    auto edge3 = get_closest_entity(mesh, 1, {(x1 + x2)/2, 0.5,  0});
    auto edge4 = get_closest_entity(mesh, 1, {x1,          0.25, 0});
    auto edge5 = get_closest_entity(mesh, 1, {x1 + 0.5,    0.75, 0});
    auto edge6 = get_closest_entity(mesh, 1, {(x1 + x2)/2, 1.0, 0});
    auto edge7 = get_closest_entity(mesh, 1, {x1,          0.75, 0});

    auto el1 = get_closest_entity(mesh, 2, {(x1 + x2)/2,  0.25, 0});
    auto el2 = get_closest_entity(mesh, 2, {(x1 + x2)/2,  0.75, 0});

    int myRank = utils::impl::comm_rank(meshComm);

    if (myRank == 0)
    {
      EXPECT_EQ(entityDests->get_num_comp(v1, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(v2, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(v3, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(v4, 0), 4);
      EXPECT_EQ(entityDests->get_num_comp(v5, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(v6, 0), 2);

      EXPECT_EQ(entityDests->get_num_comp(edge1, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(edge2, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(edge3, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(edge4, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(edge5, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(edge6, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(edge7, 0), 1);

      EXPECT_EQ(entityDests->get_num_comp(el1, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(el2, 0), 1);

      EXPECT_EQ(get_remote(entityDests, v1, 2), mesh::RemoteSharedEntity(2, 0));

      EXPECT_EQ(get_remote(entityDests, v2, 2), mesh::RemoteSharedEntity(2, 1));
      EXPECT_EQ(get_remote(entityDests, v2, 4), mesh::RemoteSharedEntity(4, 0));

      EXPECT_EQ(get_remote(entityDests, v3, 2), mesh::RemoteSharedEntity(2, 2));
      EXPECT_EQ(get_remote(entityDests, v3, 3), mesh::RemoteSharedEntity(3, 0));

      EXPECT_EQ(get_remote(entityDests, v4, 2), mesh::RemoteSharedEntity(2, 3));
      EXPECT_EQ(get_remote(entityDests, v4, 3), mesh::RemoteSharedEntity(3, 1));
      EXPECT_EQ(get_remote(entityDests, v4, 4), mesh::RemoteSharedEntity(4, 1));
      EXPECT_EQ(get_remote(entityDests, v4, 5), mesh::RemoteSharedEntity(5, 0));

      EXPECT_EQ(get_remote(entityDests, v5, 3), mesh::RemoteSharedEntity(3, 2));

      EXPECT_EQ(get_remote(entityDests, v6, 3), mesh::RemoteSharedEntity(3, 3));
      EXPECT_EQ(get_remote(entityDests, v6, 5), mesh::RemoteSharedEntity(5, 1));

      EXPECT_EQ(get_remote(entityDests, edge1, 2), mesh::RemoteSharedEntity(2, 0));

      EXPECT_EQ(get_remote(entityDests, edge2, 2), mesh::RemoteSharedEntity(2, 1));
      EXPECT_EQ(get_remote(entityDests, edge2, 4), mesh::RemoteSharedEntity(4, 0));

      EXPECT_EQ(get_remote(entityDests, edge3, 2), mesh::RemoteSharedEntity(2, 2));
      EXPECT_EQ(get_remote(entityDests, edge3, 3), mesh::RemoteSharedEntity(3, 0));

      EXPECT_EQ(get_remote(entityDests, edge4, 2), mesh::RemoteSharedEntity(2, 3));

      EXPECT_EQ(get_remote(entityDests, edge5, 3), mesh::RemoteSharedEntity(3, 1));
      EXPECT_EQ(get_remote(entityDests, edge5, 5), mesh::RemoteSharedEntity(5, 0));

      EXPECT_EQ(get_remote(entityDests, edge6, 3), mesh::RemoteSharedEntity(3, 2));

      EXPECT_EQ(get_remote(entityDests, edge7, 3), mesh::RemoteSharedEntity(3, 3));

      EXPECT_EQ(get_remote(entityDests, el1, 2), mesh::RemoteSharedEntity(2, 0));
      EXPECT_EQ(get_remote(entityDests, el2, 3), mesh::RemoteSharedEntity(3, 0));


    } else if (myRank == 1)
    {
      EXPECT_EQ(entityDests->get_num_comp(v1, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(v2, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(v3, 0), 4);
      EXPECT_EQ(entityDests->get_num_comp(v4, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(v5, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(v6, 0), 1); 

      EXPECT_EQ(entityDests->get_num_comp(edge1, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(edge2, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(edge3, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(edge4, 0), 2);
      EXPECT_EQ(entityDests->get_num_comp(edge5, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(edge6, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(edge7, 0), 2);      

      EXPECT_EQ(entityDests->get_num_comp(el1, 0), 1);
      EXPECT_EQ(entityDests->get_num_comp(el2, 0), 1);

      EXPECT_EQ(get_remote(entityDests, v1, 2), mesh::RemoteSharedEntity(2, 1));
      EXPECT_EQ(get_remote(entityDests, v1, 4), mesh::RemoteSharedEntity(4, 0));

      EXPECT_EQ(get_remote(entityDests, v2, 4), mesh::RemoteSharedEntity(4, 2));

      EXPECT_EQ(get_remote(entityDests, v3, 2), mesh::RemoteSharedEntity(2, 3));
      EXPECT_EQ(get_remote(entityDests, v3, 3), mesh::RemoteSharedEntity(3, 1));
      EXPECT_EQ(get_remote(entityDests, v3, 4), mesh::RemoteSharedEntity(4, 1));
      EXPECT_EQ(get_remote(entityDests, v3, 5), mesh::RemoteSharedEntity(5, 0));

      EXPECT_EQ(get_remote(entityDests, v4, 4), mesh::RemoteSharedEntity(4, 3));
      EXPECT_EQ(get_remote(entityDests, v4, 5), mesh::RemoteSharedEntity(5, 2));

      EXPECT_EQ(get_remote(entityDests, v5, 3), mesh::RemoteSharedEntity(3, 3));
      EXPECT_EQ(get_remote(entityDests, v5, 5), mesh::RemoteSharedEntity(5, 1));

      EXPECT_EQ(get_remote(entityDests, v6, 5), mesh::RemoteSharedEntity(5, 3));

      EXPECT_EQ(get_remote(entityDests, edge1, 4), mesh::RemoteSharedEntity(4, 1));
            
      EXPECT_EQ(get_remote(entityDests, edge2, 4), mesh::RemoteSharedEntity(4, 2));

      EXPECT_EQ(get_remote(entityDests, edge3, 4), mesh::RemoteSharedEntity(4, 3));
      EXPECT_EQ(get_remote(entityDests, edge3, 5), mesh::RemoteSharedEntity(5, 1));

      EXPECT_EQ(get_remote(entityDests, edge4, 2), mesh::RemoteSharedEntity(2, 1));
      EXPECT_EQ(get_remote(entityDests, edge4, 4), mesh::RemoteSharedEntity(4, 0));

      EXPECT_EQ(get_remote(entityDests, edge5, 5), mesh::RemoteSharedEntity(5, 2));

      EXPECT_EQ(get_remote(entityDests, edge6, 5), mesh::RemoteSharedEntity(5, 3));

      EXPECT_EQ(get_remote(entityDests, edge7, 3), mesh::RemoteSharedEntity(3, 1));
      EXPECT_EQ(get_remote(entityDests, edge7, 5), mesh::RemoteSharedEntity(5, 0));

      EXPECT_EQ(get_remote(entityDests, el1, 4), mesh::RemoteSharedEntity(4, 0));
      EXPECT_EQ(get_remote(entityDests, el2, 5), mesh::RemoteSharedEntity(5, 0));
    }
  } else if (color == 1)
  {
    EXPECT_EQ(mesh::count_valid(meshScattered->get_vertices()), 4);
    EXPECT_EQ(mesh::count_valid(meshScattered->get_edges()), 4);
    EXPECT_EQ(mesh::count_valid(meshScattered->get_elements()), 1);

    EXPECT_NE(entityOrigins, nullptr);
    EXPECT_EQ(entityDests, nullptr);

    int myRank = utils::impl::comm_rank(meshComm);
    utils::Point lowerLeftCorner;
    switch (myRank)
    {
      case 0: {lowerLeftCorner = {0,   0,   0}; break; }
      case 1: {lowerLeftCorner = {0.0, 0.5, 0}; break; }
      case 2: {lowerLeftCorner = {0.5, 0.0, 0}; break; }
      case 3: {lowerLeftCorner = {0.5, 0.5, 0}; break; }
    }

    auto el = meshScattered->get_elements()[0];
    std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
    int nverts = mesh::get_downward(el, 0, verts.data());
    EXPECT_EQ(nverts, 4);
    expect_near(verts[0]->get_point_orig(0), lowerLeftCorner, 1e-13);
    expect_near(verts[1]->get_point_orig(0), lowerLeftCorner + utils::Point(0.5, 0,   0), 1e-13);
    expect_near(verts[2]->get_point_orig(0), lowerLeftCorner + utils::Point(0.5, 0.5, 0), 1e-13);
    expect_near(verts[3]->get_point_orig(0), lowerLeftCorner + utils::Point(0,   0.5, 0), 1e-13);

   if (myRank == 0)
    {
      EXPECT_EQ(verts[0]->count_remote_shared_entities(), 0);
      EXPECT_EQ(verts[1]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[2]->count_remote_shared_entities(), 3);
      EXPECT_EQ(verts[3]->count_remote_shared_entities(), 1);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 3));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 1));      

      EXPECT_EQ(el->get_down(0)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(1)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(2)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(3)->count_remote_shared_entities(), 0);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(1), 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(2), 1));

      EXPECT_EQ(entityOrigins->get_num_comp(verts[0], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[1], 0), 2);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[2], 0), 2);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[3], 0), 1);

      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(0), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(1), 0), 2);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(2), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(3), 0), 1);      

      EXPECT_EQ(entityOrigins->get_num_comp(el, 0), 1);      

      EXPECT_EQ(get_remote(entityOrigins, verts[0], 0), mesh::RemoteSharedEntity(0, 0));

      EXPECT_EQ(get_remote(entityOrigins, verts[1], 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ(get_remote(entityOrigins, verts[1], 1), mesh::RemoteSharedEntity(1, 0));

      EXPECT_EQ(get_remote(entityOrigins, verts[2], 0), mesh::RemoteSharedEntity(0, 3));
      EXPECT_EQ(get_remote(entityOrigins, verts[2], 1), mesh::RemoteSharedEntity(1, 2));

      EXPECT_EQ(get_remote(entityOrigins, verts[3], 0), mesh::RemoteSharedEntity(0, 2));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(0), 0), mesh::RemoteSharedEntity(0, 0));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 1), mesh::RemoteSharedEntity(1, 3));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(2), 0), mesh::RemoteSharedEntity(0, 2));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 0), mesh::RemoteSharedEntity(0, 3));

      EXPECT_EQ(get_remote(entityOrigins, el, 0), mesh::RemoteSharedEntity(0, 0));
    } else if (myRank == 1)
    {
      EXPECT_EQ(verts[0]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[1]->count_remote_shared_entities(), 3);
      EXPECT_EQ(verts[2]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[3]->count_remote_shared_entities(), 0);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 3));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 3));

      EXPECT_EQ(el->get_down(0)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(1)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(2)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(3)->count_remote_shared_entities(), 0);   

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(0), 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(1), 3));

      EXPECT_EQ(entityOrigins->get_num_comp(verts[0], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[1], 0), 2);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[2], 0), 2);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[3], 0), 1);

      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(0), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(1), 0), 2);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(2), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(3), 0), 1);      

      EXPECT_EQ(entityOrigins->get_num_comp(el, 0), 1);      

      EXPECT_EQ(get_remote(entityOrigins, verts[0], 0), mesh::RemoteSharedEntity(0, 2));

      EXPECT_EQ(get_remote(entityOrigins, verts[1], 0), mesh::RemoteSharedEntity(0, 3));
      EXPECT_EQ(get_remote(entityOrigins, verts[1], 1), mesh::RemoteSharedEntity(1, 2));

      EXPECT_EQ(get_remote(entityOrigins, verts[2], 0), mesh::RemoteSharedEntity(0, 5));
      EXPECT_EQ(get_remote(entityOrigins, verts[2], 1), mesh::RemoteSharedEntity(1, 4));

      EXPECT_EQ(get_remote(entityOrigins, verts[3], 0), mesh::RemoteSharedEntity(0, 4));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(0), 0), mesh::RemoteSharedEntity(0, 2));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 0), mesh::RemoteSharedEntity(0, 4));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 1), mesh::RemoteSharedEntity(1, 6));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(2), 0), mesh::RemoteSharedEntity(0, 5));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 0), mesh::RemoteSharedEntity(0, 6));

      EXPECT_EQ(get_remote(entityOrigins, el, 0), mesh::RemoteSharedEntity(0, 1));      
    } else if (myRank == 2)
    {
      EXPECT_EQ(verts[0]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[1]->count_remote_shared_entities(), 0);
      EXPECT_EQ(verts[2]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[3]->count_remote_shared_entities(), 3);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[2], 3));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 3));

      EXPECT_EQ(el->get_down(0)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(1)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(2)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(3)->count_remote_shared_entities(), 1);   

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(2), 3));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(3), 0));   

      EXPECT_EQ(entityOrigins->get_num_comp(verts[0], 0), 2);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[1], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[2], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[3], 0), 2);

      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(0), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(1), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(2), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(3), 0), 2);      

      EXPECT_EQ(entityOrigins->get_num_comp(el, 0), 1);      

      EXPECT_EQ(get_remote(entityOrigins, verts[0], 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ(get_remote(entityOrigins, verts[0], 1), mesh::RemoteSharedEntity(1, 0));

      EXPECT_EQ(get_remote(entityOrigins, verts[1], 1), mesh::RemoteSharedEntity(1, 1));

      EXPECT_EQ(get_remote(entityOrigins, verts[2], 1), mesh::RemoteSharedEntity(1, 3));

      EXPECT_EQ(get_remote(entityOrigins, verts[3], 0), mesh::RemoteSharedEntity(0, 3));
      EXPECT_EQ(get_remote(entityOrigins, verts[3], 1), mesh::RemoteSharedEntity(1, 2));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(0), 1), mesh::RemoteSharedEntity(1, 0));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 1), mesh::RemoteSharedEntity(1, 1));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(2), 1), mesh::RemoteSharedEntity(1, 2));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 1), mesh::RemoteSharedEntity(1, 3));

      EXPECT_EQ(get_remote(entityOrigins, el, 1), mesh::RemoteSharedEntity(1, 0));
    } else if (myRank == 3)
    {
      EXPECT_EQ(verts[0]->count_remote_shared_entities(), 3);
      EXPECT_EQ(verts[1]->count_remote_shared_entities(), 1);
      EXPECT_EQ(verts[2]->count_remote_shared_entities(), 0);
      EXPECT_EQ(verts[3]->count_remote_shared_entities(), 1);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 0));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 1));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[0], 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[1], 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(verts[3], 1));

      EXPECT_EQ(el->get_down(0)->count_remote_shared_entities(), 1);
      EXPECT_EQ(el->get_down(1)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(2)->count_remote_shared_entities(), 0);
      EXPECT_EQ(el->get_down(3)->count_remote_shared_entities(), 1);

      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(0), 2));
      EXPECT_NO_THROW(mesh::get_remote_shared_entity(el->get_down(3), 1));  

      EXPECT_EQ(entityOrigins->get_num_comp(verts[0], 0), 2);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[1], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[2], 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(verts[3], 0), 2);

      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(0), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(1), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(2), 0), 1);
      EXPECT_EQ(entityOrigins->get_num_comp(el->get_down(3), 0), 2);

      EXPECT_EQ(entityOrigins->get_num_comp(el, 0), 1);      

      EXPECT_EQ(get_remote(entityOrigins, verts[0], 0), mesh::RemoteSharedEntity(0, 3));
      EXPECT_EQ(get_remote(entityOrigins, verts[0], 1), mesh::RemoteSharedEntity(1, 2));

      EXPECT_EQ(get_remote(entityOrigins, verts[1], 1), mesh::RemoteSharedEntity(1, 3));

      EXPECT_EQ(get_remote(entityOrigins, verts[2], 1), mesh::RemoteSharedEntity(1, 5));

      EXPECT_EQ(get_remote(entityOrigins, verts[3], 0), mesh::RemoteSharedEntity(0, 5));
      EXPECT_EQ(get_remote(entityOrigins, verts[3], 1), mesh::RemoteSharedEntity(1, 4));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(0), 1), mesh::RemoteSharedEntity(1, 2));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(1), 1), mesh::RemoteSharedEntity(1, 4));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(2), 1), mesh::RemoteSharedEntity(1, 5));

      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 0), mesh::RemoteSharedEntity(0, 4));
      EXPECT_EQ(get_remote(entityOrigins, el->get_down(3), 1), mesh::RemoteSharedEntity(1, 6));

      EXPECT_EQ(get_remote(entityOrigins, el, 1), mesh::RemoteSharedEntity(1, 1));         
    }

    check_remotes_centroids(meshScattered);
    mesh::check_topology(meshScattered);

    if (myRank == 0) {
      EXPECT_EQ((*elementOrigins)(el, 0, 0), mesh::RemoteSharedEntity(0, 0));
    } else if (myRank == 1) {
      EXPECT_EQ((*elementOrigins)(el, 0, 0), mesh::RemoteSharedEntity(0, 1));
    } else if (myRank == 2) {
      EXPECT_EQ((*elementOrigins)(el, 0, 0), mesh::RemoteSharedEntity(1, 0));
    } else {
      EXPECT_EQ((*elementOrigins)(el, 0, 0), mesh::RemoteSharedEntity(1, 1));
    }
  }

  MPI_Comm_free(&meshComm);
}


TEST(MeshScatter, 2x2FromTwoProcsToTwoProcs)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 3)
    GTEST_SKIP();

  std::shared_ptr<mesh::Mesh> mesh;
  std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, mesh);
  MPI_Comm unionComm = MPI_COMM_WORLD;
  MPI_Comm inputMeshComm, outputMeshComm;
  int color1 = utils::impl::comm_rank(unionComm) <= 1 ? 0 : MPI_UNDEFINED;
  int color2 = utils::impl::comm_rank(unionComm) >= 1 ? 0 : MPI_UNDEFINED;
  MPI_Comm_split(unionComm, color1, 0, &inputMeshComm);
  MPI_Comm_split(unionComm, color2, 0, &outputMeshComm);


  if (color1 == 0)
  {
    mesh = mesh::make_empty_mesh(inputMeshComm);
    int myRank = utils::impl::comm_rank(inputMeshComm);
    double xStart = myRank * 0.5;
    auto v1 = mesh->create_vertex(xStart + 0,   0,   0);
    auto v2 = mesh->create_vertex(xStart + 0.5, 0,   0);
    auto v3 = mesh->create_vertex(xStart + 0,   0.5, 0);
    auto v4 = mesh->create_vertex(xStart + 0.5, 0.5, 0);
    auto v5 = mesh->create_vertex(xStart + 0,   1,   0);
    auto v6 = mesh->create_vertex(xStart + 0.5, 1,   0);
    
    auto el1 = mesh->create_quad_from_verts(v1, v2, v4, v3);
    auto el2 = mesh->create_quad_from_verts(v3, v4, v6, v5);

    if (myRank == 0)
    {
      v2->add_remote_shared_entity({1, 0});
      v4->add_remote_shared_entity({1, 2});
      v6->add_remote_shared_entity({1, 4});

      el1->get_down(1)->add_remote_shared_entity({1, 3});
      el2->get_down(1)->add_remote_shared_entity({1, 6});

      scatterSpec->add_destination(el1, 1);
      scatterSpec->add_destination(el2, 1);
    } else
    {
      v1->add_remote_shared_entity({0, 1});
      v3->add_remote_shared_entity({0, 3});
      v5->add_remote_shared_entity({0, 5});

      el1->get_down(3)->add_remote_shared_entity({0, 1});
      el2->get_down(3)->add_remote_shared_entity({0, 4});  

      scatterSpec->add_destination(el1, 1);
      scatterSpec->add_destination(el2, 2);          
    }
  }

  mesh::impl::MeshScatter scatter(scatterSpec, mesh, outputMeshComm);
  auto meshScattered  = scatter.scatter();
  auto elementOrigins = scatter.get_element_origins();

  if (color2 == 0)
  {
    int myRank = utils::impl::comm_rank(outputMeshComm);

    if (myRank == 0)
    {
      EXPECT_EQ(mesh::count_valid(meshScattered->get_vertices()), 8);
      EXPECT_EQ(mesh::count_valid(meshScattered->get_edges()), 10);
      EXPECT_EQ(mesh::count_valid(meshScattered->get_elements()), 3);

      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0,   0,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 0,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0,   0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0,   1,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 1,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1,   0,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1,   0.5, 0}));

      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0,    0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0,    0.75, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.5,  0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.5,  0.75, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {1.0,  0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.25, 0,    0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.25, 0.5,  0}));      
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.25, 1.0,  0}));      
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 0,    0}));      
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 0.5,  0}));

      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.25, 0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.25, 0.75, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.75, 0.25, 0}));
    } else
    {
      EXPECT_EQ(mesh::count_valid(meshScattered->get_vertices()), 4);
      EXPECT_EQ(mesh::count_valid(meshScattered->get_edges()), 4);
      EXPECT_EQ(mesh::count_valid(meshScattered->get_elements()), 1);

      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1.0, 0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 1.0, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1.0, 1.0, 0}));
 
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 1.0, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.5,  0.75, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {1.0,  0.75, 0}));
 
      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.75, 0.75, 0}));
    }

    auto v1 = get_closest_entity(meshScattered, 0, {0.5, 0.5, 0});
    auto v2 = get_closest_entity(meshScattered, 0, {1.0, 0.5, 0});
    auto v3 = get_closest_entity(meshScattered, 0, {0.5, 1.0, 0});
    std::set<mesh::MeshEntityPtr, mesh::MeshEntityCompare> shared_verts = {v1, v2, v3};

    auto edge1 = get_closest_entity(meshScattered, 1, {0.75, 0.50, 0});
    auto edge2 = get_closest_entity(meshScattered, 1, {0.50, 0.75, 0});
    std::set<mesh::MeshEntityPtr, mesh::MeshEntityCompare> shared_edges = {edge1, edge2};

    for (auto& v : meshScattered->get_vertices())
      if (v)
      {
        if (shared_verts.count(v) == 0)
        {
          EXPECT_EQ(v->count_remote_shared_entities(), 0);
        } else
        {
          EXPECT_EQ(v->count_remote_shared_entities(), 1);
        }
      }

    for (auto& edge : meshScattered->get_edges())
      if (edge)
      {
        if (shared_edges.count(edge) == 0)
        {
          EXPECT_EQ(edge->count_remote_shared_entities(), 0);
        } else
        {
          EXPECT_EQ(edge->count_remote_shared_entities(), 1);
        }
      }


    check_remotes_centroids(meshScattered);
    mesh::check_topology(meshScattered);

    if (myRank == 0) {
      auto el1 = get_closest_entity(meshScattered, 2, {0.25, 0.25, 0});
      auto el2 = get_closest_entity(meshScattered, 2, {0.25, 0.75, 0});
      auto el3 = get_closest_entity(meshScattered, 2, {0.75, 0.25, 0});
      EXPECT_EQ((*elementOrigins)(el1, 0, 0), mesh::RemoteSharedEntity(0, 0));
      EXPECT_EQ((*elementOrigins)(el2, 0, 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ((*elementOrigins)(el3, 0, 0), mesh::RemoteSharedEntity(1, 0));
    } else {
      auto el1 = get_closest_entity(meshScattered, 2, {0.75, 0.75, 0});
      EXPECT_EQ((*elementOrigins)(el1, 0, 0), mesh::RemoteSharedEntity(1, 1));
    }
  }

  if (inputMeshComm != MPI_COMM_NULL)
    MPI_Comm_free(&inputMeshComm);

  if (outputMeshComm != MPI_COMM_NULL)
    MPI_Comm_free(&outputMeshComm);    
}

TEST(MeshScatter, MultipleDestinations)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  std::shared_ptr<mesh::Mesh> mesh;
  std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, mesh);
  MPI_Comm unionComm = MPI_COMM_WORLD;
  MPI_Comm inputMeshComm, outputMeshComm;
  int myrankOnUnionComm = utils::impl::comm_rank(unionComm);
  int color1 = myrankOnUnionComm <= 1 ? 0 : MPI_UNDEFINED;
  int color2 = myrankOnUnionComm >= 1 && myrankOnUnionComm <= 2 ? 0 : MPI_UNDEFINED;
  MPI_Comm_split(unionComm, color1, 0, &inputMeshComm);
  MPI_Comm_split(unionComm, color2, 0, &outputMeshComm);

  if (color1 == 0)
  {
    mesh = mesh::make_empty_mesh(inputMeshComm);
    int myRank = utils::impl::comm_rank(inputMeshComm);
    double xStart = myRank * 0.5;
    auto v1 = mesh->create_vertex(xStart + 0,   0,   0);
    auto v2 = mesh->create_vertex(xStart + 0.5, 0,   0);
    auto v3 = mesh->create_vertex(xStart + 0,   0.5, 0);
    auto v4 = mesh->create_vertex(xStart + 0.5, 0.5, 0);
    auto v5 = mesh->create_vertex(xStart + 0,   1,   0);
    auto v6 = mesh->create_vertex(xStart + 0.5, 1,   0);
    
    auto el1 = mesh->create_quad_from_verts(v1, v2, v4, v3);
    auto el2 = mesh->create_quad_from_verts(v3, v4, v6, v5);

    if (myRank == 0)
    {
      v2->add_remote_shared_entity({1, 0});
      v4->add_remote_shared_entity({1, 2});
      v6->add_remote_shared_entity({1, 4});

      el1->get_down(1)->add_remote_shared_entity({1, 3});
      el2->get_down(1)->add_remote_shared_entity({1, 6});

      scatterSpec->add_destination(el1, 1);
      scatterSpec->add_destination(el2, 1);
    } else
    {
      v1->add_remote_shared_entity({0, 1});
      v3->add_remote_shared_entity({0, 3});
      v5->add_remote_shared_entity({0, 5});

      el1->get_down(3)->add_remote_shared_entity({0, 1});
      el2->get_down(3)->add_remote_shared_entity({0, 4});  

      scatterSpec->add_destination(el1, 1);
      scatterSpec->add_destination(el1, 2);
      scatterSpec->add_destination(el2, 2);          
    }
  }

  mesh::impl::MeshScatter scatter(scatterSpec, mesh, outputMeshComm);
  auto meshScattered  = scatter.scatter();
  auto elementOrigins = scatter.get_element_origins();

  if (color2 == 0)
  {
    int myRank = utils::impl::comm_rank(outputMeshComm);

    if (myRank == 0)
    {
      EXPECT_EQ(mesh::count_valid(meshScattered->get_vertices()), 8);
      EXPECT_EQ(mesh::count_valid(meshScattered->get_edges()),    10);
      EXPECT_EQ(mesh::count_valid(meshScattered->get_elements()), 3);

      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0,   0,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 0,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0,   0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0,   1,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 1,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1,   0,   0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1,   0.5, 0}));

      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0,    0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0,    0.75, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.5,  0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.5,  0.75, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {1.0,  0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.25, 0,    0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.25, 0.5,  0}));      
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.25, 1.0,  0}));      
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 0,    0}));      
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 0.5,  0}));

      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.25, 0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.25, 0.75, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.75, 0.25, 0}));
    } else
    {
      EXPECT_EQ(mesh::count_valid(meshScattered->get_vertices()), 6);
      EXPECT_EQ(mesh::count_valid(meshScattered->get_edges()), 7);
      EXPECT_EQ(mesh::count_valid(meshScattered->get_elements()), 2);

      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 0.0, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1.0, 0.0, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1.0, 0.5, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {0.5, 1.0, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 0, {1.0, 1.0, 0}));
 
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 0.0,  0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {1.00, 0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.5,  0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 0.5,  0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.75, 1.0,  0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {0.5,  0.75, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 1, {1.0,  0.75,  0}));
 
      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.75, 0.25, 0}));
      EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.75, 0.75, 0}));
    }

    EXPECT_EQ(get_closest_entity(meshScattered, 0, {0.5, 0,   0})->count_remote_shared_entities(), 1);
    EXPECT_EQ(get_closest_entity(meshScattered, 0, {1.0, 0,   0})->count_remote_shared_entities(), 1);
    EXPECT_EQ(get_closest_entity(meshScattered, 0, {0.5, 0.5, 0})->count_remote_shared_entities(), 1);
    EXPECT_EQ(get_closest_entity(meshScattered, 0, {0.5, 1.0, 0})->count_remote_shared_entities(), 1);
    EXPECT_EQ(get_closest_entity(meshScattered, 0, {1.0, 0.5, 0})->count_remote_shared_entities(), 1);

    EXPECT_EQ(get_closest_entity(meshScattered, 1, {0.75, 0.5,  0})->count_remote_shared_entities(), 1);
    EXPECT_EQ(get_closest_entity(meshScattered, 1, {0.5,  0.75, 0})->count_remote_shared_entities(), 1);

    check_remotes_centroids(meshScattered);
    mesh::check_topology(meshScattered);

    if (myRank == 0) {
      auto el1 = get_closest_entity(meshScattered, 2, {0.25, 0.25, 0});
      auto el2 = get_closest_entity(meshScattered, 2, {0.25, 0.75, 0});
      auto el3 = get_closest_entity(meshScattered, 2, {0.75, 0.25, 0});
      EXPECT_EQ((*elementOrigins)(el1, 0, 0), mesh::RemoteSharedEntity(0, 0));
      EXPECT_EQ((*elementOrigins)(el2, 0, 0), mesh::RemoteSharedEntity(0, 1));
      EXPECT_EQ((*elementOrigins)(el3, 0, 0), mesh::RemoteSharedEntity(1, 0));
    } else {
      auto el1 = get_closest_entity(meshScattered, 2, {0.75, 0.25, 0});
      auto el2 = get_closest_entity(meshScattered, 2, {0.75, 0.75, 0});
      EXPECT_EQ((*elementOrigins)(el1, 0, 0), mesh::RemoteSharedEntity(1, 0));
      EXPECT_EQ((*elementOrigins)(el2, 0, 0), mesh::RemoteSharedEntity(1, 1));
    }
  }

  if (inputMeshComm != MPI_COMM_NULL)
    MPI_Comm_free(&inputMeshComm);

  if (outputMeshComm != MPI_COMM_NULL)
    MPI_Comm_free(&outputMeshComm); 
}

TEST(MeshScatter, ThreeDestinations)
{
  if (utils::impl::comm_size(MPI_COMM_WORLD) != 4)
    GTEST_SKIP();

  std::shared_ptr<mesh::Mesh> mesh;
  std::shared_ptr<mesh::impl::MeshScatterSpec> scatterSpec = std::make_shared<mesh::impl::MeshScatterSpec>(MPI_COMM_WORLD, mesh);
  MPI_Comm unionComm = MPI_COMM_WORLD;
  MPI_Comm inputMeshComm, outputMeshComm;
  int myrankOnUnionComm = utils::impl::comm_rank(unionComm);
  int color1 = myrankOnUnionComm <= 1 ? 0 : MPI_UNDEFINED;
  int color2 = myrankOnUnionComm >= 1 ? 0 : MPI_UNDEFINED;
  MPI_Comm_split(unionComm, color1, 0, &inputMeshComm);
  MPI_Comm_split(unionComm, color2, 0, &outputMeshComm);

  if (color1 == 0)
  {
    mesh = mesh::make_empty_mesh(inputMeshComm);
    int myRank = utils::impl::comm_rank(inputMeshComm);
    double xStart = myRank * 0.5;
    auto v1 = mesh->create_vertex(xStart + 0,   0,   0);
    auto v2 = mesh->create_vertex(xStart + 0.5, 0,   0);
    auto v3 = mesh->create_vertex(xStart + 0,   0.5, 0);
    auto v4 = mesh->create_vertex(xStart + 0.5, 0.5, 0);
    auto v5 = mesh->create_vertex(xStart + 0,   1,   0);
    auto v6 = mesh->create_vertex(xStart + 0.5, 1,   0);
    
    auto el1 = mesh->create_quad_from_verts(v1, v2, v4, v3);
    auto el2 = mesh->create_quad_from_verts(v3, v4, v6, v5);

    if (myRank == 0)
    {
      v2->add_remote_shared_entity({1, 0});
      v4->add_remote_shared_entity({1, 2});
      v6->add_remote_shared_entity({1, 4});

      el1->get_down(1)->add_remote_shared_entity({1, 3});
      el2->get_down(1)->add_remote_shared_entity({1, 6});
    } else
    {
      v1->add_remote_shared_entity({0, 1});
      v3->add_remote_shared_entity({0, 3});
      v5->add_remote_shared_entity({0, 5});

      el1->get_down(3)->add_remote_shared_entity({0, 1});
      el2->get_down(3)->add_remote_shared_entity({0, 4});          
    }

    for (int i=1; i < 4; ++i)
    { 
      scatterSpec->add_destination(el1, i);
      scatterSpec->add_destination(el2, i);
    }
  }

  mesh::impl::MeshScatter scatter(scatterSpec, mesh, outputMeshComm);
  auto meshScattered  = scatter.scatter();
  auto elementOrigins = scatter.get_element_origins();

  if (color2 == 0)
  {
    EXPECT_EQ(mesh::count_valid(meshScattered->get_vertices()), 9);
    EXPECT_EQ(mesh::count_valid(meshScattered->get_edges()),    12);
    EXPECT_EQ(mesh::count_valid(meshScattered->get_elements()), 4);

    double delta = 0.5;
    for (int i=0; i < 3; ++i)
      for (int j=0; j < 3; ++j)
      {
        auto vert = get_closest_entity(meshScattered, 0, {i*delta, j*delta,   0});
        EXPECT_TRUE(vert);
        EXPECT_EQ(vert->count_remote_shared_entities(), 2);
      }

    for (int i=0; i < 3; ++i)
    {
      auto edge1 = get_closest_entity(meshScattered, 1, {0.25, i*delta, 0});
      EXPECT_TRUE(edge1);
      EXPECT_EQ(edge1->count_remote_shared_entities(), 2);

      auto edge2 = get_closest_entity(meshScattered, 1, {0.75, i*delta, 0});     
      EXPECT_TRUE(edge2);
      EXPECT_EQ(edge2->count_remote_shared_entities(), 2);
    }

    for (int j=0; j < 3; ++j)
    {
      auto edge1 = get_closest_entity(meshScattered, 1, {j*delta, 0.25, 0});
      EXPECT_TRUE(edge1);
      EXPECT_EQ(edge1->count_remote_shared_entities(), 2);

      auto edge2 = get_closest_entity(meshScattered, 1, {j*delta, 0.75, 0});     
      EXPECT_TRUE(edge2);
      EXPECT_EQ(edge2->count_remote_shared_entities(), 2);
    }

    for (int i=0; i < 2; ++i)
      for (int j=0; j < 2; ++j)
        EXPECT_TRUE(get_closest_entity(meshScattered, 2, {0.25 + i*delta, 0.25 + j*delta, 0}));
            

    check_remotes_centroids(meshScattered);
    mesh::check_topology(meshScattered);

    auto el1 = get_closest_entity(meshScattered, 2, {0.25, 0.25, 0});
    auto el2 = get_closest_entity(meshScattered, 2, {0.25, 0.75, 0});
    auto el3 = get_closest_entity(meshScattered, 2, {0.75, 0.25, 0});
    auto el4 = get_closest_entity(meshScattered, 2, {0.75, 0.75, 0});
    EXPECT_EQ((*elementOrigins)(el1, 0, 0), mesh::RemoteSharedEntity(0, 0));
    EXPECT_EQ((*elementOrigins)(el2, 0, 0), mesh::RemoteSharedEntity(0, 1));
    EXPECT_EQ((*elementOrigins)(el3, 0, 0), mesh::RemoteSharedEntity(1, 0));
    EXPECT_EQ((*elementOrigins)(el4, 0, 0), mesh::RemoteSharedEntity(1, 1));
  }

  if (inputMeshComm != MPI_COMM_NULL)
    MPI_Comm_free(&inputMeshComm);

  if (outputMeshComm != MPI_COMM_NULL)
    MPI_Comm_free(&outputMeshComm); 
}


// test when a single edge is on 2 or more procs
