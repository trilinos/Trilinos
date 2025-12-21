#include <stk_transfer/ReducedDependencyGeometricTransfer.hpp>
#include <gtest/gtest.h>
#include <stk_search/BoundingBox.hpp>
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_topology/topology.hpp"
#include "stk_transfer/ReducedDependencyCommData.hpp"
#include "stk_transfer/ReducedDependencyGeometricTransfer.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_io/FillMesh.hpp"
#include "stk_coupling/SplitComms.hpp"

namespace {

// a simplified reference element.  Only works for axis-aligned elements.
class ReferenceElementQuad
{
  public:
    using Point2 = std::array<double, 2>;

    explicit ReferenceElementQuad(const std::array<Point2, 4>& coords) :
      m_coords(coords)
    {}

    Point2 compute_parametric_coords(const Point2& pt) const
    {
      double xi  = (pt[0] - m_coords[0][0]) / (m_coords[1][0] - m_coords[0][0]);
      double eta = (pt[1] - m_coords[0][1]) / (m_coords[3][1] - m_coords[0][1]);

      return {xi, eta};
    }

    double compute_distance(double xi, double eta) const
    {
      double delta_xi  = std::abs(xi  - m_centroid[0]);
      double delta_eta = std::abs(eta - m_centroid[1]);

      return std::min(delta_xi, delta_eta);
    }

    double interpolate(const std::array<double, 4>& nodeValues, const Point2& pt) const
    {
      double xi = pt[0];
      double eta = pt[1];
      std::array<double, 4> basisValues{(1-xi)*(1-eta), xi*(1-eta), xi*eta, (1-xi)*eta};
      double valInterp = 0;
      for (size_t i=0; i < 4; ++i)
      {
        valInterp += basisValues[i]*nodeValues[i];
      }

      return valInterp;
    }

  private:
    // corner coordinates, starting from bottom left, proceeding counter-clockwise
    std::array<Point2, 4> m_coords;
    static constexpr Point2 m_centroid{0.5, 0.5};
};

struct MeshAndFieldSpec
{
  int nx;
  int ny;
  int nz;
  double field_value_offset = 0;
};

double get_field_value(const MeshAndFieldSpec& spec, double y, double z)
{
  int i = std::floor(y / spec.ny);
  int j = std::floor(z / spec.nz);
  double delta_y = 1.0/spec.ny;
  double delta_z = 1.0/spec.nz;

  double centroid_y = (i + 0.5)*delta_y;
  double centroid_z = (j + 0.5)*delta_z;

  return y + z + 1 + centroid_y + 2*centroid_z + spec.field_value_offset;
}

void set_face_field(stk::mesh::BulkData& bulk, stk::mesh::FieldBase* field, const MeshAndFieldSpec& spec)
{
  auto coordData = bulk.mesh_meta_data().coordinate_field()->data<double>();
  auto fieldData = field->data<double, stk::mesh::ReadWrite>();
  for (stk::mesh::Bucket* bucket : bulk.get_buckets(stk::topology::FACE_RANK, *field))
  {
    for (stk::mesh::Entity face : *bucket)
    {
      auto faceValues = fieldData.entity_values(face);
      stk::mesh::ConnectedEntities nodes = bulk.get_connected_entities(face, stk::topology::NODE_RANK);
      for (stk::mesh::CopyIdx j=0_copy; j < int(nodes.size()); ++j)
      {
        auto nodeCoords = coordData.entity_values(nodes[j]);
        double y = nodeCoords(1_comp);
        double z = nodeCoords(2_comp);
        double val = get_field_value(spec, y, z);
        faceValues(j, 0_comp) = val;
        faceValues(j, 1_comp) = val + 1;
      }
    }
  }
}

void check_nodal_field(stk::mesh::BulkData& bulk, stk::mesh::FieldBase* field, const MeshAndFieldSpec& spec)
{
  auto coordData = bulk.mesh_meta_data().coordinate_field()->data<double>();
  auto fieldData = field->data<double, stk::mesh::ReadOnly>();
  for (stk::mesh::Bucket* bucket : bulk.get_buckets(stk::topology::NODE_RANK, *field))
  {
    for (stk::mesh::Entity node : *bucket)
    {
      auto nodeValues = fieldData.entity_values(node);
      auto nodeCoords = coordData.entity_values(node);
      double y = nodeCoords(1_comp);
      double z = nodeCoords(2_comp);
      double val = get_field_value(spec, y, z);
      EXPECT_DOUBLE_EQ(nodeValues(0_comp), val);
      EXPECT_DOUBLE_EQ(nodeValues(1_comp), val + 1);
    }
  }
}

class FromMeshInterface
{
  public:
    // These are required by the ReducedDependencyGeometricTransfer interface
    using EntityKey     = stk::mesh::EntityKey; // must be global IDs
    using EntityProc    = stk::search::IdentProc<EntityKey, int>;
    using BoundingBox   = std::pair<stk::search::Box<double>, EntityProc>;
    using EntityProcVec = std::vector<EntityProc>;

    virtual ~FromMeshInterface() = default;

    virtual void update_values() = 0;

    virtual void bounding_boxes(std::vector<BoundingBox>& boxes) = 0;
};

class FromMesh : public FromMeshInterface
{
  public:
    using BoundingShape = FromMeshInterface::BoundingBox::first_type;
    using Point2        = ReferenceElementQuad::Point2;

    using FromMeshInterface::EntityKey;
    using FromMeshInterface::EntityProc;
    using FromMeshInterface::BoundingBox;
    using FromMeshInterface::EntityProcVec;

    FromMesh(std::shared_ptr<stk::mesh::BulkData>& mesh, stk::mesh::FieldBase* field,
             stk::ParallelMachine sharedComm) :
      m_mesh(mesh),
      m_field(field),
      m_coordField(m_mesh->mesh_meta_data().coordinate_field()->data<double, stk::mesh::ReadOnly>()),
      m_sharedComm(sharedComm)
    {
      STK_ThrowRequireMsg(field->entity_rank() == stk::topology::FACE_RANK,
                          "This implementation of FromMesh requires a Discontinuous Galerkin type face rank field");
    } 
    
    void update_values() override {};

    void bounding_boxes(std::vector<BoundingBox>& boxes) override
    {
      int myRank = stk::parallel_machine_rank(m_sharedComm);
      boxes.clear();
      for (stk::mesh::Bucket* bucket : m_mesh->get_buckets(stk::topology::FACE_RANK, *m_field))
      {
        for (stk::mesh::Entity face : *bucket)
        {
          std::array<double, 3> min{std::numeric_limits<double>::max(), 
                                    std::numeric_limits<double>::max(),
                                    std::numeric_limits<double>::max()};
          std::array<double, 3> max{std::numeric_limits<double>::min(), 
                                    std::numeric_limits<double>::min(),
                                    std::numeric_limits<double>::min()};
          for (stk::mesh::Entity node : m_mesh->get_connected_entities(face, stk::topology::NODE_RANK))
          {
            auto nodeCoords = m_coordField.entity_values(node);
            for (stk::mesh::ComponentIdx d : nodeCoords.components())
            {
              min[d] = std::min(min[d], nodeCoords(d));
              max[d] = std::max(max[d], nodeCoords(d));
            }
          }

          BoundingShape box(min[0], min[1], min[2], max[0], max[1], max[2]);
          EntityProc identProc(m_mesh->entity_key(face), myRank);
          boxes.emplace_back(box, identProc);
        }
      }
    }

    // These are not required by the ReducedDependencyGeometricTransfer,
    // but are used in the implementation of the Interpolate class

    std::array<double, 3> get_parametric_coords_and_distance(EntityKey faceKey, const std::array<double, 3>& point) const
    {
      stk::mesh::Entity face = m_mesh->get_entity(faceKey);
      std::array<Point2, 4> coords = get_face_coords(face);

      ReferenceElementQuad refElement(coords);
      ReferenceElementQuad::Point2 ptInYZPlane{point[1], point[2]};
      auto [xi, eta] = refElement.compute_parametric_coords(ptInYZPlane);
      double dist = refElement.compute_distance(xi, eta);

      return {xi, eta, dist};
    }

    std::vector<double> interpolate_to_points(const EntityProcVec& fromEntityKeysMasked, const std::vector<Point2>& toPointParametricCoordsOnFromMesh)
    {
      constexpr size_t NodesPerFace = 4;

      std::vector<double> valuesInterp;
      if (fromEntityKeysMasked.size() > 0)
      {
        size_t valuesPerEntity = stk::mesh::field_scalars_per_entity(*m_field, m_mesh->get_entity(fromEntityKeysMasked[0].id()))/NodesPerFace;
        valuesInterp.reserve(fromEntityKeysMasked.size()*valuesPerEntity);
        auto fieldData = m_field->data<double>();

        for (size_t i=0; i < fromEntityKeysMasked.size(); ++i)
        {
          stk::mesh::Entity face = m_mesh->get_entity(fromEntityKeysMasked[i].id());
          auto faceValues = fieldData.entity_values(face);
          ReferenceElementQuad refElem(get_face_coords(face));

          assert(faceValues.num_components() == int(valuesPerEntity));
          assert(faceValues.num_copies()     == int(NodesPerFace));

          for (stk::mesh::ComponentIdx j : faceValues.components())
          {
            std::array<double, NodesPerFace> nodeValues;
            for (stk::mesh::CopyIdx k : faceValues.copies())
            {
              nodeValues[k] = faceValues(k, j);
            }

            valuesInterp.push_back(refElem.interpolate(nodeValues, toPointParametricCoordsOnFromMesh[i]));
          }
        }
      }

      return valuesInterp;
    }

  protected:

    std::array<Point2, 4> get_face_coords(stk::mesh::Entity face) const
    {
      std::array<Point2, 4> coords;
      int idx = 0;
      for (stk::mesh::Entity node : m_mesh->get_connected_entities(face, stk::topology::NODE_RANK))
      {
        auto nodeCoords = m_coordField.entity_values(node);

        // this only works when the surface is in the YZ plane
        coords[idx][0] = nodeCoords(1_comp);
        coords[idx][1] = nodeCoords(2_comp);
        idx++;
      } 
      
      return coords;
    }

    std::shared_ptr<stk::mesh::BulkData> m_mesh;
    stk::mesh::FieldBase* m_field;
    stk::mesh::ConstFieldData<double> m_coordField;
    stk::ParallelMachine m_sharedComm;
    bool m_meshModifiedSinceLastSearch = false;
    int m_numSearches = 0;
};

class FromMeshWithRepeatSearch : public FromMesh
{
  public:
    using FromMeshInterface::BoundingBox;

    FromMeshWithRepeatSearch(std::shared_ptr<stk::mesh::BulkData>& mesh, stk::mesh::FieldBase* field,
             stk::ParallelMachine sharedComm) :
      FromMesh(mesh, field, sharedComm)
    {}

    // These functions are required by the ReducedDependencyGeometricTransfer API

    bool need_repeat_search() const { return m_meshModifiedSinceLastSearch; }

    void bounding_boxes(std::vector<BoundingBox>& boxes) override
    {
      FromMesh::bounding_boxes(boxes);

      m_meshModifiedSinceLastSearch = false;
      m_numSearches++;
    }

    // These functions are used for testing

    void modify_mesh(std::shared_ptr<stk::mesh::BulkData>& mesh, stk::mesh::FieldBase* field)
    {
      m_mesh    = mesh;
      m_field   = field;
      m_coordField = m_mesh->mesh_meta_data().coordinate_field()->data<double, stk::mesh::ReadOnly>();
      m_meshModifiedSinceLastSearch = true;
    }

    int get_num_searches() const { return m_numSearches; }    

  private:
    bool m_meshModifiedSinceLastSearch = false;
    int m_numSearches = 0;
};  
            

class ToMeshInterface
{
  public:

    virtual ~ToMeshInterface() = default;

    using EntityKey                 = stk::mesh::EntityKey;  // must be global IDs
    using EntityProc                = stk::search::IdentProc<EntityKey, int>;
    using BoundingBox               = std::pair<stk::search::Point<double>, EntityProc>;
    using Point                     = stk::search::Point<double>;
    using ToPointsContainer         = std::vector<Point>;
    using ToPointsDistanceContainer = std::vector<double>;
    using EntityProcVec             = std::vector<EntityProc>;

    virtual void update_values()  = 0;

    virtual void get_to_points_coordinates(const EntityProcVec& toEntityKeys, ToPointsContainer& toPointsOnToMesh) = 0;

    virtual void bounding_boxes(std::vector<BoundingBox>& boxes) = 0;
};

class ToMesh : public ToMeshInterface
{
  public:
    using BoundingShape = ToMeshInterface::BoundingBox::first_type;

    using ToMeshInterface::EntityKey;
    using ToMeshInterface::EntityProc;
    using ToMeshInterface::BoundingBox;
    using ToMeshInterface::Point;
    using ToMeshInterface::ToPointsContainer;
    using ToMeshInterface::ToPointsDistanceContainer;
    using ToMeshInterface::EntityProcVec;

    ToMesh(std::shared_ptr<stk::mesh::BulkData>& mesh, stk::mesh::FieldBase* field,
           stk::ParallelMachine sharedComm) :
      m_mesh(mesh),
      m_field(field),
      m_sharedComm(sharedComm)
    {}

    void update_values() override {};

    void get_to_points_coordinates(const EntityProcVec& toEntityKeys, ToPointsContainer& toPointsOnToMesh) override
    {
      auto coordField = m_mesh->mesh_meta_data().coordinate_field()->data<double>();

      toPointsOnToMesh.clear();
      for (EntityProc entityAndProc : toEntityKeys)
      {
        stk::mesh::Entity node = m_mesh->get_entity(entityAndProc.id());
        auto nodeCoords = coordField.entity_values(node);
        toPointsOnToMesh.emplace_back(nodeCoords(0_comp), nodeCoords(1_comp), nodeCoords(2_comp));
      }
    }

    void bounding_boxes(std::vector<BoundingBox>& boxes) override
    {
      int myrank = stk::parallel_machine_rank(m_sharedComm);
      std::vector<stk::mesh::Entity> nodes;
      stk::mesh::get_entities(*m_mesh, stk::topology::NODE_RANK, *m_field, nodes);
      auto coord_field = m_mesh->mesh_meta_data().coordinate_field()->data<double>();

      boxes.clear();
      for (stk::mesh::Entity node : nodes)
      {
        auto node_coords = coord_field.entity_values(node);
        BoundingShape shape(node_coords(0_comp), node_coords(1_comp), node_coords(2_comp));
        EntityProc identProc(m_mesh->entity_key(node), myrank);
        boxes.emplace_back(shape, identProc);
      }
    }

    // These are not required by the ReducedDependencyGeometricTransfer,
    // but are used in the implementation of the Interpolate class    

    int get_values_per_entity() const
    {
      std::vector<stk::mesh::Bucket*> buckets = m_mesh->get_buckets(stk::topology::NODE_RANK, *m_field);
      STK_ThrowRequireMsg(buckets.size() > 0, "proc does not have sideset entities");
      return stk::mesh::field_scalars_per_entity(*m_field, *(buckets[0]));
    }

    void write_values_to_mesh(const EntityProcVec& toEntityKeysMasked, const std::vector<double>& recvData)
    {
      auto fieldData = m_field->data<double, stk::mesh::ReadWrite>();
      size_t idx = 0;
      for (size_t i=0; i < toEntityKeysMasked.size(); ++i)
      {
        stk::mesh::Entity node = m_mesh->get_entity(toEntityKeysMasked[i].id());
        auto entityValues = fieldData.entity_values(node);
        for (stk::mesh::ScalarIdx j : entityValues.scalars())
        {
          entityValues(j) = recvData[idx++];
        }
      }
    }

  protected:
    std::shared_ptr<stk::mesh::BulkData> m_mesh;
    stk::mesh::FieldBase* m_field;
    stk::ParallelMachine m_sharedComm;
    bool m_meshModifiedSinceLastSearch = false;
    int m_numSearches = 0;
};

class ToMeshWithRepeatSearch : public ToMesh
{
  public:
    ToMeshWithRepeatSearch(std::shared_ptr<stk::mesh::BulkData>& mesh, stk::mesh::FieldBase* field,
           stk::ParallelMachine sharedComm) :
      ToMesh(mesh, field, sharedComm)
    {}

    // Required ReducedDependencyGeometricTransfer API

    bool need_repeat_search() { return m_meshModifiedSinceLastSearch; }

    void bounding_boxes(std::vector<BoundingBox>& boxes) override
    {
      ToMesh::bounding_boxes(boxes);

      m_meshModifiedSinceLastSearch = false;
      m_numSearches++;      
    }

    // used for testing

    void modify_mesh(std::shared_ptr<stk::mesh::BulkData>& mesh, stk::mesh::FieldBase* field)
    {
      m_mesh    = mesh;
      m_field   = field;
      m_meshModifiedSinceLastSearch = true;
    }

    int get_num_searches() const { return m_numSearches; }    

  private:
    bool m_meshModifiedSinceLastSearch = false;
    int m_numSearches = 0;  
};


template <typename FromMeshType, typename ToMeshType>
class InterpolateInterface
{
  public:
    // these typedefs are required to exist by ReducedDependencyGeometricTransfer
    using MeshA                 = FromMeshType;
    using MeshB                 = ToMeshType;
    using EntityKeyA            = typename MeshA::EntityKey;
    using EntityKeyB            = typename MeshB::EntityKey;
    using EntityProcA           = typename MeshA::EntityProcVec::value_type;
    using EntityProcB           = typename MeshB::EntityProcVec::value_type;
    using EntityProcRelation    = std::pair<EntityProcB, EntityProcA>;  // Note: B is first
    using EntityProcRelationVec = std::vector<EntityProcRelation>;
    using BoundingBoxA          = typename MeshA::BoundingBox;
    using BoundingBoxB          = typename MeshB::BoundingBox;

    // these typedefs are not required by the ReducedDependencyGeometricTransfer, but they
    // are used in function signatures below, so its easier to have them here
    using EntityProcVecA            = typename MeshA::EntityProcVec;
    using EntityProcVecB            = typename MeshB::EntityProcVec;
    using ToPointsContainer         = typename MeshB::ToPointsContainer;
    using ToPointsDistanceContainer = typename MeshB::ToPointsDistanceContainer;

    virtual ~InterpolateInterface() = default;

    virtual void obtain_parametric_coords(const EntityProcVecA& fromEntities, const MeshA& meshA,
                                          const ToPointsContainer& toPointsOnFromMesh,
                                          ToPointsDistanceContainer& toPointsDistanceOnFromMesh) = 0;

    virtual void mask_parametric_coords(const std::vector<int>& filterMaskFrom, int fromCount) = 0;

    virtual void apply(ToMesh* toMesh, FromMesh* fromMesh,
                       const EntityProcVecB& toEntityKeysMasked,
                       const EntityProcVecA& fromEntityKeysMasked,
                       const stk::transfer::ReducedDependencyCommData& commData) = 0;

    virtual void post_coarse_search_filter(EntityProcRelationVec& BtoA, MeshA* mesha, MeshB* meshb) = 0;
};

template <typename FromMeshType, typename ToMeshType>
class Interpolate : public InterpolateInterface<FromMeshType, ToMeshType>
{
  public:
    using Base = InterpolateInterface<FromMeshType, ToMeshType>;
    using typename Base::MeshA;
    using typename Base::MeshB;
    using typename Base::EntityProcRelationVec;

    using EntityProcVecA            = typename MeshA::EntityProcVec;
    using EntityProcVecB            = typename MeshB::EntityProcVec;
    using ToPointsContainer         = typename MeshB::ToPointsContainer;
    using ToPointsDistanceContainer = typename MeshB::ToPointsDistanceContainer;

    void obtain_parametric_coords(const EntityProcVecA& fromEntities, const MeshA& meshA,
                                  const ToPointsContainer& toPointsOnFromMesh,
                                  ToPointsDistanceContainer& toPointsDistanceOnFromMesh) override
    {
      STK_ThrowRequireMsg(fromEntities.size() == toPointsOnFromMesh.size(),
                          "size mismatch between fromEntities and toPointsOnFromMesh");
      STK_ThrowRequireMsg(fromEntities.size() == toPointsDistanceOnFromMesh.size(),
                          "size mismatch between fromEntities and toPointsDistanceOnFromMesh");

      m_toPointParametricCoordsOnFromMesh.resize(fromEntities.size());
      toPointsDistanceOnFromMesh.resize(fromEntities.size());
      for (size_t i=0; i < fromEntities.size(); ++i)
      {
        std::array<double, 3> pt = {toPointsOnFromMesh[i][0], toPointsOnFromMesh[i][1], toPointsOnFromMesh[i][2]};
        auto [xi, eta, dist] = meshA.get_parametric_coords_and_distance(fromEntities[i].id(), pt);
        m_toPointParametricCoordsOnFromMesh[i] = {xi, eta};
        toPointsDistanceOnFromMesh[i] = dist;
      }
    }

    void mask_parametric_coords(const std::vector<int>& filterMaskFrom, int fromCount) override
    {
      STK_ThrowRequireMsg(filterMaskFrom.size() == m_toPointParametricCoordsOnFromMesh.size(),
                          "filterMaskFrom size does not match size of parametric coords previously cached");

      std::vector<Point2> parametricCoordsMasked(fromCount);
      size_t outputIdx = 0;
      for (size_t i=0; i < filterMaskFrom.size(); ++i)
      {
        if (filterMaskFrom[i])
        {
          parametricCoordsMasked[outputIdx++] = m_toPointParametricCoordsOnFromMesh[i];
        }
      }

      m_toPointParametricCoordsOnFromMesh = parametricCoordsMasked;
    }

    void apply(ToMesh* toMesh, FromMesh* fromMesh,
               const EntityProcVecB& toEntityKeysMasked,
               const EntityProcVecA& fromEntityKeysMasked,
               const stk::transfer::ReducedDependencyCommData& commData) override
    {
      std::vector<double> sendData;
      std::vector<double> recvData;
      int valuesPerEntity = -1;
      if (fromMesh)
      {
        sendData = fromMesh->interpolate_to_points(fromEntityKeysMasked, m_toPointParametricCoordsOnFromMesh);
        valuesPerEntity = sendData.size() / fromEntityKeysMasked.size();
      }

      if (toMesh)
      {
        valuesPerEntity = toMesh->get_values_per_entity();
      }

      STK_ThrowRequireMsg(valuesPerEntity >= 0, "Every process must have at least one of the meshes");
      recvData.resize(toEntityKeysMasked.size()*valuesPerEntity);

      stk::transfer::do_communication(commData, sendData, recvData, valuesPerEntity);

      if (toMesh)
      {
        toMesh->write_values_to_mesh(toEntityKeysMasked, recvData);
      }
    }

    void post_coarse_search_filter(EntityProcRelationVec& BtoA, MeshA* mesha, MeshB* meshb) override {}

  private:
    using Point2 = typename ReferenceElementQuad::Point2;
    std::vector<Point2> m_toPointParametricCoordsOnFromMesh;
};


template <typename FromMeshType, typename ToMeshType>
class ReducedDependencyTransferFixtureT : public ::testing::Test
{
  public:

    void setup_spmd(const MeshAndFieldSpec& fromMeshSpec, const MeshAndFieldSpec& toMeshSpec)
    {

      localComm = stk::parallel_machine_world();
      std::tie(fromBulk, fromField) = create_from_mesh(fromMeshSpec, localComm);
      fromMesh = std::make_shared<FromMeshType>(fromBulk, fromField, localComm);  

      std::tie(toBulk, toField) = create_to_mesh(toMeshSpec, localComm);
      toMesh = std::make_shared<ToMeshType>(toBulk, toField, localComm);      

      create_transfer(localComm);
    }

    void setup_mpmd(int numFromProcs, int numToProcs, const MeshAndFieldSpec& fromMeshSpec, const MeshAndFieldSpec& toMeshSpec)
    {
      stk::ParallelMachine commWorld = stk::parallel_machine_world();
      STK_ThrowRequireMsg(stk::parallel_machine_size(commWorld) == numFromProcs + numToProcs, "Incorrect number of procs specified");
      int color = stk::parallel_machine_rank(commWorld) < numFromProcs ? 0 : 1;

      comms = std::make_shared<stk::coupling::SplitComms>(commWorld, color);
      comms->set_free_comms_in_destructor(true);
      localComm  = comms->get_split_comm();
      stk::ParallelMachine sharedComm = color == 0 ? comms->get_pairwise_comm(1) : comms->get_pairwise_comm(0);

      if (color == 0)
      {
        std::tie(fromBulk, fromField) = create_from_mesh(fromMeshSpec, localComm);
        fromMesh = std::make_shared<FromMeshType>(fromBulk, fromField, sharedComm);  
      }

      if (color == 1)
      {
        std::tie(toBulk, toField) = create_to_mesh(toMeshSpec, localComm);
        toMesh = std::make_shared<ToMeshType>(toBulk, toField, sharedComm);      
      }

      create_transfer(sharedComm);
    }

    void modify_from_mesh(const MeshAndFieldSpec& spec)
    {
      if constexpr (std::is_same_v<FromMeshType, FromMeshWithRepeatSearch>)
      {
        if (fromMesh)
        {
          create_from_mesh(spec, localComm);
          fromMesh->modify_mesh(fromBulk, fromField);
        }
      } else
      {
        STK_ThrowRequireMsg(false, "FromMeshType does not support mesh modification");
      }
    }

    void modify_to_mesh(const MeshAndFieldSpec& spec)
    {
      if constexpr (std::is_same_v<ToMeshType, ToMeshWithRepeatSearch>)
      {
        if (toMesh)
        {
          create_to_mesh(spec, localComm);
          toMesh->modify_mesh(toBulk, toField);
        }
      } else
      {
        STK_ThrowRequireMsg(false, "ToMeshType does not support mesh modification");
      }
    }

    void run_and_check_transfer(const MeshAndFieldSpec& spec)
    {
      if (fromBulk)
      {
        set_face_field(*fromBulk, fromField, spec);
      }

      transfer->apply();

      if (toBulk)
      {
        check_nodal_field(*toBulk, toField, spec);
      }
    }

    std::shared_ptr<stk::mesh::BulkData> fromBulk;
    stk::mesh::FieldBase* fromField;

    std::shared_ptr<stk::mesh::BulkData> toBulk;
    stk::mesh::FieldBase* toField;    

    std::shared_ptr<FromMeshType> fromMesh;
    std::shared_ptr<ToMeshType> toMesh;
    std::shared_ptr<Interpolate<FromMeshType, ToMeshType>> interpolator;
    std::shared_ptr<stk::transfer::ReducedDependencyGeometricTransfer<Interpolate<FromMeshType, ToMeshType>>> transfer;
    std::shared_ptr<stk::coupling::SplitComms> comms;
    stk::ParallelMachine localComm;

  private:

    std::pair<std::shared_ptr<stk::mesh::BulkData>, stk::mesh::FieldBase*>
    create_from_mesh(const MeshAndFieldSpec& spec, stk::ParallelMachine localCommArg)
    {
      std::array<double, 8> initValue = {0, 0, 0, 0, 0, 0, 0, 0};
      std::string meshString = "generated:" + std::to_string(spec.nx) + "x" + std::to_string(spec.ny) +
                                        "x" + std::to_string(spec.nz) + "|bbox:-1,0,0,0,1,1|sideset:X";

      std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(localCommArg).set_spatial_dimension(3).set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA).create();
      stk::mesh::FieldBase* field = &(bulk->mesh_meta_data().declare_field<double>(stk::topology::FACE_RANK, "fromField"));
      stk::io::fill_mesh(meshString, *bulk);
      bulk->mesh_meta_data().enable_late_fields();
      stk::mesh::Part* fromSideset = bulk->mesh_meta_data().get_part("surface_1");
      stk::mesh::put_field_on_mesh(*field, *fromSideset, 2, 4, initValue.data());

      return {bulk, field};
    }

    std::pair<std::shared_ptr<stk::mesh::BulkData>, stk::mesh::FieldBase*>
    create_to_mesh(const MeshAndFieldSpec& spec, stk::ParallelMachine localCommArg)
    {
      std::array<double, 2> initValue = {0, 0};
      std::string meshString = "generated:" + std::to_string(spec.nx) + "x" + std::to_string(spec.ny) +
                                        "x" + std::to_string(spec.nz) + "|bbox:0,0,0,1,1,1|sideset:x";

      std::shared_ptr<stk::mesh::BulkData> bulk  = stk::mesh::MeshBuilder(localCommArg).set_spatial_dimension(3).set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA).create();
      stk::mesh::FieldBase* field = &(bulk->mesh_meta_data().declare_field<double>(stk::topology::NODE_RANK, "toField"));
      stk::io::fill_mesh(meshString, *bulk);
      bulk->mesh_meta_data().enable_late_fields();
      stk::mesh::Part* toSideset = bulk->mesh_meta_data().get_part("surface_1");
      stk::mesh::put_field_on_mesh(*field, *toSideset, 2, initValue.data()); 

      return {bulk, field};
    }  

    void create_transfer(stk::ParallelMachine sharedComm)
    {
      interpolator = std::make_shared<Interpolate<FromMeshType, ToMeshType>>();
      transfer = std::make_shared<stk::transfer::ReducedDependencyGeometricTransfer<Interpolate<FromMeshType, ToMeshType>>>(fromMesh, toMesh, "my_transfer", sharedComm, *interpolator);
      transfer->initialize();
    }
};

using ReducedDependencyTransferFixture = ReducedDependencyTransferFixtureT<FromMeshWithRepeatSearch, ToMeshWithRepeatSearch>;
using ReducedDependencyTransferFixtureFromMeshMod = ReducedDependencyTransferFixtureT<FromMeshWithRepeatSearch, ToMesh>;
using ReducedDependencyTransferFixtureToMeshMod = ReducedDependencyTransferFixtureT<FromMesh, ToMeshWithRepeatSearch>;

}

TEST_F(ReducedDependencyTransferFixture, SPMD)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }

  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_spmd(fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);
}


TEST_F(ReducedDependencyTransferFixture, SPMDMultipleTransfers)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }  
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_spmd(fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  fromMeshSpec.field_value_offset += 1;
  run_and_check_transfer(fromMeshSpec);
}


TEST_F(ReducedDependencyTransferFixture, 4ProcMPMD)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4)
  {
    GTEST_SKIP();
  }
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_mpmd(2, 2, fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);
}

TEST_F(ReducedDependencyTransferFixture, 3To2ProcMPMD)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 5)
  {
    GTEST_SKIP();
  }
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_mpmd(3, 2, fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);
}

namespace {

class ClassWithNeedRepeatSearch
{
  public:
    bool need_repeat_search() { return true; }
};

class ClassWithoutNeedRepeatSearch
{
  public:
    bool foo() { return true; }
};

class BaseClassWithNeedRepeatSearch
{
  public:
    virtual ~BaseClassWithNeedRepeatSearch() = default;

    bool need_repeat_search() { return true; }
};

class ClassWithInheritedNeedRepeatSearch : public BaseClassWithNeedRepeatSearch
{
  public:
    void foo() {}
};
}

TEST(ReducedDependencyGeometricTransfer, NeedRepeatSearchTrait)
{
  static_assert(stk::transfer::impl::HasNeedRepeatSearch<ClassWithNeedRepeatSearch>::value, "class must have need_repeat_search member function");
  static_assert(!stk::transfer::impl::HasNeedRepeatSearch<ClassWithoutNeedRepeatSearch>::value, "class must not have need_repeat_search member function");
  static_assert(stk::transfer::impl::HasNeedRepeatSearch<ClassWithInheritedNeedRepeatSearch>::value, "class must have need_repeat_search member function");
}


TEST_F(ReducedDependencyTransferFixture, SPMDWithFromMeshMod)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3 ||
      stk::util::get_common_coupling_version() < 19)
  {
    GTEST_SKIP();
  }

  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_spmd(fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);
  EXPECT_EQ(fromMesh->get_num_searches(), 1);
  EXPECT_EQ(toMesh->get_num_searches(),   1);

  MeshAndFieldSpec fromMeshSpecMod{5, 5, 5, 2.0};
  modify_from_mesh(fromMeshSpecMod);
  run_and_check_transfer(fromMeshSpecMod);
  EXPECT_EQ(fromMesh->get_num_searches(), 2);
  EXPECT_EQ(toMesh->get_num_searches(),   2);

  run_and_check_transfer(fromMeshSpecMod);
  EXPECT_EQ(fromMesh->get_num_searches(), 2);
  EXPECT_EQ(toMesh->get_num_searches(),   2);
}

TEST_F(ReducedDependencyTransferFixture, SPMDWithToMeshMod)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3 ||
      stk::util::get_common_coupling_version() < 19)
  {
    GTEST_SKIP();
  }

  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_spmd(fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);
  EXPECT_EQ(fromMesh->get_num_searches(), 1);
  EXPECT_EQ(toMesh->get_num_searches(),   1);

  MeshAndFieldSpec toMeshSpecMod{5, 5, 5};
  modify_to_mesh(toMeshSpecMod);
  run_and_check_transfer(fromMeshSpec);
  EXPECT_EQ(fromMesh->get_num_searches(), 2);
  EXPECT_EQ(toMesh->get_num_searches(),   2);

  run_and_check_transfer(fromMeshSpec);
  EXPECT_EQ(fromMesh->get_num_searches(), 2);
  EXPECT_EQ(toMesh->get_num_searches(),   2);
}

TEST_F(ReducedDependencyTransferFixture, 4ProcMPMDWithFromMeshMod)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4||
      stk::util::get_common_coupling_version() < 19)
  {
    GTEST_SKIP();
  }
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_mpmd(2, 2, fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);
  if (fromMesh)
  {
    EXPECT_EQ(fromMesh->get_num_searches(), 1);
  }
  if (toMesh)
  {
    EXPECT_EQ(toMesh->get_num_searches(),   1);
  }

  MeshAndFieldSpec fromMeshSpecMod{5, 5, 5, 2.0};
  modify_from_mesh(fromMeshSpecMod);
  run_and_check_transfer(fromMeshSpecMod);
  if (fromMesh)
  {
    EXPECT_EQ(fromMesh->get_num_searches(), 2);
  }
  if (toMesh)
  {
    EXPECT_EQ(toMesh->get_num_searches(),   2);
  }

  run_and_check_transfer(fromMeshSpecMod);
  if (fromMesh)
  {
    EXPECT_EQ(fromMesh->get_num_searches(), 2);
  }
  if (toMesh)
  {
    EXPECT_EQ(toMesh->get_num_searches(),   2);
  }
}

TEST_F(ReducedDependencyTransferFixture, 4ProcMPMDWithToMeshMod)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4||
      stk::util::get_common_coupling_version() < 19)
  {
    GTEST_SKIP();
  }
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_mpmd(2, 2, fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);
  if (fromMesh)
  {
    EXPECT_EQ(fromMesh->get_num_searches(), 1);
  }
  if (toMesh)
  {
    EXPECT_EQ(toMesh->get_num_searches(),   1);
  }

  MeshAndFieldSpec toMeshSpecMod{5, 5, 5};
  modify_to_mesh(toMeshSpecMod);
  run_and_check_transfer(fromMeshSpec);
  if (fromMesh)
  {
    EXPECT_EQ(fromMesh->get_num_searches(), 2);
  }
  if (toMesh)
  {
    EXPECT_EQ(toMesh->get_num_searches(),   2);
  }

  run_and_check_transfer(fromMeshSpec);
  if (fromMesh)
  {
    EXPECT_EQ(fromMesh->get_num_searches(), 2);
  }
  if (toMesh)
  {
    EXPECT_EQ(toMesh->get_num_searches(),   2);
  }
}


//-----------------------------------------------------------------------------
// Tests where only the FromMesh supports mesh modification

TEST_F(ReducedDependencyTransferFixtureFromMeshMod, SPMDMultipleTransfers)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }  
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_spmd(fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  fromMeshSpec.field_value_offset += 1;
  run_and_check_transfer(fromMeshSpec);
}

TEST_F(ReducedDependencyTransferFixtureFromMeshMod, SPMDWithMeshModThrow)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }  
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_spmd(fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  MeshAndFieldSpec fromMeshSpecMod{5, 5, 5, 2.0};
  modify_from_mesh(fromMeshSpecMod);
  EXPECT_ANY_THROW(run_and_check_transfer(fromMeshSpecMod));
}


TEST_F(ReducedDependencyTransferFixtureFromMeshMod, 4ProcMPMD)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4)
  {
    GTEST_SKIP();
  }
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_mpmd(2, 2, fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  fromMeshSpec.field_value_offset += 1;
  run_and_check_transfer(fromMeshSpec);  
}

TEST_F(ReducedDependencyTransferFixtureFromMeshMod, 4ProcMPMDMeshModThrow)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4)
  {
    GTEST_SKIP();
  }
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_mpmd(2, 2, fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  MeshAndFieldSpec fromMeshSpecMod{5, 5, 5, 2.0};
  modify_from_mesh(fromMeshSpecMod);
  EXPECT_ANY_THROW(run_and_check_transfer(fromMeshSpec)); 
}


//-----------------------------------------------------------------------------
// Tests where only the ToMesh supports mesh modification

TEST_F(ReducedDependencyTransferFixtureToMeshMod, SPMDMultipleTransfers)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }  
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_spmd(fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  fromMeshSpec.field_value_offset += 1;
  run_and_check_transfer(fromMeshSpec);
}

TEST_F(ReducedDependencyTransferFixtureToMeshMod, SPMDWithMeshModThrow)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) > 3)
  {
    GTEST_SKIP();
  }  
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_spmd(fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  MeshAndFieldSpec toMeshSpecMod{5, 5, 5};
  modify_to_mesh(toMeshSpecMod);
  EXPECT_ANY_THROW(run_and_check_transfer(fromMeshSpec));
}


TEST_F(ReducedDependencyTransferFixtureToMeshMod, 4ProcMPMD)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4)
  {
    GTEST_SKIP();
  }
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_mpmd(2, 2, fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  fromMeshSpec.field_value_offset += 1;
  run_and_check_transfer(fromMeshSpec);  
}

TEST_F(ReducedDependencyTransferFixtureToMeshMod, 4ProcMPMDMeshModThrow)
{
  if (stk::parallel_machine_size(stk::parallel_machine_world()) != 4)
  {
    GTEST_SKIP();
  }
  MeshAndFieldSpec fromMeshSpec{3, 3, 3, 1.0};
  MeshAndFieldSpec toMeshSpec{4, 4, 4};

  setup_mpmd(2, 2, fromMeshSpec, toMeshSpec);
  run_and_check_transfer(fromMeshSpec);

  MeshAndFieldSpec toMeshSpecMod{5, 5, 5,};
  modify_to_mesh(toMeshSpecMod);
  EXPECT_ANY_THROW(run_and_check_transfer(fromMeshSpec)); 
}
