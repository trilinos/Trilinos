/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef MOCK_APPS_STK_MESH_HPP
#define MOCK_APPS_STK_MESH_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <memory>
#include <utility>

namespace mock {

class StkMesh
{
public:
  using Point = stk::search::Point<double>;
  using ToPointsContainer = std::vector<Point>;
  using EntityKey = stk::mesh::EntityKey;
  using EntityProc = stk::search::IdentProc<EntityKey>;
  using EntityProcVec = std::vector<EntityProc>;
  using BoundingBox = std::pair<stk::search::Box<double>, EntityProc>;
  using BoundingSphere = std::pair<stk::search::Sphere<double>, EntityProc>;

  StkMesh(std::shared_ptr<stk::mesh::BulkData> bulk, const std::string& surfaceName)
  : m_comm(bulk->parallel()),
    m_bulk(bulk),
    m_surface(bulk->mesh_meta_data().get_part(surfaceName))
  {
    STK_ThrowRequireMsg(m_surface != nullptr,
                    "StkMesh error, no part found for surface-name '"<<surfaceName<<"'");
  }

  MPI_Comm comm() const {return m_comm;}

  std::shared_ptr<stk::mesh::BulkData> get_stk_mesh() { return m_bulk; }

  double get_stk_field_value(const EntityKey & entityKey,
                         const std::string & fieldName)
  {
    stk::mesh::EntityRank rank = entityKey.rank();
    stk::mesh::Entity entity = m_bulk->get_entity(entityKey);
    const stk::mesh::FieldBase* field = m_bulk->mesh_meta_data().get_field(rank, fieldName);
    if (field->defined_on(entity)) {
      auto fieldData = field->data<double>();
      auto fieldValue = fieldData.entity_values(entity);
      return fieldValue(0_comp);
    }
    return -99.9;
  }

  void set_stk_field_value(const EntityKey & entityKey,
                       const std::string & fieldName,
                       const double & fieldValue)
  {
    stk::mesh::EntityRank rank = entityKey.rank();
    stk::mesh::Entity entity = m_bulk->get_entity(entityKey);
    const stk::mesh::FieldBase* field = m_bulk->mesh_meta_data().get_field(rank, fieldName);
    if (field->defined_on(entity)) {
      auto fieldData = field->data<double, stk::mesh::ReadWrite>();
      auto entityFieldData = fieldData.entity_values(entity);
      entityFieldData(0_comp) = fieldValue;
    }
  }

  void set_stk_field_values(const std::string & fieldName,
                            const double & fieldValue)
  {
    const stk::mesh::FieldBase* field = stk::mesh::get_field_by_name(fieldName, m_bulk->mesh_meta_data());
    STK_ThrowRequireMsg(field != nullptr, "StkMesh::set_stk_field_values failed to find field with name "<<fieldName);
    stk::mesh::field_fill(fieldValue, *field);
  }

  bool verify_stk_field_values(const std::string & fieldName,
                               const double expectedFieldValue)
  {
    const stk::mesh::FieldBase* field = stk::mesh::get_field_by_name(fieldName, m_bulk->mesh_meta_data());
    STK_ThrowRequireMsg(field != nullptr, "StkMesh::verify_stk_field_values failed to find field with name "<<fieldName);
    stk::mesh::Selector fieldSelector(*field);
    stk::mesh::EntityVector sides;
    stk::mesh::get_entities(*m_bulk, stk::topology::FACE_RANK, fieldSelector, sides);
    auto fieldData = field->data<double>();
    for(stk::mesh::Entity side : sides) {
      auto sideFieldData = fieldData.entity_values(side);
      if (std::abs(sideFieldData(0_comp) - expectedFieldValue) > 1.e-6) {
        return false;
      }
    }
    return true;
  }

  unsigned get_field_size() const { return 1; }

  stk::search::Box<double> get_box(stk::mesh::Entity side) const
  {
    constexpr double maxDouble = std::numeric_limits<double>::max();
    double minXYZ[3] = {maxDouble, maxDouble, maxDouble};
    double maxXYZ[3] = {0.0, 0.0, 0.0};
    const stk::mesh::FieldBase* coordField = m_bulk->mesh_meta_data().coordinate_field();
    auto coordFieldData = coordField->data<double>();

    const stk::mesh::Entity* nodes = m_bulk->begin_nodes(side);
    const unsigned numNodes = m_bulk->num_nodes(side);
    for(unsigned i=0; i<numNodes; ++i) {
      auto coords = coordFieldData.entity_values(nodes[i]);
      minXYZ[0] = std::min(minXYZ[0], coords(0_comp));
      minXYZ[1] = std::min(minXYZ[1], coords(1_comp));
      minXYZ[2] = std::min(minXYZ[2], coords(2_comp));
      maxXYZ[0] = std::max(maxXYZ[0], coords(0_comp));
      maxXYZ[1] = std::max(maxXYZ[1], coords(1_comp));
      maxXYZ[2] = std::max(maxXYZ[2], coords(2_comp));
    }

    constexpr double tol = 1.e-5;
    return stk::search::Box<double>(minXYZ[0]-tol, minXYZ[1]-tol, minXYZ[2]-tol,
                                    maxXYZ[0]+tol, maxXYZ[1]+tol, maxXYZ[2]+tol);
  }

  void stk_source_bounding_boxes(std::vector<BoundingBox> & domain_vector) const
  {
    stk::mesh::EntityVector sides;
    stk::mesh::Selector ownedSurface = m_bulk->mesh_meta_data().locally_owned_part() & *m_surface;
    stk::mesh::get_entities(*m_bulk, m_bulk->mesh_meta_data().side_rank(), ownedSurface, sides);
    domain_vector.clear();
    const int thisProc = m_bulk->parallel_rank();
    for(stk::mesh::Entity side : sides) {
      EntityProc entityProc(m_bulk->entity_key(side), thisProc);
      domain_vector.emplace_back(get_box(side), entityProc);
    }
  }

  Point get_centroid(stk::mesh::Entity side) const
  {
    double sumXYZ[3] = {0.0, 0.0, 0.0};
    const stk::mesh::FieldBase* coordField = m_bulk->mesh_meta_data().coordinate_field();
    auto coordFieldData = coordField->data<double>();

    const stk::mesh::Entity* nodes = m_bulk->begin_nodes(side);
    const unsigned numNodes = m_bulk->num_nodes(side);
    for(unsigned i=0; i<numNodes; ++i) {
      auto coords = coordFieldData.entity_values(nodes[i]);
      sumXYZ[0] += coords(0_comp);
      sumXYZ[1] += coords(1_comp);
      sumXYZ[2] += coords(2_comp);
    }

    return Point(sumXYZ[0]/numNodes, sumXYZ[1]/numNodes, sumXYZ[2]/numNodes);
  }

  void stk_dest_bounding_boxes(std::vector<BoundingSphere>& range_vector) const
  {
    stk::mesh::EntityVector sides;
    stk::mesh::Selector ownedSurface = m_bulk->mesh_meta_data().locally_owned_part() & *m_surface;
    stk::mesh::get_entities(*m_bulk, m_bulk->mesh_meta_data().side_rank(), ownedSurface, sides);
    constexpr double radius = 1.e-6;
    range_vector.clear();
    const int thisProc = m_bulk->parallel_rank();
    for(stk::mesh::Entity side : sides) {
      EntityProc entityProc(m_bulk->entity_key(side), thisProc);
      range_vector.emplace_back(stk::search::Sphere<double>(get_centroid(side), radius), entityProc);
    }
  }

  void get_to_points_coordinates(const EntityProcVec & /*to_entity_keys*/, ToPointsContainer &to_points)
  {
    stk::mesh::EntityVector sides;
    stk::mesh::Selector ownedSurface = m_bulk->mesh_meta_data().locally_owned_part() & *m_surface;
    stk::mesh::get_entities(*m_bulk, m_bulk->mesh_meta_data().side_rank(), ownedSurface, sides);
    to_points.clear();
    for(stk::mesh::Entity side : sides) {
      to_points.push_back(get_centroid(side));
    }
  }

private:
  MPI_Comm m_comm;
  std::shared_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::Part* m_surface;
};

}

#endif // MOCK_APPS_STK_MESH_HPP
