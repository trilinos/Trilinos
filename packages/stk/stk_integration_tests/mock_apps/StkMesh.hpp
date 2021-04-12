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
#include <stk_util/parallel/Parallel.hpp>
#include <memory>
#include <utility>

namespace mock {

class StkMesh
{
public:
  using Point = stk::search::Point<double>;
  using ToPointsContainer = std::vector<Point>;
  using EntityKey = uint64_t;
  using EntityProc = stk::search::IdentProc<EntityKey>;
  using EntityProcVec = std::vector<EntityProc>;
  using BoundingBox = std::pair<stk::search::Box<double>, EntityProc>;
  using BoundingSphere = std::pair<stk::search::Sphere<double>, EntityProc>;

  StkMesh(MPI_Comm mpiComm)
  : m_comm(mpiComm), m_owning_rank(0), m_sourceEntityKey(3), m_destEntityKey(5)
  {
  }

  MPI_Comm comm() const {return m_comm;}

  int owning_rank() const { return m_owning_rank; }

  double get_stk_field_value(const EntityKey & entityKey,
                         const std::string & fieldName)
  {
    if (entityKey == m_sourceEntityKey) {
      return get_value_from_map(fieldName, m_sourceEntityFieldValues, -99.9);
    }
    if (entityKey == m_destEntityKey) {
      return get_value_from_map(fieldName, m_destEntityFieldValues, -99.9);
    }
    return -99.9;
  }

  void set_stk_field_value(const EntityKey & entityKey,
                       const std::string & fieldName,
                       const double & fieldValue)
  {
    if (entityKey == m_sourceEntityKey) {
      m_sourceEntityFieldValues[fieldName] = fieldValue;
    }
    if (entityKey == m_destEntityKey) {
      m_destEntityFieldValues[fieldName] = fieldValue;
    }
  }

  unsigned get_field_size() const { return 1; }

  EntityKey get_stk_source_entity_key() const { return m_sourceEntityKey; }
  EntityKey get_stk_dest_entity_key() const { return m_destEntityKey; }

  stk::search::Box<double> get_source_box() const { return stk::search::Box<double>(0., 0., 0., 1., 1., 1.); }
  Point get_dest_point() const {return Point(0.5, 0.5, 0.5);}

  void stk_source_bounding_boxes(std::vector<BoundingBox> & domain_vector) const
  {
    if (stk::parallel_machine_rank(m_comm) == owning_rank()) {
      EntityProc entityProc(m_sourceEntityKey, owning_rank());
      domain_vector.emplace_back(get_source_box(), entityProc);
    }
  }

  void stk_dest_bounding_boxes(std::vector<BoundingSphere>& range_vector) const
  {
    if (stk::parallel_machine_rank(m_comm) == owning_rank()) {
      EntityProc entityProc(m_destEntityKey, owning_rank());
      const double radius = 1.e-6;
      range_vector.emplace_back(stk::search::Sphere<double>(get_dest_point(), radius), entityProc);
    }
  }

  void get_to_points_coordinates(const EntityProcVec &to_entity_keys, ToPointsContainer &to_points)
  {
    if (stk::parallel_machine_rank(m_comm) == owning_rank()) {
      to_points.push_back(get_dest_point());
    }
  }

private:
  using FieldValues = std::map<std::string,double>;
  double get_value_from_map(const std::string& name, const FieldValues& fieldValues, const double& defaultIfNotFound)
  {
    FieldValues::const_iterator iter = fieldValues.find(name);
    if (iter != fieldValues.end()) {
      return iter->second;
    }
    return defaultIfNotFound;
  }

  MPI_Comm m_comm;
  int m_owning_rank = 0;
  EntityKey m_sourceEntityKey;
  EntityKey m_destEntityKey;
  FieldValues m_sourceEntityFieldValues;
  FieldValues m_destEntityFieldValues;
};

}

#endif // MOCK_APPS_STK_MESH_HPP
