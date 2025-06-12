/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef MOCK_APPS_SPARC_MESH_HPP
#define MOCK_APPS_SPARC_MESH_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <memory>
#include <utility>

namespace mock {

constexpr static unsigned maxNodesPerSide = 4;

struct SparcSide
{
  uint64_t key;
  unsigned numNodes;
  double nodeCoords[maxNodesPerSide][3];
};

class SparcMesh
{
public:
  using Point = stk::search::Point<double>;
  using ToPointsContainer = std::vector<Point>;
  using EntityKey = uint64_t;
  using EntityProc = stk::search::IdentProc<EntityKey>;
  using EntityProcVec = std::vector<EntityProc>;
  using BoundingBox = std::pair<stk::search::Box<double>, EntityProc>;
  using BoundingSphere = std::pair<stk::search::Sphere<double>, EntityProc>;

  SparcMesh(MPI_Comm mpiComm, const std::vector<SparcSide>& inputSides)
  : m_comm(mpiComm), m_sides(inputSides)
  {
  }

  MPI_Comm comm() const {return m_comm;}

  double get_sparc_field_value(const EntityKey & entityKey,
                         const std::string & fieldName)
  {
    KeyAndField keyAndField(entityKey, fieldName);
    return get_value_from_map(keyAndField, m_entityFieldValues, -99.9);
  }

  void set_sparc_field_value(const EntityKey & entityKey,
                       const std::string & fieldName,
                       const double & fieldValue)
  {
    KeyAndField keyAndField(entityKey, fieldName);
    m_entityFieldValues[keyAndField] = fieldValue;
  }

  void set_sparc_field_values(const std::string & fieldName,
                              const double & fieldValue)
  {
    for(const SparcSide& side : m_sides) {
      KeyAndField keyAndField(side.key, fieldName);
      m_entityFieldValues[keyAndField] = fieldValue;
    }
  }

  unsigned get_field_size() const { return 1; }

  stk::search::Box<double> get_box(const SparcSide& side) const
  {
    constexpr double maxDouble = std::numeric_limits<double>::max();
    double minXYZ[3] = {maxDouble, maxDouble, maxDouble};
    double maxXYZ[3] = {0.0, 0.0, 0.0};

    const unsigned numNodes = side.numNodes;
    for(unsigned i=0; i<numNodes; ++i) {
      for(unsigned d=0; d<3; ++d) {
        minXYZ[d] = std::min(minXYZ[d], side.nodeCoords[i][d]);
        maxXYZ[d] = std::max(maxXYZ[d], side.nodeCoords[i][d]);
      }
    }

    constexpr double tol = 1.e-5;
    return stk::search::Box<double>(minXYZ[0]-tol, minXYZ[1]-tol, minXYZ[2]-tol,
                                    maxXYZ[0]+tol, maxXYZ[1]+tol, maxXYZ[2]+tol);
  }

  Point get_centroid(const SparcSide& side) const
  {
    double sumXYZ[3] = {0.0, 0.0, 0.0};

    const unsigned numNodes = side.numNodes;
    for(unsigned i=0; i<numNodes; ++i) {
      for(unsigned d=0; d<3; ++d) {
        sumXYZ[d] += side.nodeCoords[i][d];
      }
    }
    return Point(sumXYZ[0]/numNodes, sumXYZ[1]/numNodes, sumXYZ[2]/numNodes);
  }

  void sparc_source_bounding_boxes(std::vector<BoundingBox> & domain_vector) const
  {
    domain_vector.clear();
    const int thisProc = stk::parallel_machine_rank(comm());
    for(const SparcSide& side : m_sides) {
      EntityProc entityProc(side.key, thisProc);
      domain_vector.emplace_back(get_box(side), entityProc);
    }
  }

  void sparc_dest_bounding_boxes(std::vector<BoundingBox>& range_vector) const
  {
    range_vector.clear();
    const int thisProc = stk::parallel_machine_rank(comm());
    for(const SparcSide& side : m_sides) {
      EntityProc entityProc(side.key, thisProc);
      range_vector.emplace_back(get_box(side), entityProc);
    }
  }

  void get_to_points_coordinates(const EntityProcVec & /*to_entity_keys*/, ToPointsContainer &to_points)
  {
    to_points.clear();
    for(const SparcSide& side : m_sides) {
      to_points.push_back(get_centroid(side));
    }
  }

private:
  using KeyAndField = std::pair<EntityKey,std::string>;
  using FieldValues = std::map<KeyAndField,double>;
  double get_value_from_map(const KeyAndField& keyAndField, const FieldValues& fieldValues, const double& defaultIfNotFound)
  {
    FieldValues::const_iterator iter = fieldValues.find(keyAndField);
    if (iter != fieldValues.end()) {
      return iter->second;
    }
    return defaultIfNotFound;
  }

  MPI_Comm m_comm;
  FieldValues m_entityFieldValues;
  std::vector<SparcSide> m_sides;
};

}

#endif // MOCK_APPS_SPARC_MESH_HPP
