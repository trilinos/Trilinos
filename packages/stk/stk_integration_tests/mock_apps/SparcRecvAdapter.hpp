/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef MOCK_APPS_SPARC_RECV_ADAPTER_HPP
#define MOCK_APPS_SPARC_RECV_ADAPTER_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include "SparcMesh.hpp"
#include <memory>
#include <utility>

namespace mock {

class SparcRecvAdapter
{
public:
  SparcRecvAdapter(MPI_Comm mpiComm, SparcMesh& mockMesh, const std::string & fieldName)
  : m_comm(mpiComm),
    m_mesh(mockMesh),
    m_fieldName(fieldName)
  {
  }
  using EntityKey = SparcMesh::EntityKey;
  using EntityProc = SparcMesh::EntityProc;
  using EntityProcVec = SparcMesh::EntityProcVec;
  using BoundingBox = SparcMesh::BoundingBox;

  //Used for Reduced dependency
  using Point = stk::search::Point<double>;
  using ToPointsContainer = std::vector<Point>;
  using ToPointsDistance = double;
  using ToPointsDistanceContainer = std::vector<ToPointsDistance>;

  MPI_Comm comm() const {return m_comm;}

  void update_values() {called_update_values = true;}

  bool called_update_values = false;

  void set_field_value(const EntityKey & entityKey, const double & value)
  {
    m_mesh.set_sparc_field_value(entityKey, m_fieldName, value);
  }

  void bounding_boxes(std::vector<BoundingBox> & range_vector) const
  {
    m_mesh.sparc_dest_bounding_boxes(range_vector);
    int procInSearchComm = stk::parallel_machine_rank(m_comm);
    for(BoundingBox& box : range_vector) {
      box.second.set_proc(procInSearchComm);
    }
  }

  void get_to_points_coordinates(const EntityProcVec &to_entity_keys, ToPointsContainer &to_points)
  {
    m_mesh.get_to_points_coordinates(to_entity_keys, to_points);
  }

private:
  MPI_Comm m_comm;
  SparcMesh& m_mesh;
  std::string m_fieldName;
};

}

#endif // MOCK_APPS_SPARC_RECV_ADAPTER_HPP
