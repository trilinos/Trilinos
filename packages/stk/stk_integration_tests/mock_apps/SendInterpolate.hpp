/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef MOCK_APPS_SEND_INTERPOLATE_HPP
#define MOCK_APPS_SEND_INTERPOLATE_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_transfer/GeometricTransfer.hpp>
#include <stk_transfer/ReducedDependencyGeometricTransfer.hpp>
#include <stk_transfer/TransferBase.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <memory>
#include <utility>

namespace mock {

template<typename SendAdapter, typename RecvAdapter>
class SendInterpolate
{
public:
  using MeshA = SendAdapter;
  using MeshB = RecvAdapter;
  using EntityKeyA = typename MeshA::EntityKey;
  using EntityKeyB = typename MeshB::EntityKey;
  using EntityProcA = typename MeshA::EntityProc;
  using EntityProcB = typename MeshB::EntityProc;
  using EntityKeyMap = std::multimap<EntityKeyB, EntityKeyA>;
  using EntityProcRelation = std::pair<EntityProcB, EntityProcA>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;
  static void filter_to_nearest(EntityKeyMap & /*local_range_to_domain*/, const MeshA & /*mesha*/, const MeshB & /*meshb*/)
  {
    //no filtering needed since map is one-to-one
  }

  ~SendInterpolate()
  {
  }

  //Specific to single point case right now
  void obtain_parametric_coords(typename MeshA::EntityProcVec entities_to_copy_from,
      MeshA &/*sendAdapter*/,
      const typename MeshB::ToPointsContainer & /*to_points_on_from_mesh*/,
      typename MeshB::ToPointsDistanceContainer & to_points_distance_on_from_mesh)
  {
    for (unsigned i = 0; i < entities_to_copy_from.size(); ++i)
    {
      to_points_distance_on_from_mesh.push_back(0.0);
    }
  }

  void mask_parametric_coords(const std::vector<int> & /*filter_mask_from*/, int /*from_count*/)
  {
  }

  void
  apply(MeshB * /*recvAdapter*/,
      MeshA * sendAdapter,
      const typename MeshB::EntityProcVec & to_entity_keys_masked,
      const typename MeshA::EntityProcVec & from_entity_keys_masked,
      const stk::transfer::ReducedDependencyCommData & comm_data)
  {
    std::vector<double> tempResults(from_entity_keys_masked.size());
    std::vector<double> tempResultsReceived(to_entity_keys_masked.size());
    size_t fieldSizePerEntity = 1; //hard-coded for mock-apps case of 1 scalar field.

    STK_ThrowRequire(sendAdapter != nullptr);
    load_from_mesh_data(tempResults, from_entity_keys_masked, *sendAdapter);

    stk::transfer::do_communication(comm_data, tempResults, tempResultsReceived, fieldSizePerEntity);
  }

  void load_from_mesh_data(std::vector<double> & tempResults,
      const typename MeshA::EntityProcVec & from_entity_keys_masked,
      MeshA & sendAdapter)
  {
    for(unsigned i=0; i<from_entity_keys_masked.size(); ++i) {
      typename MeshA::EntityKey key = from_entity_keys_masked[i].id();
      tempResults[i] = sendAdapter.get_field_value(key);
    }
  }
};

} // namespace mock

#endif // MOCK_APPS_SEND_INTERPOLATE_HPP
