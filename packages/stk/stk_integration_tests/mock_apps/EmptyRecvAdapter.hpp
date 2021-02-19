/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef MOCK_APPS_EMPTY_RECV_ADAPTER_HPP
#define MOCK_APPS_EMPTY_RECV_ADAPTER_HPP

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <memory>
#include <utility>

namespace mock {

class EmptyRecvAdapter
{
public:
  using EntityKey = uint64_t;
  using EntityProc = stk::search::IdentProc<EntityKey>;
  using EntityProcVec = std::vector<EntityProc>;
  using BoundingBox = std::pair<stk::search::Sphere<double>, EntityProc>;

  using Point = stk::search::Point<double>;
  using ToPointsContainer = std::vector<Point>;
  using ToPointsDistance = double;
  using ToPointsDistanceContainer = std::vector<ToPointsDistance>;

  void update_values() { ThrowErrorMsg("EmptyRecvAdapter shouldn't be called."); }

  void bounding_boxes(std::vector<BoundingBox> & /* range_vector */) const
  {
    ThrowErrorMsg("EmptyRecvAdapter shouldn't be called.");
  }

  void get_to_points_coordinates(const EntityProcVec & /* to_entity_keys */, ToPointsContainer & /* to_points */)
  {
    ThrowErrorMsg("EmptyRecvAdapter shouldn't be called.");
  }
};

} // namespace mock

#endif // MOCK_APPS_EMPTY_RECV_ADAPTER_HPP
