/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICINTERP_HPP_
#define STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICINTERP_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <vector>
#include <map>

#include <stk_search/FilterCoarseSearch.hpp>
#include <stk_transfer/TransferTypes.hpp>
#include <stk_util/util/ReportHandler.hpp>
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace transfer {
namespace spmd {

template <class SENDMESH, class RECVMESH>
class GeometricInterp {
 public:
  using SendMesh = SENDMESH;
  using RecvMesh = RECVMESH;
  using SendEntityKey = typename SENDMESH::EntityKey;
  using RecvEntityKey = typename RECVMESH::EntityKey;
  using SendEntityProc = typename SENDMESH::EntityProc;
  using RecvEntityProc = typename RECVMESH::EntityProc;

  using MeshA = SENDMESH;
  using MeshB = RECVMESH;
  using EntityKeyA = typename SENDMESH::EntityKey;
  using EntityKeyB = typename RECVMESH::EntityKey;
  using EntityProcA = typename SENDMESH::EntityProc;
  using EntityProcB = typename RECVMESH::EntityProc;

  using EntityProcRelation = std::pair<RecvEntityProc, SendEntityProc>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;

  using EntityKeyMap = std::multimap<EntityKeyB, EntityKeyA> ;

  static void filter_to_nearest(EntityKeyMap& /*RangeToDomain*/, MeshA& /*FromMesh*/, MeshB& /*ToMesh*/)
  {
    STK_ThrowRequireMsg(false, "GeometricInterp::filter_to_nearest() not implemented");
  }

  static void apply(MeshB& /*ToMesh*/, MeshA& /*FromMesh*/, const EntityKeyMap& /*RangeToDomain*/)
  {
    STK_ThrowRequireMsg(false, "GeometricInterp::apply() wrong overload, can't take a map any more.");
  }

  static void apply(RecvMesh& recvMesh, SendMesh& sendMesh,
                    const EntityProcRelationVec& RangeToDomain,
                    const stk::search::FilterCoarseSearchResult<RecvMesh>& filterResult)
  {
    const int maxNumFields = recvMesh.max_num_values();

    InterpolationData data(maxNumFields);

    const unsigned spatialDimension = recvMesh.spatial_dimension();
    std::vector<double> coordVec(spatialDimension, 0.0);
    std::vector<double> parametricCoords;

    for(auto&& ii : RangeToDomain) {
      const RecvEntityKey recvEntity = ii.first.id();
      const SendEntityKey sendEntity = ii.second.id();

      filterResult.get_parametric_coordinates(recvEntity, parametricCoords);

      recvMesh.coordinates(recvEntity, coordVec);

      data.nFields = recvMesh.num_values(recvEntity);
      for(unsigned n = 0; n < data.nFields; ++n) {
        data.fieldPtr[n]       = recvMesh.value(recvEntity, n);
        data.fieldSize[n]      = recvMesh.value_size(recvEntity, n);
        data.fieldKey[n]       = recvMesh.value_key(recvEntity, n);
        data.fieldDataIndex[n] = recvMesh.get_index(n);
      }

      data.debug = false;
      sendMesh.interpolate_fields(sendEntity, coordVec, parametricCoords, data);
    }
  }
};

} // namespace spmd
} // namespace transfer
} // namespace stk


#endif /* STK_STK_TRANSFER_UTIL_STK_TRANSFER_UTIL_SPMD_GEOMETRICINTERP_HPP_ */

