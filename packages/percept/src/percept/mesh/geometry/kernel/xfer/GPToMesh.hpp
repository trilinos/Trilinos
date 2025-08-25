// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef GPToMesh_h
#define GPToMesh_h

#include <percept/PerceptMesh.hpp>
#include <percept/xfer/ToMesh.hpp>
#include "GPFromMesh.hpp"

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

namespace percept {

  class GPToMeshAdaptor;

  class GPToMesh : public ToMesh<GPToMeshAdaptor> {
  public :

    typedef ToMesh<GPToMeshAdaptor> Base;

    GPToMesh(PerceptMesh& eMesh, const stk::mesh::PartVector& parts,
             const std::vector<stk::mesh::FieldBase *>& fields,
             TransferType transferType,
             SrcFieldType srcFieldType=SRC_FIELD,
             const double                 radius=.0001) :
      Base(*eMesh.get_bulk_data(),
           eMesh.get_coordinates_field(),
           0,
           eMesh.get_bulk_data()->parallel(),
           transferType,
           srcFieldType,
           radius),
      m_eMesh(eMesh),
      toParts_(parts)
    {
      Base::toFields_ = fields;
    }

    ~GPToMesh(){};

    PerceptMesh& m_eMesh;
    stk::mesh::PartVector         toParts_;

    typedef std::set<stk::mesh::Entity> EntitySet;
    typedef std::map<stk::mesh::Entity, EntitySet > EntityToEntitySetMap;
    EntityToEntitySetMap m_searchMap;
  };

  class GPToMeshAdaptor {
  public:

    typedef std::multimap<stk::mesh::EntityKey, stk::mesh::EntityKey> EntityKeyMap;

    // return true to skip the standard processing
    template<class FromMesh, class ToMesh>
    bool
    apply(ToMesh &ToPointsIn,
          const FromMesh  &FromElemIn,
          const EntityKeyMap &RangeToDomain)
    {
      const GPFromMesh& FromElem = static_cast<const GPFromMesh& > (FromElemIn);
      GPToMesh& ToPoints = static_cast<GPToMesh& > (ToPointsIn);

      const stk::mesh::BulkData &fromBulkData = FromElem.fromBulkData_;
      stk::mesh::BulkData         &toBulkData = ToPoints.toBulkData_;

      //stk::mesh::EntityRank fromRank = FromElem.fromField_->entity_rank();

      ToPoints.m_searchMap.clear();

      typename EntityKeyMap::const_iterator ii;
      for(ii=RangeToDomain.begin(); ii!=RangeToDomain.end(); ++ii )
        {
          const stk::mesh::EntityKey thePt  = ii->first;
          const stk::mesh::EntityKey theBox = ii->second;

          stk::mesh::Entity node =   toBulkData.get_entity(thePt);
          stk::mesh::Entity face = fromBulkData.get_entity(theBox);

          //if (GPSTK_DEBUG)
          //      std::cout << "relation[ " << i << "]= {" << relation[i].first << "} --> { " << relation[i].second << "}" << std::endl;
          VERIFY_OP_ON(toBulkData.is_valid(node), ==, true, "bad node");
          VERIFY_OP_ON(fromBulkData.is_valid(face), == , true, "bad face");

          std::set<stk::mesh::Entity>& set = ToPoints.m_searchMap[node];
          set.insert(face);
        }
      //std::cout << "ToPoints.m_searchMap.size= " << ToPoints.m_searchMap.size() << std::endl;

      return true;
    }

    // return true to skip the standard processing
    template<class FromMesh, class ToMesh>
    bool
    filter_to_nearest(const EntityKeyMap &/*RangeToDomain*/,
                      const FromMesh  &/*FromElem*/,
                      ToMesh &/*ToPoints*/)
    {
      // no-op for this specialization
      return true;
    }

    template<class ToMesh>
    void
    modify_selector(const ToMesh& ToPointsIn, stk::mesh::Selector& sel)
    {
      const GPToMesh& ToPoints = static_cast<const GPToMesh& > (ToPointsIn);
      // allow shared and ghost nodes too
      sel = stk::mesh::selectUnion(ToPoints.toParts_);
      //sel = sel & stk::mesh::selectUnion(ToPoints.toParts_);
    }

  };


} // namespace percept

#endif
