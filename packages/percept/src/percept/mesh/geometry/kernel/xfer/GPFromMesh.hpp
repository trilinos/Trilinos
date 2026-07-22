// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef GPFromMesh_h
#define GPFromMesh_h

#include <percept/PerceptMesh.hpp>
#include <percept/xfer/FromMesh.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>

namespace percept{

  //--------------------------------------------------------------------
  /**
   * modified from Nalu:
   *
   *  @class  FromMesh
   *  @author Stefan P. Domino
   *  @date   July, 2013
   *  @brief  Class that holds source mesh for transfer.
   *
   */
  //--------------------------------------------------------------------

  class GPFromMeshAdaptor;

  class GPFromMesh : public FromMesh<GPFromMeshAdaptor> {
  public :
    typedef FromMesh<GPFromMeshAdaptor> Base;

    GPFromMesh(PerceptMesh& eMesh, const stk::mesh::PartVector& parts,
               const std::vector<stk::mesh::FieldBase *>& fields) :
      Base(*eMesh.get_bulk_data(),
           eMesh.get_coordinates_field(),
           0,
           eMesh.get_bulk_data()->parallel()),
      m_eMesh(eMesh),
      fromParts_(parts),
      m_inflationFactor(2.0)
    {
      Base::fromFields_ = fields;
    }

    ~GPFromMesh(){}

    PerceptMesh                  &m_eMesh;
    stk::mesh::PartVector         fromParts_;
    double m_inflationFactor;
  };

  class GPFromMeshAdaptor {
  public:

    template<class FM>
    static void
    get_entities(const FM& fromMeshIn, stk::mesh::EntityRank /*rank*/, std::vector<stk::mesh::Entity>& vecFaces)
    {
      // Maybe use a ptr here before casting? and then verify that the cast
      // was successful, followed by a VERIFY_OP_ON check.
      const GPFromMesh& fromMesh = static_cast<const GPFromMesh& >(fromMeshIn);
      //VERIFY_OP_ON(&fromMesh, !=, 0, "bad cast");

      //stk::mesh::Selector sel = fromMesh.fromMetaData_.locally_owned_part();
      stk::mesh::Selector sel = fromMesh.fromMetaData_.universal_part();
      sel = sel & stk::mesh::selectUnion(fromMesh.fromParts_);

      std::vector<stk::mesh::Entity>  vecShells;
      stk::mesh::get_selected_entities(sel , fromMesh.fromBulkData_.buckets(fromMesh.fromMetaData_.side_rank()), vecFaces);
      stk::mesh::get_selected_entities(sel , fromMesh.fromBulkData_.buckets(stk::topology::ELEMENT_RANK), vecShells);
      vecFaces.insert(vecFaces.end(), vecShells.begin(), vecShells.end());
    }

    template<class FM>
    static void
    modify_bounding_box(const FM& fromMeshIn,
                        stk::search::Point<double>& min_corner,
                        stk::search::Point<double>& max_corner)
    {
      const GPFromMesh& fromMesh = static_cast<const GPFromMesh& >(fromMeshIn);
      //VERIFY_OP_ON(&fromMesh, !=, 0, "bad cast");

      unsigned nDim = fromMesh.fromMetaData_.spatial_dimension();
      double centroid[3] = {0,0,0};
      double maxDiam = 0.0;
      for (unsigned iDim = 0; iDim < nDim; iDim++)
        {
          centroid[iDim] = 0.5*(min_corner[iDim] + max_corner[iDim]);
          maxDiam = std::max(maxDiam, max_corner[iDim] - min_corner[iDim]);
        }
      maxDiam *= fromMesh.m_inflationFactor;
      for (unsigned iDim = 0; iDim < nDim; iDim++)
        {
          min_corner[iDim] = centroid[iDim] - 0.5*maxDiam;
          max_corner[iDim] = centroid[iDim] + 0.5*maxDiam;
        }
    }

  };

} // namespace percept

#endif
