// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_BuildBoundingBoxes_hpp
#define percept_BuildBoundingBoxes_hpp

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <percept/function/ElementOp.hpp>
#include <percept/FieldTypes.hpp>

#define EXTRA_PRINT 0

  namespace percept
  {

    class BuildBoundingBoxes : public ElementOp
    {

    public:
      typedef stk::search::IdentProc<uint64_t,unsigned> IdentProc;
      typedef stk::search::Point<double> Point;
      typedef stk::search::Box<double> Box;
      typedef std::pair<Point,IdentProc> BoundingPoint;
      typedef std::pair<Box,IdentProc> AABoundingBox;

      std::vector<AABoundingBox>& m_boxes;
      CoordinatesFieldType *m_coords_field;
      bool m_notInitialized;
    public:
      BuildBoundingBoxes(std::vector<AABoundingBox>& boxes, CoordinatesFieldType *coords_field) :  m_boxes(boxes), m_coords_field(coords_field),
                                                                                              m_notInitialized(false)
      {
      }

      virtual ~BuildBoundingBoxes() {
    	  m_coords_field = NULL;
    	  m_notInitialized = true;
      }

      void init_elementOp()
      {
      }
      void fini_elementOp()
      {
        m_notInitialized=true;  // force this object to be used only once
      }
      bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase */*field*/,  const stk::mesh::BulkData& bulkData)
      {
        if (m_notInitialized)
          throw std::runtime_error("BuildBoundingBoxes::operator(): you must re-construct this object before reusing it");

        AABoundingBox bb = getBoundingBox(element, bulkData);
        if (EXTRA_PRINT) std::cout << "bb = (" << bb.first << "," << bb.second << ")" << std::endl;
        m_boxes.push_back(bb);
        return false;  // never break out of the enclosing loop
      }

      AABoundingBox getBoundingBox(const stk::mesh::Entity element, const stk::mesh::BulkData& bulkData)
      {
        const unsigned SpatialDim = bulkData.mesh_meta_data().spatial_dimension();
        Point min_corner, max_corner;
        const MyPairIterRelation elem_nodes(bulkData, element, stk::topology::NODE_RANK );
        unsigned numNodes = elem_nodes.size();
        for (unsigned iNode = 0; iNode < numNodes; iNode++)
        {
          stk::mesh::Entity node = elem_nodes[iNode].entity();
          double * coord_data = stk::mesh::field_data( *m_coords_field, node);
          if (iNode == 0)
          {
            for (unsigned iDim = 0; iDim < SpatialDim; iDim++)
            {
              min_corner[iDim] = coord_data[iDim];
              max_corner[iDim] = coord_data[iDim];
            }
          }
          else
          {
            for (unsigned iDim = 0; iDim < SpatialDim; iDim++)
            {
              min_corner[iDim] = std::min(min_corner[iDim], coord_data[iDim]);
              max_corner[iDim] = std::max(max_corner[iDim], coord_data[iDim]);
            }
          }
        }
        //DJS: should proc be set to bulk_data.parallel_rank()?
        int proc = 0;
        return AABoundingBox( Box(min_corner,max_corner), IdentProc(bulkData.identifier(element),proc));
      }
    };

  } //namespace percept

#endif
