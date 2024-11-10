// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_STKSearcher_hpp
#define percept_STKSearcher_hpp

#include <percept/Percept.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>



#include <Intrepid2_FunctionSpaceTools.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <percept/function/FieldFunction.hpp>

#include <percept/function/internal/Searcher.hpp>
#include <percept/function/internal/BuildBoundingBoxes.hpp>

#include <percept/FieldTypes.hpp>

#define EXTRA_PRINT 0

  namespace percept
  {
    class STKSearcher : public Searcher
    {
      typedef BuildBoundingBoxes BBB;
      typedef BBB::AABoundingBox BBox;
      typedef BBB::BoundingPoint BPoint;
      typedef BBB::IdentProc IdProc;
      typedef std::vector<std::pair<IdProc, IdProc> > IdentProcRelation;

      stk::mesh::BulkData *m_bulk;
      std::vector< BBox > m_boxes;

    public:

      STKSearcher(stk::mesh::BulkData *bulk);

      virtual ~STKSearcher();

      void setupSearch();

      void tearDownSearch();

      /**
       *  Dimensions of input_phy_points = ([P]=1, [D]) 
       *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
       */

      virtual const stk::mesh::Entity findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, 
                                                   unsigned& found_it, const stk::mesh::Entity hint_element );

    private:


    };

  }

#endif
