// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_SimpleSearcher_hpp
#define stk_encr_SimpleSearcher_hpp

#include <cmath>
#include <math.h>
#include <vector>
#include <utility>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/function/Function.hpp>

#include <percept/function/internal/Searcher.hpp>
#include <percept/function/internal/IsInElement.hpp>

#include <percept/FieldTypes.hpp>

  namespace percept
  {
    class FieldFunction;

    class SimpleSearcher : public Searcher
    {
      stk::mesh::BulkData *m_bulk;
    public:
      SimpleSearcher(stk::mesh::BulkData *bulk);

      virtual ~SimpleSearcher() {}

      /**
       *  Dimensions of input_phy_points = ([P]=1, [D])
       *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
       */

      virtual const stk::mesh::Entity findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates,
                                                   unsigned& found_it, const stk::mesh::Entity hint_element ) override;
    };

  }

#endif
