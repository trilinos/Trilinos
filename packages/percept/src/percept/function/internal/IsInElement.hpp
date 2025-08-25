// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_IsInElement_hpp
#define stk_encr_IsInElement_hpp

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
#include <percept/function/ElementOp.hpp>

  namespace percept
  {

    class IsInElement : public ElementOp
    {

    public:
      bool m_found_it;
      MDArray& m_input_phy_points;
      MDArray& m_found_parametric_coordinates;
      stk::mesh::Entity m_foundElement;
    public:
      IsInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates);

      bool operator()(const stk::mesh::Entity element, const stk::mesh::BulkData& bulkData);
      bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase* field, const stk::mesh::BulkData& bulkData) override;
      void init_elementOp() override;
      void fini_elementOp() override;

    private:


      /**
       *  Dimensions of input_phy_points = ([P]=1, [D])
       *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
       */
      void isInElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates, unsigned& found_it, const stk::mesh::Entity element,
                       const stk::mesh::BulkData& bulkData);

    };

  }
#endif
