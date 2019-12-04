// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_Searcher_hpp
#define stk_encr_Searcher_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <percept/function/Function.hpp>

  namespace percept
  {

    class Searcher
    {
    public:

      /** Find the element containing this physical point and return if found (also set the found_it flag to 1, else 0).
       *  If hint_element is non-null, use it to check first if it contains the point to potentially avoid a more costly search.
       *
       *  Dimensions of input_phy_points = ([P]=1, [D])
       *  Dimensions of found_parametric_coordinates = ([P]=1, [D])
       */
      virtual const stk::mesh::Entity findElement(MDArray& input_phy_points, MDArray& found_parametric_coordinates,
						  unsigned& found_it, const stk::mesh::Entity hint_element )=0;
      virtual void setupSearch() {}
      virtual void tearDownSearch() {}
      virtual ~Searcher() {}
    };
  }

#endif
