// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_Colorer_hpp
#define adapt_Colorer_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>
#include <percept/stk_mesh.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <unordered_set>

  namespace percept {

    typedef stk::mesh::Entity ColorerStoredEntity;

    typedef std::vector<ColorerStoredEntity> ColorerSetType;

    typedef std::unordered_set<stk::mesh::EntityId> ColorerNodeSetType;
    typedef std::unordered_set<stk::mesh::EntityId> ColorerElementSetType;

    class Colorer
    {
    public:

      Colorer(std::vector< ColorerSetType >& element_colors, std::vector<stk::mesh::EntityRank> ranks);
      Colorer(std::vector<stk::mesh::EntityRank> ranks );
      void color(percept::PerceptMesh& eMesh, unsigned *elementType = 0,  stk::mesh::PartVector* fromParts = 0, stk::mesh::FieldBase *element_color_field=0);

      /// Set to true to avoid the coloring step and just return elements in a single color
      void setNoColoring(bool no_coloring);
      bool getNoColoring();

      std::vector< ColorerSetType >& getElementColors();

    private:
      std::vector< ColorerSetType >& m_element_colors;
      std::vector< ColorerSetType > m_element_colors_internal;
      std::vector<stk::mesh::EntityRank> m_entityRanks;
      bool m_noColoring;
    };

  }

#endif
