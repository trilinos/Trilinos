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

#define STK_ADAPT_COLORER_STORED_ENTITYID 0

#define STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR 1

/// define only one of these to be 1
#define STK_ADAPT_COLORER_SET_TYPE_BOOST 1

#if STK_ADAPT_COLORER_SET_TYPE_BOOST
#include <boost/unordered_set.hpp>
#endif


  namespace percept {

#if STK_ADAPT_COLORER_STORED_ENTITYID
    typedef stk::mesh::EntityId ColorerStoredEntity;
#else
    typedef stk::mesh::Entity ColorerStoredEntity;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef std::vector<ColorerStoredEntity> ColorerSetType;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_BOOST
#  if !STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef boost::unordered_set<ColorerStoredEntity> ColorerSetType;
#  endif
    typedef boost::unordered_set<stk::mesh::EntityId> ColorerNodeSetType;
    typedef boost::unordered_set<stk::mesh::EntityId> ColorerElementSetType;
#endif

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
