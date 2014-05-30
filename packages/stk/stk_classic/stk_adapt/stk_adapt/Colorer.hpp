#ifndef stk_adapt_Colorer_hpp
#define stk_adapt_Colorer_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>
#include <stk_percept/stk_mesh.hpp>

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>

#include <stk_util/parallel/Parallel.hpp>

#define STK_ADAPT_COLORER_STORED_ENTITYID 0

#define STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR 1

/// define only one of these to be 1
#define STK_ADAPT_COLORER_SET_TYPE_BOOST 1
#define STK_ADAPT_COLORER_SET_TYPE_TR1 0
#define STK_ADAPT_COLORER_SET_TYPE_STD 0


#if STK_ADAPT_COLORER_SET_TYPE_BOOST
#include <boost/unordered_set.hpp>
#endif

#if STK_ADAPT_COLORER_SET_TYPE_STD
#include <set>
#endif

#if STK_ADAPT_COLORER_SET_TYPE_TR1
#include <tr1/unordered_set>
#endif

namespace stk_classic {
  namespace adapt {

#if STK_ADAPT_COLORER_STORED_ENTITYID
    typedef mesh::EntityId ColorerStoredEntity;
#else
    typedef mesh::Entity *ColorerStoredEntity;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef std::vector<ColorerStoredEntity> ColorerSetType;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_BOOST
#  if !STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef boost::unordered_set<ColorerStoredEntity> ColorerSetType;
#  endif
    typedef boost::unordered_set<mesh::EntityId> ColorerNodeSetType;
    typedef boost::unordered_set<mesh::EntityId> ColorerElementSetType;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_TR1
#  if !STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef tr1::unordered_set<ColorerStoredEntity> ColorerSetType;
#  endif
    typedef tr1::unordered_set<mesh::EntityId> ColorerNodeSetType;
    typedef tr1::unordered_set<mesh::EntityId> ColorerElementSetType;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_STD
#  if !STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef std::set<ColorerStoredEntity> ColorerSetType;
#  endif
    typedef std::set<mesh::EntityId> ColorerNodeSetType;
    typedef std::set<mesh::EntityId> ColorerElementSetType;
#endif


    class Colorer
    {
    public:

      Colorer(std::vector< ColorerSetType >& element_colors, std::vector<mesh::EntityRank> ranks);
      Colorer(std::vector<mesh::EntityRank> ranks );
      void color(percept::PerceptMesh& eMesh, unsigned *elementType = 0,  mesh::PartVector* fromParts = 0, mesh::FieldBase *element_color_field=0);

      /// Set to true to avoid the coloring step and just return elements in a single color
      void setNoColoring(bool no_coloring);
      bool getNoColoring();

      std::vector< ColorerSetType >& getElementColors();

    private:
      std::vector< ColorerSetType >& m_element_colors;
      std::vector< ColorerSetType > m_element_colors_internal;
      std::vector<mesh::EntityRank> m_entityRanks;
      bool m_noColoring;
    };

  }
}

#endif
