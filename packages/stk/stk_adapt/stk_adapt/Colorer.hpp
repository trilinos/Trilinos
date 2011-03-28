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

namespace stk {
  namespace adapt {

#if STK_ADAPT_COLORER_STORED_ENTITYID
    typedef EntityId ColorerStoredEntity;
#else
    typedef Entity *ColorerStoredEntity;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef std::vector<ColorerStoredEntity> ColorerSetType;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_BOOST
#  if !STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef boost::unordered_set<ColorerStoredEntity> ColorerSetType;
#  endif
    typedef boost::unordered_set<EntityId> ColorerNodeSetType;
    typedef boost::unordered_set<EntityId> ColorerElementSetType;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_TR1
#  if !STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef tr1::unordered_set<ColorerStoredEntity> ColorerSetType;
#  endif
    typedef tr1::unordered_set<EntityId> ColorerNodeSetType;
    typedef tr1::unordered_set<EntityId> ColorerElementSetType;
#endif

#if STK_ADAPT_COLORER_SET_TYPE_STD
#  if !STK_ADAPT_COLORER_SET_TYPE_USE_VECTOR
    typedef std::set<ColorerStoredEntity> ColorerSetType;
#  endif
    typedef std::set<EntityId> ColorerNodeSetType;
    typedef std::set<EntityId> ColorerElementSetType;
#endif


    class Colorer
    {
    public:


      //Colorer(std::vector< ColorerSetType >& element_colors, std::vector<EntityRank> ranks = std::vector<EntityRank>());
      //Colorer(std::vector<EntityRank> ranks = std::vector<EntityRank>());
      Colorer(std::vector< ColorerSetType >& element_colors, std::vector<EntityRank> ranks);
      Colorer(std::vector<EntityRank> ranks );
      void color(percept::PerceptMesh& eMesh, unsigned *elementType = 0,  PartVector* fromParts = 0, FieldBase *element_color_field=0);

      std::vector< ColorerSetType >& getElementColors();

    private:
      std::vector< ColorerSetType >& m_element_colors;
      std::vector< ColorerSetType > m_element_colors_internal;
      std::vector<EntityRank> m_entityRanks;
    };

  }
}

#endif
