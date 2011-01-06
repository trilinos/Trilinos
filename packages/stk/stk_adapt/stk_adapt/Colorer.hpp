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
#include <stk_util/environment/ProgramOptions.hpp>

/// define only one of these to be 1
#define SET_TYPE_BOOST 1
#define SET_TYPE_TR1 0
#define SET_TYPE_STD 0


#if SET_TYPE_BOOST
#include <boost/unordered_set.hpp>
#endif

#if SET_TYPE_STD
#include <set>
#endif

#if SET_TYPE_TR1
#include <tr1/unordered_set>
#endif

namespace stk {
  namespace adapt {


#if SET_TYPE_BOOST
    typedef boost::unordered_set<EntityId> ColorerSetType;
#endif

#if SET_TYPE_TR1
    typedef tr1::unordered_set<EntityId> ColorerSetType;
#endif

#if SET_TYPE_STD
    typedef std::set<EntityId> ColorerSetType;
#endif


    class Colorer
    {
    public:


      Colorer(std::vector<EntityRank> ranks = std::vector<EntityRank>()) : m_entityRanks()
      {
        //EntityRank ranks[2] = {Face, Element};
        if (ranks.size())
          {
            m_entityRanks = ranks;
          }
        else
          {
            m_entityRanks.push_back(mesh::Face);
            m_entityRanks.push_back(mesh::Element);
          }
      }
      void color(percept::PerceptMesh& eMesh, unsigned *elementType = 0, FieldBase *element_color_field=0);

      std::vector< ColorerSetType >& getElementColors();

    private:
      std::vector< ColorerSetType > m_element_colors;
      std::vector<EntityRank> m_entityRanks;
    };

  }
}

#endif
