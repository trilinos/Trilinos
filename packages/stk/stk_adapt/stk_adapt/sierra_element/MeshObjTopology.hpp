/**-------------------------------------------------------------------*
 *    Copyright 1999 - 2009 Sandia Corporation.                       *
 *    Under the terms of Contract DE-AC04-94AL85000, there is a       *
 *    non-exclusive license for use of this work by or on behalf      *
 *    of the U.S. Government.  Export of this program may require     *
 *    a license from the United States Government.                    *
 *--------------------------------------------------------------------*/

#ifndef stk_adapt_sierra_element_MeshObjTopology_hpp
#define stk_adapt_sierra_element_MeshObjTopology_hpp

#include <stk_adapt/sierra_element/CellTopology.hpp>

namespace stk {
  namespace adapt {
    namespace Elem {

      class MeshObjTopology
      {
      public:
        MeshObjTopology(
                        const CellTopologyData *      cell_topology_data)
          : m_cellTopology(cell_topology_data)
        {}
  
      public:
        ~MeshObjTopology()
        {}
      
        Elem::CellTopology getCellTopology() const {
          return m_cellTopology;
        }
  
      protected:
        Elem::CellTopology                    m_cellTopology;
  
      private:
        MeshObjTopology();
        MeshObjTopology(const MeshObjTopology &);
        MeshObjTopology & operator = (const MeshObjTopology &);
      };

    } // namespace Elem
  } // namespace adapt
} // namespace stk

#endif // stk_adapt_sierra_element_MeshObjTopology_hpp
